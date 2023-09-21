#SCFDE implementation: paper: Broadband wireless Using Single Carrier and Frequency Domain Equalization by David D Falconer

from .msk import *
from scipy.linalg import  cho_factor, cho_solve
import numba
from pyfftw.interfaces.numpy_fft import fft, ifft

#Expand the frequency domain vector X's last dimension, taking CIR and padding zero to frame_len
def cir_interp(X, frame_len):
	ax = len(X.shape) - 1
	return fft(ifft(X, axis=ax), axis=ax, n=frame_len)

def make_dfe(v_f, fb_len, frame_len):
	if fb_len != 0:
		fb_select = range (1, fb_len + 1)
		v_t = fft(v_f)
		cir_index = np.array([ [ i - j if i >= j else frame_len + i - j  for j in fb_select] for i in fb_select])
		return np.conj(cho_solve(cho_factor(v_t[cir_index]), -v_t[fb_select]))
	else:
		return np.zeros((0,))


@numba.njit
def hard_decision(y):
	o1 = np.ones(y.shape)
	return 2 * o1 * (np.real(y) >= 0) - o1

@numba.njit
def backward_equalize(frame_len, prev_decision, data_frame_fwd, backward_eq):
	equalized = np.zeros((frame_len,), dtype=numba.complex128)
	for j in range(frame_len): # time-domain feedback equalization
		equalized[j] = eqval = data_frame_fwd[j] - prev_decision.astype(numba.complex128) @ backward_eq
		prev_decision = np.roll(prev_decision, 1)
		prev_decision[0] = 1 if np.real(eqval) >= 0 else -1
	return equalized


def frame_equalize(FORWARD_EQ, backward_eq, DATA, pilot):
	''' ! equalize a single frame
	:param FORWARD_EQ: complex (frame_len, SPBCH),
	:param backward_eq: complex (fb_len,)
	:param DATA: frequency domain data, already prescaled and handled the SPB stuffs
	:param pilot: pilot text for decision feedback initially. Will flip so  lower index means more recent sample
	:returns:  equalized complex (frame_len,) ready to take hard decision
	'''
	frame_len, SPBCH = FORWARD_EQ.shape
	assert DATA.shape[0] == frame_len and DATA.shape[1] == SPBCH
	data_frame_fwd = ifft((FORWARD_EQ * DATA).sum(axis=1), axis=0)
	if backward_eq.size > 0:
		return backward_equalize(
			frame_len = frame_len,
			prev_decision =np.flip(pilot[pilot.size-backward_eq.size:]),
			data_frame_fwd = data_frame_fwd,
			backward_eq=backward_eq )
	else: return  data_frame_fwd
def scfde_packetize(pilot, prpt, symbol, payload_len, spb):
	pilot_len, CH, L = pilot.shape[0], symbol.shape[1], symbol.shape[0]
	assert 1 == len(pilot.shape)
	frame_len = payload_len + pilot_len + 1; # FFT length: +1 means the tail symbol
	data_part = symbol[(pilot_len*prpt*spb):, :]
	
	#discard first pilot_len (non-cyclic ISI)
	training_part = rearrange(symbol[pilot_len*spb:(pilot_len*prpt)*spb, :],
			'(pilot_rpt pilot_len_spb) CH -> pilot_rpt CH pilot_len_spb', pilot_len_spb=pilot_len * spb) 
	return training_part, data_part, frame_len

from functools import lru_cache

def indirect_equalizer(h, frame_len, fb_len):
	''' ! for each frequency do MRC: Assumes the prescaler so s2 is implicitly =1
	paper: Broadband, formula (7) (8)
	:param frame_len: desired fft length for final fft equalizer
	:param fb_len: desired dfe taps.
	:returns: (FORWARD_EQ, backward_eq) '''
	
	H = fft(h, n=frame_len, axis=0)
	Hs = np.sum(H * np.conj(H), axis=1)
	eq_gain = 1 / (1 + Hs)# s2 = 1
	backward_eq = make_dfe(eq_gain, fb_len, frame_len)
	BACK_EQ = fft(np.hstack([1, backward_eq]), n=frame_len)
	FORWARD_EQ = np.conj(H) * (BACK_EQ * eq_gain).reshape([-1, 1])
	return FORWARD_EQ,  backward_eq
	


def indirect_channel_training(training_part, pilot):
	''' ! Channel training : H, A -> h, sigma
	:param training_part (pilot_rpt, CH, pilot_len * spb)
	:param pilot (pilot_len,) the time_domain pilot waveform
	:returns: (h, prescaler)
		prescaler: float (1, SPBCH): scaling each channel such that the s2 equals to 1 ,
			which is essential for MRC combining, leveling each channel's noise floor
		s2 is average variance, calibrate by pilot_rpt/(pilot_rpt-1), must=1 after prescaling
		h (pilot_len, SPBCH) the averaged CIR, operating AFTER the prescaler ''' 

	pilot_rpt, CH, pilot_len_spb = training_part.shape
	pilot_len = pilot.shape[0]
	assert pilot_len_spb % pilot_len == 0
	spb = pilot_len_spb // pilot_len
	SPBCH = CH * spb


	A = fft(pilot)
	X = rearrange(fft(training_part, axis=2), 'pilot_rpt CH (spb pilot_len) -> pilot_rpt (spb CH) pilot_len', spb=spb)
	Xm = np.sum(X, axis=0) / pilot_rpt; # (SPBCH, pilot_len)
	div = X - Xm.reshape([1, SPBCH, pilot_len])
	s2 = np.sum((div * np.conj(div)).real, axis=(0, 2)) / (pilot_len * pilot_len * (pilot_rpt - 1) )
	# one pilot_len from averaging, another from the FFT gain.

	prescaler = (1 / np.sqrt(s2)).reshape([1, SPBCH])
	h = ifft(Xm / A.reshape([1, -1]), axis=1) * prescaler.T; # (SPBCH, pilot_len)
	return h.T, prescaler

def indirect_frame_grad(h,  DATA, decision):
	''' ! LMS update of one packet, returns gradiant to the CIR
	:param h: the original CIR shape=[pilot_len, SPBCH]
	:param pilot_len desired dof to limit to, integer
	:param DATA: current scaled data frame to take grad shape=[frame_len, SPBCH]
	:param decision: decision of data (1 or -1) to reference shape=[frame_len,]
	:returns: shape=[pilot_len, CH], h+=grad*learn_rate; to make indirect_equalizer() for next frame
	'''
	frame_len, SPBCH = DATA.shape
	pilot_len = h.shape[0]
	assert len(h.shape) == 2 and h.shape[1] == SPBCH
	assert len(decision.shape) == 1 and decision.shape[0] == frame_len 

	C = fft(decision).reshape([frame_len, 1])
	H = fft(h, axis=0, n=frame_len) # [frame_len, SPBCH]
	err = np.conj(C) * (DATA - H * C)
	return ifft(err, axis=0)[:pilot_len, :]

@lru_cache
def dftmat(M, icols, period = None):
	if period is None:
		period = M
	mn = (np.arange(icols) % period).reshape([1, -1]) * (np.arange(M) % period).reshape([-1,1])
	return np.exp(-2j*np.pi * mn / period)
def RLS_R(decision, frame_len, pilot_len):
	C = fft(decision)
	C2 = (np.conj(C) * C).reshape([-1, 1]); # must match the grad's coefficient about h
	return ifft(C2 * dftmat(frame_len, pilot_len), axis=0)[:pilot_len, :]

class Equalizer:
	def calc(self):
		self.FORWARD_EQ, self.backward_eq = indirect_equalizer(self.h, self.frame_len, self.fb_len)
	def __init__(self, h, frame_len,  pilot, fb_len, LR):
		self.backward_eq= self.FORWARD_EQ = None
		self.frame_len = frame_len; self.pilot = pilot;  self.fb_len = fb_len; self.LR = LR
		self.pilot_len =  pilot.size
		self.h= h
		if LR < 0:
			self.R = RLS_R(pilot, self.pilot_len, self.pilot_len)
		self.backward_eq= self.FORWARD_EQ = None; self.calc()


	def update(self, DATA, payload_decision):
		decision = np.hstack([payload_decision, [1], self.pilot])
		grad = indirect_frame_grad(self.h, DATA, decision)
		if self.LR >= 0:
			self.h += grad * self.LR / self.frame_len
		else:
			self.R = self.R * (-self.LR) + RLS_R(decision, self.frame_len, self.pilot_len)
			self.h += cho_solve(cho_factor(self.R), grad)
		self.calc()
		return self
	def __call__(self, DATA):
		return frame_equalize(self.FORWARD_EQ, self.backward_eq, DATA, self.pilot)


def indirect_equalization_main(pilot, prpt, symbol,
		payload_len, frame_cnt, fb_len, frame_ecc, LR, spb, eq_iter):
	'''! indirection equalization mian
	:param pilot: pilot sequence {+-1}
	:param prpt: pilot repetition encoded, will use (prpt-1) for channel estimation
	:param symbol the decimated symbol, one symbol per bit, fractional sampling is implemented now as multiple channels
		instead, it is not the right way, instead it should channelize AFTER the FFT
	:param payload_len
	:param frame_cnt
	:param fb_len: number of dfe feedback taps
	:param frame_ecc: callback function, pass the decoded bits and return the corrected bits (both in +-1 metric)
		to correct the payload that will be used for calculating the gradient of channel estimation. For preliminary testing we can feed ground-truth data bits
	:param LR: Positive value: learning rate using the LMS algorithm
	Negative value: the forgetting factor using the RLS algorithm
	:param spb: sample per bit
	:param eq_iter: number of iterations (uses the decision to refine channel estimation then decide again)
	'''

	training_part, data_part, frame_len = scfde_packetize(pilot, prpt, symbol, payload_len, spb)
	h0, prescaler = indirect_channel_training(training_part, pilot)
	print(f'prescaler={prescaler}, hrms={cnorm(h0.ravel())/h0.size}')
	E = Equalizer(h0, frame_len, pilot, fb_len, LR)
	#return h0, s2, data_part
	def make_cfo_target(h):
		return np.conj(fft(h,axis=0))

	cfo_target = make_cfo_target(E.h)
	equalized = np.zeros((eq_iter, frame_cnt, frame_len), dtype=complex)
	ce = np.zeros((frame_cnt, frame_len), dtype=complex)
	packet_begin = 0

	pilot_len =  pilot.size
	packet_begin_all = np.zeros((frame_cnt,))
	for i in range(frame_cnt) :
		packet_begin_all[i] = packet_begin
		data = data_part[spb*packet_begin:(packet_begin + frame_len)*spb, :]
		DATA = rearrange(fft(data, axis = 0), '(spb frame_len) CH -> frame_len (spb CH)', spb=spb) * prescaler
		payload_decision = None
		for j in range(eq_iter):
			if payload_decision is None:
				Ec = E
			else:
				Ec = copy.deepcopy(E).update(DATA, payload_decision)
			equalized[j, i,:] = Ec(DATA)
			payload_decision = hard_decision(equalized[j, i, :payload_len])
		decision_ecc =frame_ecc(payload_decision)
		Ec = copy.deepcopy(E).update(DATA, decision_ecc)
		ce[i,:] = Ec(DATA)

		corr = np.abs(ifft(cfo_target * fft(Ec.h, axis=0), axis=0)).sum(axis=1); # sum over all channels,
		k = (np.argmax(corr) + pilot_len//2) % pilot_len - pilot_len//2
		if k == 1 or k == -1:
			print(f'new cfo target k={k} on{i}-th frame')
			packet_begin += k
			E.h = np.roll(E.h, -k, axis=0)
			data = data_part[spb*packet_begin:(packet_begin + frame_len)*spb, :]
			DATA = rearrange(fft(data, axis = 0), '(spb frame_len) CH -> frame_len (spb CH)', spb=spb) * prescaler
			E.update(DATA, decision_ecc)
			cfo_target = make_cfo_target(E.h)
			packet_begin += frame_len
		else:
			if abs(k) != 0:
				print(f"WARNING: Channel continuity fail, k={k} on {i}-th frame")
			packet_begin += frame_len
			E = Ec

	return equalized, ce, packet_begin_all
def scfde_unpacketize(payload_len, pilot, prpt, data):
	# used in simulations, generates tail bit 1 and put together the pilots, ready to convolve with time varying channel
	assert data.size % payload_len == 0
	packets = data.size // payload_len
	data_part = np.hstack([data.reshape([packets, payload_len]), np.tile([1, *pilot], (packets, 1))]).ravel()
	training_part = np.tile(pilot, (prpt,))
	return np.hstack([training_part, data_part])
