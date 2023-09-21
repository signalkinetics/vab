# Paper: Stojanovic, Milica, Josko Catipovic, and John G. Proakis. "Adaptive multichannel combining and equalization for underwater acoustic communications." The Journal of the Acoustical Society of America 94.3 (1993): 1621-1631.
# Use the RLS rule and direct weight update


from msk import *
from linregs import awgn
from scipy.linalg import cho_factor, cho_solve, pinvh
from plots import *

# For samping offset fast convolution of forward taps
from pyfftw.interfaces.numpy_fft import fft, ifft
from pyfftw import next_fast_len



import numba


class Basic_RLS():
	'''
	:param lam: forgetting factor
	:param P: Current Kalman uncertainty matrix, the inverse of Ruu
	:param c: Current equalizer coefficient (conjugate), generating output =c.conjugate() @ u
	'''
	def __init__(self, u_all, d_all, lam, prev_RLS):
		''' Constructor (prev=None) or batch-update (when choosing new CFO targets)
		:param lam: forgetting factor
		:param u_all: All input data for training, shape = (step, ch_all), already wrapped by DFE_Equalizer
		:param d_all: ground truth decisions for training symbols, shape = (step,), {-1,+1} for BPSK
		'''
		step, ch = u_all.shape
		assert d_all.size == step
		self.lam = lam
		loglam = np.log(lam)
		atten = np.flip(np.exp(loglam * np.arange(step))).reshape([step, 1])
		u_all_atten = u_all.astype(np.complex128) * atten
		Rud = d_all.conjugate() @ u_all_atten  # (step,) @ (step, ch) -> (ch,)
		Ruu = u_all.T @ u_all_atten.conjugate()  # (ch, step) @ (step, ch) -> (ch, ch)
		atten_d = np.exp(loglam * step)
		if prev_RLS is not None:
			prev_Ruu = pinvh(prev_RLS.P) * atten_d
			Ruu += prev_Ruu
			Rud += prev_Ruu @ prev_RLS.c
		else:
			Ruu += np.eye(ch) * atten_d
		self.P, rank = pinvh(Ruu, return_rank=True)
		if rank < ch:
			raise AssertionError(f'expect rank = ch = {ch}, but rank is {rank}, d_all.size={d_all.size}')
		self.c = self.P @ Rud  # (ch, ch) @ (ch,) -> (ch,)

	@staticmethod
	@numba.njit
	def __step(c, P, lam, u):
		# Matlab-style RLS inverse  See the math documented in rls.md
		ds = c.conjugate() @ u
		d = 1 if ds.real > 0.0 else -1
		e = ds - d
		pu = P @ u
		upu = u.conjugate() @ pu
		k = pu * (1 / (lam + upu))
		c1 = c -  e.conjugate() * k
		kup = k.reshape((-1, 1)) * pu.conjugate()
		P1 = (P - kup) * (1 / lam)
		return c1, P1, d, e

	def step(self, u):
		self.c, self.P, d, e = self.__step(self.c, self.P, self.lam, u)
		return d, e


import ctypes as ct
import os
dn = os.path.dirname(os.path.realpath(__file__))
accel = ct.CDLL(os.path.join(dn, 'linreg_mc.so'))
accel.new_ldlt_rls.argtypes = (ct.c_int, ct.c_int, ct.POINTER(ct.c_double),  ct.POINTER(ct.c_double), ct.c_double, ct.c_void_p)
accel.new_ldlt_rls.restype = ct.c_void_p
accel.free_ldlt_rls.argtypes = (ct.c_void_p, )
accel.free_ldlt_rls.restype = None
accel.deepcopy_ldlt_rls.argtypes = (ct.c_void_p,)
accel.deepcopy_ldlt_rls.restype= ct.c_void_p
accel.ldlt_rls_get_c.argtypes = (ct.c_void_p, ct.POINTER(ct.c_double))
accel.ldlt_rls_get_c.restype = None
accel.ldlt_rls_update.argtypes = (ct.c_void_p, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double))
accel.ldlt_rls_update.restype = None

def dptr(array):
	return array.ctypes.data_as(ct.POINTER(ct.c_double))

class LDLT_RLS(): # binding of xxx_ldlt_rls functions in linreg_mc.cpp
	def __deepcopy__(self, memo):  
		return LDLT_RLS(None, None, None, self)
	def __copy__(self):
		raise ValueError("Copy isn't allowed, which would cause double free of the pointer. Use deepcopy instead")
	def __del__(self):
		accel.free_ldlt_rls(self.ptr)
	def __init__(self, u_all, d_all, lam, prev_RLS):
		if u_all is None:
			assert prev_RLS is not None
			self.ptr = accel.deepcopy_ldlt_rls(prev_RLS.ptr)
			self.ch = prev_RLS.ch
		else:
			self.ch = u_all.shape[1]
			d_all_c1 = d_all.astype(np.complex128)
			self.ptr = accel.new_ldlt_rls(self.ch, u_all.shape[0], dptr(u_all), dptr(d_all_c1), lam, None if prev_RLS is None else prev_RLS.ptr)
	@property
	def c(self):
		c = np.zeros((self.ch,), dtype=np.complex128)
		accel.ldlt_rls_get_c(self.ptr, dptr(c))
		return c

	def step(self, u):
		assert u.dtype == np.complex128 and u.size == self.ch
		ds = self.c.conjugate() @ u
		d = 1 if ds.real > 0.0 else -1
		e = np.array([ds - d])
		accel.ldlt_rls_update(self.ptr, dptr(u), dptr(e))
		return d, e.item()

class DFE_Equalizer():
	''' RLS based DFE_Equalizer, essentially a wrapper of Basic_RLS that generates the concatenated inputs. 
			Make Toeplitz shaped (formula 24) and concatenate across channels (formula 32).
	:param past_d: shape(bw_len,)
	:param past_sp: shape(fw_len, spb, ch) '''
	def shift_sp(self, sp):
		''' :param sp: To avoid numerical issue, caller should scale globally the input samples to RMS of 1!
		:return: concatenated u vector, windowed sample of (fw_len, spb, ch) first, followed by decision (bw_len)
		'''
		self.past_sp = np.concatenate((self.past_sp[1:, :, :],  sp[np.newaxis]), axis=0)
		shifted_sp =np.hstack((self.past_sp.ravel(),  self.past_d))
		#return shifted_sp;
		return shifted_sp
	def get_fwdeq(self):
		''' :return: time domain of forward taps [fw_len * spb, ch] to inspecst correlation between cfo targets.
		'''
		return self.RLS.c[:self.fw_len * self.spb * self.ch].reshape([self.fw_len * self.spb, self.ch]);
	def shift_d(self, d):
		''' Update the decision part of state with new decision d '''
		self.past_d = np.hstack((self.past_d[1:], d))

	def __init__(self, spb, lam, fw_len, bw_len, d_all, _train_sp, prev = None, impl = LDLT_RLS):
		""" Constructor (prev=None), will append zeros assume the first symbols has no ISI.
		Or an updater when there is a 'prev' structure that tells about the past_sp and past_d, when changing the
		sampling offset target and retrain from the middle.

		:param spb: sample per bit
		:param lam: forgetting factor
		:param fw_len: forward taps
		:param bw_len: feedback taps
		:param d_all: BPSK {-1,+1} bits of training part shape=(train_steps,)
		:param _train_sp: Training data samples of shape((train_steps) * spb, ch) """
		train_sp= rearrange(_train_sp, '(s spb) ch -> s spb ch', spb = spb)
		train_steps, self.spb, self.ch = train_sp.shape
		self.fw_len = fw_len
		assert d_all.size == train_steps
		if prev is None:
			self.past_d = np.zeros((bw_len,), dtype=np.complex128)
			self.past_sp = np.zeros((fw_len, spb, self.ch), dtype=np.complex128)
		else:
			self.past_d = prev.past_d.copy()
			self.past_sp = prev.past_sp.copy()

		u_all = np.zeros((train_steps, fw_len * spb * self.ch + bw_len), dtype=np.complex128)
		for i in range(train_steps):
			u_all[i, :] = self.shift_sp(train_sp[i, :, :])
			self.shift_d(d_all[i])
		self.RLS = impl(u_all, d_all, lam, None if prev is None else prev.RLS)

	def step(self, new_sp):
		u = self.shift_sp(new_sp.reshape([self.spb, -1]))
		d, e = self.RLS.step(u)
		self.shift_d(d)
		return d, e

def dfe_sfo_main(sample, loop, spb, lam, fw_len, bw_len, train_gt, impl, enable_sfo = False):
	'''
	Working interface is all sample segment + spb, doesn't have the step function because it read the samples twice, just as the scfde_main
	Wrapper of DFE_Equalizer, save a spare deepcopy instance of it, as the CFO target.
	:param sample (bits*spb, CH) train+decide bpsk symbols
	:param spb lam fw_len bw_len impl: same as DFE_Equalizer input
	:param loop: desired steps (bits) after the training that needs to be recovered.


	Then detect the sample offset at some interval (now: every bit) by perform a correlation (convolution of the conjugate flipped)
	of the current  forward tap fw[t2,:] with the backup fw[t1,:]. Ignore the >2 sample offset. If find 1-sample offset then,
	Discard the current equalizer at t2 with incorrect sampling offset:
		Make new equalizer from the EQ[t1]  train data with new offset till t2. We can't time-shift the coefficients directly because 1) there is missing 
		taps 2) we can't shift the covariance matrix. So retrain is the easiest.

	Init the DFE outside and pass the instance here? No.  It should destroy the unused DFE instance so it should keep the referene count=1

	'''
	train_len = train_gt.size

	target_len = next_fast_len(fw_len * spb * 2 - 1)
	fft1 = lambda _fw:  fft(_fw, axis=0, n=target_len)
	if not np.allclose([-1, 1], np.unique(train_gt)):
		raise ValueError('train_gt must be bpsk format but we find value other than [-1, +1]');

	def make_cfo_target(train_gt, train_samples, prev):
		eq = DFE_Equalizer(spb, lam, fw_len, bw_len, train_gt, train_samples, prev, impl)
		eq_b = copy.deepcopy(eq)
		cfo_target = fft1(np.flip(eq_b.get_fwdeq().conjugate(), axis=0))
		return eq, eq_b, cfo_target
	train_samples = sample[:train_len * spb, : ] # reference
	#conditioner: make train_samples to have rms of 1
	sample *= np.sqrt(train_samples.size / cnorm(train_samples.ravel()))
	eq, eq_b, cfo_target = make_cfo_target(train_gt=train_gt, train_samples=train_samples, prev=None)
	next_sample = train_len * spb
	last_sf_i = 0
	d_all = np.zeros((loop,))
	e_all = np.zeros((loop,), np.complex128)
	for i in range(loop):
		d_all[i], e_all[i] = eq.step(sample[next_sample:next_sample+spb, :]);
		cfo_correlation = np.abs(ifft(fft1(eq.get_fwdeq()) * cfo_target, axis=0)).sum(axis=1);
		k = np.argmax(cfo_correlation) - spb * fw_len + 1

		if i > last_sf_i + train_len and (k == 1 or k == -1) and enable_sfo:
			print(f'new cfo target k={k} on i={i} bit' 
					+ 'rx data is %s than reader clock'%('faster' if k < 0 else 'slower'));
			next_sample += spb + k;
			eq, eq_b, cfo_target = make_cfo_target(
				train_samples =  sample[next_sample - spb * (i - last_sf_i):next_sample, :],
				train_gt = d_all[last_sf_i + 1:i + 1],
				prev= eq_b);
			
			last_sf_i = i
		else:
			next_sample += spb;
			if i > last_sf_i + train_len and abs(k) != 0 and False:
				print(f'WARNING: cfo continuity fail, k={k} on i={i}-th, last_sf_i={last_sf_i}')

	return d_all, e_all, eq


def equalization_test():
	train_len = 400
	fw_len=4
	bw_len=4

	lam = 0.997
	spb=1
	ch = 8
	source = np.random.randint(2, size=6666)* 2 - 1;
	source_isi = awgn(np.tile(convolve(source, [1,0.8-0.3j]).reshape((-1,1)), (1,ch)),20)
	eq = DFE_Equalizer(spb, lam, fw_len, bw_len, source[:train_len], source_isi[:train_len*spb,:])
	rms = 0; loop = source.size - train_len
	for i in range(loop):
		d, e = eq.step(source_isi[(i + train_len)*spb:(i + train_len+1)*spb:, :])
		rms += e.real**2+e.imag**2
	return rms/loop*ch
