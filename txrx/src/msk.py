import numpy as np
import copy
from einops import rearrange
from scipy.signal import convolve, remez, lfilter

def cnorm(vec):
	if (vec.size == max(vec.shape)):
		vec = vec.ravel()
		return np.real(vec @ np.conj(vec.T))
	raise ValueError(f'must input a vector, but input has shape {vec.shape}')
def half_sine(spb):
	return np.sin(np.arange(spb * 2) / (2 * spb) * np.pi).reshape([-1,1])

_c0, _c1 = -np.roots([1, np.pi, 1])
_square_offset=0.633-0.263j; # from linreg_mc
_square_gain = 1.2546+0.3871j
#using Laurent's linear representation.
def fm0_laurent_samples(bits):
	bits = np.hstack((bits, np.logical_xor.reduce(bits))); # padding the end condition
	n = bits.shape[0]
	bits_c = np.logical_xor(bits, np.arange(n) % 2)
	bits_r = np.hstack([1, bits_c]).cumsum() % 2
	mod = (bits_r * 2 - 1) * np.exp(0.5j * np.pi * (np.arange(n + 1) % 2))
	prs = convolve(mod, [_c1, 1])[1:-1]
	return prs

def fm0_encode(bits, spb): #  spb means sample per bit
	bits = np.hstack((bits, np.logical_xor.reduce(bits))); # padding the end condition
	n = bits.shape[0]
	is_balanced = n % 2 == 1
	nc = spb * (n + 2)
	subcar = np.exp(np.arange(nc) / spb * np.pi * 1.5j)# starting from -T

	mask = half_sine(spb).reshape((1, -1))

	bits_c = np.logical_xor(bits, np.arange(n) % 2)
	bits_r = np.hstack([1, bits_c]).cumsum() % 2
	bpsk = (bits_r * 2 - 1).reshape((-1, 1))

	pad = ((0, is_balanced * spb),(spb, (not is_balanced) * spb))
	real_subcar = np.real(subcar)[pad[0][0]:nc - pad[0][1]].reshape((-1, 2 * spb))
	imag_subcar = np.imag(subcar)[pad[1][0]:nc - pad[1][1]].reshape((-1, 2 * spb))
	real_part =  np.pad((bpsk[0::2, :] * mask * real_subcar).reshape((-1,)), (pad[0],))
	imag_part =  np.pad((bpsk[1::2, :] * mask * imag_subcar).reshape((-1,)), (pad[1],))
	modulated = real_part - 1j * imag_part
	return np.real(modulated) + np.imag(modulated)

def wrap_square(is_square, wf, spb):
	def branch_metric(samples, i, scale,  ref):
		deviation = (samples[i, :] @ scale.T).item()  - ref
		scalegrad = deviation * np.conj(samples[i, :])
		return np.real(deviation * np.conj(deviation)), scalegrad

	def branch_metric_square(samples, i, scale,  ref): #return bm, and d(err)/d(scale)
		return branch_metric(samples, i, scale, ref * _square_gain - _square_offset * np.power(-_c1, i))
	if (is_square):
		return branch_metric_square,  np.pad(wf, ((spb,spb), (0,0)))
	else:
		return branch_metric, wf
def fm0_ideal_pm_lms(wf, spb, bits, LR0, scale, is_square=False): 
	branch_metric, wf = wrap_square(is_square, wf, spb)
	samples = fm0_bb_decimate(wf, spb)
	scale=copy.deepcopy(scale)
	s_ref = fm0_laurent_samples(bits)
	assert len(scale.shape) == 1 # scale is a row vector
	assert samples.shape[1] == scale.shape[0] 

	n = samples.shape[0]
	LR = LR0 * cnorm(scale) / scale.shape[0]
	bm = np.zeros((n,))
	for i in range(n):
		bm[i], grad = branch_metric(samples, i, scale, s_ref[i])
		#correct learning rate for different signal level, potentially unstable
		scale += -LR * grad
	pm_grow = np.cumsum(bm)

	return {'pm_grow' : pm_grow, 'snr' : cnorm(s_ref) / pm_grow[-1] }


## wf = (..complex_preamble_matching...())
# fm0_decode(wf) is upper side band
# fm0_decode(conj(wf)) is lower side band,
# fm0_decode(real(wf)) is double side band
def fm0_decode(wf, spb, LR0, scale, is_square=False, scale_override=None): 
	branch_metric, wf = wrap_square(is_square, wf, spb)
	samples = fm0_bb_decimate(wf, spb)
	assert len(scale.shape) == 1 # scale is a row vector
	assert samples.shape[1] == scale.shape[0] 
	n = samples.shape[0] + 1
	pm = np.zeros((n + 1, 2))
	tr = np.zeros((n + 1, 2), dtype=np.int8)
	dp_scale = np.zeros((n + 1, 2, scale.shape[0]), dtype=scale.dtype)
	LR = LR0 *  cnorm(scale) / scale.shape[0]

	ref = np.array([[-1-1j*_c1, -1 + 1j*_c1, 1 - 1j * _c1, 1 + 1j * _c1]])
	#mapping rule:(symbol,bits_r) =  ((1,1),(1j,1),(-1,0), (-1j,0), )
	refj = np.vstack((ref, 1j * np.conj(ref))); # [imag, bit_prev << 1 | bit_now]
	
	# i -> i + 1
	def update(i, j):  #new_pm, j>>1, new_scale
		bm, grad = branch_metric(samples, i, dp_scale[i, j >> 1, :],  refj[i & 1, j])
		return pm[i][j >> 1] + bm, j >> 1, dp_scale[i, j >> 1, :] + -LR * grad

	pm[0, :] = [ np.inf, 0 ]
	for i in range(2):
		dp_scale[0, i, :] = scale
	for i in range(0, n - 1): #forward
		for j in range(2):
			v0, v2 = update(i, j | 0), update(i, j | 2)
			pm[i + 1, j], tr[i + 1, j], dp_scale[i + 1, j, :] = v0 if v0[0] < v2[0] else v2

	bits_r = np.zeros(n, dtype='bool')
	j = 0;   fb = int(n % 4 == 1 or n % 4 == 2)
	tr[n, j] = fb
	for i in range(n - 1, -1, -1): #backward
		j = tr[i + 1, j]
		bits_r[i] = j
	#undo bits_r, bits_c mapping
	bits_c = np.logical_xor(bits_r[1:], bits_r[:-1])
	return {'path_metric': pm[n - 1, fb], 'decoded_bits': np.logical_xor(bits_c, np.arange(n - 1) % 2)[:-1] }; #, 'bits_r' : bits_r, 'pm' : pm, 'tr' : tr }


def fm0_lowpass_filt(w1, spb):
	lp_dl = spb * 96
	lowpass = remez(2*lp_dl+1, [0, 0.75 * 0.85, 0.75 * 0.90, 0.5 * spb], [1, 0], Hz=spb).reshape([-1,1])
	return convolve(w1, lowpass, 'same')

def fm0_subcar_filt(wf, spb):
	if (len(wf.shape) == 1):
		wf = wf.reshape([-1, 1])
	return fm0_lowpass_filt(wf * np.exp(-1.5j * np.pi / spb *  np.arange(wf.shape[0]).reshape([-1, 1])), spb)



def fm0_bb_decimate(w2, spb, zpad=False, fractional = False): # matched filter, decimate, whiten, do_fm0_mix
	# undo the gain matched filter (spb==max(conv(half_sine,half_sine))), and ssb mixing halves amplitude
	w3 = convolve(w2, half_sine(spb)/  (0.5 * spb))
	if fractional:
		samples= rearrange(w3[0 : w3.shape[0] // spb * spb, :], '(l spb) ch -> l (spb ch)', spb=spb)
	else:
		samples = w3[ :: spb, :]
	whitened = lfilter([1], [1, _c1], samples, axis=0) * np.pi * _c1
	whitened = whitened[1:-1, :] if zpad else  whitened[2:-2, :]
	L = whitened.shape[0]
	if fractional:
		whitened_f = rearrange(whitened, 'l (spb ch) -> (l spb)  ch', spb = spb)
		return  np.exp(np.pi/(2j * spb) * np.arange(L * spb)).reshape([L * spb, -1]) * whitened_f
	else:
		return np.exp(np.pi/2j * np.arange(L)).reshape([L, -1]) * whitened

def fm0_encode_square(bits, spb, invert = False, tail = True):
	assert(spb % 2 == 0)
	if (tail):
		bits = np.hstack((bits, np.logical_xor.reduce(bits)))
	n = bits.shape[0]
	last_enc = invert
	out = np.zeros((n * 2, 1))
	box = np.ones((1, spb // 2))
	for i in range(n):
		out[i * 2, 0] = not last_enc
		bit_out = not last_enc if bits[i] else last_enc
		last_enc = out[i * 2 + 1, 0] = bit_out
	out = ((out*2-1) * box).reshape((-1,))
	return out

def fm0_test():
	raise NotImplementedError() # fm0 decimate still uses old interface and dimensions
	spb = 8
	for length in range(5,9):
		for bits in range(2 ** length):
			source = np.array([bits >> i & 1 for i in range(length)], dtype=bool)
			fm0 = fm0_encode(source, spb)
			decode= fm0_decode(fm0, spb, 0)
			try:
				assert(decode['path_metric'] < 0.1)
				assert (not any(decode['decoded_bits'] ^ source))
			except AssertionError as e:
				print(f'source={source}')
				raise e
	print(f'test passed')
from linregs import make_pn

def run_fm0_simulation(spb): # should return Q function: twice BER as QPSK 
	raise NotImplementedError() # fm0 decimate still uses old interface and dimensions
