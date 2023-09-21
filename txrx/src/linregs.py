import numpy as np
import os
import ctypes as ct
import plotly.graph_objects as go
from numpy.testing import assert_allclose
from scipy.linalg import cho_factor, cho_solve
from scipy.signal import  convolve, remez, find_peaks
import re
from scipy.signal.windows import chebwin
from pyfftw.interfaces.numpy_fft import fft, ifft
from soxr import resample
from numba import njit,stencil
from plots import *

abssq = lambda y: np.real( np.conj(y)  * y)  # absolute value squared
inner = lambda x, y : np.sum(np.conj(x) * y)
nom = lambda x : np.real(inner(x, x))
dbvol = lambda x: np.log10(abs(x)) * 20
dbpwr = lambda x: np.log10(abs(x)) * 10


def linreg_offset(rx,tpl):
	assert rx.shape[0] == tpl.shape[0] 
	assert 1 == len(rx.shape) 
	assert 1 == len(tpl.shape)
	m = tpl.shape[0]; res = {}
	res['sx'] = sx = rx.sum(); sy = tpl.sum()
	res['sxxb'] = sxxb = nom(rx) - nom(sx) / m
	res['syyb'] = syyb = nom(tpl) - nom(sy) / m
	res['sxyb'] = sxyb = inner(rx, tpl) - inner(sx, sy) / m

	res['scale'] = sxyb / sxxb
	res['bias'] = (sy - sx * res['scale']) / m
	res['snr'] = -1 + syyb / (syyb - sxxb * nom(res['scale']))
	return res


_herm = lambda m : np.conj(m).transpose()


def linreg_mc(rx, tpl): # uses normal equation
	if len(tpl.shape) != 2:
		tpl = tpl.reshape((-1, 1))
	assert len(rx.shape) == 2 
	assert rx.shape[0] == tpl.shape[0] 
	m = tpl.shape[0]; res = {}
	#print(f'linreg_mc: rx has shape{rx.shape}')
	sxy = _herm(rx) @ tpl
	res['scale'] = cho_solve(cho_factor(_herm(rx)@rx), sxy)
	syy = np.real(_herm(tpl) @ tpl).item()
	loss = syy - (_herm(sxy) @ res['scale']).item()
	res['snr'] = -1 + syy / loss
	return res


def measure_power_complex(signal):
	return np.real(np.sum(np.conj(signal) * signal)) / (signal.size)


def measure_power_real(signal):
	return np.sum(signal * signal) / ( signal.size)


def awgn(signal, snr): # this equates the ENERGY ratio of signal and noise
	# if has long zeros then the per-power SNR should be compensated by duty ratio.
	def noise(variance):
		return np.random.normal(0, np.sqrt(variance), signal.shape).astype(signal.dtype)
	if (np.iscomplexobj(signal)):
		variance = measure_power_complex(signal) *  10**(-0.1 * snr) / 2
		return signal +  noise(variance) + noise(variance)*1j
	else:
		variance =  measure_power_real(signal) * 10**(-0.1 * snr)
		return signal + noise(variance)


def _od_pad(rx, m):
	N = 1 << (2 + m.bit_length())
	bn = N - m + 1
	rxpwr = measure_power_complex(rx)
	n = rx.shape[0]
	rx = np.concatenate([np.random.normal(0, rxpwr**0.5, [bn + (-n) % bn, *rx.shape[1:]]), rx], axis=0).astype('complex64')
	columns = rx.shape[0] // bn

	def unpad(x): #unpad:  remove the padding induced zeros in the result so it has the same shape as input.
		xsh = x.shape
		return x.reshape([xsh[0] * xsh[1], -1]) [ (-n) % bn:, :].reshape([-1, *xsh[2:]])
	return N, bn, rx, columns, unpad

dn = os.path.dirname(os.path.realpath(__file__))
accel = ct.CDLL(os.path.join(dn, 'linreg_mc.so'))
accel.ldlt_linreg.argtypes = (ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.c_int, ct.c_int, ct.c_int, ct.POINTER(ct.c_float))
accel.ldlt_linreg.restype = None
@stencil
def re_a_h_b(a, b):
	return a[0].real * b[0].real + a[0].imag * b[0].imag


@stencil
def re_a_h_b_2d(a, b):
	return a[0,0].real * b[0,0].real + a[0,0].imag * b[0,0].imag


def linreg_mc_od(rx, tpl): #sliding window, adaptive combination,
	m = tpl.shape[0]
	assert len(rx.shape) == 2 
	if len(tpl.shape) != 2:
		tpl = tpl.reshape((-1, 1))
	tpl = tpl.astype('complex64')
	N, bn, rx, columns, unpad = _od_pad(rx, m)
	k =  rx.shape[1]
	h_tpl = fft(np.flip(tpl, axis=0), n = N, axis=0)
	syy = np.real(_herm(tpl) @ tpl).item()

	scale_all = np.zeros((columns - 1, bn, k), dtype='complex64')
	loss_all = np.zeros((columns - 1, bn), dtype='float32')
# 	print(f'columns={columns}')

	for i in range(1, columns):
		view = rx[(i + 1)* bn - N: (i + 1) * bn, :]
		h_view =  fft(np.conj(view), n = N, axis=0)
		sxy_seg = ifft(h_view * h_tpl, n=N, axis=0)[m - 1 :, :]
		if 0: # backup in case of numerical 
			sxx_cumsum =  np.cumsum(np.conj(view).reshape((-1, k, 1)) * view.reshape((-1, 1, k)), axis = 0)
			sxx_all = sxx_cumsum[m - 1:, :, :] - np.concatenate((np.zeros((1,k,k)),sxx_cumsum[:-m, :, :]), axis=0)
			if i == 4:
				# print(sxx_all[0, :, :])
				# print()
				# print(sxy_seg[0, :])
				pass

			for j in range(0, bn):
				sxy = sxy_seg[j, :].T
				scale_all[i - 1, j, :] = cho_solve(cho_factor(sxx_all[j, :, :]), sxy)
		else:
			_wp = lambda x : x.ctypes.data_as(ct.POINTER(ct.c_float))
			accel.ldlt_linreg(_wp(view), _wp(sxy_seg), bn, m, k,  _wp(scale_all[i - 1, :, :]))
		if i == 4:
			# print(scale_all[i - 1, 0, :])
			pass

		loss_all[i - 1,: ] = syy - re_a_h_b_2d(sxy_seg[:, :],  scale_all[i - 1, :, :]).sum(axis=1)
	snr = -1 + syy / loss_all
	return (unpad(scale_all), unpad(snr))


#TODO need a streaming version: no padding, no offset
#view is a callback that returns one [bn] of data, may call syscall to listen on fifos, perform some downconversion and feed 
def _segment_sum(w, m): #sum of a sliding window of length m; we have $bn$ such results.
	c = np.cumsum(w)
	return c[m - 1:] - np.hstack(([0],c[:-m]))


def linreg_od(rx, tpl):
	m = tpl.shape[0]
	assert 1 == len(rx.shape) and 1 == len(tpl.shape)
	N, bn, rx, columns, unpad = _od_pad(rx, m)
	h_tpl = fft(np.flip(tpl), N)
	syyb = nom(tpl)
	scale_all = np.zeros((columns - 1, bn), dtype='complex64')
	bias_all = np.zeros_like(scale_all)
	snr_all = np.zeros((columns - 1, bn), dtype='float32')
	for i in range(1, columns):
		view = rx[(i + 1)* bn - N: (i + 1) * bn]
		sxy = ifft(fft(np.conj(view), N) * h_tpl)[m - 1 :]
		sxxb = _segment_sum(abssq(view), m)
		sxyb = sxy
		scale_all[i - 1, :] = scale = sxyb / sxxb
		residual=syyb - sxxb * abssq(scale)
		snr_all[i-1,:] = snr = -1 + syyb / residual
	return unpad(scale_all), unpad(snr_all)


def linreg_offset_od(rx, tpl):
	m = tpl.shape[0]
	assert 1 == len(rx.shape) and 1 == len(tpl.shape)

	N, bn, rx, columns, unpad = _od_pad(rx, m)
	h_tpl = fft(np.flip(tpl), N)
	sy = tpl.sum()

	syyb = nom(tpl) - nom(sy) / m

	#print(f'bn = {bn}, rx.shape[0]={rx.shape[0]}, columns={columns}, mod = {rx.shape[0]%bn}')
	
	scale_all = np.zeros((columns - 1, bn), dtype='complex64')
	bias_all = np.zeros_like(scale_all)
	snr_all = np.zeros((columns - 1, bn), dtype='float32')
	for i in range(1, columns):
		view = rx[(i + 1)* bn - N: (i + 1) * bn]
		sxy = ifft(fft(np.conj(view), N) * h_tpl)[m - 1 :]
		sx = _segment_sum(view, m)
		sxxb = _segment_sum(abssq(view), m) - abssq(sx) / m
		sxyb = sxy - np.conj(sx) * sy / m
		scale_all[i - 1, :] = scale = sxyb / sxxb
		bias_all[i - 1, :] = bias = (sy - sx * scale) / m
		residual=syyb - sxxb * abssq(scale)
		snr_all[i-1,:] = snr = -1 + syyb / residual
	return unpad(scale_all), unpad(bias_all), unpad(snr_all)


from pylfsr import LFSR


def make_pn(fpoly, string=False):
	order = max(fpoly)
	#capatible with firmware pnseq.h
	seq = LFSR(fpoly=fpoly).runKCycle(order + 2**(order) - 1)[order:] ^ 1
	if (string):
		return ''.join(map(str,seq))
	else: return seq;

pn5_str = make_pn([5, 2], True)
pn10_seq = make_pn([10, 3], False)
pn15_seq = make_pn([15, 1], False)
pn15_str = make_pn([15, 1], True)
pn7_str = make_pn([7, 1], True)
pn8_str = make_pn([8, 6, 5,1], True)
pn9_str = make_pn([9, 4,], True)


def make_seq(string, spb):
	prv = '0'; rt = []
	assert spb % 2 == 0;
	string=re.sub('\s','',string)
	for i in string:
		rt.append(np.hstack([[prv == i], np.ones((spb//2 - 1,))]) * (1 if i=='1' else -1))
		prv = i
	return np.hstack(rt).astype('float32')


def test_linreg():
	test_seq1 = '00011101'
	tpl_wf = make_seq(test_seq1, 8)
	#print(f'tpl_wf has shape{tpl_wf.shape}')
	delay = np.zeros(170,dtype='complex64')
	occurrence = [20,80,160]

	delay[occurrence] = 1+1j
	sig = convolve(tpl_wf.astype('complex64'), delay)
	duty_dB = dbpwr(len(occurrence) * tpl_wf.shape[0] / sig.shape[0])
	sig_n = awgn(sig, 20 + duty_dB)
	scale, bias, snr = linreg_offset_od(sig_n, tpl_wf)
	sig_mc = np.tile(sig, (4,1)).T
	sig_mc_n = awgn(sig_mc, 20 + duty_dB) 
	scale_mc, snr_mc = linreg_mc_od(sig_mc_n, tpl_wf)

	for i in occurrence:
		ilt = i + tpl_wf.size - 1
		res = linreg_offset(sig_n[i : ilt + 1],tpl_wf)
		assert_allclose(res['scale'], scale[ilt], rtol=1e-5)
		assert_allclose(res['bias'], bias[ilt], rtol=1e-5)
		assert_allclose(res['snr'], snr[ilt], rtol=3e-5)
		res_mc = linreg_mc(sig_mc_n[i : ilt + 1, :], tpl_wf)
		assert_allclose(res_mc['scale'], scale_mc[ilt, :].reshape((-1,1)), rtol=1e-4)
		assert_allclose(res_mc['snr'], snr_mc[ilt], rtol=3e-4)
	pcs(dbpwr(snr))
	pcs(dbpwr(snr_mc))
	print('test passed!')


def arange_mod_fs_int(n, fs):
	a = np.arange(fs)
	return np.hstack((np.tile(a, (n // fs,)), a[:n % fs]))


def mixer(Fc, Fs, segment, Fs_new):
	# 1. Fc might be an array, do te range extension.

	Fc = np.atleast_1d(Fc)[np.newaxis, np.newaxis, :]

	#print(f'Fc is {Fc}')
	if len(segment.shape) == 1:
		segment = segment[:, np.newaxis, np.newaxis]
	elif len(segment.shape) == 2:
		segment = segment[:, :, np.newaxis]
	nsamp = segment.shape[0]
	car_ph = (np.mod(np.arange(nsamp).reshape((-1, 1, 1)) * Fc, Fs) / Fs * 2 * np.pi).astype('float32')
	mixed = np.hstack((segment * np.cos(car_ph), 
		segment * np.sin(car_ph))).reshape((nsamp, -1))  #  n m k -> n (2 m, k), pysoxr doesn't like complex input.

	n_ch = segment.shape[1] * Fc.size
	resampled = resample(mixed, Fs, Fs_new)
	return resampled[:, 0 : n_ch] + 1j * resampled[:, n_ch: n_ch * 2]
def mixer_cfo(Fc_reader, Fc_tag, tag_dev, Fs_nom, segment, Fs_new):

	assert 0.99 < tag_dev and tag_dev < 1.01;

	#2: Cartesian product of Fc_reader (multiple carrier frequency and Fc_tag (plus and minus sideband)
	Fc_reader_array= np.atleast_1d(Fc_reader).reshape((-1, 1)) / tag_dev
	Fc_tag_array = np.atleast_1d(Fc_tag).reshape((1, -1))
	#3:  Decoder process **negative** Fc_tag (lower sideband) normally and upper (positive) needs a complex conjugate
	return mixer(
			Fc=  (np.sign(-Fc_tag_array) * (Fc_reader_array + Fc_tag_array)).ravel(),
			Fs= Fs_nom / tag_dev,
			segment = segment, 
			Fs_new = Fs_new)

def peak_snr(snr, n):
	p,v= find_peaks(snr, distance=7e3)
	s=[(snr[k], k) for k  in p]; s.sort(reverse=True)
	return sorted([i[1] for i in s[:n]])
