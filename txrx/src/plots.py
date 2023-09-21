
import numpy as np;
import plotly.graph_objects as go
from scipy.signal import  convolve, remez, freqz;
from scipy.signal.windows import chebwin
from pyfftw.interfaces.numpy_fft import fft, ifft, fftshift, fftfreq;

def visfilter(Num, Fs, worN = 4096):
	w, h = freqz(Num, worN = worN, fs=Fs);
	fig=go.Figure();
	db = lambda x : 20 * np.log(np.abs(x))
	fig.add_trace(go.Scatter(x=w,y=db(h)))
	fig.show(render='notebook');
	return fig;

def pc(y, x=None, fig=None, name=None, Fs=None, xlabel=None, ylabel=None):
	if (y.size == max(y.shape)):
		y = y.ravel();
	if (fig is None):
		fig=go.Figure(); 
	name = [f'{name}-R', f'{name}-I'] if name else ['Re', 'Im'];
	if Fs is not None and x is None:
		x = np.arange(y.size) / Fs
	fig.add_trace(go.Scatter(x=x, y=np.real(y), name=name[0]));
	if (np.iscomplexobj(y)):
		fig.add_trace(go.Scatter(x=x, y=np.imag(y), name=name[1])); 
	return fig;

def pxx(sig, Fs, avg = 1, fig=None, name = None):
	
	sz = sig.shape[0] // avg;
	print(f'sz = {sz}')
	pxx_1 = np.zeros((avg, sz))
	def fftpwr(wave):
		window = chebwin(wave.shape[0], 150);
		window = window * wave.shape[0] / window.sum();
		spec = fftshift(fft(wave * window));
		return np.abs(spec * np.conj(spec)) / wave.shape[0]**2;
	for i in range(avg):
		pxx_1[i,:] = fftpwr(sig[i * sz : (i + 1) * sz])
	f = fftshift(fftfreq(sz, 1 / Fs));
	pxxdb = 10 * np.log10(pxx_1.sum((0,))/ avg); # power average
	#hide negative if signal is real.
	if not np.iscomplexobj(sig):
		f = f[sz //2 :];
		pxxdb = pxxdb[sz //2 :];
		annotation = [0, sz // 4];
	else:
		annotation =  [sz // 2, sz // 4, sz * 3 // 4]
	#pxxdb = 10 * np.log10(pxx_1).sum((0,))/ avg; # video average
	
	if fig is None:
		fig = go.Figure();
	#print(f'p0={pxxdb[sz//2]}, plm={pxxdb[sz//4]},prm={pxxdb[sz*3//4]}')

	fig.add_trace(go.Scatter(x=f,y=pxxdb, name=name))
	fig.update_layout(margin={'b':0,'l':0,'r':0,'t':0})
	for x in annotation:
		fig.add_annotation(x=f[x], y=pxxdb[x],
				text=f"{pxxdb[x]:.2f}",
				arrowhead=1)
	return fig;

def pxxdict(sigdict, Fs, avg, fig=None):
	for k,v in sigdict.items():
		fig = pxx(sig=v,Fs=Fs,avg=avg,fig=fig,name=k);
	return fig;

def pcs(*args, **kwargs):
	fig = pc(*args, **kwargs);
	fig.show(renderer='notebook')
	return fig;
