# from plotly.subplots import make_subplots
import numpy as np;
# import plotly.express as px
# import plotly.graph_objects as go

from scipy.signal.windows import chebwin;
from pyfftw.interfaces.numpy_fft import fft, ifft;
from scipy.signal import periodogram, convolve, remez, freqz;

from linregs import *;
from msk import *;
# from dbutils import *;

# pcs = lambda x: pc(x.ravel()).show(renderer='notebook');
# from inspect import getsource

# !clang++-12 test_fpga.cpp -O2 -lpthread -o test_fpga
# commit array mapping:
# mapping: ch0 : 2; ch1: 0; ch2:  3; ch3: 1; ck4:7 ch5:4 ch6:6 ch7: 5

# import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument('rx-file',
#                     help="name for the rx file")
# parser.add_argument('--n',default=7,type=int,help="number of channels")
# parser.add_argument('--nsamps', default=7000000,type=int,help='number of samples to receive')

# args = parser.parse_args()

# @save_src_decor
def test_lite(myarray, enable_gdb, rxch = tuple(range(8))):
    chkval = np.abs(myarray.reshape((-1,))).max();
    if (0.9999 < chkval):
        raise OverflowError;
    
    writearray =(myarray[:, [1, 3, 0, 2, 5, 7, 6, 4]] * (2**31 - 1)).astype('int32');
    writearray.tofile('/dev/shm/dac.bin');
    import subprocess;
    proc = subprocess.Popen(['/home/prwang/pyaudio_test/test_fpga', '/dev/shm/dac.bin', '/dev/shm/adc.bin']);
    # if (enable_gdb):
    #     !(sleep 1; ncat localhost 23333 < /dev/null) # message to gdb
    proc.wait();
    raw = np.fromfile('/dev/shm/adc.bin',dtype='int32').reshape((-1,2,8)).astype('float64') / (2**23 - 1); 
    #signal back is 24 bit sign extended;
    merged = (raw[:,0,rxch] - raw[:,1,rxch]).astype('float32'); #two adcs sampling opposite polarity, 3db increase snr
    print(f'rx data has shape {merged.shape}')
    return merged;

Fs=192000; # on 49152 OCXO
Fc=18500;

# @save_src_decor
def siso_test(seconds, amp_all,enable_gdb=False):
    Phase = np.arange(8).reshape((1,-1))*0;
    Amplitude = np.array([[0,0,0,1, 0,0,0,0]]) * amp_all; 
    myarray = np.sin(np.arange(Fs * seconds).reshape([-1, 1])*Fc*2*np.pi/Fs + Phase)  *Amplitude
    return test_lite(myarray, enable_gdb,rxch=(0,))

def record(filename,nch,nsamps,amp_all=0.07):
    Phase = np.arange(8).reshape((1,-1))*0;
    Amplitude = np.array([[1,0,0,0, 0,0,0,0]]) * amp_all; 
    myarray = np.sin(np.arange(nsamps).reshape([-1, 1])*Fc*2*np.pi/Fs + Phase)*Amplitude

    rx = test_lite(myarray,False,tuple(range(0,nch)))
    rx.tofile(filename)
