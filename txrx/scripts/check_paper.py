import numpy as np;
import plotly.graph_objects as go
from numpy.linalg import *
from scipy.signal import periodogram, convolve, remez, freqz, find_peaks, lfilter, welch;
from scipy.signal.windows import chebwin;
import pandas as pds;
from datetime import datetime
from linregs import *;
from plots import *;
from msk import *;
# from dbutils import *;
from equalization import *; 
from jack_record import record;
import os
from datetime import datetime
import subprocess

import plotly.io as pio
import profile;

import argparse

parser = argparse.ArgumentParser()
# parser.add_argument('filename',type=str)
parser.add_argument('dist',type=int)
parser.add_argument('--volt',type=int)
parser.add_argument('--bps',type=int)
args = parser.parse_args()

dist = args.dist
volt = args.volt
bps = args.bps

print(f"Distance={dist}m, Volt={volt}rms, DR={bps}")

# File options ###
home="/home/jradema/Documents/sk/oceans/"

# tx_amp = 0.08 # 60Vrms
tx_amp = 0.05 # 40Vrms
# tx_amp = 0.015 # 10Vrms
# tx_amp = 0.028 # 20Vrms
# tx_amp = 0
# trial=1
num_trials = 24
trial_start = 0
angles = range(0,1)
nch=7

msk_ord = 0
data_rate = bps
fs = 2e5
nbits = 10000

rx_dir=os.path.join(home,"rx_outputs");
tx_dir=os.path.join(home,"tx_outputs");
expname='River_PAB2_Van_Atta';
now = datetime.now()
datestr=now.strftime("%m-%d-%Y");
datestr="Range_0124_0127_0130"
expset=expname + "_" + datestr;
# expset="./"

try:
    os.mkdir(os.path.join(rx_dir,expset));
    # print(f"Created directory {expset} in {rx_dir} for experiment set.")
except FileExistsError as error:
    pass
    # print(f"Directory {expset} already exists")
    
# rx_root="noise_test2"
# rx_root = "test_data_resamp_int_tp.bin"

tx_filename = f"msk_data_dr={data_rate}_ord={msk_ord}_nbits=10k_fs=2e5.dat";
tx_path = os.path.join(tx_dir,tx_filename);

def get_wf( filename, spb, Fc_c, data_rate, msk_ord, Fs, stretch=1):
    Fc_t = -data_rate * (msk_ord*2+3) /4;
    Fs_new = data_rate * spb;
    segment=np.fromfile(filename,dtype=np.float32).reshape([-1, 7])[:, :]
#     pcs(segment[:,0])
#     print(f'segment has shape {segment.shape}, Fs={Fs}, Fs_new={Fs_new}')
    return mixer_cfo(Fc_c, [Fc_t, -Fc_t], stretch, Fs, segment, Fs_new), segment

def fm0_bpsk1(source):
    last_sym = 1;
    polarity = -1;
    output = np.zeros((source.size * 2));
    for i in range(source.size):
        if last_sym != source[i]: # xor = 1;
            output[i*2:i*2+2] = -polarity, -polarity
            polarity = -polarity;
        else: # xor = 0
            output[i*2:i*2+2] = -polarity, polarity
        last_sym = source[i];
    return output;

### EQUALIZE vs. VARIABLE (angle, range, etc.) ### 
spb=2
source = pn15_seq[:nbits] * 2 - 1;    

fw_len=15
bw_len=4
train_len = 1000;
lam = 0.999
Fc_c = 18500
Fs=192000

max_step = source.size-train_len

preamble_bits = source[:train_len];
template_fm0 = fm0_bpsk1(preamble_bits);
template = fm0_subcar_filt(template_fm0, spb)
var_list = angles
exp_bits = source[train_len:]

snr_arr = np.zeros((len(var_list),))
ber_arr = np.zeros((len(var_list),))
# volt = 10

# spectrums = np.zeros((len(var_list),Nfft.astype("uint32")))

for i in range(0,len(var_list)):
    snrs = []
    bers = []
    
    angle = var_list[i]
    # print(f"Angle: {angle}")

    rx_root = f"fixed_006A_dr={data_rate}bps_ord={msk_ord}_Vrms={volt}_{dist}m_1m_single_foam_sep_purui_rx";
    
    for trial in range(trial_start,num_trials+trial_start):
        rx_filename = f"{rx_root}_{trial}.dat"
        rx_path = os.path.join(rx_dir,expset,rx_filename); 

        if not os.path.isfile(rx_path):
            print(f'Could not find {rx_filename}, skipping...')
            continue

        print(f"Trial: {trial}")

        # if trial >= 0 and trial < 4:
        #     continue

        dat1, segment = get_wf(rx_path, spb, Fc_c, data_rate, msk_ord, Fs,)

        # sig = passband[:, 6];
        # window_size = np.floor(len(sig)/100);
        # window = chebwin(window_size.astype("uint32"), att);

        # spectrums[i,:] = welch(sig,Fs,window=window,nfft=Nfft);

        # print(f'dat has shape{dat1.shape}')
        rxfilt = fm0_lowpass_filt(dat1, spb);
        # print(f"rxfilt has shape {rxfilt.shape}")
        scale,corr = linreg_mc_od(rxfilt, template);
        
        offset = fw_len//3
        #offset  =0;

        pos= np.argmax(corr[:]) + (offset - preamble_bits.size) * spb + 1;
        # print(f'pos={pos}')
        bpsk_sample = fm0_bb_decimate(rxfilt[pos:, :], spb, True, True);

        try:
            d_all, e_all, eq = dfe_sfo_main(bpsk_sample, max_step, spb, lam, fw_len, bw_len, preamble_bits, LDLT_RLS);

            ber = np.sum(source[train_len:train_len+d_all.size] != d_all) / max_step;
            bplot = np.cumsum(source[train_len:train_len+d_all.size] != d_all) ;
            avg_len = 20
            mse = 10*np.log10(convolve(np.abs(e_all)**2,np.ones((avg_len,)))/avg_len)
            snr = d_all.size/np.sum(np.abs(e_all+d_all-exp_bits)**2);

        except:
            snr = 1
            ber = 0.5

        print('snr is %.2f dB' % (10*np.log10(snr)))
        print('ber is %.2E' % ber)
        
        bers.append(ber)
        snrs.append(snr)

    snr_arr[i] = 10*np.log10(np.median(snrs));
    ber_arr[i] = np.median(bers) if np.median(bers) >= 1e-4 else 1e-4;

print("Median SNR is %.2f dB" % snr_arr)
print("Median BER is %.2E" % ber_arr)

# pcs(snr_arr)
# pcs(ber_arr)
    