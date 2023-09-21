#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from IPython.core.display import HTML, display
display(HTML('<script>Jupyter.notebook.clear_all_output(); '+
             'Jupyter.notebook.kernel.restart();</script>'))


# In[1]:


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
from itertools import combinations
import argparse
    

import plotly.io as pio
import profile; 


# In[20]:

parser = argparse.ArgumentParser()
# parser.add_argument('filename',type=str)
parser.add_argument('dist',type=int)
parser.add_argument('--volt',type=int)
args = parser.parse_args()

dist = args.dist
volt = args.volt

# save_filename = args.filename
# distance = args.dist
# print(distance)

### File options ###
home="/home/jradema/Documents/sk/oceans/"

# dist=150
# tx_amp = 0.08 # 60Vrms
tx_amp = 0.06 # 40Vrms
# tx_amp = 0.015 # 10Vrms
# tx_amp = 0.025 # 20Vrms
# volt = 40 #np.round(tx_amp*1000).astype("uint8")
# trial=1
num_trials = 24
trial_start = 0
angles = range(-90,91,5)
nch=7

msk_ord = 0
data_rate = 500
fs = 2e5
nbits = 10000

rx_dir=os.path.join(home,"rx_outputs");
tx_dir=os.path.join(home,"tx_outputs");
expname='River_PAB2_Van_Atta';
now = datetime.now()
datestr=now.strftime("%m-%d-%Y");
datestr="Range_0124_0127_0130"
expset=expname + "_" + datestr;

try:
    os.mkdir(os.path.join(rx_dir,expset));
    print(f"Created directory {expset} in {rx_dir} for experiment set.")
except FileExistsError as error:
    print(f"Directory {expset} already exists")
    
# rx_root = f"fixed_vanatta4x2_dr={data_rate}bps_ord={msk_ord}_Vrms={volt}_{dist}m_1m_single_foam_sep_purui_rx";

#rx_filename = f"vanatta4x2_{nch}ch_rx_nostag_006B_006F_006A_006C_x_001A_004A_004B_004D_18,5kfc_0,0deg_0,5kbps_usrp_2,5m_depth_005B_{volt}Vrms_{dist}m_{dist}m_1m_2foam_sep_{trial}.dat";

tx_filename = f"msk_data_dr={data_rate}_ord={msk_ord}_nbits=10k_fs=2e5.dat";
tx_path = os.path.join(tx_dir,tx_filename);


# In[21]:


def get_wf( filename, spb, Fc_c, data_rate, msk_ord, Fs, hyd_num, stretch=1):
    Fc_t = -data_rate * (msk_ord*2+3) /4;
    Fs_new = data_rate * spb;
    segment=np.fromfile(filename,dtype=np.float32).reshape([-1, 7])[:, hyd_num]
#     pcs(segment[:,0])
#     print(f'segment has shape {segment.shape}, Fs={Fs}, Fs_new={Fs_new}')
    return mixer_cfo(Fc_c, [Fc_t, -Fc_t], stretch, Fs, segment, Fs_new)

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


# In[25]:


### EQUALIZE vs. VARIABLE (angle, range, etc.) ### 

# Nfft = 2**(np.ceil(np.log2(Fs*10)));
# att = 300; # attenuation of chebwin in dB

# var_list = [22,23,24,25,27,28,29,30,32,34,36,38,40,42,44,47]
# var_list = [40,47,50,53,56,60,63,66,70,78,80,81,83]
# var_list = [81,83,90,100]
# var_list = [10,20,30,40,50,60]
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
exp_bits = source[train_len:]

# var_list = [70,80,90,100,110,120,130,140,150]
#snr_arr = np.zeros((len(var_list),))
#ber_arr = np.zeros((len(var_list),))
# volt = 40
# snr_dis = np.zeros((1,7))
# ber_dis = np.zeros((1,7))
# spectrums = np.zeros((len(var_list),Nfft.astype("uint32")))

# for i in range(0,len(var_list)):   
#     dist = var_list[i]
#     print(f"Distance: {dist}")
    
rx_root = f"fixed_vanatta4x2_dr={data_rate}bps_ord={msk_ord}_Vrms={volt}_{dist}m_1m_single_foam_sep_purui_rx";

snr_arr = np.zeros((127,num_trials))
ber_arr = np.zeros((127,num_trials))

hyd_index = 0
import gc
for hyd in range(1,nch+1):

    # snr_h = []
    # ber_h = []
    hyd_comb = combinations([0,1,2,3,4,5,6],hyd)
    print(f"Num_Hyd: {hyd}")

    for h in hyd_comb:
        
        print(f"hyd_index {hyd_index}")
        for trial in range(0,num_trials):
            if (dist == 40 or dist == 50) and trial >= 0 and trial < 4 and data_rate == 500:
                continue
                
            if dist == 150 and trial >=4 and trial < 8 and data_rate == 1000:
                continue

            rx_filename = f"{rx_root}_{trial}.dat"
            rx_path = os.path.join(rx_dir,expset,rx_filename); 

            if not os.path.isfile(rx_path):
                #print(f'Could not find trial {trial}, skipping...')
                continue
            
            print(f"Trial: {trial}")

            #print(h)
            dat1 = get_wf(rx_path, spb, Fc_c, data_rate, msk_ord, Fs, h,)

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
            #print(f'pos={pos}')
            bpsk_sample = fm0_bb_decimate(rxfilt[pos:, :], spb, True, True);

            try:
                d_all, e_all, eq = dfe_sfo_main(bpsk_sample, max_step, spb, lam, fw_len, bw_len, preamble_bits, LDLT_RLS);

                ber = np.sum(exp_bits != d_all) / max_step;
                bplot = np.cumsum(exp_bits != d_all) ;
                avg_len = 20
                mse = 10*np.log10(convolve(np.abs(e_all)**2,np.ones((avg_len,)))/avg_len)
                snr = d_all.size/np.sum(np.abs(e_all+d_all-exp_bits)**2);

                del eq
            except:
                snr = 1
                ber = 0.5

            #print('snr is %.2f dB' % (10*np.log10(snr)))
            #print('ber is %.2E' % ber)
            if ber < 1e-4:
                ber = 1e-4
                
            snr_arr[hyd_index, trial] = snr
            ber_arr[hyd_index, trial] = ber

        hyd_index = hyd_index + 1

    # snr_hyd[hyd] = np.median(snr_h)   
    # ber_hyd[hyd] = np.median(ber_h)
    
# snr_dis[i,:] = snr_hyd
# ber_dis[i,:] = ber_hyd

snr_arr.astype("float32").tofile(f"matplot/range_results/snr_multihyd_single_foam_{dist}m_{data_rate}bps_{volt}Vrms_4x2va_{datestr}.bin")
ber_arr.astype("float32").tofile(f"matplot/range_results/ber_multihyd_single_foam_{dist}m_{data_rate}bps_{volt}Vrms_4x2va_{datestr}.bin")
    
    #snr_arr[i] = 10*np.log10(np.median(snrs));
    #ber_arr[i] = np.median(bers) if np.median(bers) >= 1e-4 else 1e-4;





# %%
