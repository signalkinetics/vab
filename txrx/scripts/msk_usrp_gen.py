from msk_modulator_usrp import make_tx_data
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('data-rate',default=500,type=int)
parser.add_argument('order',default=0,type=int)
parser.add_argument('nbits',default=10000,type=int)

args = parser.parse_args()

dr = vars(args)['data-rate']
order = vars(args)['order']
nbits = vars(args)['nbits']
fs = 2e5

# order = 0
# dr = 250
# nbits = 10000
# fs = 2e5

sig = make_tx_data(order,fs,dr,nbits)
filename = f'../tx_outputs/msk_data_dr={dr}_ord={order}_nbits={nbits // 1000}k_fs=2e5.dat';
sig.tofile(filename)
print(f"Shape is {sig.shape}")
print("done")