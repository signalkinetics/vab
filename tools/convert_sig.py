import numpy as np
from soxr import resample
import argparse
import os

parser = argparse.ArgumentParser()
# parser.add_argument('filename',type=str)
parser.add_argument('filename',type=str)
args = parser.parse_args()

fname = args.filename
print(f"Converting {os.path.basename(fname)}...")

sig = np.fromfile(fname, dtype="float32").reshape((-1,7))
conv = resample(sig, 192000, 44100)
conv_int = (conv*(2**23-1)).astype("int32")

np.transpose(conv_int).tofile(f"{os.path.splitext(fname)[0]}_44k1_int32.bin")
