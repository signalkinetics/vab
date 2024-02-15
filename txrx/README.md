
# Concepts
## MSK
Continuous phase, progress by $+\pi/2$, or $-\pi/2$ per bit
## FM0
For each bit,  flip the output polarity once (bit 1, $\pm\pi$, ) or twice (bit 0, $\pm2\pi$)

## MSK+Subcarrier
For each bit, flip N times to represent 1, or flip N+1 times to represent 0.

## OQPSK representation of MSK
OQPSK modulation: shift every other BPSK symbol by $\pi/2$ phase ($\times1j$). The  phase change between symbols is then $\pi/2$ or $-\pi/2$.

 Convolve an OQPSK impulse train with a two-symbol half-sine pulse will result in a MSK waveform, with the new  i-th MSK bits as the XOR of $(i-1)$-th and $i$-th OQPSK bit, because OQPSK is in absolute phase but MSK is in phase difference (i.e. frequency).

That is to say, if the transmitter does XOR precoding to modify the source before sending to MSK(FM0 or subcarrier) encoder, the information will be in the OQPSK format and information is encoded in **absolute phase**. The benefit of use such precoding is that it drops the phase continuity constraint in FM0 that doubles the BER at the same SNR compared to OQPSK.

## OQPSK decode pipeline (SPB : #sample per bit):
1. Matched filter of the half-sine pulse
2. Noise-whitening of the previous match filter step: IIR H(z)=1/(1+0.359z^(-SPB))
3. Modulate with exp(-1j / SPB * pi / 2), this converts the OQPSK back to BPSK
4. Feed to a vanilla BPSK DFE equalizer.

# How to Compile
The equalization operates on the square root of the covariance matrics, which is supported via C++/Eigen's Cholesky module. Therefore, the equalization needs to compile the underlying `linreg_mc.cpp` in order to function properly.

Currently, the supported platform are Ubuntu 20.04 and 22.04, X86-64. (ARM machines, especially the Apple silicon Macs are not supported) On those machines you may
```bash
(cd src; make)
```
and the `linreg_mc.so` will be compiled, so the `equalization.py` can call it properly.

# Code structure
`` scfde.py`` single carrier frequency domain equalization

`` equalization.py `` : Time domain equalization, accepts BPSK data or the preprocessed single-sideband MSK data from `msk.py`


`` msk.py `` MSK half-sine match filter and whiten filter. The subcarrier is not maintained and should instead be handled by ``mixer_cfo`` in linregs.py. After the pre-processing it should become BPSK samples and handled by `equalization.py`. The Viterbi decoder was intended for single-tap channel and isn't useful/maintained. 


`` linregs.py `` mixer function, linear regression for single and multiple channel single-tap preamble alignment
