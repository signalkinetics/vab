# SK Piezoelectric Backscatter Van Atta

The Van Atta is a specialized antenna array first created for the RF domain in 1955. It is a retroreflective structure, in that it reflects an incoming wave back in the same direction that it came from. You can view the original patent [here](https://patents.google.com/patent/US2908002A/en).

This project is an extension of the Van Atta array to piezoelectric devices in the underwater acoustic domain. The method is effectively the same, but there are unique challenges specific to piezoelectric devices, specifically maximum power transfer matching. Our group's derivations for far field beam patterns and maximum power transfer matching are given in [this document](https://drive.google.com/file/d/1DULdc8UHSOuENSTepJ27DRMAoALPb1Ii/view?usp=sharing).

## Installation

Visit our [How to Set-Up USRP Tutorial](https://docs.google.com/document/d/1ePhpxb1y15XBV4cvTJPKnmYtWSkvokSHL-_EcKC4X9U/edit) for UHD build, installation, and getting started instructions. Installing this repo is a prerequisite for compiling the C++ sources that control the USRP. 

Once UHD has been `make`'d and `sudo make install`'d in a separate directory, if installing for the first time, run: <br>
`mkdir build && cd build && cmake ../ && make -j` <br>

The binaries will be built into executables in the build folder.  

## Use

#### `/src/matlab/get_sparams.m`

Converts a node RLC model into the equivalent 2-port network S-parameters. The derivation of this code is shown in the theory document above. 

#### `/src/matlab/impedance_read.m`

A small snippet to read impedance CSVs and plot them. 

#### `/src/matlab/rlc_modeler.m`

Uses a least-squares algorithm to model the impedance of a node as an equivalent RLC, and provides the R, L, and C values that best fit the impedance. Needs to be run a couple times before convergence.

#### `/src/matlab/spectrogram_plot.m`

An unfinished piece of spectrogram code. Mostly unused.

#### `/src/matlab/subcarrier_extractor.m`

Accepts datasets with backscatter peaks next to the carrier. Extracts their location based on a provided modulation frequency and averages their power across multiple trials and angles of rotation.

#### `/src/matlab/subcarrier_extractor_baseband.m`

Similar to the subcarrier extractor, except this code provides filtering and downconversion, and then computes an SNR-like metric for received backscatter power. 

#### `/src/matlab/van_atta_beam_patterns.m`

Computes a theoretical van atta beam pattern based on the theory given in the document above. Incorporates computation of S-parameters and their effect on the pattern. Accepts node RLC models and array parameters such as element spacing, frequency, element number, and angle of incidence. Plots the received pattern for two different modes of array pattern and direction pattern. 

#### `/src/c/rx_multi_samples.cpp`

A modified UHD source that connects to multiple USRPs and places each channel into the real and imaginary part of the .dat file

## Datasets

Datasets for measured node [impedances](https://drive.google.com/drive/u/1/folders/1Lwp_WQQzYhoznev65niFUcd4ucXDke_k) and [collected backscatter measurements](https://drive.google.com/drive/u/1/folders/1og-sW82Xx33VXqlI1ZAbbL-HCPHXB5KW) are located in the google drive. 