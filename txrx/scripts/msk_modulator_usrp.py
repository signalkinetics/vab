from linregs import *
def make_tx_data(msk_order, output_sr, msk_dr, src_len):
	dup = 10
	vol = 0.0
	voh = 0.8
	source = pn15_seq;
	raw_spb = ((msk_order + 1) * (msk_order + 2)) * dup
	osc = np.zeros((2, raw_spb), bool) # all start with state of 1
	for i in range(raw_spb):
		osc[1, i] = i % ((msk_order + 2) * dup * 2)  < ((msk_order + 2) * dup)
		osc[0, i] = i % ((msk_order + 1) * dup * 2)  < ((msk_order + 1) * dup)
	
	last_bit = 1;
	polarity = 1;
	generated_raw = np.zeros((src_len, raw_spb), np.double);
	for i in range(src_len):
		xorbit = source[i] ^ last_bit
		last_bit = source[i];
		generated_raw[i, :]  =  \
			vol  + (voh - vol) * np.logical_xor(polarity, osc[xorbit, :])
		polarity = polarity ^ xorbit ^  (msk_order & 1)
	a = resample(generated_raw.ravel(), raw_spb * msk_dr, output_sr)
	return a.astype(np.complex64)
