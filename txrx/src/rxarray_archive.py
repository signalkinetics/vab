Fc = 20000;
Fs = 192000;
spb = 8;
source=pn10_seq

@save_src_decor
def array_rx(ii, data_rate, sub_ch, sideband, LRs):
	
	Fs_new = spb * data_rate;
	segment=get_wf(ii).reshape([-1,8])[:,sub_ch];
	dat1 = mixer(Fc, Fs, segment, Fs_new);
	if (sideband == 'lsb' or sideband == 'usb'):
		if (sideband == 'lsb'):
			dat1 = np.conj(dat1);	
		rxfilt = fm0_subcar_filt(dat1, spb);
	else:
		rxfilt = np.hstack([fm0_subcar_filt(dat1, spb),fm0_subcar_filt(np.conj(dat1), spb)])
	pn5_wf=make_seq(pn5_str*3, spb);
	template = fm0_subcar_filt(pn5_wf, spb)
	scale,snr = linreg_mc_od(rxfilt, template);

	def peak_snr(snr, n):
		p,v= find_peaks(snr, distance=7e3)
		s=[(snr[k], k) for k  in p]; s.sort(reverse=True)
		return sorted([i[1] for i in s[:n]])
	pos4 = peak_snr(snr, 2 if data_rate == 500 else 4)

	def do_ber(bits):
		if (source.size != bits.size):
			return 1;
		return np.logical_xor(bits, source).sum() / source.size
	@save_result_decor
	def array_packet_handle(packetnum, ii, data_rate, LR0, sub_ch, spb):
		pos = pos4[packetnum];
		packet = rxfilt[pos + 1:pos+1+(source.size + 1)*spb, :]
		lm1 = fm0_ideal_pm_lms(packet, spb, source, LR0, scale[pos, :], True)
		dd = fm0_decode(packet, spb, LR0, scale[pos, :], True);
		ber = do_ber(dd['decoded_bits'])
		return {'datapos' : pos, 'BER' : do_ber(dd['decoded_bits']), 'SNR' : 10*np.log10(lm1['snr']), 'sideband' : sideband }
	for LR0 in LRs:
		for i in range(len(pos4)):
			try:
				array_packet_handle(
					packetnum = i,
					ii = ii,
					data_rate = data_rate, 
					LR0 = LR0,
					sub_ch = sub_ch,
					spb = spb
				);
			except:
				pass

for data_rate in [1000]:
	for combination  in [[0,1,2,3,4,6,7]]: # note 5 is broken;
		for sideband in [ 'lsb', 'usb']:
			for i in db_exp.find({'rxarray' : True}):
				array_rx(i['_id'], data_rate, combination, sideband, [0.05]);
from pprint import pprint
for i in db_exp.find({'rxarray' : True, 'bsdist' : 10}):
	a = list(db_exp.find({'ii' : i['_id'], 'data_rate':500, }, {'BER', 'sideband', 'packetnum', 'SNR'}).sort([('packetnum', 1), ('BER', 1)]))
	for b in a: print(b)
