function[data_out, prbs_seq] = prbs_gen(POLY_LENGTH, POLY_TAP, NBITS, NUM)

prbs = ones(1,POLY_LENGTH);

for j = 1:NUM
	for i = 1:NBITS
		prbs_xor = xor(prbs(POLY_LENGTH),prbs(POLY_TAP));
		prbs = [prbs_xor, prbs(1:POLY_LENGTH-1)];
		data_serial(i) = prbs_xor;
	end
	data_out(:,j) = data_serial;
end
prbs_seq = reshape(data_out,1,NBITS*NUM);

