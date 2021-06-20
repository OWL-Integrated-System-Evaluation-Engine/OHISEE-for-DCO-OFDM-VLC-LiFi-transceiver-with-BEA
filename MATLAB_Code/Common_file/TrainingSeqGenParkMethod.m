function[LongTrainingSeqOutput,SeqAfterIFFTOutput,SeqAfterHermitianSymOutput] = TrainingSeqGenParkMethod(N,Ncp,RepeatTime,GainSubcarrierLinear,Channel_Gain)

if (mod(log2(N/2),1) ~= 0 | N==2)
	error('N is invalid');
end

% S_N sequence generation
for i = 1:log2(N/2)
	intial_seq = [1,1];
	if i == 1
		SeqOneOFDMSymTimeDom = intial_seq;
	else
		SeqOneOFDMSymTimeDom = [SeqOneOFDMSymTimeDom,SeqOneOFDMSymTimeDom(1:2^(i-1)/2),-1.*SeqOneOFDMSymTimeDom(2^(i-1)/2+1:end)];
	end
end

% Hermitian Symmery
SeqAfterHermitianSym = [0,SeqOneOFDMSymTimeDom(1:N/2-1),0,fliplr(conj(SeqOneOFDMSymTimeDom(1:N/2-1)))];
SeqAfterHermitianSymOutput = SeqAfterHermitianSym;
% IFFT calculation
SeqAfterIFFT = ifft(SeqAfterHermitianSym,N,2);
SeqAfterIFFTOutput = SeqAfterIFFT;
if Channel_Gain
	SeqAfterHermitianSym = [0,SeqOneOFDMSymTimeDom(1:N/2-1).*GainSubcarrierLinear(2:N/2),0,fliplr(conj(SeqOneOFDMSymTimeDom(1:N/2-1).*GainSubcarrierLinear(2:N/2)))];
	SeqAfterIFFT = ifft(SeqAfterHermitianSym,N,2);	% IFFT calculation
end
% Repeat
SeqAfterIFFTRepeat = repmat(SeqAfterIFFT,1,RepeatTime);
% Adding CP
LongTrainingSeqOutput = [SeqAfterIFFT(N-Ncp+1:end),SeqAfterIFFTRepeat];
