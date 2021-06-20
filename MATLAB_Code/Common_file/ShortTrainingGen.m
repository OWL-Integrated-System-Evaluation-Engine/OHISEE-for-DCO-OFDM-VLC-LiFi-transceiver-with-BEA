function [TrainingSeqOutput] = ShortTrainingGen(GainSubcarrierLinear, Channel_Gain)
%clc;
%clear all;
%close all;
%% Parameter list:
LowPaddingNum = 5; % padding number of low frequency
HighPaddingNum = 5; % padding number of high frequency
N = 128;
Ncp = 32;
%WORD_LENGTH = load('..\Channel_model\Parameter_Cal_Result\DAC_RES.txt');
%FRACTION_LENGTH = load('..\Channel_model\Parameter_Cal_Result\FRAC_LENGTH.txt');

WORD_LENGTH = load('DAC_RES.txt');
FRACTION_LENGTH = load('FRAC_LENGTH.txt');

TrainingSeq = sqrt(13/6) .* [0,0,1+1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,0,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0];
TrainingSeqAddPadding = [zeros(1,LowPaddingNum),TrainingSeq,zeros(1,HighPaddingNum)];
% Hermitian Symmtry
if Channel_Gain
	TrainingSeqHermitianSymmetry = [0,TrainingSeqAddPadding.*GainSubcarrierLinear(2:N/2),0,fliplr(conj(TrainingSeqAddPadding.*GainSubcarrierLinear(2:N/2)))];
else
	TrainingSeqHermitianSymmetry = [0,TrainingSeqAddPadding,0,fliplr(conj(TrainingSeqAddPadding))];
end
% IFFT transformation
TrainingSeqAfterIFFT = ifft(TrainingSeqHermitianSymmetry,N);
TrainingSeqAddingCP = [TrainingSeqAfterIFFT(1,N-Ncp+1:N),TrainingSeqAfterIFFT(1,N-Ncp+1:N),TrainingSeqAfterIFFT,TrainingSeqAfterIFFT];

QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
TrainingSeqOutput = quantize(QuantizerInst,TrainingSeqAddingCP);

%% Used to verify the design of FPGA
%TrainingSeqAddingCPReal = real(TrainingSeqAddingCP);
%TrainingSeqAddingCPImag = imag(TrainingSeqAddingCP);
%[TrainingSeqAddingCPRealQuan,~,TrainingSeqAddingCPRealBinary] = DACInputGen(TrainingSeqAddingCPReal, WORD_LENGTH, FRACTION_LENGTH);
%[TrainingSeqAddingCPImagQuan,~,TrainingSeqAddingCPImagBinary] = DACInputGen(TrainingSeqAddingCPImag, WORD_LENGTH, FRACTION_LENGTH);
 