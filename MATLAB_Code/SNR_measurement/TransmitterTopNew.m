function [Output2DAC,DataSymbolQuanReshape,DataSerial]= TransmitterTopNew(N,M,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,WORD_LENGTH,FRACTION_LENGTH)
%clear all;
%close all;
%clc;

% % Parameters used for test
% N = 128;
% M = 4;
% Ncp = 32;
% RepeatTime = 2;
% FrameNum = 100;
% LowPaddingNum = 0;
% HighPaddingNum = 0;
% SumcarrierNum = 59;
% WORD_LENGTH = 14;
% FRACTION_LENGTH = 12; 
% or it is not added. including short training sequence, long training sequence and effective message.
%% Parameter Aera
POLY_LENGTH = 7;
POLY_TAP = 1;
NBITS = 1; %Data width of parallel output
NUM = FrameNum*SubCarrierNum*M; %Byte number
GainSubcarrierLinear =[];
%% Parameters Checking
if N ~= (SubCarrierNum + LowPaddingNum + HighPaddingNum + 1) * 2
	error('The number of subcarrier is invalid');
else
	disp('Parameter checed, simulation begins ...');
end

[DataLogic,~] = prbs_gen(POLY_LENGTH, POLY_TAP, NBITS, NUM);
DataSerial = double(DataLogic);

%% Symbol Mapping
DataSymbol =  qammod(DataSerial',2^M,'gray','InputType','bit','UnitAveragePower',true);
DataSymbolReal = real(DataSymbol);
DataSymbolImag = imag(DataSymbol);
[DataSymbolRealQuan,~,DataSymbolRealBinary] = DACInputGen(DataSymbolReal, WORD_LENGTH, FRACTION_LENGTH);
[DataSymbolImagQuan,~,DataSymbolImagBinary] = DACInputGen(DataSymbolImag, WORD_LENGTH, FRACTION_LENGTH);
QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
DataSymbolQuan = quantize(QuantizerInst,DataSymbol);

DataSymbolQuanReshape = reshape(DataSymbolQuan,SubCarrierNum,FrameNum);
DataSymbolAddingCP = [];

% register used in scrambler
ScramblerRegister = ones(1,7);

for nFrame = 1:FrameNum
	% scrambler data generation
	ScramblerData = xor(ScramblerRegister(4),ScramblerRegister(7));
	ScramblerRegister = [ScramblerData,ScramblerRegister(2:7)];
	% pilot generation
	if ScramblerData == 0
		PilotSeq = [1,-1,1,1];
    elseif ScramblerData ==1
		PilotSeq = [-1,1,-1,-1];
	else
		error('ScramblerData is invalid');
    end
	% Pick up one symbol to process
	DataSymbolOneFrame =  DataSymbolQuanReshape(:,nFrame);
	% Hermitian Symmetry and IFFT transformation
	% Adding Padding and pilot to the symbol
	DataSymbolOneFrameAddPadiing = [zeros(LowPaddingNum,1);DataSymbolOneFrame;zeros(HighPaddingNum,1)];
	DataSymbolOneFrameHermitianSymmetry = [0;DataSymbolOneFrameAddPadiing;0;flipud(conj(DataSymbolOneFrameAddPadiing))];
	%% IFFT transformation
	DataSymbolOneFrameAfterIFFT = ifft(DataSymbolOneFrameHermitianSymmetry,N,1); %Every coloum one IFFT symbol

	%% Adding CP
	DataSymbolOneFrameAddingCP = [DataSymbolOneFrameAfterIFFT(N-Ncp+1:N);DataSymbolOneFrameAfterIFFT];
	DataSymbolAddingCP = [DataSymbolAddingCP,DataSymbolOneFrameAddingCP];
    end
QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
DataSymbolAddingCPQuan = quantize(QuantizerInst,reshape(DataSymbolAddingCP,1,size(DataSymbolAddingCP,1)*size(DataSymbolAddingCP,2)));
%% Short Training Sequence Generation
[ShortTrainingSeqOutput] = ShortTrainingGen(GainSubcarrierLinear, 0);
%% Long Training Sequence Generation
[LongTrainingSeqOutput,~,~] = TrainingSeqGenParkMethod(N,Ncp,RepeatTime,GainSubcarrierLinear,0);
%% Adding Training Sequence to output signal
DataPackage = [ShortTrainingSeqOutput,LongTrainingSeqOutput,DataSymbolAddingCPQuan];
%% Transmit the data package to DAC code
BIT_OFFSET = 2; % the final data need to right shift BIT_OFFSET to match the fully utilize the digital to analog converter.
Output2DAC = (DataPackage *2^BIT_OFFSET + 2^(WORD_LENGTH - FRACTION_LENGTH - 1)) .* (2^(FRACTION_LENGTH));