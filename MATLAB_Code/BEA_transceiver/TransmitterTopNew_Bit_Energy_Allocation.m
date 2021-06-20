function [Output2DAC,DataSymbolQuan,DataSerial]= TransmitterTopNew_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,WORD_LENGTH,FRACTION_LENGTH)

%% Parameter Aera
POLY_LENGTH = 7;
POLY_TAP = 1;
%Pilot_pos = [7,21,43,57]; % this is the position of pilot
NBITS = 1; %Data width of parallel output
NUM = FrameNum*RSum; %Byte number

GainSubcarrierLinear =[];
%% Parameters Checking
if N ~= (SubCarrierNum + LowPaddingNum + HighPaddingNum + 1) * 2 
	error('The number of subcarrier is invalid');
else
	disp('SubCarrierNum, LowPaddingNum, HighPaddingNum and PILOT_NUM checked');
end

if sum(RArray) ~= RSum
	error('The bit allcation result is not consist with the design');
else 
	disp('RArray and RSum checked');
end
disp('Transmitter Simulation begins ...');
[DataLogic,~] = prbs_gen(POLY_LENGTH, POLY_TAP, NBITS, NUM);
DataSerial = double(DataLogic);
DataSerialReshape = reshape(DataSerial, RSum, FrameNum);
%% Symbol Mapping, each loop processes one subcarrier.
IndexLastTime = 0;
DataSymbol = [];
for nSubcarrier = 1:SubCarrierNum
    if RArray(nSubcarrier)==0
        DataSymbolOneSubcarrier = zeros(1,FrameNum);
    else    
        DataSerialOneSubcarrier = DataSerialReshape(IndexLastTime+1:IndexLastTime+RArray(nSubcarrier),:);
        IndexLastTime =  sum(RArray(1:nSubcarrier));
        DataSymbolOneSubcarrier =  qammod(DataSerialOneSubcarrier,2^RArray(nSubcarrier),'gray','InputType','bit','UnitAveragePower',true) .* (sqrt(SArray(nSubcarrier)/2));
    end
	DataSymbol = [DataSymbol;DataSymbolOneSubcarrier];
end

QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
DataSymbolQuan = quantize(QuantizerInst,DataSymbol);

% draw the constellation of the tranmistted signal.
for fig_index = 1:SubCarrierNum
    fig = figure(fig_index);
    scatter(real(DataSymbolQuan(fig_index,:)),imag(DataSymbolQuan(fig_index,:)),'.');
    xlim([min(-1,min(real(DataSymbolQuan(fig_index,:)))),max(1,max(real(DataSymbolQuan(fig_index,:))))]);
    ylim([min(-1,min(imag(DataSymbolQuan(fig_index,:)))),max(1,max(imag(DataSymbolQuan(fig_index,:))))]);
    set(fig,'visible','off');
    saveas(fig,fullfile('.\Generated_pic\Transmitted_Constellation\',[num2str(fig_index),'th','Subcarrier.jpg']));
end

DataSymbolAddingCP = [];
for nFrame = 1:FrameNum
	% Pick up one symbol to process
	DataSymbolOneFrame =  DataSymbolQuan(:,nFrame);
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
[LongTrainingSeqOutput,~,~] = TrainingSeqGenParkMethod(N,Ncp,RepeatTime,GainSubcarrierLinear, 0);
%% Adding Training Sequence to output signal
DataPackage = [ShortTrainingSeqOutput,LongTrainingSeqOutput,DataSymbolAddingCPQuan];
%% Transmit the data package to DAC code
BIT_OFFSET = 2; % the final data need to right shift BIT_OFFSET to match the fully utilize the digital to analog converter.
Output2DAC = (DataPackage *2^BIT_OFFSET + 2^(WORD_LENGTH - FRACTION_LENGTH - 1)) .* (2^(FRACTION_LENGTH));