function [SNR_dB] = Receiver_SNR_GUI
%% Receiving Part Code
%% Add files path
addpath('.\Common_file');
addpath('.\Channel_model');
addpath('.\Channel_model\Parameter_Cal_Result');
%% Global Parameter
N= load('Parameter_Cal_Result\N.txt');
M = load('Parameter_Cal_Result\M.txt'); %16QAM modulation scheme
Ncp = load('Parameter_Cal_Result\Ncp.txt');
RepeatTime = 2;
FrameNum = load('Parameter_Cal_Result\FrameNum.txt');  % ofdm symbol number of each package transmitted
TotalFrameNum = load('Parameter_Cal_Result\FrameNum.txt'); % total ofdm symbol number you want to get
SubCarrierNum = load('Parameter_Cal_Result\SubCarrierNum.txt'); % every OFDM frame has 63 subcarriers
LowPaddingNum = load('Parameter_Cal_Result\LowPaddingNum.txt'); % padding number of low frequency
HighPaddingNum = load('Parameter_Cal_Result\HighPaddingNum.txt'); % padding number of high frequency
ThresholdValue = 35;

WORD_LENGTH = load('Parameter_Cal_Result\DAC_RES.txt');
FRACTION_LENGTH = load('Parameter_Cal_Result\FRAC_LENGTH.txt');
%% read the data from file
fid_1=fopen('..\Matlab_ADS_Data\BasebandInput_SNR.tim','r'); % Original output signal from ADS. 
TIME_AND_VAR_rece =[]; 
while 1     
	tline=fgetl(fid_1);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	TIME_AND_VAR_rece = [TIME_AND_VAR_rece;tline];
end
fclose(fid_1);
DigitalReceived = TIME_AND_VAR_rece(1:end,2).';% Save received 'input_SNR' data only into one row.

% downsampling
span = load('Parameter_Cal_Result\span.txt'); % span point of raised cosine FIR filter
OverSampling = load('.\Channel_model\Parameter_Cal_Result\OverSampling.txt'); % Samples in one symbol.
rolloff = load('.\Channel_model\Parameter_Cal_Result\rolloff.txt');
filter_h = rcosdesign(rolloff, span, OverSampling);
DigitalReceivedDownSampling = upfirdn(DigitalReceived,filter_h,1,OverSampling);

% add noise here
% DigitalReceivedDownSampling = awgn(DigitalReceivedDownSampling,20,'measured'); % 20dB SNR.

fid_2=fopen('..\Matlab_ADS_Data\DataSymbolQuanReshapeReal_SNR.txt','r'); % Save real part of QAM signal.
DataSymbolQuanReshapeReal =[]; 
while 1     
	tline=fgetl(fid_2);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	DataSymbolQuanReshapeReal = [DataSymbolQuanReshapeReal;tline];
end
fclose(fid_2);

fid_3=fopen('..\Matlab_ADS_Data\DataSymbolQuanReshapeImag_SNR.txt','r'); % Save image part of QAM signal.
DataSymbolQuanReshapeImag =[]; 
while 1     
	tline=fgetl(fid_3);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	DataSymbolQuanReshapeImag = [DataSymbolQuanReshapeImag;tline];
end
fclose(fid_3);

DataSymbolQuanReshape = complex(DataSymbolQuanReshapeReal,DataSymbolQuanReshapeImag); % Combine Real data with Imag data.

fid_4=fopen('..\Matlab_ADS_Data\DataSerialOri_SNR.txt','r');
DataSerialOri =[]; 
while 1     
	tline=fgetl(fid_4);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	DataSerialOri = [DataSerialOri;tline];
end
fclose(fid_4);
[SNR_dB] = Receiving_Realtime_SNR(N,M,Ncp,RepeatTime,FrameNum,TotalFrameNum,SubCarrierNum,LowPaddingNum,HighPaddingNum,ThresholdValue,DigitalReceivedDownSampling,DataSymbolQuanReshape,DataSerialOri)