clear all;
close all;
clc;
%% Add files path
addpath('.\Common_file');
addpath('.\Channel_model');
addpath('.\Channel_model\Parameter_Cal_Result');
%% Global Parameter
N= load('Parameter_Cal_Result\N.txt');
SAMPLING_RATE = load('Parameter_Cal_Result\SAMPLING_RATE.txt');
T_S = 1/SAMPLING_RATE;
M = load('Parameter_Cal_Result\M.txt'); %16QAM modulation scheme
Ncp = load('Parameter_Cal_Result\Ncp.txt');
RepeatTime = 2;
FrameNum = load('Parameter_Cal_Result\FrameNum.txt'); % ofdm symbol number of each package transmitted
TotalFrameNum = load('Parameter_Cal_Result\FrameNum.txt'); % total ofdm symbol number you want to get
SubCarrierNum = load('Parameter_Cal_Result\SubCarrierNum.txt'); % every OFDM frame has 63 subcarriers
LowPaddingNum = load('Parameter_Cal_Result\LowPaddingNum.txt'); % padding number of low frequency
HighPaddingNum = load('Parameter_Cal_Result\HighPaddingNum.txt'); % padding number of high frequency
ThresholdValue = 55;

WORD_LENGTH = load('Parameter_Cal_Result\DAC_RES.txt');
FRACTION_LENGTH = load('Parameter_Cal_Result\FRAC_LENGTH.txt');

% intermediate variables
Output2DAC = [];
DataSymbolQuanReshape = [];
DataSerialOri = [];
[Output2DAC,DataSymbolQuanReshape,DataSerialOri] = TransmitterTopNew(N,M,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,WORD_LENGTH,FRACTION_LENGTH);
TimeArray = 0:T_S:(length(Output2DAC)-1)*T_S;
%% save .dat file
% save DigitalOutputRoll_64QAM
TIME_AND_VAR = [TimeArray.', Output2DAC.'];
fid_w_1 = fopen('..\Matlab_ADS_Data\BasebandOutput_SNR.tim','w');
fprintf(fid_w_1,'%s\n','BEGIN TIMEDATA');
fprintf(fid_w_1,'%s\n','% time voltage');
[row,col] = size(TIME_AND_VAR);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_1,'%e\n',TIME_AND_VAR(i,j)); %using scientific notation
		else
			fprintf(fid_w_1,'%e\t',TIME_AND_VAR(i,j)); %using scientific notation 
		end
	end
end
fprintf(fid_w_1,'%s','END');
fclose(fid_w_1);

%% Receiving Part Code
% read the data from file
fid_1=fopen('..\Matlab_ADS_Data\BasebandInput_SNR.tim','r');
TIME_AND_VAR_rece =[]; 
while 1     
	tline=fgetl(fid_1);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	TIME_AND_VAR_rece = [TIME_AND_VAR_rece;tline];
end
DigitalReceived = TIME_AND_VAR_rece(1:end,2).';
Receiving_Realtime_SNR(N,M,Ncp,RepeatTime,FrameNum,TotalFrameNum,SubCarrierNum,LowPaddingNum,HighPaddingNum,ThresholdValue,DigitalReceivedDownSampling,DataSymbolQuanReshape,DataSerialOri)
