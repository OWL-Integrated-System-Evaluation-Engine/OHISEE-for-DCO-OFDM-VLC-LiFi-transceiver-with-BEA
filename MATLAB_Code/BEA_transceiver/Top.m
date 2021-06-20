clear all;
close all;
clc;

%% Add files path
addpath('./Common_file');
addpath('./Channel_model');
addpath('./Channel_model/Parameter_Cal_Result');

%% Global Parameter
N= load('Parameter_Cal_Result\N.txt');
SAMPLING_RATE = load('Parameter_Cal_Result\SAMPLING_RATE.txt');
T_S = 1/SAMPLING_RATE;
Ncp = load('Parameter_Cal_Result\Ncp.txt');
RepeatTime = load('Parameter_Cal_Result\RepeatTime.txt');
FrameNum = load('Parameter_Cal_Result\FrameNum.txt'); % ofdm symbol number of each package transmitted
TotalFrameNum = load('Parameter_Cal_Result\TotalFrameNum.txt'); % total ofdm symbol number you want to get
SubCarrierNum = load('Parameter_Cal_Result\SubCarrierNum.txt'); % every OFDM frame has 59 subcarriers
LowPaddingNum = load('Parameter_Cal_Result\LowPaddingNum.txt'); % padding number of low frequency
HighPaddingNum = load('Parameter_Cal_Result\HighPaddingNum.txt'); % padding number of high frequency
ThresholdValue = 55;
WORD_LENGTH = load('Parameter_Cal_Result\DAC_RES.txt');
FRACTION_LENGTH = load('Parameter_Cal_Result\FRAC_LENGTH.txt');

% Calcuate RArray and SArray
RArray = [];
SArray = [];
% load gn value from SNR measurement result
SER = 1e-7;
gn =[];
fid_r_0=fopen('Parameter_Cal_Result/SNRLinearResult.txt','r'); 
while 1     
    tline=fgetl(fid_r_0);     
    if ~ischar(tline),break;
    end     
    tline=str2num(tline);     
    gn = [gn,tline];
end 
fclose(fid_r_0);
%[SArray,RArray,RSum,ESum,SubCarrierNumUsed,BER]=DMTRA(gn,SubCarrierNum,SER);
[SArray,RArray,RSum,ESum,SubCarrierNumUsed,BER]=LCRA(gn,SubCarrierNum,SER);
disp('Target Bit Error Rate is:')
disp(BER);

% intermedia varibles
Output2DAC = [];
DataSerialOri = [];
[Output2DAC,~,DataSerialOri]= TransmitterTopNew_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,WORD_LENGTH,FRACTION_LENGTH);
[DataRate,BitErrorRate] = Receiving_Realtime_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,TotalFrameNum,SubCarrierNum,LowPaddingNum,HighPaddingNum,ThresholdValue,Output2DAC,DataSerialOri);
