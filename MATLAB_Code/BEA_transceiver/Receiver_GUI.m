function[OFDMSymbolWithoutPaddingAcc,TimeReceived, DigitalReceived, DataRate, BitErrorRate] = Receiver_GUI
%% Add files path
addpath('.\Common_file');
addpath('.\Channel_model');
addpath('.\Channel_model\Parameter_Cal_Result');

N= load('Parameter_Cal_Result\N.txt');
Ncp = load('Parameter_Cal_Result\Ncp.txt');
RSum= load('Parameter_Cal_Result\RSum.txt');
RepeatTime = load('Parameter_Cal_Result\RepeatTime.txt');
FrameNum = load('Parameter_Cal_Result\FrameNum.txt'); % ofdm symbol number of each package transmitted
TotalFrameNum = load('Parameter_Cal_Result\TotalFrameNum.txt'); % total ofdm symbol number you want to get
SubCarrierNum = load('Parameter_Cal_Result\SubCarrierNum.txt'); % every OFDM frame has 63 subcarriers
LowPaddingNum = load('Parameter_Cal_Result\LowPaddingNum.txt'); % padding number of low frequency
HighPaddingNum = load('Parameter_Cal_Result\HighPaddingNum.txt'); % padding number of high frequency
WORD_LENGTH = load('Parameter_Cal_Result\DAC_RES.txt');
FRACTION_LENGTH = load('Parameter_Cal_Result\FRAC_LENGTH.txt');
ThresholdValue = 35;
% read the data from file
RArray = [];
SArray = [];
fid_r_0=fopen('Parameter_Cal_Result\RArray.txt','r'); 
while 1     
    tline=fgetl(fid_r_0);     
    if ~ischar(tline),break;
    end     
    tline=str2num(tline);     
    RArray = [RArray,tline];
end 
fclose(fid_r_0);
fid_r_1=fopen('Parameter_Cal_Result\SArray.txt','r'); 
while 1     
    tline=fgetl(fid_r_1);     
    if ~ischar(tline),break;
    end     
    tline=str2num(tline);     
    SArray = [SArray,tline];
end 
fclose(fid_r_1);

%% read the data from file
fid_1=fopen('..\Matlab_ADS_Data\BasebandInput.tim','r');
TIME_AND_VAR_rece =[]; 
while 1     
    tline=fgetl(fid_1);     
    if ~ischar(tline),break;
    end     
    tline=str2num(tline);     
    TIME_AND_VAR_rece = [TIME_AND_VAR_rece;tline];
end
fclose(fid_1);
TimeReceived = TIME_AND_VAR_rece(1:end,1).';
DigitalReceived = TIME_AND_VAR_rece(1:end,2).';
% downsampling
span = load('Parameter_Cal_Result\span.txt'); % span point of raised cosine FIR filter
OverSampling = load('.\Channel_model\Parameter_Cal_Result\OverSampling.txt');
rolloff =load('.\Channel_model\Parameter_Cal_Result\rolloff.txt');
filter_h = rcosdesign(rolloff, span, OverSampling);
DigitalReceivedDownSampling = upfirdn(DigitalReceived,filter_h,1,OverSampling);
% DigitalReceivedDownSampling = awgn(DigitalReceivedDownSampling,20,'measured');
fid_2=fopen('..\Matlab_ADS_Data\DataSerialOri.txt','r');
DataSerialOri =[]; 
while 1     
    tline=fgetl(fid_2);     
    if ~ischar(tline),break;
    end     
    tline=str2num(tline);     
    DataSerialOri = [DataSerialOri;tline];
end
fclose(fid_2);
[OFDMSymbolWithoutPaddingAcc,DataRate,BitErrorRate] = Receiving_Realtime_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,TotalFrameNum,SubCarrierNum,LowPaddingNum,HighPaddingNum,ThresholdValue,DigitalReceivedDownSampling,DataSerialOri);
