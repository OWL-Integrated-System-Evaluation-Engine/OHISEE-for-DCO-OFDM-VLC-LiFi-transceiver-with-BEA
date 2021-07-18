 function [] = Parameter_Load(N,Ncp,SubCarrierNum,LowPaddingNum,HighPaddingNum,FrameNum,SAMPLING_RATE,OverSampling,target_ser,DAC_RES,FRAC_LENGTH)

%% Add file path
addpath('.\Channel_model\Measurement_Result');
addpath('.\Channel_model\Parameter_Cal_Result');

% N = 128;
save '.\Channel_model\Parameter_Cal_Result\N.txt' -ascii N;
M = 2;
save '.\Channel_model\Parameter_Cal_Result\M.txt' -ascii M;
%Ncp = 32;
save '.\Channel_model\Parameter_Cal_Result\Ncp.txt' -ascii Ncp;
%SubCarrierNum = 63;% every OFDM frame has 63 subcarriers
save '.\Channel_model\Parameter_Cal_Result\SubCarrierNum.txt' -ascii SubCarrierNum;
%LowPaddingNum = 0; % padding number of low frequency
save '.\Channel_model\Parameter_Cal_Result\LowPaddingNum.txt' -ascii LowPaddingNum;
%HighPaddingNum = 0; % padding number of high frequency
save '.\Channel_model\Parameter_Cal_Result\HighPaddingNum.txt' -ascii HighPaddingNum;
RepeatTime = 2;
save '.\Channel_model\Parameter_Cal_Result\RepeatTime.txt' -ascii RepeatTime;
%FrameNum = 400;
save '.\Channel_model\Parameter_Cal_Result\FrameNum.txt' -ascii FrameNum;
TotalFrameNum = FrameNum;
save '.\Channel_model\Parameter_Cal_Result\TotalFrameNum.txt' -ascii TotalFrameNum;
%SAMPLING_RATE = 50e6; % padding number of high frequency
save '.\Channel_model\Parameter_Cal_Result\SAMPLING_RATE.txt' -ascii SAMPLING_RATE;
%OverSampling = 4;
save '.\Channel_model\Parameter_Cal_Result\OverSampling.txt' -ascii OverSampling;
span = 12;
save '.\Channel_model\Parameter_Cal_Result\span.txt' -ascii span;
rolloff = 0.25;
save '.\Channel_model\Parameter_Cal_Result\rolloff.txt' -ascii rolloff;
%target_ser = 1e-4;
save '.\Channel_model\Parameter_Cal_Result\target_ser.txt' -ascii target_ser;
disp('Parameter load done')
%DAC_RES = 16;
%save '.\Channel_model\Parameter_Cal_Result\DAC_RES.txt' -ascii DAC_RES;
%FRAC_LENGTH = 14;
%save '.\Channel_model\Parameter_Cal_Result\FRAC_LENGTH.txt' -ascii FRAC_LENGTH;
