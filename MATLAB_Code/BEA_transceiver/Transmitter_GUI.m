function[DataSymbolQuan,TimeArray,Output2DAC,V_LED_DRIVER_UP_SAMPLING,PAPR_VALUE] = Transmitter_GUI
%% Add files path
addpath('.\Common_file');
addpath('.\Channel_model');
addpath('.\Channel_model\Parameter_Cal_Result');

N= load('Parameter_Cal_Result\N.txt');
SAMPLING_RATE = load('Parameter_Cal_Result\SAMPLING_RATE.txt');
T_S = 1/SAMPLING_RATE;
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
% intermedia varibles
Output2DAC = [];
DataSerialOri = [];
[Output2DAC,DataSymbolQuan,DataSerialOri]= TransmitterTopNew_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,WORD_LENGTH,FRACTION_LENGTH)
%%[Output2DAC,~,DataSerialOri]= TransmitterTopNew_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,LowPaddingNum,HighPaddingNum,SubCarrierNum,WORD_LENGTH,FRACTION_LENGTH);
% write the result in the file
TimeArray = 0:T_S:(length(Output2DAC)-1)*T_S;
% transform digital to the DAC to the dac analog current output
DAC_MAX_ABS = load('Parameter_Cal_Result\DAC_MAX_ABS.txt'); % Maximum absolute current of DAC output, Unit: A
I_OUTA = Output2DAC/(2^WORD_LENGTH) * DAC_MAX_ABS;
I_OUTB = ((2^WORD_LENGTH-1-Output2DAC)/(2^WORD_LENGTH)) * DAC_MAX_ABS;
% transform current to voltage output of amplifier and LED driver
LED_DRIVER_GAIN = load('Parameter_Cal_Result\LED_DRIVER_GAIN.txt');
V_LED_DRIVER = LED_DRIVER_GAIN * (I_OUTA - I_OUTB);
TIME_AND_VAR = [TimeArray.', V_LED_DRIVER.'];
% upsampling
span = load('Parameter_Cal_Result\span.txt'); % span point of raised cosine FIR filter 
OverSampling = load('.\Channel_model\Parameter_Cal_Result\OverSampling.txt');
rolloff = load('.\Channel_model\Parameter_Cal_Result\rolloff.txt');
% shape filter generation
filter_h = rcosdesign(rolloff, span, OverSampling);
V_LED_DRIVER_UP_SAMPLING = upfirdn(V_LED_DRIVER,filter_h,OverSampling);
TimeArrayUpsampling = 0:T_S/OverSampling:(length(V_LED_DRIVER_UP_SAMPLING)-1)*(T_S/OverSampling);
TIME_AND_VAR_UP_SAMPLING = [TimeArrayUpsampling.', V_LED_DRIVER_UP_SAMPLING.'];
PAPR_VALUE = PAPR_cal(V_LED_DRIVER_UP_SAMPLING);
if(exist('..\Matlab_ADS_Data\BasebandOutput.tim'))
    delete('..\Matlab_ADS_Data\BasebandOutput.tim');
end
fid_w_0 = fopen('.\..\Matlab_ADS_Data\BasebandOutput.tim','w');
fprintf(fid_w_0,'%s\n','BEGIN TIMEDATA');
fprintf(fid_w_0,'%s\n','% time voltage');
[row,col] = size(TIME_AND_VAR_UP_SAMPLING);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_0,'%e\n',TIME_AND_VAR_UP_SAMPLING(i,j)); %using scientific notation
		else
			fprintf(fid_w_0,'%e\t',TIME_AND_VAR_UP_SAMPLING(i,j)); %using scientific notation 
		end
	end
end
fprintf(fid_w_0,'%s','END');
fclose(fid_w_0);

if(exist('..\Matlab_ADS_Data\DataSerialOri.txt'))
    delete('..\Matlab_ADS_Data\DataSerialOri.txt');
end
fid_w_1 = fopen('.\..\Matlab_ADS_Data\DataSerialOri.txt','w');
[row,col] = size(DataSerialOri);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_1,'%e\n',DataSerialOri(i,j)); %using scientific notation
		else
			fprintf(fid_w_1,'%e\t',DataSerialOri(i,j)); %using scientific notation 
		end
	end
end
fclose(fid_w_1);
disp("Transmitter Simulation done.");
delete('.\Channel_model\Parameter_Cal_Result\ChipRate.txt');

% update the parameters file to ADS
ChipRate = SAMPLING_RATE;
Num = length(TIME_AND_VAR_UP_SAMPLING);
ADSVAR=[1 OverSampling ChipRate Num];
fp=fopen('.\Channel_model\Parameter_Cal_Result\ChipRate.txt', 'a');
fprintf(fp,'%s\n','BEGIN ADSVAR');
fprintf(fp,'%s\n','% index(real) OverSampling(real) ChipRate(real) Num(real)');
fprintf(fp,'%g\t',ADSVAR);
fprintf(fp,'\n');
fprintf(fp,'%s','END');
fclose(fp);
