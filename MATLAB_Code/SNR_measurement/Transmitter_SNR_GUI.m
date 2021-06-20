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
FrameNum = load('Parameter_Cal_Result\FrameNum.txt');  % ofdm symbol number of each package transmitted
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
% transform digital to the DAC to the dac analog current output
DAC_MAX_ABS = load('Parameter_Cal_Result\DAC_MAX_ABS.txt'); % Maximum absolute current of DAC output, Unit: A
I_OUTA = Output2DAC/(2^WORD_LENGTH) * DAC_MAX_ABS;
I_OUTB = ((2^WORD_LENGTH-1-Output2DAC)/(2^WORD_LENGTH)) * DAC_MAX_ABS;
% transform current to voltage output of amplifier and LED driver
LED_DRIVER_GAIN = load('Parameter_Cal_Result\LED_DRIVER_GAIN.txt');
V_LED_DRIVER = LED_DRIVER_GAIN * (I_OUTA - I_OUTB);
%% save the generated data to the files.
% save DigitalOutputRoll_64QAM
TIME_AND_VAR = [TimeArray.', V_LED_DRIVER.'];
% upsampling through raised cosine FIR filter
span = load('Parameter_Cal_Result\span.txt'); % span point of raised cosine FIR filter 
OverSampling = load('.\Channel_model\Parameter_Cal_Result\OverSampling.txt');
rolloff = load('.\Channel_model\Parameter_Cal_Result\rolloff.txt');
% shape filter generation
filter_h = rcosdesign(rolloff, span, OverSampling);
V_LED_DRIVER_UP_SAMPLING = upfirdn(V_LED_DRIVER,filter_h,OverSampling);
TimeArrayUpsampling = 0:T_S/OverSampling:(length(V_LED_DRIVER_UP_SAMPLING)-1)*(T_S/OverSampling);
TIME_AND_VAR_UP_SAMPLING = [TimeArrayUpsampling.', V_LED_DRIVER_UP_SAMPLING.'];
if(exist('..\Matlab_ADS_Data\BasebandOutput_SNR.tim'))
    delete('..\Matlab_ADS_Data\BasebandOutput_SNR.tim');
end
fid_w_1 = fopen('.\..\Matlab_ADS_Data\BasebandOutput_SNR.tim','w');
fprintf(fid_w_1,'%s\n','BEGIN TIMEDATA');
fprintf(fid_w_1,'%s\n','% time voltage');
[row,col] = size(TIME_AND_VAR_UP_SAMPLING);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_1,'%e\n',TIME_AND_VAR_UP_SAMPLING(i,j)); %using scientific notation
		else
			fprintf(fid_w_1,'%e\t',TIME_AND_VAR_UP_SAMPLING(i,j)); %using scientific notation 
		end
	end
end
fprintf(fid_w_1,'%s','END');
fclose(fid_w_1);

if(exist('..\Matlab_ADS_Data\DataSymbolQuanReshapeReal_SNR.txt'))
    delete('..\Matlab_ADS_Data\DataSymbolQuanReshapeReal_SNR.txt');
end
DataSymbolQuanReshapeReal = real(DataSymbolQuanReshape);
DataSymbolQuanReshapeImag = imag(DataSymbolQuanReshape);
fid_w_2 = fopen('.\..\Matlab_ADS_Data\DataSymbolQuanReshapeReal_SNR.txt','w');
[row,col] = size(DataSymbolQuanReshapeReal);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_2,'%e\n',DataSymbolQuanReshapeReal(i,j)); %using scientific notation
		else
			fprintf(fid_w_2,'%e\t',DataSymbolQuanReshapeReal(i,j)); %using scientific notation 
		end
	end
end
fclose(fid_w_2);

if(exist('..\Matlab_ADS_Data\DataSymbolQuanReshapeImag_SNR.txt.txt'))
    delete('..\Matlab_ADS_Data\DataSymbolQuanReshapeImag_SNR.txt.txt');
end
fid_w_3 = fopen('.\..\Matlab_ADS_Data\DataSymbolQuanReshapeImag_SNR.txt','w');
[row,col] = size(DataSymbolQuanReshapeImag);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_3,'%e\n',DataSymbolQuanReshapeImag(i,j)); %using scientific notation
		else
			fprintf(fid_w_3,'%e\t',DataSymbolQuanReshapeImag(i,j)); %using scientific notation 
		end
	end
end
fclose(fid_w_3);

if(exist('..\Matlab_ADS_Data\DataSerialOri_SNR.txt'))
    delete('..\Matlab_ADS_Data\DataSerialOri_SNR.txt');
end
fid_w_4 = fopen('.\..\Matlab_ADS_Data\DataSerialOri_SNR.txt','w');
[row,col] = size(DataSerialOri);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_4,'%e\n',DataSerialOri(i,j)); %using scientific notation
		else
			fprintf(fid_w_4,'%e\t',DataSerialOri(i,j)); %using scientific notation 
		end
	end
end
fclose(fid_w_4);

delete('.\Channel_model\Parameter_Cal_Result\ChipRate.txt');
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
disp("Signal generation for SNR measurement done.");

