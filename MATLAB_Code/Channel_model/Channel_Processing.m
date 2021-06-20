function [FrameOutput] = Channel_Processing(DAC_Code)
% this file is the channel model
% this part only keep the data processing
% the parameter calculation is moved into other function to save the simulation time
% Four parts are included
% DC bias, LED, Channel Model, Photo detector are included
% Input of this function is one OFDM frame 
% Output of this function is the OFDM frame received by photo detector

%% Variable List: 
% FrameOutput: OFDM frame received by photodetector
% DAC_Code: Digital to Analog input.
% DISTANCE: distance between transmitter and receiver
% TRANSMITTER_ANGLE: Radius is used, not degree! using pi to represent 180
% RECEIVER_ANGLE: Radius is used, not degree! using pi to represent 180
%% Parameter List: These parameters can be set as the input of the function if need
DAC_MAX_ABS = load('Parameter_Cal_Result\DAC_MAX_ABS.txt'); % Maximum absolute current of DAC output, Unit: A
DC_BIAS = load('Parameter_Cal_Result\DC_BIAS.txt'); % DC bias, decided by the property of LED, Unit: A
REFER_FLUX = load('Parameter_Cal_Result\REFER_FLUX.txt'); % typical luminous flux of LED is 425lm at 350mA, Unit: lm
WORD_LENGTH = load('Parameter_Cal_Result\DAC_RES.txt');
FRACTION_LENGTH = load('Parameter_Cal_Result\FRAC_LENGTH.txt');
%% First part: DAC model and DC bias model
% transform digital to the DAC to the dac analog current output
I_OUTA = DAC_Code/(2^WORD_LENGTH) * DAC_MAX_ABS;
I_OUTB = ((2^WORD_LENGTH-1-DAC_Code)/(2^WORD_LENGTH)) * DAC_MAX_ABS;
% transform current to voltage output of amplifier and LED driver
LED_DRIVER_GAIN = load('Parameter_Cal_Result\LED_DRIVER_GAIN.txt');
V_LED_DRIVER = LED_DRIVER_GAIN * (I_OUTA - I_OUTB) + DC_BIAS;
%% Second part: LED model

% First part: non-linear effect
% transform voltage to LED current
fid_r_0=fopen('Parameter_Cal_Result\Volt2CurrPoly.txt','r');
Volt2CurrPoly =[]; 
while 1     
	tline=fgetl(fid_r_0);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	Volt2CurrPoly = [Volt2CurrPoly,tline];
end
fclose(fid_r_0);
CurrentIntoLed = polyval(Volt2CurrPoly,V_LED_DRIVER);
% transform current into LED to the output power
% read the parameter of channel model from txt file
fid_r_1=fopen('Parameter_Cal_Result\Curr2RelaFluxPoly.txt','r');
Curr2RelaFluxPoly =[]; 
while 1     
	tline=fgetl(fid_r_1);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	Curr2RelaFluxPoly = [Curr2RelaFluxPoly,tline];
end
fclose(fid_r_1);
LedFluxRelative = polyval(Curr2RelaFluxPoly,CurrentIntoLed);
LedFluxAbs = REFER_FLUX*LedFluxRelative;
 
% Second part: bandwidth limitation effect
% apply filter on the relavant and absolute output flux
%read the parameter of channel model from txt file
fid_r_2=fopen('Parameter_Cal_Result\FilterCof.txt','r');
%fid_1=fopen('baseband_output_imag.tim','r'); 
FilterCof =[]; 
while 1     
	tline=fgetl(fid_r_2);     
	if ~ischar(tline),break;
	end     
	tline=str2num(tline);     
	FilterCof = [FilterCof,tline];
end
fclose(fid_r_2);
LedFluxRelativeBandLim = filtfilt(FilterCof,1,CurrentIntoLed);

%% Third part: LED effect, path loss and PD effect
PdCurrentRefer = load('Parameter_Cal_Result\PdCurrentRefer.txt'); 
PdCurrent = 9 * PdCurrentRefer.*LedFluxRelative; % using bandwidth limitation version

% Forth part: PD noise model
% Noise adding processing
% Quantum shot noise
%NOISE_VARIANCE_DARK_CURRENT = load('Parameter_Cal_Result\NOISE_VARIANCE_DARK_CURRENT.txt');
%NOISE_VARIANCE_THERMAL = load('Parameter_Cal_Result\NOISE_VARIANCE_THERMAL.txt');
q = load('Parameter_Cal_Result\q.txt');
BANDWIDTH = load('Parameter_Cal_Result\BANDWIDTH.txt');
AVERAGE_RECEIVED_CURRENT = sum(PdCurrent)/length(PdCurrent);
NOISE_VARIANCE_QUANTUM = 2*q*AVERAGE_RECEIVED_CURRENT*BANDWIDTH;
PD_NOISE_DATASHEET = load('Parameter_Cal_Result\PD_NOISE_DATASHEET.txt');
TIA_NOISE_DATASHEET = load('Parameter_Cal_Result\TIA_NOISE_DATASHEET.txt');
% Total noise variance calculation
NOISE_VARIANCE_TOTAL = NOISE_VARIANCE_QUANTUM + PD_NOISE_DATASHEET^2 + TIA_NOISE_DATASHEET^2;
NoiseCurrent = wgn(1,length(PdCurrent),NOISE_VARIANCE_TOTAL,'linear','real');
FrameOutput = PdCurrent + NoiseCurrent - AVERAGE_RECEIVED_CURRENT;
% FrameOutput = PdCurrent - AVERAGE_RECEIVED_CURRENT;% Used for testing, received current without average current