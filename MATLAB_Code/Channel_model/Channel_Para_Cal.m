clear all;
close all;
clc
% this part is used to calcalate the parameter of channel model
%% Add file path
addpath('.\Measurement_Result');
addpath('.\Parameter_Cal_Result');

%% Parameter list
N = 128; % the number of IFFT number
save 'Parameter_Cal_Result\N.txt' -ascii N;
WORD_LENGTH = 14;
save 'Parameter_Cal_Result\DAC_RES.txt' -ascii WORD_LENGTH;
FRACTION_LENGTH = 12;
save 'Parameter_Cal_Result\FRAC_LENGTH.txt' -ascii FRACTION_LENGTH;
DISTANCE = 0.5;
TRANSMITTER_ANGLE = 0;
RECEIVER_ANGLE = 0;
SAMPLING_RATE = 2e7; % Sampling Rate is 20Mhz
save 'Parameter_Cal_Result\SAMPLING_RATE.txt' -ascii SAMPLING_RATE;

DAC_MAX_ABS = 0.02; % Maximum absolute current of DAC output, Unit: A
save 'Parameter_Cal_Result\DAC_MAX_ABS.txt' -ascii DAC_MAX_ABS;
LED_DRIVER_GAIN = 100; % Unit is Ohm, from current to voltage
save 'Parameter_Cal_Result\LED_DRIVER_GAIN.txt' -ascii LED_DRIVER_GAIN;
DC_BIAS = 27.45; % DC bias, decided by the property of LED, Unit: V
save 'Parameter_Cal_Result\DC_BIAS.txt' -ascii DC_BIAS;
REFER_FLUX = 30.6; % typical luminous flux of LED is 425lm at 350mA, Unit: lm
save 'Parameter_Cal_Result\REFER_FLUX.txt' -ascii REFER_FLUX;
HALF_POWER_ANGLE = 48/180*pi; % Radius is used. not degree.
TRANSMISSION_OP_BPFILTER = 1;
GAIN_NON_IMG_CONCENTRATOR = 1;
PD_AREA = 7e-6; % unit m^2
q = 1.6e-19; % electron charge
save 'Parameter_Cal_Result\q.txt' -ascii q;
PD_NPE = 3.5e-15; % noise equivalent power of photodetector, unit: W/Hz^1/2
TIA_NPE = 2.1e-12; % noise equivalent power of tranimpedence amplifier, unit: A/Hz^1/2
% DARK_CURRENT = 1e-9; % parameter of photo detector,unit is A
% k = 1.38e-23;
% R_SH_DIODE = 1e9; %the internal resistance of PD(This is an assumption number)
% TEMPERATURE = 300; % enviroment temperature, Unit: K
%% Second part: LED model
% First Part: no-linear effect
% Calculate the Polynomial from voltage to current
LedVolt9 = xlsread('Measurement_Result\Volt_Current.xlsx','A3:A17');%the measurement voltage of 9 LEDs
LedCurrent9 = xlsread('Measurement_Result\Volt_Current.xlsx','B3:B17');%the measurement current of 9 LEDs
PLOY_ORDER_Volt2Current = 3; % the order of ploynomial to discribe the relationship between voltage and current of LEDs, the order is larger, the model is more accurate
Volt2CurrPoly = polyfit(LedVolt9, LedCurrent9, PLOY_ORDER_Volt2Current);% Calculate the coefficient of polynomial which transform voltage to the current
fid_w_0 = fopen('Parameter_Cal_Result/Volt2CurrPoly.txt','w');
[row,col] = size(Volt2CurrPoly);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_0,'%g\n',Volt2CurrPoly(i,j));
		else
			fprintf(fid_w_0,'%g\t',Volt2CurrPoly(i,j));
		end
	end
end
fclose(fid_w_0);
% Calculate the Polynomial from current to relative output flux
LedSampleCurrent = xlsread('Measurement_Result\Current_RelativeFlux.xlsx','A2:A8'); % This is from datasheet, Unit: mA
RelativeFlux = xlsread('Measurement_Result\Current_RelativeFlux.xlsx','B2:B8'); % This is from datasheet, this is relative flux value
PLOY_ORDER = 3; % the order of ploynomial to discribe the relationship between output current to relative flux, the order is larger, the model is more accurate
Curr2RelaFluxPoly = polyfit(LedSampleCurrent, RelativeFlux, PLOY_ORDER);% Calculate the coefficient of polynomial which transform current to relative flux
% Write the result in the file
fid_w_1 = fopen('Parameter_Cal_Result/Curr2RelaFluxPoly.txt','w');
[row,col] = size(Curr2RelaFluxPoly);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_1,'%g\n',Curr2RelaFluxPoly(i,j));
		else
			fprintf(fid_w_1,'%g\t',Curr2RelaFluxPoly(i,j));
		end
	end
end
fclose(fid_w_1);

% Second part: bandwidth limitation effect
% Below is the measurement result
FreqSamp = xlsread('Measurement_Result\Bandwidth.xlsx','A55:A87');% This sampling rate in our model is 20Mbps, the maximum frequency in model can reach to 10MHz (Nyquist Frequency)
GainDbSamp = xlsread('Measurement_Result\Bandwidth.xlsx','D55:D87');
GainLinearSamp = 10.^(GainDbSamp./10);

% Data fitting to estimate accurate data
NYQUIST_FREQ = SAMPLING_RATE/2; % Half of Sampling Rate
Freq = linspace(0,NYQUIST_FREQ,NYQUIST_FREQ);

[population,~] = fit(FreqSamp.*1e6,GainDbSamp,'pchipinterp');% the unit of FreqSamp in excel is MHz
GainDb = population(Freq);
[population,~] = fit(FreqSamp.*1e6,GainLinearSamp,'pchipinterp');
GainLinear = population(Freq);
% % frequency normalization to design filter
% NormalizedFreq = Freq./NYQUIST_FREQ;
% % filter design, the result is coefficient of the filter
% FILTER_ORDER = 30; % order of bandwidth limitation effect simulation filter
% FilterCof = firls(FILTER_ORDER,NormalizedFreq,GainLinear);
% fvtool(FilterCof);
% % Write the result in the file
% fid_w_1 = fopen('Parameter_Cal_Result/FilterCof.txt','w');
% [row,col] = size(FilterCof);
% for i = 1:row
% 	for j = 1:col
% 		if(j == col)
% 			fprintf(fid_w_1,'%g\n',FilterCof(i,j));
% 		else
% 			fprintf(fid_w_1,'%g\t',FilterCof(i,j));
% 		end
% 	end
% end
% fclose(fid_w_1);

%% calculate the gain in the frequency domain.
plot(Freq,GainLinear);
hold on;
% Energy loss in frequency domain
NormlizedFreq = 0:N/2-1;
CenterFreqSubcarrier = NYQUIST_FREQ*NormlizedFreq./(N/2-1);
GainSubcarrierLinear = population(CenterFreqSubcarrier);
scatter(CenterFreqSubcarrier,GainSubcarrierLinear);
hold off;
% Write the result in the file
fid_w_2 = fopen('Parameter_Cal_Result/GainSubcarrierLinear.txt','w');
[row,col] = size(GainSubcarrierLinear);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_2,'%g\n',GainSubcarrierLinear(i,j));
		else
			fprintf(fid_w_2,'%g\t',GainSubcarrierLinear(i,j));
		end
	end
end
fclose(fid_w_2);

% Third part: LED effect, path loss and PD effect
% below data is from data sheet
LambdaSample = xlsread('Measurement_Result\Lambda_RelativeFlux.xlsx','A3:A26'); % Wavelength of light, Unit: nm
NormalizedFluxDensitySample = xlsread('Measurement_Result\Lambda_RelativeFlux.xlsx','B3:B26');

Lambda = linspace(380,780,401);
[population, ~] = fit(LambdaSample,NormalizedFluxDensitySample,'pchipinterp');
NormalizedFluxDensity = zeros(1,401);
NormalizedFluxDensity = population(Lambda);

VLambda = 1.019*exp(-285.4*(Lambda./1e3 - 0.559).^2);
NormlizedFlux = 683 * trapz(Lambda, NormalizedFluxDensity.*(VLambda)');
figure(2)
plot(Lambda,NormalizedFluxDensity);

% PD response curve calculation,it will be used calculate PD current
PdLambdaSample = xlsread('Measurement_Result\Lambda_Responsivity.xlsx','A3:A11');
PdResponsivitySample = xlsread('Measurement_Result\Lambda_Responsivity.xlsx','B3:B11');
PdLambda = linspace(380,780,401);
[population,gof] = fit(PdLambdaSample,PdResponsivitySample,'pchipinterp');
PdResponsivity = population(PdLambda);
figure(3)
plot(Lambda,PdResponsivity);

% Pass loss calculation, it will be used calculate PD current
ORDER_LAMBERTIAN_EMISSION = -log(2)/log(cos(HALF_POWER_ANGLE));
PATH_LOSS = PD_AREA*(ORDER_LAMBERTIAN_EMISSION+1)/(2*pi*DISTANCE^2)*((cos(TRANSMITTER_ANGLE))^ORDER_LAMBERTIAN_EMISSION)*cos(RECEIVER_ANGLE)*TRANSMISSION_OP_BPFILTER*GAIN_NON_IMG_CONCENTRATOR;
PdCurrentRefer = PATH_LOSS*REFER_FLUX/NormlizedFlux*trapz(Lambda,NormalizedFluxDensity.*PdResponsivity);% To make this value have physical meaning, the PD current when the LED forward current is 350mA is calculated. 
save 'Parameter_Cal_Result\PdCurrentRefer.txt' -ascii PdCurrentRefer;
% Forth part: PD noise model
% adding the noise of photo detector and transimpedence amplifier
BANDWIDTH = SAMPLING_RATE/2; % Bandwidth is half of sampling rate
PD_TYP_RESPONSIVITY = 0.45; % typical responsivity of photo detector, unit A/W
PD_NOISE_DATASHEET = PD_TYP_RESPONSIVITY*PD_NPE*sqrt(BANDWIDTH);
save 'Parameter_Cal_Result\PD_NOISE_DATASHEET.txt' -ascii PD_NOISE_DATASHEET;
TIA_NOISE_DATASHEET = TIA_NPE*sqrt(BANDWIDTH);
save 'Parameter_Cal_Result\TIA_NOISE_DATASHEET.txt' -ascii TIA_NOISE_DATASHEET;
% save 'Parameter_Cal_Result\BANDWIDTH.txt' -ascii BANDWIDTH;
% % Dark current noise
% NOISE_VARIANCE_DARK_CURRENT = 2*q*DARK_CURRENT*BANDWIDTH;
% save 'Parameter_Cal_Result\NOISE_VARIANCE_DARK_CURRENT.txt' -ascii NOISE_VARIANCE_DARK_CURRENT;
% % Thermal noise
% NOISE_VARIANCE_THERMAL = 4*k*TEMPERATURE/R_SH_DIODE*BANDWIDTH;
% save 'Parameter_Cal_Result\NOISE_VARIANCE_THERMAL.txt' -ascii NOISE_VARIANCE_THERMAL;

