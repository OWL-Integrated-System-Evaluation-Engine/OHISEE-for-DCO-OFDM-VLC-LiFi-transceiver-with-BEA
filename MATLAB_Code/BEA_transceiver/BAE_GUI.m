function [SArray,RArray] = BAE_GUI();
% Add files path
addpath('.\Common_file');
addpath('.\Channel_model');
addpath('.\Channel_model\Parameter_Cal_Result');
% Calcuate RArray and SArray
RArray = [];
SArray = [];
% load gn value from SNR measurement result
SER = load('Parameter_Cal_Result\target_ser.txt');
SubCarrierNum = load('Parameter_Cal_Result\SubCarrierNum.txt'); % every OFDM frame has 59 subcarriers
gn =[];
fid_r_0=fopen('Parameter_Cal_Result\SNRLinearResult.txt','r'); 
while 1     
    tline=fgetl(fid_r_0);     
    if ~ischar(tline),break;
    end     
    tline=str2num(tline);     
    gn = [gn,tline];
end 
fclose(fid_r_0);
%% Result of LCRA function call mismatch
%[en_bar,bn_bar,bit_sum,energy_sum,N_used,BER]=LCRA(gn,SubCarrierNum,SER)
[SArray,RArray,RSum,ESum,SubCarrierNumUsed,BER]=LCRA(gn,SubCarrierNum,SER);
% write allocation result in the files
save '.\Channel_model\Parameter_Cal_Result\RSum.txt' -ascii RSum;
save '.\Channel_model\Parameter_Cal_Result\ESum.txt' -ascii ESum;
fid_w_0 = fopen('.\Channel_model\Parameter_Cal_Result\SArray.txt','w');
[row,col] = size(SArray);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_0,'%e\n',SArray(i,j)); %using scientific notation
		else
			fprintf(fid_w_0,'%e\t',SArray(i,j)); %using scientific notation 
		end
	end
end
fclose(fid_w_0);

fid_w_1 = fopen('.\Channel_model\Parameter_Cal_Result\RArray.txt','w');
[row,col] = size(RArray);
for i = 1:row
	for j = 1:col
		if(j == col)
			fprintf(fid_w_1,'%e\n',RArray(i,j)); %using scientific notation
		else
			fprintf(fid_w_1,'%e\t',RArray(i,j)); %using scientific notation 
		end
	end
end
fclose(fid_w_1);
