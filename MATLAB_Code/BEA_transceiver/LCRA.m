function [en_bar,bn_bar,bit_sum,energy_sum,N_used,BER]=LCRA(gn,SubCarrierNum,SER)
% gn: SNR of nth subcarrier with unit input energy, linear input rather than db input
% SubCarrierNum: the number of subcarrier could loading data
% SER: target symbol error rate, used to calculate the gap value. 
% en_bar: energy scaling factor allocated to each subcarrier.
% bn_bar: bit number allocated to each subcarrier.
% N_used: used subcarrier.
% BER: theory bit error rate.
Ebar = 1; % normalized energy

% parameter checked:
if length(gn) ~= SubCarrierNum
	error('the length of gn is not consistent with SubCarrierNum.');
end

% from SER to gap calculation
snrGap = 1/3*((qfuncinv(SER/4))^2);
snrGapDb = 10*log10(snrGap);

% initialization
en=zeros(1,SubCarrierNum);
bn=zeros(1,SubCarrierNum);

% Total energy so far
E_so_far=0;
% decision table for QAM
decision_table(1:SubCarrierNum) = 2 * snrGap ./ gn(1:SubCarrierNum);
while(1)
	[y,index] = min(decision_table);
	E_so_far = E_so_far + y;
	if E_so_far > Ebar * SubCarrierNum
		break;
	else
		en(index) = en(index) + y;
		bn(index) = bn(index) + 1; 
		decision_table(index) = 2 * decision_table(index);
	end
end
en_bar = en;
bn_bar = bn;
bit_sum = sum(bn_bar);
energy_sum = sum(en_bar);
N_used = length(find(bn_bar~=0));
BER = SER/(bit_sum/N_used);

fig_1 = figure(1);
scatter(1:SubCarrierNum,bn_bar,'LineWidth',2);
xlim([0,60]);
ylim([0,6]);
for subcarrierindex = 1:SubCarrierNum
    line([subcarrierindex,subcarrierindex],[bn(subcarrierindex),0]);
end
title('Number of Bit on Each Subcarrier');
xlabel('Index of Subcarrier');
ylabel('Number of Bit');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(fig_1,"Generated_pic/Bit_Allocation",'jpg');

fig_2 = figure(2);
scatter(1:SubCarrierNum,en_bar,'LineWidth',2);
xlim([0,60]);
ylim([0,3]);
for subcarrierindex = 1:SubCarrierNum
    line([subcarrierindex,subcarrierindex],[en_bar(subcarrierindex),0]);
end
title('Energy Allocation Scaling Factor on Each Subcarrier');
xlabel('Index of Subcarrier');
ylabel('Energy Scaling Factor');
set(gca, 'fontsize', 16);
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'LineWidth', 1.5);
saveas(fig_2,"Generated_pic/Energy_Allocation",'jpg');