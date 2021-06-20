function[PAPR_VALUE] = PAPR_cal(OFDMSequenceTimeDomain)
% This function is used to calculated peak to average power ratio (PAPR) of ofdm signal in time domain

Peak = max(abs(OFDMSequenceTimeDomain).^2); % maximum power of output signal
Ave = mean(abs(OFDMSequenceTimeDomain).^2); % mean power of output signal
PAPR_VALUE = 10*log10(Peak/Ave);