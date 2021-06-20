function [OFDMSymbolWithoutPaddingAcc,DataRate,BitErrorRate] = Receiving_Realtime_Bit_Energy_Allocation(N,RArray,RSum,SArray,Ncp,RepeatTime,FrameNum,TotalFrameNum,SubCarrierNum,LowPaddingNum,HighPaddingNum,ThresholdValue,DigitalReceivedDownSampling,DataSerialOri)

% Starting data acquisition and BER calculation
%% Synchronization
i = 1; % index of input data
% Some registers used in Frame Syn
L1 = 32;
L2 = 70;
L3 = N;
L4 = N*(RepeatTime-1) + Ncp;
SynchronizationOffset = 0;% this parameter is the offset position of the beginning of the sequence. 
                                  % if this value is 0, Then the sequence
                                  % beginning posision is the first sample after the repeaded time peak
% = 40;
% Parameter Check
if SynchronizationOffset > Ncp
    error('The SynchronizationOffset is invalid, Please adjust the parameters.');
end

% globle registers and varibles
OFDMSymbolWithoutPaddingAcc = [];
OFDMSymbolNumArray = [];
ErrorBitNumArray = [];
ErrorBitRatioArry = [];
while(1)
     % go through the channel model
     waveformASC = DigitalReceivedDownSampling;
     LengthReceivedSeq = length(waveformASC);
     DataIn = 0;
     FirstOrderShiftRam = zeros(1,L1);
     SecondOrderShiftRam = zeros(1,L2);
     CorrelatedSTSRam = zeros(1,L1);
     CorrelatedOriRam = zeros(1,L1);
     CorrelatedSum = 0;
     CorrelatedOriSum = 0;
     
     % Some registers used in Symbol Syn
     DataShiftRam = zeros(1,L3);
     SignDataShiftRam = zeros(1,L3);
     LongTSTRam = zeros(1,L3);
     SignLongTSTRam = zeros(1,L3);
     DataShiftForChannelEstimationRam = zeros(1,L4);
     SymbolDetected = 0;
     TimeCounter = 0;%detect the range between the peak.
     BitSumPro = 0;
     Synchronized = 0;
     
     % Some register used in removing CP module

     CpDataIn = 0;
     FrameDetcted = 0; % every time ratio is larger than threshold, this value will be added to 1.
     SymbolCounter = 0; % counting the symbol after the synchronization
     FrameCounter = 0;% counting the OFDM frame after the synchronization
     SymbolSynEnable = 0;
     
     % Some registers used in Channel estimation
     OneLSTafterFFT = zeros(1,N);
     AccumulateLTSafterFFT = zeros(1,N);
     
     % Some registers used in demodulation
     OFDMSymbolWithoutPadding = [];
     plot_enable = 0;
     i = 1;
     RatioOutputArray = [];
     MproArray = [];
     BitSumProAarry = [];
     BitSumAarry = [];


     while(1)
         % the interpolation will be added later
         DataInMSB = DataIn;
         if i <= LengthReceivedSeq
             DataIn = waveformASC(i);
         elseif i > LengthReceivedSeq && i <= LengthReceivedSeq + 1 + L1 + L2 + L3 +L4
             DataIn = 0; %When the data transmission is done, the 0 or noise is transmitted in the system
             %i = i-1;
         else
             disp("Symbol synchronization cannot be found in the whole seq ...");
             break;
         end
         i = i + 1;
         FirstOrderShiftRamMSB = FirstOrderShiftRam(L1);
         FirstOrderShiftRam = [DataInMSB,FirstOrderShiftRam(1:L1-1)];
         SecondOrderShiftRamMSB = SecondOrderShiftRam(L2);
         SecondOrderShiftRam = [FirstOrderShiftRamMSB,SecondOrderShiftRam(1:L2-1)];
         CorrelatedSTS = sign(DataInMSB) * sign(FirstOrderShiftRamMSB);
         CorrelatedSTSRamMSB = CorrelatedSTSRam(L1);
         CorrelatedSTSRam = [CorrelatedSTS,CorrelatedSTSRam(1:L1-1)];
         CorrelatedSum = CorrelatedSum + CorrelatedSTS - CorrelatedSTSRamMSB;
         CorrelatedOri = sign(DataInMSB)*sign(DataInMSB);
         CorrelatedOriRamMSB = CorrelatedOriRam(L1);
         CorrelatedOriRam = [CorrelatedOri,CorrelatedOriRam(1:L1-1)];
         CorrelatedOriSum = CorrelatedOriSum + CorrelatedOri - CorrelatedOriRamMSB;
         % This part is used to plot the frame synchronization result
         if CorrelatedOriSum == 0
             RatioOutput = 0;
         else
             RatioOutput = abs(CorrelatedSum)/CorrelatedOriSum;
         end
         RatioOutputArray(i) = RatioOutput;
         % Below is Frame detection module
         if CorrelatedSum > CorrelatedOriSum/2 %threshold is 0.5
             FrameDetcted = FrameDetcted  + 1;
         end
         if FrameDetcted == 32
             disp('Frame is detected ...');
             SymbolSynEnable = 1;
         end
         % Symbol Synchronization
         if(SymbolSynEnable)
             DataShiftRamMSB = DataShiftRam(L3);
             DataShiftRam = [SecondOrderShiftRamMSB,DataShiftRam(1:L3-1)];
             DataShiftForChannelEstimationRamMSB = DataShiftForChannelEstimationRam(L4);
             DataShiftForChannelEstimationRam = [DataShiftRamMSB,DataShiftForChannelEstimationRam(1:L4-1)];
             [~,LongTSTRam,~]= TrainingSeqGenParkMethod(N,Ncp,1,[],0);
             LongTSTRam = fliplr(LongTSTRam);
             SignDataShiftRam(find(DataShiftRam>=0)) = 0;
             SignDataShiftRam(find(DataShiftRam<0)) = 1;
             SignLongTSTRam(find(LongTSTRam>=0)) = 0;
             SignLongTSTRam(find(LongTSTRam<0)) = 1;
             BitSum = sum(~xor(SignDataShiftRam,SignLongTSTRam));
             BitSumAarry(i) = BitSum;
             BitSumLastPro = BitSumPro;
             BitSumPro = 2*BitSum-(N);
             BitSumProAarry(i) = BitSumPro;
             if BitSumPro > 60 % first threshold
                 Mpro = BitSumPro - BitSumLastPro;
             else
                 Mpro = BitSumPro;
             end
             MproArray(i) = Mpro;
             if Mpro > ThresholdValue %Threshold is 120
                 if TimeCounter == 0 
                     SymbolDetected = SymbolDetected + 1;
                 elseif TimeCounter > N/2 % time range is larger than 
                     SymbolDetected = SymbolDetected + 1;
                     TimeCounter = 0;
                 end
             end
             if SymbolDetected ~= 0 && SymbolDetected ~= RepeatTime
                 TimeCounter = TimeCounter + 1;
             end
             if SymbolDetected == RepeatTime
                 Synchronized = 1;
                 SymbolDetected = 0;
                 disp('The Symbol is synchronized ...');
                 PayLoadLength = LengthReceivedSeq + L1 + L2 - i + 2;
             end
             if Synchronized
                 if(FrameCounter <= RepeatTime-1)
                     SymbolCounter = SymbolCounter + 1;
                     CpDataIn(SymbolCounter) = DataShiftForChannelEstimationRam((RepeatTime-1)*N+SynchronizationOffset);
                     if SymbolCounter == N
                         FrameCounter = FrameCounter + 1;
                         SymbolCounter = 0;
                         OneLTSafterFFT = fft(CpDataIn); %FFT transformation
                         AccumulateLTSafterFFT = AccumulateLTSafterFFT + OneLTSafterFFT;
                         if FrameCounter == RepeatTime
                             AverageLTSafterFFT = AccumulateLTSafterFFT./RepeatTime;
                             ReceivedLTSFreqDom = AverageLTSafterFFT(2:N/2);
                             [~,~,LTSFreqDom] = TrainingSeqGenParkMethod(N,Ncp,RepeatTime,[],0);
                             H = ReceivedLTSFreqDom./LTSFreqDom(2:N/2); % Channel information 
                         end
                     end
                 else
                     if PayLoadLength-(FrameCounter-RepeatTime+1)*(N+Ncp) >= 0 && FrameCounter - RepeatTime < FrameNum
                         SymbolCounter = SymbolCounter + 1;
                         CpDataIn(SymbolCounter) = DataShiftForChannelEstimationRam((RepeatTime-1)*N+SynchronizationOffset);
                     else
                         disp('Received Done, Demodulation processing begin ...');
                         plot_enable = 1;
                         if FrameCounter-RepeatTime == 0
                             disp("the sequence is too short to include any effective payload.")
                             break;
                         end
                         OFDMSymbolWithoutPaddingAcc = [OFDMSymbolWithoutPaddingAcc,OFDMSymbolWithoutPadding]; %OFDM symbols smapled in several times
                         OFDMSymbolNumArray = [OFDMSymbolNumArray,FrameCounter-RepeatTime]; % the OFDM number getted each sampling time
                         % demodulation and calculate the bit error rate
                         DataSerialReshape = [];
                         for nSubcarrier = 1:SubCarrierNum
                             if RArray(nSubcarrier) == 0
                                 DataSerialReshape = DataSerialReshape;
                             else
                                 OFDMSymbolWithoutPaddingOneSubcarrier = OFDMSymbolWithoutPadding(nSubcarrier,:);
                                 DataSerialOneSubcarrier =  qamdemod(OFDMSymbolWithoutPaddingOneSubcarrier ./ (sqrt(SArray(nSubcarrier)/2)),2^RArray(nSubcarrier),'gray','OutputType','bit','UnitAveragePower',true);
                                 DataSerialReshape = [DataSerialReshape;DataSerialOneSubcarrier];
                             end
                         end
                         DataSerial = reshape(DataSerialReshape,1,RSum*(FrameCounter-RepeatTime));
                         [ErrorNum, ErrorRate] = biterr(DataSerialOri(:,1:RSum*(FrameCounter-RepeatTime)),DataSerial);
                         ErrorBitNumArray = [ErrorBitNumArray,ErrorNum];
                         ErrorBitRatioArry = [ErrorBitRatioArry,ErrorRate];
                         disp('Bit Error Rate in this time is');
                         disp(ErrorRate);
                         break;
                     end
                     if SymbolCounter == N+Ncp
                         FrameCounter = FrameCounter + 1;
                         SymbolCounter = 0;
                         % removing CP
                         OneOFDMSymbolWithoutCP = CpDataIn(Ncp+1:end);
                         % FFT transformation
                         OneOFDMSymbolFreqDom = fft(OneOFDMSymbolWithoutCP);
                         OneOFDMSymbolWithoutHermitian = OneOFDMSymbolFreqDom(2:N/2);
                         OneOFDMSymbolAfterChannelEstimation = OneOFDMSymbolWithoutHermitian./H;
                         % Pilot extraction, no channel estimation
                         OneOFDMSymbolWithoutPadding = OneOFDMSymbolAfterChannelEstimation(LowPaddingNum+1:LowPaddingNum+SubCarrierNum);
                         OFDMSymbolWithoutPadding = [OFDMSymbolWithoutPadding,OneOFDMSymbolWithoutPadding.'];
                     end
                 end
             end 
         end        
     end
     disp('Drawing the Symbol Synchronization ...');
     figure(3)
     %plot(1:i,RatioOutputArray);
     plot(1:1000,RatioOutputArray(1:1000));
     %title('Frame Detection');
     xlabel('Sampling Point Index');
     ylabel('Ratio of the Delay Correlation and Energy');
     set(gca, 'fontsize', 16);
     set(gca, 'XMinorTick', 'on');
     set(gca, 'YMinorTick', 'on');
     set(gca, 'XGrid', 'on');
     set(gca, 'YGrid', 'on');
     set(gca, 'LineWidth', 1.5);
     if(SymbolSynEnable == 1) %when the frame is detected, the symbol detection is started, or it doesn's work
         figure(4)
         %plot(1:i,MproArray,'r');
        % plot(1:5000,MproArray(1:5000));
         %title('Symbol Synchronization');
         xlabel('Index of Subcarrier');
         ylabel('Sum of the Correlation');
         set(gca, 'fontsize', 16);
         set(gca, 'XMinorTick', 'on');
         set(gca, 'YMinorTick', 'on');
         set(gca, 'XGrid', 'on');
         set(gca, 'YGrid', 'on');
         set(gca, 'LineWidth', 1.5);
     end
	% Drawing the receiveing constellation
    if(sum(OFDMSymbolNumArray)~=0)
        for fig_index = 1:SubCarrierNum
            fig = figure(fig_index+4);
            scatter(real(OFDMSymbolWithoutPaddingAcc(fig_index,:)),imag(OFDMSymbolWithoutPaddingAcc(fig_index,:)),'.');
           xlim([min(-1,min(real(OFDMSymbolWithoutPaddingAcc(fig_index,:)))),max(1,max(real(OFDMSymbolWithoutPaddingAcc(fig_index,:))))]);
           ylim([min(-1,min(imag(OFDMSymbolWithoutPaddingAcc(fig_index,:)))),max(1,max(imag(OFDMSymbolWithoutPaddingAcc(fig_index,:))))]);
            set(fig,'visible','off');
            saveas(fig,fullfile('.\Generated_pic\Received_Constellation\',[num2str(fig_index),'th','Subcarrier.jpg']));
        end
    end
	if (sum(OFDMSymbolNumArray) >= TotalFrameNum)
	    disp('Receiving Done, Begin to calculate ....')
        % calculate the final bit error rate
        TotalBitReceived = sum(OFDMSymbolNumArray)*RSum;
        TotalErrorBitReceived = sum(ErrorBitNumArray);
        BitErrorRate = TotalErrorBitReceived/TotalBitReceived;
        disp('Bit Error Rate is: ')
        disp(BitErrorRate);
        % calculate the data rate
        SAMPLING_RATE = load('Parameter_Cal_Result\SAMPLING_RATE.txt');
        DataRate = (RSum)/(1/SAMPLING_RATE*(N+Ncp));
        disp('Data Rate is: ');
        disp(DataRate/1e6);
        disp('Mbps');
	    break;
	end
end

