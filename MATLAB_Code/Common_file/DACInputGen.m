function[FLOAT_OUTPUT,BINARY_OUTPUT,BINARY_OUTPUT_COMPLE] = DACInputGen(INPUT, WORD_LENGTH, FRACTION_LENGTH)
% the BINARY_OUTPUT includes three parts: sign part (1 bit)
% + integer part(WORD_LENGTH-1-FTRACTION_LENGTH bit)
% + fraction part(FRACTION_LENGTH bit)

%%script form for testing
%INPUT = [-1.725;1.75];
%WORD_LENGTH = 8;
%FRACTION_LENGTH = 6;

[row, col] = size(INPUT);
if (col ~= 1 && row ==1)
	Num = col;
elseif (row ~= 1 && col ==1)
	Num = row;
elseif (row ==1 && col ==1)
	Num = 1;
else
	error('matrix input is not allowed')
end

% quantize
QuantizerInst = quantizer('fixed','Nearest','saturate',[WORD_LENGTH,FRACTION_LENGTH]);
FLOAT_OUTPUT = quantize(QuantizerInst,INPUT);
FIX_INPUT = FLOAT_OUTPUT;
% sign bit generation
SignPart = sign(FIX_INPUT);
SignPart(find(SignPart == 1)) = 0;
SignPart(find(SignPart == -1)) = 1;
FIX_INPUT = abs(FIX_INPUT);
if (col ~= 1 && row ==1)
	SignPart = SignPart';
end
% integer part binary generation
IntergerPart = dec2bin(FIX_INPUT); %now the data is char array
IntergerBinary = double(IntergerPart)-48; % char -> ASCII -> double
if size(IntergerBinary,2) < (WORD_LENGTH - 1- FRACTION_LENGTH)
	IntergerBinary = [zeros(Num,WORD_LENGTH - 1- FRACTION_LENGTH-size(IntergerBinary,2)),IntergerBinary];
elseif size(IntergerBinary,2) > (WORD_LENGTH - 1- FRACTION_LENGTH)
	IntergerBinary = IntergerBinary(:,((1+size(IntergerBinary,2)-(WORD_LENGTH-1-FRACTION_LENGTH)):end));
end
% fraction part binary generation
FractionPart = mod(FIX_INPUT,1);
count = 0;
if(FRACTION_LENGTH == 0)
	flag = 0;
	FractionBinary = [];
else
	flag = 1;
	FractionBinary = zeros(Num,FRACTION_LENGTH);
end
while(flag)
	count = count + 1;
	if (count > FRACTION_LENGTH)
		flag = 0;
		break;
	end
	FractionPart = FractionPart .* 2;
	IndexLessThanOne = find(FractionPart < 1);
	FractionBinary(IndexLessThanOne,count) = 0;
	IndexGreatThanOne = find(FractionPart >= 1);
	FractionPart(IndexGreatThanOne) = FractionPart(IndexGreatThanOne) - 1;
	FractionBinary(IndexGreatThanOne,count) = 1;
end

BINARY_OUTPUT = [SignPart, IntergerBinary,FractionBinary];


% Using for to calculate 2 complementary code one by one
BINARY_OUTPUT_COMPLE = [SignPart,zeros(Num,WORD_LENGTH-1)];
for num = 1:Num
	OriCodeTemp = BINARY_OUTPUT(num,2:end);
	if(~SignPart(num))
		BINARY_OUTPUT_COMPLE(num,2:end) = OriCodeTemp;
	else
		l = length(OriCodeTemp);
		for i = 1:l % except sign bit, every bit nagates
			if OriCodeTemp(i) == 1
				OriCodeTemp(i) = 0;
			else
				OriCodeTemp(i) = 1;
			end
		end
		% add one on last bit
		temp_bit = l;
		while(temp_bit ~= 0)
			if (OriCodeTemp(temp_bit) == 0) 
				OriCodeTemp(temp_bit) = 1;
				temp_bit = 0;
			else
				OriCodeTemp(temp_bit) = 0;
				temp_bit = temp_bit - 1;
			end
		end
	end
		BINARY_OUTPUT_COMPLE(num,2:end) = OriCodeTemp;
end
	

