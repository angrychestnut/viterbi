close all;
clear;
clc;
N = 1e5;   % simulate N bits  each  transmission (one  block)
maxNumErrs = 100; % get at least  100  bit  errors (more is  better)
maxNum = 1e7; % OR stop if  maxNum  bits  have  been  simulated
EbN0 =  -1:0.5:12; % power  efficiency  range -2:0.5:2.5

constraint_length = 5;
trellis = poly2trellis(constraint_length,[23 33]);
spect = distspec(trellis,20);
BER_UB = bercoding(EbN0,'conv','soft',1/2,spect); % BER bound


% These are the constellations to choose from. All gray-coded.
BPSKconst = [-1,1];
QPSKconst = [1+1j, 1-1j, -1+1j, -1-1j];
AMPMconst = [-1j+1, -3+3j, 3j+1, -3-1j, -3j+3, 1j-1, 3+1j, -1-3j];

% Here you choose the constellation to use:
currentCons = BPSKconst;
bitsPerSymbol = log2(length(currentCons)); 

Es = mean(abs(currentCons).^2); %Average symbol energy
coderate = 1/2;
N0 = Es./(bitsPerSymbol * coderate * 10.^(EbN0/10)); %Power spectrum density




% Other  Options  ...
% ======================================================================= %
% Simulation  Chain

% ======================================================================= %

BER = zeros(1, length(EbN0 )); % pre -allocate a vector  for the  BER  results
for i = 1: length(EbN0) % use  parfor ('help  parfor ') to  parallelize
totErr = 0;   % running  number  of  errors  observed
num = 0; % running  number  of bits  processed





while (( totErr  < maxNumErrs) && (num < maxNum ) )
% ===================================================================== %
% Begin  processing  one  block  of  information
% ===================================================================== %
% [SRC] generate N information  bits



variance = N0(i)/2; % Noise variance
%variance = N0(i);

%Generate N number of random bits. Make sure N is a multiple of
%bitsPerSymbol WARNING: This might be changed for the coded case.
N = N-mod(N,bitsPerSymbol);
bits = randi([0,1],1,N);

% [ENC] convolutional  encoder

% todo

encodedbits=convenc(bits,trellis);


% [MOD] symbol  mapper

messages = buffer(encodedbits,bitsPerSymbol); %Change bit array to matrix of messages
messageIndexes = bi2de(messages','left-msb'); %Convert bit message to message index
symbols = currentCons(1+messageIndexes); %Map all messages to symbols.

% [CHA] add  Gaussian  noise

%noise = sqrt(variance)*(randn(1,length(symbols))+1j*randn(1,length(symbols)));
noise = sqrt(variance)*(randn(1,length(symbols)));
y = noise+symbols;


% % [HR] Hard  Receiver
% 
% % Allocate arrays for storing symbols and corresponding messageindexes
% estimated_symbols = zeros(1,length(y));
% estimated_messageindex = zeros(1,length(y));
% 
% for mi = 1:length(y) %For all the received symbols pick the symbol closest in the constellation
%     
%   [M, I] = min(abs(y(mi)-currentCons));
%   estimated_symbols(mi) = currentCons(I);
%   estimated_messageindex(mi) = I;
%   
% end
% 
% %Map message indexes to bits. Lowest message index is 1 which is mapped to
% %[0 0]
% estimated_messages = de2bi(estimated_messageindex-1,'left-msb');
% 
% estimated_bits = reshape(estimated_messages',length(encodedbits),1)';
% 
% decodedbits_hard=vitdec(estimated_bits,trellis,50,'trunc','hard');

% [SR] Soft  Receiver
% rey=real(y);
% imy=imag(y);
% Y = upsample(rey,2,0)+upsample(imy,2,1);

decodedbits_soft=vitdec(y,trellis,5*constraint_length,'trunc','unquant');

% ===================================================================== %
% End  processing  one  block  of  information
% ===================================================================== %
%bits=[bits,0 0 0 0]
BitErrs = sum(decodedbits_soft ~= bits); % count  the bit  errors  and  evaluate  the bit  error  rate
totErr = totErr + BitErrs;
num = num + N;
disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
num2str(num) '/' num2str(maxNum) ' bits. Projected  error  rate = '...
num2str(totErr/num , '%10.1e') '. +++']);
end
BER(i) = totErr/num;
end

% % ======================================================================= %
% % End
% ===================================================================== %