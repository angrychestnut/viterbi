close all;clear all; clc;
N = 1e5;   % simulate N bits  each  transmission (one  block)
maxNumErrs = 100; % get at least  100  bit  errors (more is  better)
maxNum = 1e7; % OR stop if  maxNum  bits  have  been  simulated
EbN0 =  -1:0.5:12; % power  efficiency  range

% These are the constellations to choose from. All gray-coded.
BPSKconst = [-1,1];
QPSKconst = [1, 1j, -1j, -1];
AMPMconst = [-1j+1, -3+3j, 3j+1, -3-1j, -3j+3, 1j-1, 3+1j, -1-3j];

% Here you choose the constellation to use:
currentCons = BPSKconst;
bitsPerSymbol = log2(length(currentCons)); 

Es = mean(abs(currentCons).^2); %Average symbol energy

N0 = Es./(10.^(EbN0/10)*bitsPerSymbol); %Power spectrum density


% BER_theory = qfunc(sqrt(Es./N0));
% 
% semilogy(EbN0, BER_theory);
pause(1)


% Other  Options  ...
% ======================================================================= %
% Simulation  Chain

% ======================================================================= %

BER = zeros(1, length(EbN0 )); % pre -allocate a vector  for the  BER  results
for i = 1: length(EbN0) % use  parfor ('help  parfor ') to  parallelize
totErr = 0;   % running  number  of  errors  observed
num = 0; % running  number  of bits  processed





while (( totErr  < maxNumErrs) && (num < maxNum ))
% ===================================================================== %
% Begin  processing  one  block  of  information
% ===================================================================== %
% [SRC] generate N information  bits



variance = N0(i)/2; % Noise variance


%Generate N number of random bits. Make sure N is a multiple of
%bitsPerSymbol WARNING: This might be changed for the coded case.
N = N-mod(N,bitsPerSymbol);
bits = randi([0,1],1,N);

% [ENC] convolutional  encoder

% todo




% [MOD] symbol  mapper

messages = buffer(bits,bitsPerSymbol); %Change bit array to matrix of messages
messageIndexes = bi2de(messages'); %Convert bit message to message index
symbols = currentCons(1+messageIndexes); %Map all messages to symbols.

% [CHA] add  Gaussian  noise

%noise = sqrt(variance)*(randn(1,length(symbols))+1j*randn(1,length(symbols)));
noise = sqrt(variance)*(randn(1,length(symbols)));
y = noise+symbols;

% [HR] Hard  Receiver

% Allocate arrays for storing symbols and corresponding messageindexes
estimated_symbols = zeros(1,length(y));
estimated_messageindex = zeros(1,length(y));

for mi = 1:length(y) %For all the received symbols pick the symbol closest in the constellation
    
  [M, I] = min(abs(y(mi)-currentCons));
  estimated_symbols(mi) = currentCons(I);
  estimated_messageindex(mi) = I;
  
end

%Map message indexes to bits. Lowest message index is 1 which is mapped to
%[0 0]
estimated_messages = de2bi(estimated_messageindex-1);

estimated_bits = reshape(estimated_messages',N,1)';

% [SR] Soft  Receiver


% ===================================================================== %
% End  processing  one  block  of  information
% ===================================================================== %
BitErrs = sum(estimated_bits ~= bits); % count  the bit  errors  and  evaluate  the bit  error  rate
totErr = totErr + BitErrs;
num = num + N;
disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
num2str(num) '/' num2str(maxNum) ' bits. Projected  error  rate = '...
num2str(totErr/num , '%10.1e') '. +++']);
end
BER(i) = totErr/num;
end

hold on
semilogy(EbN0, BER)
hold off
% ======================================================================= %
% End
% ======================================================================= %