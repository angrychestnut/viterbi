function decodedBits = viterbi(inputbits)

% The viterbi decoding algorithm for X1 generator polynomial as given in project PM.

pathsToState = [0,2; %From here we get to state 0. We get to state 0 from either state 0 or 2
    0,2; % To state 1 etc...
    1,3;
    1,3];

outputAtTransistion = [0,0; %The bit output corresponding to the transition from pathsToState(x) to x
    1,1; %There are two paths to each state. Therefore two bit output alternative. First and second in the same way as pathsToState.
    
    1,1;
    0,0;
    
    0,1;
    1,0;
    
    1,0;
    0,1];


%To run the algorithm in linear time we cannot copy any arrays at
%each iteration. (MATLAB do not support passing by reference a.k.a pointers). We Therefore simulate a single linked list
%by a value array and a index array which gives index to previous entry in the array.

previousbitIndex = zeros(1,length(inputbits+5)*2); %The array index of the previous bit in bit sequence
bitvector = zeros(1,length(inputbits+5)*2); %Here are the bits corresponding to indexes in PreviousStepArrayIndex
lastentry = 1; %Pointer to keep track of where to put new bits into the bitvector

lastBitIndexForState = [0 0 0 0]; %Indexes to keep track of the last bit for each state
tmplastBitIndexForState = [0 0 0 0]; %Temporary storing location during calculations of new last bits.

weightToState = [0, 100, 100, 100]; %The hamming weight to get to the corresponding state.There is no chance that we started in a state other than 0. I therefore give these higher weight
weightToStateNext = [0,0,0,0]; %Temporary array for weight during calculation


for i = 1:2:length(inputbits) %For all bits in sequence
    currentbits = inputbits(i:i+1); %Looked at two and two
    
    for currentState = 0:3 %We calculate the shortest path to all 4 states.
        
        bitsIfFirst = outputAtTransistion(2*(currentState+1)-1,:); %There are two paths to every state. These are the output bits for the two alternative paths
        bitsIfSecond = outputAtTransistion(2*(currentState+1),:);
        
        %The weight to the state given that we come at the first
        %resp. second path.
        weightFirst = weightToState(pathsToState(currentState+1,1)+1) + sum(currentbits ~= bitsIfFirst);
        weightSecond = weightToState(pathsToState(currentState+1,2)+1) + sum(currentbits ~= bitsIfSecond);
        
        
        if weightFirst < weightSecond %If there is a total shorter hamming distance for the first alternative this is more likely.
            weightToStateNext(currentState+1) = weightFirst;
            shortestAlternative = 1;
        else
            weightToStateNext(currentState+1) = weightSecond;
            shortestAlternative = 2;
        end
        
        
        bitvector(lastentry) = mod(currentState,2); %Save the input bit corresponding to the path in the 'linked list'- To get to even state nr a zero is needed, otherwise a 1
        previousbitIndex(lastentry) = lastBitIndexForState(pathsToState(currentState+1,shortestAlternative)+1); %Link to the previous bits in sequence
        
        tmplastBitIndexForState(currentState+1) = lastentry; %Save a pointer to the most likely bit
        
        lastentry = lastentry + 1; %Move pointer to empty location.
        
    end
    
    
    lastBitIndexForState = tmplastBitIndexForState;
    weightToState = weightToStateNext;
    
end

N = length(inputbits)/2-2; %Number of outputbits

%Follow path to get bit sequence reverse
decodedBits = zeros(1,N+2);

prev = lastBitIndexForState(1); %Pointer to the previous bit in sequence

for j = 0:N-1+2 %We need to find two extra bits. Because the machine generates 2 more bits/dimension then needed.
    decodedBits(N+2-j) = bitvector(prev);
    prev = previousbitIndex(prev); %Update pointer
end

decodedBits = decodedBits(1:end-2); %Remove last bits
