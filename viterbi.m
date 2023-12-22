clear all
warning off
%%% QPSK over AWGN
%%%%% SIGNAL CONSTELLATION %%%%%%%%%
symbolBook=[1 -1];
bitBook=[0;1];
nBitPerSym=size(bitBook,2);
M=length(symbolBook);
%%%%%%%%%%%%%%CHANNEL%%%%%%%%%%%%%%%%%%
channel = [0.74 -0.514 0.37 0.216 0.062];
N = length(channel)-1; %We have 4 past values.
%%%%%%%%%%%%%%%%%%%States%%%%%%%%%%%%%%%%%
number_of_states = M^(N);
array = linspace(0,number_of_states-1,number_of_states);
bit_array = de2bi(array,N,'left-msb'); 
states = zeros(16,4); %This array represents the states s_{k-1}
for i = 1:number_of_states
    for j = 1:N
        if bit_array(i,j) == 0
            states(i,j) = 1;
        else
            states(i,j) = -1;
        end
    end
end
states_i_0 = zeros(number_of_states,N); 
states_i_0(:,2:end) = states(:,1:3);
states_i_0(:,1) = 1;%These are the s_k states when input is 1

states_i_1 = zeros(number_of_states,N); 
states_i_1(:,2:end) = states(:,1:3);
states_i_1(:,1) = -1;%These are the s_k states when input is -1

diagram = zeros(16,4,2);
diagram(:,:,1) = states_i_0;
diagram(:,:,2) = states_i_1; %The 3-D matrix which stores all s_k's.
%%%%%%%%%%%%%%%%%%%%%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputs = zeros(number_of_states,M);
for i = 1:M
    for j = 1:number_of_states
        outputs(j,i) = dot([diagram(j,:,i) states(j,N)],channel);
    end
end
%%%%%%%%%%%%%%%%%%%%%MAP%%%%%%%%%%%%%%%%%%%%%%%%
zero_mapping = zeros(1,number_of_states); %The enumeration of the next state for input 1.
for i = 1:number_of_states
    buffer_zero = diagram(i,:,1);
    buffer_zero(buffer_zero>0) = 0;
    buffer_zero = buffer_zero * (-1);
    zero_mapping(1,i) = bin2dec(sprintf('%d',buffer_zero))+1;
end

one_mapping = zeros(1,number_of_states); %The enumeration of the next state for input -1.
for i = 1:number_of_states
    buffer_one = diagram(i,:,2);
    buffer_one(buffer_one>0) = 0;
    buffer_one = buffer_one * (-1);
    one_mapping(1,i) = bin2dec(sprintf('%d',buffer_one))+1;
end

%%%%%%%%%%%%%% MONTE CARLO PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
nSymPerFrame=3000;
nBitsPerFrame=nSymPerFrame*nBitPerSym;
max_nFrame=500;
snr_db=0:20;

nBitErrors=zeros(length(snr_db), 1);
nTransmittedFrames=zeros(length(snr_db), 1);
nErroneusFrames=zeros(length(snr_db), 1);
SYMBOLBOOK = repmat(transpose(symbolBook),1,nSymPerFrame);
%%%%%%%%%%%%%%%%%%%Viterbi Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    this_snr=snr_db(nEN);
    sigma_noise = 1/sqrt(10^(this_snr/10));
    while (nTransmittedFrames(nEN)<max_nFrame)
        nTransmittedFrames(nEN) = nTransmittedFrames(nEN) + 1;
        %%%%%%%%%% INFORMATION GENERATION %%%%%%%%%%
        trSymIndices=randi(M,[1,nSymPerFrame]);
        trSymVec=symbolBook(trSymIndices);
        trBitsMat=bitBook(trSymIndices,:)';
        %%%%%%%%%%%%%CHANNEL%%%%%%%%%%%%%%%%%%%%%
        intersymbol_interference = conv(trSymVec,channel); %The channel effect
        noise=1/sqrt(2)* [randn(1, length(intersymbol_interference))];
        recSigVec = intersymbol_interference + noise*sigma_noise;
        recSigVec = recSigVec(1,1:end-N);% Normalization to prevent overflowing.
        %%%%%%%%%%%%% Viterbi Algorithm %%%%%%%%%%%%
        path_Length = zeros(number_of_states,nSymPerFrame+1);
        path_Length(2:end,1) = inf; %Initialization
        for i = 1:nSymPerFrame
            for j = 1:number_of_states
                zero_index = zero_mapping(1,j);
                branch_Length = (recSigVec(1,i)-outputs(j,1))^(2);
                if path_Length(zero_index,i+1) == 0
                    path_Length(zero_index,i+1) = path_Length(j,i) + branch_Length;
                else
                    if path_Length(zero_index,i+1) > path_Length(j,i) + branch_Length
                        path_Length(zero_index,i+1) = path_Length(j,i) + branch_Length;
                    end
                end
                one_index = one_mapping(1,j);
                branch_Length_1 = (recSigVec(1,i)-outputs(j,2))^(2);
                if path_Length(one_index,i+1) == 0
                    path_Length(one_index,i+1) = path_Length(j,i) + branch_Length_1;
                else
                    if path_Length(one_index,i+1) > path_Length(j,i) + branch_Length_1
                        path_Length(one_index,i+1) = path_Length(j,i) + branch_Length_1;
                    end
                end
            end
        end
        [V,I] = min(path_Length);  %Decision
        I = I(2:end);
        I(I<=8) = 1;
        I(I>8) = -1;
        %%%%%%%%%%%%%%%%%%%DETECTION%%%%%%%%%%%%%%%%%%%%%%
        RECSIGVEC=repmat(I,length(symbolBook),1);
        distance_mat=abs(SYMBOLBOOK-RECSIGVEC);
        [~, det_sym_ind]=min(distance_mat,[],1);
        detected_bits=[bitBook(det_sym_ind, :)]';
        err = sum(sum(abs(trBitsMat-detected_bits)));
        nBitErrors(nEN)=nBitErrors(nEN)+err;
    end % End of while loop
    sim_res=[nBitErrors nTransmittedFrames]
end %end for (SNR points)
sim_res=[nBitErrors nTransmittedFrames]


semilogy(snr_db, nBitErrors./nTransmittedFrames/nBitsPerFrame, 'r-x')
grid on
legend("Viterbi Algorithm","Interpreter","Latex")
xlabel('$\frac{E_{b}}{N_{0}}(SNR)$','Interpreter','Latex')
ylabel('BER')
