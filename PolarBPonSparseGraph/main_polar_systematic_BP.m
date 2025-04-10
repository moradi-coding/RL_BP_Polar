close all; clc; clear; format short;

% for parfor
CoreNum = 2;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(CoreNum);
else
    disp('matlab pool already started');
end

% output folder
dirname = 'output';
if ~exist(dirname, 'dir')
    mkdir(dirname)
end

% parameters
N = 64;  K = 32;  Rc = K/N;
% load('H256.mat')
load('H64.mat')
maxBP_iter = 200;  

% rate profile
SNR_Cons = 2;
type = 1; % 1--> GA, 3--> Tse_RMpolar, 4 --> RM
RP = RM_Polar_Profile(N, K, SNR_Cons, type);
RP = bitrevorder(RP);

EbNo_vec = 1:0.5:5;
maxRun = 1E8;
maxFE = [1000*ones(1,length(EbNo_vec)-3) 500 500 500];
ClustNum = [1e3 1e3 1e3 1e3 2e3 1e4*ones(1,length(EbNo_vec)-5)];

FER = zeros(1,length(EbNo_vec));
FE = zeros(1,length(EbNo_vec));  % #of frame errors
BER = zeros(1,length(EbNo_vec));
Nruns = zeros(1,length(EbNo_vec)); % #of actual runs

decoder = comm.LDPCDecoder('ParityCheckMatrix',H ,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount', maxBP_iter,'NumIterationsOutputPort',1,'IterationTerminationCondition', 'Parity check satisfied');

fprintf('-------------------------------------\n');
for EbNo_count = 1:length(EbNo_vec)
    tic;
    EbNo_dB = EbNo_vec(EbNo_count);
    EbNo = 10^(EbNo_dB/10);
    sigma = 1/sqrt(2*Rc*EbNo);

    Nblkerrs = 0;
    Nbiterrs = 0;
    fprintf('[%02d:%02d] Starting! SNR = %.1f\n',0,0, EbNo_dB);
    for i = 1:maxRun/ClustNum(EbNo_count)
        parfor j = 1:ClustNum(EbNo_count)
            %Generate random message
            msg = randi([0 1], K, 1);
            code = zeros(N, 1);
            code(RP) = msg;
            % Polar transformation
            x = polarTransform(code, RP);
            x(~RP) = 0; % systematic
            x = polarTransform(x, RP);

            % BPSK modulation
            modulated = 1 - 2*x;
            % AWGN
            r = modulated + randn(N, 1)*sigma;
            % Decoding
            Lch = 2*r./sigma.^2;  % calc LLRs
        
            % extend Lch vector by 0 positions
            Lch_ext = zeros(size(H, 2), 1);
            Lch_ext((end - N+1) : end) = Lch; % assuming the last positions are channel positions
            [LLRxhat, iter] = step(decoder, Lch_ext);     %and decode
%             numiter(i) = iter;    % save iteration count

            xhat = (-sign(LLRxhat))/2 + 0.5;
            c = xhat((end - N+1) : end);
            mas_cap = c(RP);

            Nblkerrs = Nblkerrs + any(mas_cap~=msg);
            Nbiterrs = Nbiterrs + sum(mas_cap~=msg);
        end
        if Nblkerrs >= maxFE(EbNo_count)
            break;
        end
        t = toc;
        elapsed_m = t/60;
        elapsed_h = floor(elapsed_m/60);
        elapsed_m = floor(elapsed_m - elapsed_h*60);
        fprintf('[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,(i*ClustNum(EbNo_count)),Nblkerrs);
    end
    t = toc;
    elapsed_m = t/60;
    elapsed_h = floor(elapsed_m/60);
    elapsed_m = floor(elapsed_m - elapsed_h*60);
    temp = (i*ClustNum(EbNo_count));
    fprintf(2,'[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,temp,Nblkerrs);
    
    FER(EbNo_count) = Nblkerrs/temp;
    FE(EbNo_count) = Nblkerrs;
    BER(EbNo_count) = Nbiterrs./ temp/K;
    Nruns(EbNo_count) = temp;
    fprintf('-------------------------------------\n');
end



% Creating sim_state
res.Nruns = Nruns;
res.FE = FE;
res.FER = FER;
res.BER = BER;
res.maxBP_iter = maxBP_iter;

res.SNR = EbNo_vec;
res.K = K;
res.N = N;
res.maxRun = maxRun;
res.maxFE = maxFE;


filename = [dirname sprintf('/polarBP_PRP_N%d_K%d_maxBP_iter%d_scaleMinSum.mat',N,K,maxBP_iter)];
save(filename,'res');
