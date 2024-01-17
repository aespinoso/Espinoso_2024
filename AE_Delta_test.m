function [d,Test_D_A,Test_D_B] = AE_Delta_test(signal,Fs,filt_idx)
% Function that computes the contribution measure and the corresponding two
% one sided surrogate test for a selected filtering setting.
%
% INPUTS:
%          - signal: multichannel EEG recordings, LxN (number of channels x
%            number of samples).
%          - Fs: sampling frequency
%          - filt_idx: filtering setting to use:
%                   -> filt_idx = 1: low frequency [4,30] Hz [1]
%                   -> filt_idx = 2: high frequency [80,150] Hz [1]
%
% OUTPUTS:
%          - d: phase-locking contribution measure L x num_win x Ss (number
%            of channels x number of windows x (num surrogates + 1 original))
%          - Test_D_A: outcome A of the phase-locking contribution test 
%            L x num_win (number of channels x number of windows). 
%            Rejection (1) or non rejection (0) in the direction of outcome
%            A.
%          - Test_D_B: outcome B of the phase-locking contribution test 
%            L x num_win (number of channels x number of windows).
%            Rejection (1) or non rejection (0) in the direction of outcome
%            B.
%          
% References
% [1] Bandarabadi, M., Gast, H., Rummel, C., Bassetti, C., Adamantidis, A.,
%     Schindler, K., & Zubler, F. (2019). Assessing epileptogenicity using
%     phase-locked high frequency oscillations: a systematic comparison of
%     methods. Frontiers in neurology, 10, 1132.
% [2] Andrzejak, R. G., Espinoso, A., García-Portugués, E., Pewsey, A.,
%     Epifanio, J., Leguia, M. G., & Schindler, K. (2023). High expectations
%     on phase locking: Better quantifying the concentration of circular data.
%     Chaos: An Interdisciplinary Journal of Nonlinear Science, 33(9).
% -------------------------------------------------------------------------

%% i) Information of signal

L = size(signal,2); % number of channels
N = size(signal,1); % number of samples
time = [0:N]/Fs; time = time(1:end-1); % time vector

%% ii) Parameter of the surrogates
ParamSurro.Number = 19;        % Number of surrogates
ParamSurro.MaxIter = 120;      % Maximal Iterations
ParamSurro.type = 1;           % 1: perfect amplitudes, 2: perfect periodogram

%% iii) Filtering setting
EEG_orig = signal;
if filt_idx == 1 % low freq
    freqLF = [4 30];
    orderLF = round(3*(Fs/freqLF(1))); % filter order for LF

    fir1CoefLF = fir1(orderLF,[freqLF(1),freqLF(2)]./(Fs/2)); % filter coeff. for LF
    EEG = filtfilt(fir1CoefLF,1,EEG_orig);
elseif filt_idx == 2 % high freq
    freqHFO = [80 150];
    orderHFO = round(3*(Fs/freqHFO(1))); % filter order for HFO

    fir1CoefHFO = fir1(orderHFO,[freqHFO(1),freqHFO(2)]./(Fs/2)); % filter coeff. for HFO
    EEG = filtfilt(fir1CoefHFO,1,EEG_orig);
end

%% iv) Windows info
time_w = 20; % length windows
time_ol = 5; % overlap = time_w - time_ol

all_win = 0:time_ol*Fs:N-time_w*Fs;
num_win = length(all_win); % total number of windows
time_vec_w = 0:time_ol:N/Fs-time_w; % time vector with windows

%% v) Analysis window-wise mean resultant length

% Number of workers to perform the parfor
pool = parpool; % Automatically determines the number of workers
numWorkers = pool.NumWorkers;

myCluster = parcluster('local');
myCluster.NumWorkers = numWorkers;
saveProfile(myCluster);

parpool('local',numWorkers)

all_chan = [1:L]; % vector of num channels

% Start vectors mean resultant length
Ss = ParamSurro.Number + 1; % number of surrogates + 1 original signal
MRL_ww_all = zeros(num_win,Ss); % windows x (num surrogates + 1)
MRL_ww_chan = zeros(L,num_win,Ss); % channels x windows x (num surrogates + 1)

% Computation mean resultant length for every channel and surrogate
parfor (ss = 1:Ss,numWorkers)

    idx_w = 1;
    MRL_kk_w = zeros(L,num_win); % channels x windows
    MRL_w = zeros(num_win,1); % 1 x windows

    for ww = 0:time_ol*Fs:N-time_w*Fs % windows

        interval = [ww + 1:time_w*Fs + ww];

        if ss == 1 % Original
            EEG_ss = EEG';
            EEG_ww = EEG_ss(:,interval);
        else % Surrogates
            EEG_ss = EEG';
            EEG_ww = ASR_SurrogateMulti(EEG_ss(:,interval),ParamSurro);
        end

        % MRL window-wise all network

        [MRL_w(idx_w)] = mean_resultant_length(EEG_ww); 

        % MRL from all the network less selected channel window-wise

        for kk = 1:L
            idx_MRL = setdiff(all_chan,kk);
            [MRL_kk_w(kk,idx_w)] = mean_resultant_length(EEG_ww(idx_MRL,:));
        end

        idx_w = idx_w + 1;

    end

    disp(['MPL Surrogate ',num2str(ss),' window-wise'])

    MRL_ww_all(:,ss) = MRL_w; % MRL for all network for num surrogates + original
    MRL_ww_chan(:,:,ss) = MRL_kk_w; % MRL for all network less selected channel for num surrogates + original

end

%% vi) Computation re-normalized mean resultant length and contribution measure

% Expected value see Ref. [2]
Z_exp_all = 0.5*sqrt(pi/L); % expected value all network
Z_exp_chan = 0.5*sqrt(pi/(L-1)); % expected value when extracting one channel from the network

% Start vector contribution measure
d = zeros(L,num_win,Ss);

for ss = 1:Ss

    % MRL all network windows-wise
    MRL_all_norm = (MRL_ww_all(:,ss)-Z_exp_all)./(1-Z_exp_all);

    % MRL all netw except chan window-wise
    MRL_chan_norm = (MRL_ww_chan(:,:,ss)-Z_exp_chan)./(1-Z_exp_chan);

    d(:,:,ss) = MRL_all_norm' - MRL_chan_norm;

end

%% vii) Computation two one-sided surrogate test

% Test D for each channel
Test_D_A = zeros(L,size(d,2)); % test outcome A
Test_D_B = zeros(L,size(d,2)); % test outcome B

for kk = 1:L % number channels
    for w = 1:size(d,2) % number windows
        if d(kk,w,1) > max(d(kk,w,2:20)) % test max
            Test_D_A(kk,w) = 1;
        end

        if d(kk,w,1) < min(d(kk,w,2:20)) % test min
            Test_D_B(kk,w) = 1;
        end
    end
end

