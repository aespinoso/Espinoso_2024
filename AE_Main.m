% This is the main script to calculate the methods described in:
% "The part and the whole: how single nodes contribute to large-scale phase locking in
%  functional EEG networks"
%
% It calls the 3 main functions:
% 1) AE_Delta_test: computes the phase-locking contribution measure and
%    surrogate-test for a LxN (channels x samples)
% 2) AE_joint_test: computes the joint test for low and high frequencies
%    from the previous function
% 3) AE_lambda: computes the relative difference from the previous results
%    between SOZ and nonSOZ channels
%
% * Note: the 3 functions can be used separately, but in this script we
% show how to call them in order to obtain the steps followed in the paper

% Author: Anaïs Espinoso, 2024
% Contact: anais.espinoso@gmail.com
%..........................................................................

%% 1) Phase-locking contribution measure and test

signal = ; % load multichannel recordings
Fs = 512; % sampling frequency (examplary used in the paper)

% i) Low frequencies [4-30] Hz
[D_L,Test_low_A,Test_low_B] = AE_Delta_test(signal,Fs,1);

% ii) Low frequencies [80-150] Hz
[D_H,Test_high_A,Test_high_B] = AE_Delta_test(signal,Fs,2);

%% 2) Joint test

[Test_joint_A,Test_joint_B,Test_low,Test_high,Test_joint] = AE_joint_test(Test_low_A,Test_low_B,Test_high_A,Test_high_B);

%% 3) Relative difference between SOZ and nonSOZ channels
% This will give a similar table as table II, with one extra column of
% "all" periods

% i) Types of channels
SOZ_chan = [1,2]; % indices from the array of channels corresponding to SOZ channels (example)
nonSOZ_chan = [3,4]; % indices from the array of channels corresponding to SOZ channels (example)

% ii) We join the results computed in 1 and 2
Tests_A_B = {Test_low_A,Test_high_A,Test_joint_A,Test_low_B,Test_high_B,Test_joint_B}; % all tests with outcomes A and B

% iii) Initialize what we want to save, in this case every lambda por the
% different tests
lambda = zeros(1,6); % (low A,high A,joint A,low B,high B,joint B)
lambda_b = lambda;
lambda_d = lambda;
lambda_a = lambda;

% iv) Obtantion of lambda in different periods
for ii = 1:6 % type test (low A,high A,joint A,low B,high B,joint B)

    Test = Tests_A_B{ii};

    for jj = 1:4 % periods (all periods,before,during,after)

        if jj == 1
            % i) All time windows
            [lambda(ii),p_res,p_nores] = AE_lambda(Test,SOZ_chan,nonSOZ_chan);
        elseif jj == 2
            % ii) Before seizure
            idx_before_1 = 1; % index of window starting before seizure period (example)
            idx_before_2 = 5; % index of window finishing before seizure period (example)
            [lambda_b(ii),p_res_b,p_nores_b] = AE_lambda(Test,SOZ_chan,nonSOZ_chan,idx_before_1,idx_before_2);
        elseif jj == 3
            % iii) During seizure
            idx_during_1 = 6; % index of window starting during seizure period (example)
            idx_during_2 = 10; % index of window finishing during seizure period (example)
            [lambda_d(ii),p_res_d,p_nores_d] = AE_lambda(Test,SOZ_chan,nonSOZ_chan,idx_during_1,idx_during_2);
        elseif jj == 4
            % iv) After seizure
            idx_after_1 = 11; % index of window starting during seizure period (example)
            idx_after_2 = 15; % index of window finishing during seizure period (example)
            [lambda_a(ii),p_res_a,p_nores_a] = AE_lambda(Test,SOZ_chan,nonSOZ_chan,idx_after_1,idx_after_2);
        end

    end
end

% v) Obtained table

% Names of rows and columns
rows = {'DL_A','DH_A','DJ_A','DL_B','DH_B','DJ_B'};
columns = {'All periods', 'Before', 'During', 'After'};

% Joint the data to show
data_table = [lambda',lambda_b',lambda_d',lambda_a'];

% Create a table
T = array2table(data_table, 'RowNames', rows, 'VariableNames', columns);
disp(T);
