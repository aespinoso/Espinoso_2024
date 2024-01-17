function [Test_joint_A,Test_joint_B,Test_low,Test_high,Test_joint] = AE_joint_test(Test_low_A,Test_low_B,Test_high_A,Test_high_B)
% Function to compute the joint test from results of
% "AE_Contribution_Main". To do so, we need the results of low and high
% frequencies from outcome A and B.
%
% INPUTS:
%          - Test_low_A: outcome A of the phase-locking contribution test
%            L x num_win (number of channels x number of windows) for low
%            freqs. Rejection (1) or non rejection (0).
%          - Test_low_B: outcome B of the phase-locking contribution test
%            L x num_win (number of channels x number of windows) for low
%            freqs. Rejection (1) or non rejection (0).
%          - Test_high_A: outcome A of the phase-locking contribution test
%            L x num_win (number of channels x number of windows) for high
%            freqs. Rejection (1) or non rejection (0).
%          - Test_high_B: outcome B of the phase-locking contribution test
%            L x num_win (number of channels x number of windows) for high
%            freqs. Rejection (1) or non rejection (0).
%
% OUTPUTS:
%          - Test_joint_A: outcome A of joint test L x num_win (number of
%            channels x number of windows). Rejection (1) or non rejection
%            (0).
%          - Test_joint_B: outcome B of joint test L x num_win (number of
%            channels x number of windows). Rejection (1) or non rejection
%            (0).
%          - Test_low: outcome A (1), outcome B (-1) and outcome C (0) for
%            phase-locking contribution test applied to low frequencies.
%          - Test_high: outcome A (1), outcome B (-1) and outcome C (0) for
%            phase-locking contribution test applied to high frequencies.
%          - Test_joint: outcome A (1), outcome B (-1) and outcome C (0) for
%            joint test.
%--------------------------------------------------------------------------

%% i) Computation joint test from low and high tests

% Outcome A
Test_joint_A = Test_low_A .* Test_high_A;
% Outcome B
Test_joint_B = Test_low_B .* Test_high_B;

%% ii) Put in a same matrix outcome A and B for graphical purposes

% Low frequency
Test_low_B = Test_low_B * -1;

Test_low = zeros(size(Test_low_A));
Test_low(Test_low_A == 1) = 1;
Test_low(Test_low_B == -1) = -1;

% High frequency
Test_high_B = Test_high_B * -1;

Test_high = zeros(size(Test_high_A));
Test_high(Test_high_A == 1) = 1;
Test_high(Test_high_B == -1) = -1;

% Joint test
Test_joint_B = Test_joint_B * -1;

Test_joint = zeros(size(Test_joint_A));
Test_joint(Test_joint_A == 1) = 1;
Test_joint(Test_joint_B == -1) = -1;

end