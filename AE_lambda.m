function [lambda,p_SOZ,p_nonSOZ] = AE_lambda(Test,SOZ_chan,nonSOZ_chan,varargin)
% Function that computes lammbda values for every window and channels from
% the database of Bern seizures from the univariate phase-based tests
% INPUTS:
%        - Test: binary matrix of results of phase-based test. The
%        dimensions of the matrix are channels (L) x windows (W)
%        - SOZ_chan: SOZ channels from the analysed patient
%        - nonSOZ_chan: nonSOZ channels, complementary of SOZ_chan
%     * The nexts are varargin 1 -> 2 *
%        - idx_1: idx of windows to start analysis (varagin 1)
%        - idx_2: idx of windows to finish analysis (varargin 2)
%
% OUTPUTS:
%        - lambda: lambda value to see contrast between SOZ and nonSOZ
%          rejections.
%        - p_SOZ: fraction of rejections SOZ channels
%        - p_nonSOZ: fraction of rejections nonSOZ channels
%--------------------------------------------------------------------------

%% i) Channel category and period obtention from test
if (isempty(varargin)) || isempty(varargin{1}) % we want all columns (all windows) (no input or empty en varargin 1)
    Test_SOZ = Test(SOZ_chan,:);
    Test_nonSOZ = Test(nonSOZ_chan,:);
else % from determined period
    idx_1 = varargin{1};
    idx_2 = varargin{2};

    Test_SOZ = Test(SOZ_chan,idx_1:idx_2);
    Test_nonSOZ = Test(nonSOZ_chan,idx_1:idx_2);
end

%% ii) Fraction of rejections SOZ channels
p_SOZ = nanmean(nanmean(Test_SOZ')); 

%% iii) Fraction of rejections nonSOZ channels
p_nonSOZ = nanmean(nanmean(Test_nonSOZ'));

%% iv) Lambda
lambda = (p_SOZ-p_nonSOZ)/(p_SOZ+p_nonSOZ);

end