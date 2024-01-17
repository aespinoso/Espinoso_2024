function MRL = mean_resultant_length(data)
% Computes the mean resultant length See Ref. [A]
% Inputs:
%         - data: matrix of which we want to compute the mean resultant
%           length (LxN  (L number signals x N number of samples))
% Outputs:
%          - MRL: mean resultant length.
%
% References:
% [A] Andrzejak, R. G., Espinoso, A., García-Portugués, E., Pewsey, A.,
%     Epifanio, J., Leguia, M. G., & Schindler, K. (2023). High expectations
%     on phase locking: Better quantifying the concentration of circular data.
%     Chaos: An Interdisciplinary Journal of Nonlinear Science, 33(9).
%
% -------------------------------------------------------------------------

N = size(data,2); % number of samples
L = size(data,1); % number of channels

% Hilbert protophase
for ii = 1:L
    y_h = hilbert(data(ii,:));
    phase_ii = (atan2(imag(y_h),real(y_h)));
    phase_ii = phase_ii(round(0.05*N:0.95*N)); % cut phase 5% beginning and end
    phase_h(ii,:) = phase_ii;
end

% Sum across channels (across N)
MPS_t = abs(mean(exp(1i*phase_h)));
% Sum accros samples (across L)
MRL = abs(mean(MPS_t)); % Mean resultant length

end

