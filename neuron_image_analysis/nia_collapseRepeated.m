function [tp,yp,ep] = nia_collapseRepeated(t, y, offset, period, nsamples)
%NIA_COLLAPSEREPEATED Collapse a repeated signal to its median
%   [tp,yp,ep] = nia_collapseRepeated(t, y, offset, period, nsamples
%   collapses a wave, by taking the median over multiple waves.
%
%   This function takes the following arguments:
%
%       t - A [1xN] array specifying the time for each y value
%
%       y - A [...xN] array specifying the values to be collapsed. The
%           array y may be a multi-dimensional array.  The last dimension
%           is always assumed to represent time, and must have an identical
%           length to the last dimension of t.
%
%       offset - A scalar double representing the time point at which the
%           waveform starts.
%
%       period - A scalar double representing the period of the waveform.
%
%       nsamples - A scalar integer representing the number of time points
%           in the output waveform. All data points are interpolated so
%           that the averaged values align to the output waveform.
%
%   This function has the following outputs:
%
%       tp - A [1x(nsamples)] array storing the interpolated time values.
%
%       yp - A [...x(nsamples)] array storing the collapsed values.
%
%       ep - A [...x(nsamples)] array storing the standard error of the
%           mean of the collapsed values.

% Check arguments
if ~isfloat(t) || ~isrow(t)
    error 'The argument 't' has an invalid type';
end

y_orig_size = size(y);

if ~isfloat(y) || y_orig_size(end) ~= size(t,2)
    error 'The argument 'y' has an invalid type';
end

% Transform t and y into standard 2D form, with time in first dimension
% and different measures in different columns
t = t';
y = shiftdim(y, length(y_orig_size)-1);
y = reshape(y, size(t,1), prod(y_orig_size(1:end-1)));

% Find the time base
t_interp = linspace(0, period, nsamples)';
num_cycles = floor((max(t)-offset) ./ period);

% Create interpolated y, with cycles in last dimension
y_interp = zeros(nsamples, size(y,2), num_cycles);

% Interpolate to common time base
for idx=1:num_cycles
    y_interp(:,:,idx) = interp1(t, y, ...
        t_interp + offset + (idx-1) * period);
end

% Take the median and mean
yp = median(y_interp, 3);
ep = std(y_interp, 1, 3)./sqrt(num_cycles);
clear y_interp;

% Reshape the outputs.
tp = t_interp';

yp = reshape(yp, [nsamples, y_orig_size(1:end-1)]);
yp = shiftdim(yp, 1);

if nargout >= 3
    ep = reshape(ep, [nsamples, y_orig_size(1:end-1)]);
    ep = shiftdim(ep, 1);
else
    ep = [];
end

end