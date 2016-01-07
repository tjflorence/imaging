function allowed = inChannelRange(obj, slice_idx, ch_idx, rnge)
%INCHANNELRANGE Check if the given channel is in the range
%   allowed = inChannelRange(slice_idx, ch_idx, rnge) checks if the passed slice
%   is in the passed range. Returns true if it is in range and false
%   otherwise. This function assumes that the slice is in range as
%   per the function inSliceRange().
%
%   This function accepts the following arguments:
%
%       slice_idx - Index of slice to test.
%
%       ch_idx - Index of channel to test.
%
%       rnge - The range of slices to to use for the test. See the 
%           function isValidRange() for a description of the format
%           of the range. This function requires a normalized range
%           as produced by the function normalizeRange().
%
%   Note that for efficiency this function performs no error checking
%   on the input arguments.

% Note that in the following comparisons we make use of the
% fact that comparisons with NaN values are always false (see
% IEEE 754) to provide NaN behavior specified in documentation.

ch = obj.slices(slice_idx).channels(ch_idx).ch;
allowed = nnz(rnge.channels == ch) > 0;
        
end
