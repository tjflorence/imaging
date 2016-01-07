function allowed = inSliceRange(obj, slice_idx, rnge)
%INSLICERANGE Check if the given slice is in the range
%   allowed = inSliceRange(slice_idx, rnge) checks if the passed slice
%   is in the passed range. Returns true if it is in range and false
%   otherwise.
%
%   This function accepts the following arguments:
%
%       slice_idx - Index of slice to test.
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

time = obj.slices(slice_idx).time;
allowed = ~(time < rnge.time(1) || time > rnge.time(2));

if ~isempty(rnge.pos_vec)
    cur_pos = obj.slices(slice_idx).pos{rnge.pos_vec};
    allowed = allowed && ~any(cur_pos < rnge.pos_lims(1,:) | ...
        cur_pos > rnge.pos_lims(2,:));
end

end
