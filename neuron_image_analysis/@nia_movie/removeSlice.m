function removeSlice(obj, slice_idx)
%REMOVESLICE Remove a slice from the movie
%   removeSlice(slice_idx) removes a slice from the given movie. Note that
%   this function may modify the slice index associated with other frames.
%   If an error occurs, then the movie is unmodified.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

% Remove entry from slices array
obj.slices(slice_idx) = [];

% Invalidate caches
obj.pos_lu = {};
obj.pos_ranges = {};
obj.ch_list = [];
obj.hist_info = [];

end