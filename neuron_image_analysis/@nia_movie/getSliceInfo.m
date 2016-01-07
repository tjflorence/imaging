function slice_info = getSliceInfo(obj, slice_idx)
%GETSLICEINFO Retrieves metadata for a single slice
%   slice_info = getSliceInfo(slice_idx) retrieves the metadata for a
%   single slice of the movie. The returned slice info is a structure with
%   the fields 'time', 'pos', and 'num_ch'. The 'time' field is a scalar
%   double that contains the timestamp, the 'pos' field is a cell array of
%   row vectors describing the slice positions, and the 'num_ch' field is
%   the number of channels in the slice.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

% Retrieve slice info
slice_info.time = obj.slices(slice_idx).time;
slice_info.pos = obj.slices(slice_idx).pos;
slice_info.num_ch = length(obj.slices(slice_idx).channels);
end