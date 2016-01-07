function removeChannelDataById(obj, slice_idx, ch_id)
%REMOVECHANNELDATABYID Removes image for a single channel of a single slice
%   removeChannelDataById(slice_idx, ch_id) removes the image data for a
%   single channel of a single slice of the movie. This function will find
%   the first channel entry with the given channel identifier. If no such
%   channel is found, then this function does nothing.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

if ~nia_isScalarInteger(ch_id)
    error 'The argument ''ch_id'' has an invalid value';
end

% Find the appropriate channel
all_ch = obj.slices(slice_idx).channels(:).ch;
ch_idx = find(all_ch == ch_id, 1, 'first');

if isempty(ch_idx)
    return;
end

% Remove channel information
obj.slices(slice_idx).channels(ch_idx) = [];

% Invalidate channel list and histograms caches
obj.ch_list = [];
obj.hist_info = [];

end