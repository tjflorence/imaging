function ch_data = getChannelDataById(obj, slice_idx, ch_id)
%GETCHANNELDATABYID Retrieves image for a single channel of a single slice
%   ch_data = getChannelDataById(slice_idx, ch_id) retrieves the image data
%   for a single channel of a single slice of the movie. This function will
%   find the first channel entry with the given channel identifier. If no
%   such channel is found, then this function returns the empty array. The
%   returned channel data is a structure with the fields 'ch' and 'image'.
%   The 'ch' field contains the channel index, and the 'image' field
%   contains the image data.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

if ~nia_isScalarInteger(ch_id)
    error 'The argument ''ch_id'' has an invalid value';
end

% Find the appropriate channel
all_ch = [obj.slices(slice_idx).channels(:).ch];
ch_idx = find(all_ch == ch_id, 1, 'first');

if isempty(ch_idx)
    ch_data = [];
    return;
end

% Retrieve channel information
ch_data.ch = ch_idx;
ch_data.image = obj.slices(slice_idx).channels(ch_idx).image;

end