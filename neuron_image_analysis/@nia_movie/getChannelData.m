function ch_data = getChannelData(obj, slice_idx, ch_idx)
%GETCHANNELDATA Retrieves image for a single channel of a single slice
%   ch_data = getChannelData(slice_idx, ch_idx) retrieves the image data
%   for a single channel of a single slice of the movie. The returned
%   channel data is a structure with the fields 'ch' and 'image'. The 'ch'
%   field contains the channel identifier, and the 'image' field contains
%   the image data.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

if ~nia_isScalarInteger(ch_idx) || ...
        ch_idx < 1 || ch_idx > length(obj.slices(slice_idx).channels)
    error 'The argument ''ch_idx'' has an invalid value';
end

% Retrieve the channel info
ch_data.ch = obj.slices(slice_idx).channels(ch_idx).ch;
ch_data.image = obj.slices(slice_idx).channels(ch_idx).image;
end