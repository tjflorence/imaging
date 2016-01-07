function ch_idx = addChannelData(obj, slice_idx, ch_data)
%ADDCHANNELDATA Add channel data entry to the slice
%   ch_idx = addChannelData(slice_idx, ch_data) adds a channel to the given
%   slice and returns the channel index for the new entry. The passed
%   channel data must be a structure with the fields 'ch' and 'image'. The
%   'ch' field contains the channel identifier, and the 'image' field
%   contains the image data. If an error occurs, then the channel
%   information is unmodified.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

if ~isstruct(ch_data) || length(ch_data) ~= 1
    error 'The argument ''ch_data'' has an invalid type';
end

ch_data_names = {'ch', 'image'};
[ch_data_ok, ch_data_msg] = nia_hasValidFieldNames(...
    ch_data, ch_data_names, ch_data_names);

if ~ch_data_ok
    error(ch_data_msg, 'ch_data');
end

if ~nia_isScalarInteger(ch_data.ch)
    error 'The argument ''ch_data.ch'' has an invalid value';
end

if ~isnumeric(ch_data.image)
    error 'The argument ''ch_data.image'' has an invalid value';
end

% Append to end of channels array
if isempty(obj.slices(slice_idx).channels)
    obj.slices(slice_idx).channels.ch = ch_data.ch;
    obj.slices(slice_idx).channels.image = ch_data.image;
else
    obj.slices(slice_idx).channels(end+1) = ch_data;
end

ch_idx = length(obj.slices(slice_idx).channels);

% Invalidate cache information
obj.ch_list = [];
obj.hist_info = [];

end