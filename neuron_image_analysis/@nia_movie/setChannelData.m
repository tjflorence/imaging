function setChannelData(obj, slice_idx, ch_idx, ch_data)
%SETCHANNELDATA Sets channel data for a single channel of a single slice
%   setChannelData(slice_idx, ch_idx, ch_data) modifies the image data for
%   a single channel of a single slice of the movie. The passed channel
%   data must be a structure with the fields 'ch' and 'image'. The 'ch'
%   field contains the channel identifier, and the 'image' field contains
%   the image data. Any of these fields may be omitted, in which case that
%   element of the channel info is unmodified. If an error occurs, then the
%   channel information is unmodified.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

if ~nia_isScalarInteger(ch_idx) || ...
        ch_idx < 1 || ch_idx > length(obj.slices(slice_idx).channels)
    error 'The argument ''ch_idx'' has an invalid value';
end

if ~isstruct(ch_data) || length(ch_data) ~= 1
    error 'The argument ''ch_data'' has an invalid type';
end

ch_data_names = {'ch', 'image'};
[ch_data_ok, ch_data_msg] = nia_hasValidFieldNames(...
    ch_data, ch_data_names, {});

if ~ch_data_ok
    error(ch_data_msg, 'ch_data');
end

if isfield(ch_data, 'ch')
    if ~nia_isScalarInteger(ch_data.ch)
        error 'The argument ''ch_data.ch'' has an invalid value';
    end
end

if isfield(ch_data, 'image')
    if ~isnumeric(ch_data.image)
        error 'The argument ''ch_data.image'' has an invalid value';
    end
end

% Retrieve the channel info
if isfield(ch_data, 'ch')
    obj.slices(slice_idx).channels(ch_idx).ch = ch_data.ch;

    obj.ch_list = [];
    obj.hist_info = [];
end

if isfield(ch_data, 'image')
    obj.slices(slice_idx).channels(ch_idx).image = ch_data.image;
    
    obj.hist_info = [];
end

end