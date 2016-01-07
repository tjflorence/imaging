function removeChannelData(obj, slice_idx, ch_idx)
%REMOVECHANNELDATA Remove the specified channel
%   removeChannelData(slice_idx, ch_idx) removes the channel data
%   associated with the channel index. Note that the channel index is not
%   necessarily equal to the channel identifier, which is retrieveable
%   using the the getChannelData() function. Also note, this function 
%   may alter the channel index associated with other channels of the
%   the slice.

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

if ~nia_isScalarInteger(ch_idx) || ...
        ch_idx < 1 || ch_idx > length(obj.slices(slice_idx).channels)
    error 'The argument ''ch_idx'' has an invalid value';
end

% Remove entry from channels array
obj.slices(slice_idx).channels(ch_idx) = [];

% Invalidate channel list and histograms caches
obj.ch_list = [];
obj.hist_info = [];

end