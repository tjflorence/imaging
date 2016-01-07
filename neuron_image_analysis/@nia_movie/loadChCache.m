function loadChCache(obj)
%LOADCHCACHE Prepare channel list for the movie
%   loadChCache() scans through the current movie contents and prepares the
%   list of channel identifiers. If a channel list cache is already
%   available, then this function returns immediately.

% Check if we already have a channel list
if ~isempty(obj.ch_list)
    return;
end

% Check for trivial cast
if isempty(obj.slices)
    return;
end

% Traverse the movie
ch_list = [];

for slice_idx=1:length(obj.slices)
    for ch_idx=1:length(obj.slices(slice_idx).channels)
        new_ch = obj.slices(slice_idx).channels(ch_idx).ch;
        ch_list = unique([ch_list, new_ch]);
    end
end

obj.ch_list = ch_list;

end

