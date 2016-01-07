function loadHistCache(obj)
%LOADHISTCACHE Prepare histogram info for the movie
%   loadHistCache() scans through the current movie contents and
%   prepares the histogram cache. If a histogram cache is already
%   available, then this function returns immediately.

cfg_num_hist_bins = 100;

% Check if we already have histogram information
if ~isempty(obj.hist_info)
    return;
end

% Check for the trivial case
if isempty(obj.slices)
    return;
end

% Check if we need a new channel list
if isempty(obj.ch_list)
    obj.loadChCache();
end

% Preallocate histogram info
num_ch = length(obj.ch_list);
hist_info = [];
hist_info(num_ch).min = [];
hist_info(num_ch).max = [];
hist_info(num_ch).dist = [];
hist_info(num_ch).edges = [];

% Traverse each slice to get ranges
for slice_idx=1:length(obj.slices)
    channels = obj.slices(slice_idx).channels;
    
    for ch_idx=1:length(channels)
        ch = find(obj.ch_list == channels(ch_idx).ch);
                
        ch_data = channels(ch_idx).image;
        ch_data = ch_data(~isnan(ch_data));
        
        hist_info(ch).min = min([hist_info(ch).min, min(ch_data(:))]); %#ok<AGROW>
        hist_info(ch).max = max([hist_info(ch).max, max(ch_data(:))]); %#ok<AGROW>
    end
end

% Set zero to one scale for channels that have no data
for ch_idx=1:length(hist_info)
    if isempty(hist_info(ch_idx).min) || isempty(hist_info(ch_idx).max)
        hist_info(ch_idx).min = 0; %#ok<AGROW>
        hist_info(ch_idx).max = 1; %#ok<AGROW>
    end
end

% Identify histogram edges
for ch_idx=1:length(hist_info)
    edges = linspace(...
        hist_info(ch_idx).min, ...
        hist_info(ch_idx).max, ...
        cfg_num_hist_bins + 1);
    
    hist_info(ch_idx).edges = edges; %#ok<AGROW>
    hist_info(ch_idx).dist = zeros(2, cfg_num_hist_bins); %#ok<AGROW>
    hist_info(ch_idx).dist(1,:) = 0.5*(edges(1:end-1) + edges(2:end)); %#ok<AGROW>
end

% Traverse each slice to get histogram
for slice_idx=1:length(obj.slices)
    channels = obj.slices(slice_idx).channels;
    
    for ch_idx=1:length(channels)
        ch = find(obj.ch_list == channels(ch_idx).ch);
 
        ch_data = channels(ch_idx).image;
        ch_data = ch_data(~isnan(ch_data));
        
        counts = histc(ch_data, hist_info(ch).edges)';
        
        % There is a bug in the MATLAB implementation of histc that causes
        % it to return an empty array if the input array is empty. The
        % correct behavior would of course be to return an array of zeros.
        % we have to catch this case and workaround.
        if isempty(counts)
            counts = zeros(1, length(hist_info(ch).edges));
        end
        
        hist_info(ch).dist(2,:) = hist_info(ch).dist(2,:) + counts(1:end-1); %#ok<AGROW>
    end
end

% Remove temporary edges field
hist_info = rmfield(hist_info, 'edges');

obj.hist_info = hist_info;

end

