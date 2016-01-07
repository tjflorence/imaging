function nia_playFlatMovie(mov)
%NIA_PLAYFLATMOVIE Play the movie in the passed matrix.
%   nia_playFlatMovie(mov) displays the passed movie. The argument mov must
%   be a 4D array, where the first two dimensions correspond to height and
%   width, the third dimension corresponds to time, and the last dimension
%   corresponds to channel.
%
%   Example:
%       mov = rand(100, 100, 100, 2);
%       nia_playFlatMovie(mov);

% Check arguments
if ~isfloat(mov) || ~isreal(mov) || ...
        ndims(mov) > 4 || isempty(mov)
    error 'The argument ''mov'' has an invalid type';
end

callback_list.getImage = @getImage;
callback_list.getChannels = @getChannels;
callback_list.getPosRanges = @getPosRanges;
callback_list.getHistInfo = @getHistInfo;
callback_list.processROIs = @processROIs;

nia_playGenericMovie(mov, callback_list);
end

function [im, time] = getImage(mov, pos, colormap_spec)
% This function retrieves a single frame from the movie.

accum = zeros(size(mov, 1), size(mov, 2), 3);

for ch_idx=1:size(mov, 4)
    spec = colormap_spec(ch_idx);
    
    ch_im = nia_applyColormap(...
        squeeze(mov(:,:,pos,ch_idx)), ...
        [spec.min, spec.max], ...
        spec.colormap, spec.nan_color);
    
    accum = accum + ch_im;
end

accum(accum < 0) = 0;
accum(accum > 1) = 1;

im = accum;
time = pos;

end

function out = getChannels(mov)
% This function retrieves the list of channel identifiers

    out = 1:size(mov,4);
end

function out = getPosRanges(mov)
% This function retieves the valid position ranges

    out = [1; size(mov,3)];
end

function out = getHistInfo(mov)
% This function retrieves information on the pixel
% intensity histgrams for the movie.

cfg_max_sample = 10000;
cfg_hist_bins_denom = 5;
cfg_hist_bins_max = 64;

% prealloc histogram info
num_ch = size(mov,4);
hist_info(num_ch).min = [];
hist_info(num_ch).max = [];
hist_info(num_ch).dist = [];

for ch_idx=1:num_ch
    ch_data = mov(:,:,:,ch_idx);
    ch_data = ch_data(~isnan(ch_data));
    
    if numel(ch_data) > cfg_max_sample
        sample_indices = randperm(cfg_max_sample);
    else
        sample_indices = 1:numel(ch_data);
    end
    
    sampled_data = ch_data(sample_indices);
    
    min_val = min(sampled_data);
    max_val = max(sampled_data);
    
    % the following are just reasonable heuristics
    % for picking histograms edges
    if min_val < 0
        set_min_val = min_val;
        set_max_val = max_val;
    else
        set_min_val = 0;
        set_max_val = 2^ceil(log2(max_val));
    end
    
    hist_bins = ceil(numel(sampled_data)/cfg_hist_bins_denom);
    if hist_bins > cfg_hist_bins_max
        hist_bins = cfg_hist_bins_max;
    end
    
    hist_info(ch_idx).min = set_min_val;
    hist_info(ch_idx).max = set_max_val;
    
    hist_edges = linspace(0, set_max_val, hist_bins+1);
    dist = zeros(2, hist_bins);
    
    dist(1,:) = 0.5*(hist_edges(1:end-1) + hist_edges(2:end));
    tmp = histc(sampled_data, hist_edges);
    dist(2,:) = tmp(1:end-1);
    
    hist_info(ch_idx).dist = dist;
end

out = hist_info;

end

function dset = processROIs(mov, roi_list)
% This function processes the passed ROIs against the
% passed movie in order to produce a cell array of time
% series

dset = cell(1, length(roi_list));

for idx=1:length(roi_list)
    mask = createMask(roi_list(idx).handle);
    
    vec = zeros(2, size(mov,3));
    vec(1,:) = 1:size(mov,3);
    
    % Note: the movie must be passed whole without any array
    % subarray operations, otherwise matlab will duplicate
    % the array prior to passing it, this will obliterate
    % performance. this is easy to workaround, we just do
    % the subarray work in the evaluation function.
    vec(2,:) = nia_scanROIFlatMovie(...
        mov, mask, roi_list(idx).channel);
    
    dset{idx} = vec;
end

end
