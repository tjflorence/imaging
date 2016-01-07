function nia_playMovie(mov, pos_vec)
%NIA_PLAYMOVIE Play the movie in the passed object.
%   nia_playMovie(mov, vec_idx) displays the passed movie. The argument mov
%   must be a nia_movie object, and may contain an arbitrary number of
%   channels and an arbitrary position length. The argument pos_vec
%   specifies the position vector to use when displaying the movie. The
%   argument pos_vec is optional, and may be omitted. If pos_vec is set to
%   the empty array, then the slice index is used as the position. Note
%   that for efficiency and to ensure correctness, when this function
%   displays the ROIs it will return a NaN value for any frame that does
%   not have the same size as the current frame.
%
%   Example:
%       mov = nia_movie;
%       mov.loadPrairieTSeriesXML('TSeries-001', 'TSeries-001.xml');
%       nia_playMovie(mov);

% TODO:
%   1. Implement ROI scan for pos_vec == []

% Assign default arguments
if nargin < 2
    pos_vec = 1;
end

% Check arguments
pos_lens = mov.getPositionLengths();

if ~(isempty(pos_vec) || (nia_isScalarInteger(pos_vec) && ...
        pos_vec >= 1 && pos_vec <= length(pos_lens)))
    error 'The argument ''pos_vec'' has an invalid value';
end

% Push functions into wrapper callback structure
mov_wrapper.mov = mov;
mov_wrapper.pos_vec = pos_vec;

callback_list.getImage = @getImage;
callback_list.getChannels = @getChannels;
callback_list.getPosRanges = @getPosRanges;
callback_list.getHistInfo = @getHistInfo;
callback_list.processROIs = @processROIs;

% Invoke the main function
nia_playGenericMovie(mov_wrapper, callback_list);
end

function [im,time] = getImage(mov_wrapper, pos, colormap_spec)
% This function retrieves a single frame from the movie.

mov = mov_wrapper.mov;
pos_vec = mov_wrapper.pos_vec;

% Check for trivial case
if mov.isEmpty()
    im = [];
    time = NaN;
    return;
end

% Get the slice index
if isempty(pos_vec)
    slice_idx = pos;
else
    slice_idx = mov.getSliceIndex(pos_vec, pos);
end

% Check for missing frame
if slice_idx == 0
    im = [];
    time = NaN;
    return;
end

% Get slice information
ch_list = mov.getChannelList();
slice_info = mov.getSliceInfo(slice_idx);

% If no data available, return imediately
if slice_info.num_ch == 0
    im = [];
    time = [];
    return;
end

% Add each channel to the result
accum = [];

for ch_idx=1:slice_info.num_ch
    ch_data = mov.getChannelData(slice_idx, ch_idx);
    
    ch_loc = find(ch_list == ch_data.ch);
    spec = colormap_spec(ch_loc); %#ok<FNDSB>
    
    ch_im = nia_applyColormap(ch_data.image, [spec.min, spec.max], ...
        spec.colormap, spec.nan_color);
    
    if isempty(accum)
        accum = zeros(size(ch_im));
    end
    
    accum = accum + ch_im;
end

% Clip any pixels exceeding range
accum(accum < 0) = 0;
accum(accum > 1) = 1;

% Push result to the output
im = accum;
time = slice_info.time;
end

function out = getChannels(mov_wrapper)
% This function retrieves the channel identifiers for the movie.

mov = mov_wrapper.mov;

if mov.isEmpty()
    out = 1;
    return;
end

out = mov.getChannelList();

end

function out = getPosRanges(mov_wrapper)
% This function retrieves the position ranges for the movie.

mov = mov_wrapper.mov;
pos_vec = mov_wrapper.pos_vec;

if isempty(pos_vec)
    if ~mov.isEmpty()
        out = [1; mov.getNumSlices()]; 
    else
        out = [1; 1];
    end
else
    pos_ranges = mov.getPositionRanges();
    out = pos_ranges{pos_vec};
end

end

function out = getHistInfo(mov_wrapper)
% This function retrieves information on the pixel
% intensity histgrams for the movie.

mov = mov_wrapper.mov;

if mov.isEmpty()
   out = [];
   out.min = 0;
   out.max = 1;
   out.dist = [0, 1; 0, 0];
   return; 
end

out = mov.getHistInfo();
end

function dset = processROIs(mov_wrapper, roi_list)
% This function processes the passed ROIs against the
% passed movie in order to produce a cell array of time
% series

mov = mov_wrapper.mov;
pos_vec = mov_wrapper.pos_vec;
dset = cell(1, length(roi_list));

if mov.isEmpty()
    for roi_idx=1:length(roi_list)
        dset{roi_idx} = [1; NaN];
    end
    
    return;
end

% Find the size of every mask that we will need
for roi_idx=1:length(roi_list)
    
    home_pos = roi_list(roi_idx).home_pos;
    scan_pos = roi_list(roi_idx).scan_pos;
    channel = roi_list(roi_idx).channel;
    
    if isempty(pos_vec)
        % TODO: this could probably be implemented
        dset{roi_idx} = [0; NaN];
        continue;
    end
        
    slice_idx = mov.getSliceIndex(pos_vec, home_pos);
    ch_data = mov.getChannelDataById(slice_idx, channel);
    
    if isempty(ch_data)
        dset{roi_idx} = [0; NaN];
        continue;
    end
    
    pos_ranges = mov.getPositionRanges();
    home_pos(1,scan_pos) = pos_ranges{pos_vec}(1,scan_pos);
    scan_extents = [scan_pos; pos_ranges{pos_vec}(2, scan_pos)];
    
    switch class(roi_list(roi_idx).handle)
        case 'imellipse'
            roi_vertices = roi_list(roi_idx).handle.getVertices;
        case 'impoly'
            roi_vertices = roi_list(roi_idx).handle.getPosition;
        otherwise
            error 'Invalid ROI handle';
    end
    
    mask = poly2mask(roi_vertices(:,1), roi_vertices(:,2), ...
        size(ch_data.image, 1), size(ch_data.image,2));
    
    mean_vals = mov.scanROI(pos_vec, home_pos, scan_extents, channel, mask);
    
    vec = [mean_vals(2,:); mean_vals(3,:)];

    dset{roi_idx} = vec;
end

end
