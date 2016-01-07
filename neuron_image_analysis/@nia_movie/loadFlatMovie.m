function loadFlatMovie(obj, fmovie, freq)
%LOADFLATMOVIE Read movie from 4D array containing movie data.
%   loadFlatMovie(folder, fmovie) populates the movie with data from the
%   passed 4D array. 
%
%   The function accepts the following arguments:
%
%       fmovie - A 4D array containing the movie data. The first and second
%           dimensions of the passed array are taken as with the height and
%           width of each frame. The third dimension is taken as time, and
%           the last dimension is taken as channel. Trailing dimensions may
%           be omitted, so a 3D array is interpreted as a movie with a
%           single channel, and a 2D array is interpreted as a movie with a
%           single frame. The position vector of the resulting movie is set
%           to the index of the frame in the time dimension.
%
%       freq - A scalar floating point value that specifies the frequency
%           of the frames. If this argumentm m  is omitted, then a value of one
%           is assumed.

if nargin < 3
    freq = 1.0;
end

% Check input arguments
if ~isfloat(fmovie) || ~isreal(fmovie) || ...
        ndims(fmovie) > 4 || isempty(fmovie)
    error 'The argument ''mov'' has an invalid type';
end

if ~isfloat(freq) || ~isreal(freq) || freq < 0
    error 'The argument ''freq'' has an invalid value';
end


% Destroy the movie
obj.slices = [];
obj.ch_list = [];
obj.pos_lens = [];
obj.pos_lu = {};
obj.pos_ranges = {};
obj.hist_info = [];

% Allocate the slices array
num_frames = size(fmovie, 3);
num_channels = size(fmovie, 4);

slices(num_frames).time = [];
slices(num_frames).channels = [];
slices(num_frames).pos = {};

% Traverse each frame
for slice_idx=1:num_frames
    slices(slice_idx).time = (slice_idx-1) * freq;

    % Allocate channel data array
    ch_data = [];
    ch_data(num_channels).ch = []; %#ok<AGROW>
    ch_data(num_channels).image = []; %#ok<AGROW>
    
    for ch_idx=1:num_channels
        ch_data(ch_idx).ch = ch_idx; %#ok<AGROW>
        ch_data(ch_idx).image = squeeze(fmovie(:,:,slice_idx,ch_idx)); %#ok<AGROW>
    end
    
    slices(slice_idx).channels = ch_data;
    slices(slice_idx).pos{1} = slice_idx;
end

% Create the channel list
ch_list = 1:num_channels;

% Create position lengths entry
pos_lens_entry = 1;

% Push results into movie object
obj.slices = slices;
obj.ch_list = ch_list;
obj.pos_lens = pos_lens_entry;

end
