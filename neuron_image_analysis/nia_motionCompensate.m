function found_displ = nia_motionCompensate(mov, pos_vec, home_pos, ...
    pos_idx, channel_id, window)
%NIA_MOTIONCOMPENSATE Compensate for image motion
%   displ = nia_motionCompensate(mov, pos_vec, home_pos, pos_idx,
%   channel_id) calculates the apparent image motion, and then transforms
%   the image to eliminate that motion. It returns the calculated
%   displacement in pixels. The movie is modified by this function, and may
%   be left in a half-transformed state if this function is interrupted.
%
%   The motion compensation algorithm must traverse the movie in a time
%   ordered sequence, but because of the user-defined nature of the
%   position vector it is invalid to assume that the slice indices or a
%   specific position vector would be time-ordered. Instead, the user must
%   specify a "home" position that specifies the position vector to use for
%   the sequence. The argument pos_idx specifies the index into the
%   position vector that should be iterated in order to produce the
%   sequence, it is varied between the lowest and highest values found in
%   the movie.
%
%   This function takes the following arguments:
%
%       mov - Movie to transform
%
%       pos_vec - A scalar integer that contains the index of the position 
%           vector to use for the scan.
%
%       home_pos - A [1xN] array where N is the number of position
%           dimensions, that specifies the 'home' position for the motion
%           compensation.
%
%       pos_idx - A scalar integer that specifies the index of the position
%           vector to be iterated.
%
%       channel_id - A scalar integer that specifies the identifier of the
%           channel to use for performing motion estimation.
%
%       window - An optional argument that specifies the window to be used
%           for performing image registration. Must be of the form
%           [first_row, last_row, first_column, last_column].  If this
%           argument is set to the empty array, then the entire image is
%           used.
%
%   This function returns the following values:
%
%       found_displ - A [3xN] array where N is the number of found
%           positions that specifies the calculated displacement. The first
%           two rows specify the x and y displacement in pixels, and the
%           third row specifies the angular rotation in radians. All
%           displacements are given as the displacement that must be
%           applied to the motion compensated image in order to reproduce
%           the original image, and the angular rotation is given with
%           positive values indicating counter-clockwise rotation about a
%           point in the center of the image.

% TODO: 
%   1. what happens if home_pos does not have all frames in pos_ranges?
%   2. add argument checking
%   3. expose different algorithms

% Configuration variables
cfg.max_disp = 10;
cfg.max_t_disp = 200e-3;
cfg.win_size = 32 + 2*cfg.max_disp;
cfg.win_overlap = 0.25;
cfg.interp_xsize = 5;
cfg.interp_ysize = 5;
cfg.bump_frac = 0.01;
cfg.quiet = false;
cfg.prog_bar_len = 40;
cfg.enable_multipass = true;
cfg.algorithm = 'trans'; % 'optim', 'trans', or 'grid'

% Provide default values
if nargin < 6
    window = [];
end

% Find the positions values
pos_ranges = mov.getPositionRanges;
start_pos = pos_ranges{pos_vec}(1,pos_idx);
last_pos = pos_ranges{pos_vec}(2,pos_idx);
num_pos = last_pos - start_pos + 1;

prev_displ = zeros(3, num_pos);
cur_displ = zeros(3, num_pos);

ref_width_expn = 3;
warped_width = 1;

% Pick out reference widths to use
if cfg.enable_multipass
    min_ref_width = 2*warped_width; % always 2X warped_width
    min_ref_width = min([min_ref_width, num_pos]);

    min_ref_width_power = ceil(log(min_ref_width)/log(ref_width_expn));
    max_ref_width_power = ceil(log(num_pos)/log(ref_width_expn));
    ref_widths_list = ref_width_expn.^(...
        min_ref_width_power:max_ref_width_power);
else
    ref_widths_list = num_pos;
end

% Pick out progress bar positions
prog_marker_locs = round(linspace(0, num_pos, cfg.prog_bar_len));
prog_marker_hdr = 'Pass %d of %d:  ';
prog_marker_hdr_len = [];

% Pick out the registration function
switch cfg.algorithm
    case 'optim'
        registerImages = @registerImagesOptim;
    case 'trans'
        registerImages = @registerImagesTrans;
    case 'grid'
        registerImages = @registerImagesGrid;
end

if ~cfg.quiet
    fprintf(1, '\nMotion compensation:\n');
end

% Traverse the list of reference image averaging window widths
for ref_width_idx=1:length(ref_widths_list)
    ref_width = ref_widths_list(ref_width_idx);
    ref_sum_im = [];
    warp_sum_im = [];

    if ~cfg.quiet
        if ref_width_idx == 1
            prog_marker_hdr_len = length(sprintf(prog_marker_hdr, ...
                length(ref_widths_list), length(ref_widths_list)));
            
            prog_str = 'Progress';
            
            % Note that in calculation below we add 2 to the length
            % of the progress bar to include the '|' end caps
            fprintf(1, repmat(' ', 1, round(prog_marker_hdr_len + ...
                (cfg.prog_bar_len+2)/2 - length(prog_str)/2)));
            fprintf(1, '%s\n', prog_str);
            
            fprintf(1, repmat(' ', 1, prog_marker_hdr_len));
            fprintf(1, '|');
            fprintf(1, repmat('-', 1, cfg.prog_bar_len));
            fprintf(1, '|\n');
        end
        
        count = fprintf(1, prog_marker_hdr, ref_width_idx, ...
            length(ref_widths_list));
        
        if count < prog_marker_hdr_len
            fprintf(1, repmat(' ', 1, prog_marker_hdr_len - count));
        end
        
        fprintf(1, '|');
    end
    
    prog_marker_pos = 1;
    
    for pos_val=start_pos:last_pos   
        if ~cfg.quiet
            while pos_val > prog_marker_locs(prog_marker_pos)
                fprintf(1, '-');
                prog_marker_pos = prog_marker_pos + 1;
            end
        end

        % Build the reference image
        if isempty(ref_sum_im)
            [ref_sum_im, ref_start_pos, ref_last_pos] = ...
                initRunningSum(mov, pos_vec, home_pos, pos_idx, ...
                pos_ranges, pos_val, ref_width, channel_id, ...
                -prev_displ, - start_pos + 1);
        else
            [ref_sum_im, ref_start_pos, ref_last_pos] = ...
                updateRunningSum(mov, pos_vec, home_pos, pos_idx, ...
                pos_ranges, pos_val, ref_width, channel_id, ...
                ref_sum_im, -prev_displ, - start_pos + 1);
        end

        ref_im = ref_sum_im ./ (ref_last_pos - ref_start_pos + 1);

        % Build the warped image
        if isempty(warp_sum_im)
            [warp_sum_im, warp_start_pos, warp_last_pos] = ...
                initRunningSum(mov, pos_vec, home_pos, pos_idx, ...
                pos_ranges, pos_val, warped_width, channel_id, ...
                0*prev_displ, - start_pos + 1);
        else
            [warp_sum_im, warp_start_pos, warp_last_pos] = ...
                updateRunningSum(mov, pos_vec, home_pos, pos_idx, ...
                pos_ranges, pos_val, warped_width, channel_id, ...
                warp_sum_im, 0*prev_displ, - start_pos + 1);
        end

        warp_im = warp_sum_im ./ (warp_last_pos - warp_start_pos + 1);

        % Register the two images
        init_displ = [0.1; 0.1; 20e-3];
        
        ref_im = log(ref_im + 1);
        warp_im = log(warp_im + 1);

        if ~isempty(window)
            win_ref_im = ref_im(window(1):window(2), window(3):window(4));
            win_warp_im = warp_im(window(1):window(2), window(3):window(4));
        else
            win_ref_im = ref_im;
            win_warp_im = warp_im;
        end
        
        cur_displ(:, pos_val-start_pos+1) = registerImages(...
            win_ref_im, win_warp_im, init_displ, cfg);
    end
    
    if ~cfg.quiet
        while prog_marker_pos <= cfg.prog_bar_len
            fprintf(1, '-');
            prog_marker_pos = prog_marker_pos + 1;
        end
        
        fprintf(1, '|\n');
    end
    
    
    prev_displ = bsxfun(@plus, cur_displ, -mean(cur_displ, 2));
    cur_displ = 0*cur_displ;
end

found_displ = prev_displ;
clear prev_displ cur_displ;

% Shift the image to compensate for the found motion
for pos_val=start_pos:last_pos
    cur_pos = home_pos;
    cur_pos(pos_idx) = pos_val;
    
    slice_idx = mov.getSliceIndex(pos_vec, cur_pos);
    slice_info = mov.getSliceInfo(slice_idx);
    
    slice_displ = found_displ(:, pos_val - start_pos + 1);
    
    for ch_idx=1:slice_info.num_ch
        cur_ch_data = mov.getChannelData(slice_idx, ch_idx);
        new_ch_data.image = displaceImage(cur_ch_data.image, ...
            -slice_displ, 'cubic');
        
        mov.setChannelData(slice_idx, ch_idx, new_ch_data);
    end
end

end

function [sum_im,first_pos,last_pos] = initRunningSum(mov, pos_vec, ...
    home_pos, pos_idx, pos_ranges, pos_val, width, channel_id, displ, ...
    displ_idx_offset)
% This function initializes a running windowed average of images by
% subtracting off the frame that passed out of the window and adding in the
% frame that moved into the window. The mov argument is the movie, the
% home_pos argument is the uniterated position vector, the pos_idx is the
% index into the position vector to iterate, the pos_val is the position
% value for the current center of the moving window, the width is the width
% of the window, the channel_id is the channel identifier, the displ is a
% matrix with displacements that should be applied to each frame prior to
% performing the average, and the displ_idx_offset is a value that should
% be added to the position value in order to obtain the corresponding index
% into the displacement vector.

    first_pos = pos_val - width;
    first_pos = max([first_pos, pos_ranges{pos_vec}(1, pos_idx)]);
        
    last_pos = pos_val + width;
    last_pos = min([last_pos, pos_ranges{pos_vec}(2, pos_idx)]);
    
    sum_im = getSumImage(mov, pos_vec, home_pos, pos_idx, first_pos, ...
        last_pos, channel_id, displ, displ_idx_offset);
end

function total = getSumImage(mov, pos_vec, home_pos, pos_idx, ...
    first_pos, last_pos, channel_id, displ, displ_idx_offset)
% This function finds the summed image for the frames between a 
% range of positions

    home_pos(pos_idx) = first_pos;
    scan_extents = [pos_idx; last_pos];
    
    total = 0;
    
    total = mov.scan(pos_vec, home_pos, scan_extents, channel_id, ...
        @(im,info,udata) getSum(im, info, udata, pos_vec, pos_idx, ...
        displ, displ_idx_offset), total);
end

function udata = getSum(im, info, udata, pos_vec, pos_idx, displ, displ_idx_offset)
% This function is a helper function used by getSumImage to sum
% multiple frames.

    displ_val = displ(:, info.pos{pos_vec}(pos_idx) + displ_idx_offset);

    im = displaceImage(im, displ_val, 'linear');

    udata = udata + im;
end

function [sum_im,first_pos,last_pos] = updateRunningSum(mov, pos_vec, ...
    home_pos, pos_idx, pos_ranges, pos_val, width, channel_id, sum_im, ...
    displ, displ_idx_offset)
% This function updates a running windowed average of images by subtracting
% off the frame that passed out of the window and adding in the frame that
% moved into the window. The mov argument is the movie, the pos_vec
% argument is the index of the position vector, the home_pos argument is
% the uniterated position vector, the pos_idx is the index into the
% position vector to iterate, the pos_val is the position value for the
% current center of the moving window, the width is the width of the
% window, the channel_id is the channel identifier, the sum_im is the
% summed image value from the previous frame, the displ is a matrix with
% displacements that should be applied to each frame prior to performing
% the average, and the displ_idx_offset is a value that should be added to
% the position value in order to obtain the corresponding index into the
% displacement vector.

    first_pos = pos_val - width;
    last_pos = pos_val + width;
    
    if first_pos > pos_ranges{pos_vec}(1, pos_idx)
        % Calculate image that just moved out of the current running
        % average window, so that we can subtract it off sum
        
        retired_pos = home_pos;
        retired_pos(pos_idx) = first_pos - 1;
        slice_idx = mov.getSliceIndex(pos_vec, retired_pos);
        if slice_idx == 0
            error 'Missing frame in given position range';
        end
        
        ch_data = mov.getChannelDataById(slice_idx, channel_id);
        retired_idx = first_pos - 1 + displ_idx_offset;
        
        retired_im = displaceImage(ch_data.image, ...
            displ(:, retired_idx), 'linear');
        
        sum_im = sum_im - retired_im;
    end
    
    if last_pos <= pos_ranges{pos_vec}(2, pos_idx)
        % Calculate the image that just moved into the current runnning
        % average window, so that we can add it to the sum
        
        new_pos = home_pos;
        new_pos(pos_idx) = last_pos;
        slice_idx = mov.getSliceIndex(pos_vec, new_pos);
        if slice_idx == 0
            error 'Missing frame in given position range';
        end
        
        ch_data = mov.getChannelDataById(slice_idx, channel_id);
        new_idx = last_pos + displ_idx_offset;
        
        new_im = displaceImage(ch_data.image, ...
            displ(:, new_idx), 'linear');
        
        sum_im = sum_im + new_im;
    end
    
    first_pos = max([first_pos, pos_ranges{pos_vec}(1, pos_idx)]);
    last_pos = min([last_pos, pos_ranges{pos_vec}(2, pos_idx)]);
end

function im_trans = displaceImage(im, displ, method)
% This function displaces the passed image using the
% given displacement (in [dx; dy; dtheta] form). The 
% method must be 'nearest', 'linear', or 'cubic.

cfg_use_imwarp = false;

% Note: the easiest way to perform this transformation is using the
% imwarp function. As is typical, matlab has implemented this function
% in a way that is stupid and slow.  As of 2014b, this function runs
% 3X slower than a naive implementation using griddedInterpolant.
% In the event that mathworks eventually fixes their implementation,
% the imwarp function should be used and the naive implementation should
% be removed.

if cfg_use_imwarp
    % Create transformation matrix
    A = [cos(displ(3)),  -sin(displ(3)),    0; ...
        sin(displ(3)),   cos(displ(3)),    0; ...
        displ(1),        displ(2),    1 ];
    
    tform = affine2d(A);
    
    % Create frame of reference
    rref = imref2d(...
        size(im), ...
        [-size(im,2)/2, size(im,2)/2], ...
        [-size(im,1)/2, size(im,1)/2]);
    
    % Perform image transformation
    im_trans = imwarp(im, rref, tform, method, ...
        'OutputView', rref, ...
        'FillValues', -1);
else
    % Create transformation matrix
    A = [cos(-displ(3)),  -sin(-displ(3)),    0; ...
        sin(-displ(3)),   cos(-displ(3)),    0; ...
        -displ(1),        -displ(2),    1 ];
    
    % Pad the image
    pad_size = 1;
    im_padded = padarray(im, [pad_size, pad_size], 'symmetric');
    
    % Create coordinates for padded image
    x_coords_padded = linspace(-size(im_padded,2)/2+0.5, ...
        size(im_padded,2)/2-0.5, size(im_padded,2));
    y_coords_padded = linspace(-size(im_padded,1)/2+0.5, ...
        size(im_padded,1)/2-0.5, size(im_padded,1));
    
    % Create an interpolator function
    im_interpolator = griddedInterpolant(...
        {y_coords_padded, x_coords_padded}, ...
        im_padded, method);
    
    % Create coordinates for original image
    x_coords = linspace(-size(im,2)/2+0.5, size(im,2)/2-0.5, size(im,2));
    y_coords = linspace(-size(im,1)/2+0.5, size(im,1)/2-0.5, size(im,1));
    
    % Transform the coordinates
    [Y,X] = ndgrid(y_coords, x_coords);
    R = [X(:), Y(:), ones(numel(X),1)];
    R = R * A;
    X = reshape(R(:,1), size(im,1), size(im,2));
    Y = reshape(R(:,2), size(im,1), size(im,2));
    
    % Perform interpolation
    im_trans = im_interpolator(Y, X);
    
    % Mask out extrapolated pixels. Note that the interpolation
    % grid is based on pixel centers, but that we allow extrapolation
    % to the pixel edges.  This prevents needless filling at the
    % edges of the image.
    bnds_mask = X < -size(im,2)/2 | X > size(im,2)/2 | ...
        Y < -size(im,1)/2 | Y > size(im,1)/2;
    im_trans(bnds_mask) = -1;
end

end

function optim_displ = registerImagesOptim(im_fixed, im_moving, ...
    init_displ, cfg)
% This function returns the displacement in the x direction,
% displacement in the y direction, and the displacement in
% the angular direction that maximum the alignment between
% the fixed and moving images.

% Set up nonlinear constraints: A * x < b
Acon = [ ...
    1,  0,  0; ...
   -1,  0,  0; ...
    0,  1,  0; ...
    0, -1,  0; ...
    0,  0,  1; ...
    0,  0, -1 ];
bcon = [ ...
    cfg.max_disp;
    cfg.max_disp;
    cfg.max_disp;
    cfg.max_disp;
    cfg.max_t_disp;
    cfg.max_t_disp ];

options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'TolX', 0.0001, ...
    'UseParallel', false);

optim_x = fmincon(...
    @(x) findError(im_fixed, im_moving, x), ...
    init_displ, Acon, bcon, [], [], [], [], [], options);

optim_displ = optim_x;
end

function err = findError(im_fixed, im_moving, displ)
% This function finds the squared difference of the two images after
% the given displacement is applied (see registerImages for definition
% of arguments).

im_trans = displaceImage(im_moving, displ, 'linear');
mask = im_trans ~= -1;

err = mean((im_fixed(mask) - im_trans(mask)).^2);
end

function optim_displ = registerImagesTrans(im_fixed, im_moving, ~, cfg)
% This function returns the displacement in the x direction,
% displacement in the y direction, and the displacement in
% the angular direction that maximum the alignment between
% the fixed and moving images.

optim_displ = zeros(3, 1);
optim_displ(1:2) = findOffset(im_fixed, im_moving, cfg);

end

function optim_displ = registerImagesGrid(im_fixed, im_moving, ~, cfg)
% This function returns the displacement in the x direction,
% displacement in the y direction, and the displacement in
% the angular direction that maximum the alignment between
% the fixed and moving images.

% Calculate the displacement vector for a grid of windows that
% are distributed evenly across the image

xpos = calculateDisplPositions(size(im_fixed,2), cfg);
ypos = calculateDisplPositions(size(im_fixed,1), cfg);

displ = zeros(2, length(xpos), length(ypos));
weight = zeros(length(xpos), length(ypos));
win_center = zeros(2, length(xpos), length(ypos));

for x_idx=1:length(xpos)
    for y_idx=1:length(ypos)
        % Calculate the displacement that maximizes the
        % cross-correlation
        im_fixed_win = im_fixed(...
            ypos(y_idx):ypos(y_idx) + cfg.win_size - 1, ...
            xpos(x_idx):xpos(x_idx) + cfg.win_size - 1);
        
        im_moving_win = im_moving(...
            ypos(y_idx):ypos(y_idx) + cfg.win_size - 1, ...
            xpos(x_idx):xpos(x_idx) + cfg.win_size - 1);
        
        [displ(:,x_idx,y_idx), weight(x_idx, y_idx)] = ...
            findOffset(im_fixed_win, im_moving_win, cfg);

        % Calculate the center of the window position
        xpos_center = xpos(x_idx) + cfg.win_size / 2;
        ypos_center = ypos(y_idx) + cfg.win_size / 2;
        
        % Calculate the center position in coordinates in which the upper
        % left edge of image has position (-width/2, -height/2) and the
        % lower right edge of image has position (width/2, height/2). This
        % is the coordinate system used by the displaceImage() function
        % above.
        xpos_center = xpos_center - size(im_fixed,2)/2 + 0.5;
        ypos_center = ypos_center - size(im_fixed,1)/2 + 0.5;
        
        win_center(:,x_idx,y_idx) = [xpos_center; ypos_center];
    end
end

% Fit the found displacements to a rotation where
%
% displ(1) = y_position * rotation
% displ(2) = - x_position * rotation
%
% Perform fit by creating a predicted value that is the 
% concatentation of the x and y displacement vectors, and the
% predictor values are the concatenation of (-y_position)
% and (x_position) vectors. We use the correlation value to
% weight the fit, so that empty windows are ignored.

% Calculate the predicted value
fit_dx = displ(1,:,:);
fit_dy = displ(2,:,:);
fit_d = [fit_dx(:); fit_dy(:)];

% Calculate the predictor variable
fit_jx = [ones(length(fit_dx(:)),1); zeros(length(fit_dx(:)),1)];
fit_jy = [zeros(length(fit_dx(:)),1); ones(length(fit_dx(:)),1)];

fit_ix = win_center(2,:,:);
fit_iy = - win_center(1,:,:);
fit_i = [fit_ix(:); fit_iy(:)];

fit_ii = [fit_jx'; fit_jy'; fit_i']';

% Calculating the weighting variable
fit_w = weight(:,:);
fit_w = [fit_w(:); fit_w(:)];

% Perform least squares fit
optim_displ = lscov(fit_ii, fit_d, fit_w);

end

function win_pos = calculateDisplPositions(dim_size, cfg)
% This function calculates the start position for each window

stride = round(cfg.win_overlap * cfg.win_size);

% total_length = (num_windows-1) * stride + cfg.win_size
% the equation below solves this for total_length == dim_size
num_windows = floor((dim_size - cfg.win_size) / stride) + 1;

total_length = (num_windows-1) * stride + cfg.win_size;
offset = floor((dim_size - total_length)/2) + 1;

win_pos = (0:num_windows-1)*stride + offset;

end

function [displ, weight] = findOffset(im_fixed, im_moving, cfg)
% This function finds the displacement of im_moving provides that best
% match to im_fixed.  Matching is performed by Fourier-based
% cross-correlation of the two images followed by Gaussian sub-pixel
% interpolation.

    % Detrend the image data
    im_moving = detrendImage(im_moving);
    im_fixed = detrendImage(im_fixed);
    
    % Calculate fft sizes
    fft_xsize = 2^ceil(log2(size(im_moving,2)));
    fft_ysize = 2^ceil(log2(size(im_moving,1)));
    
    % Zero all of the pixels on edge of im1, so that correlation
    % exposes the same number of pixels regardless of displacement.
    % This is necessary to prevent zero displacement bias.
    im_moving_wind = zeros(size(im_moving));
    im_moving_wind(...
        cfg.max_disp+1:end-cfg.max_disp, cfg.max_disp+1:end-cfg.max_disp) = ...
        im_moving(cfg.max_disp+1:end-cfg.max_disp, cfg.max_disp+1:end-cfg.max_disp);

    % Take the Fourier transform of im_moving_wind and im2, and use
    % these to calculate the cross-correlation.
    im_moving_ft = fft2(im_moving_wind, fft_ysize, fft_xsize);
    im_fixed_ft = fft2(im_fixed, fft_ysize, fft_xsize);
    imX_ft = im_moving_ft .* conj(im_fixed_ft);
    imX = ifft2(imX_ft);
    
    % Rearrange the cross-correlation matrix so that rows hold
    % increasing y displacements, and the columns hold increasing
    % x displacements.
    imX = fftshift(imX);
    
    % Zero points outside of the area of interest
    [dispx,dispy] = meshgrid(-fft_xsize/2:1:fft_xsize/2-1, ...
      -fft_ysize/2:1:fft_ysize/2-1);
    imX(abs(dispx) > cfg.max_disp | abs(dispy) > cfg.max_disp) = 0.0;
    
    % Identify the largest point in the cross-correlation
    [max_i, max_j] = find(imX == max(imX(:)), 1, 'first');
    
    % Convert the row/column values to a displacement
    intg_disp_x = max_j - (fft_xsize/2+1);
    intg_disp_y = max_i - (fft_ysize/2+1);
    
    % If the displacement exceeds maximum, then return immediately
    if abs(intg_disp_x) > cfg.max_disp || abs(intg_disp_y) > cfg.max_disp
        displ(1) = min(max(intg_disp_x, -cfg.max_disp), cfg.max_disp);
        displ(2) = min(max(intg_disp_y, -cfg.max_disp), cfg.max_disp);
        weight = 0;
        return;
    end
    
    % Sample cross-correlation around the center point
    imX_center = ...
        imX(max_i-(cfg.interp_ysize-1)/2 : max_i+(cfg.interp_ysize-1)/2, ...
        max_j-(cfg.interp_xsize-1)/2 : max_j+(cfg.interp_xsize-1)/2);
    
    % If a cross-correlation value is exactly zero, then the logarithm
    % used for interpolation is ill-defined.  So we check for any values
    % less than bump_frac times the maximum value, and then reset these
    % to bump_frac times the maximum value.  Provided that bump_frac is
    % much less than one, then this has a negligible effect on the 
    % interpolation but stabilizes the calculation.
    imX_center_max = imX(max_i, max_j);
    imX_center(imX_center < cfg.bump_frac * imX_center_max) = ...
        cfg.bump_frac * imX_center_max;
    
    % Take the logarithm of the correlation, provided that the 
    % correlation is Gaussian, imX_center_log is paraboloidal.
    imX_center_log = log(imX_center);
    
    % Calculate the coordinates associated with the values in the
    % samples cross-correlation map
    x_coords = linspace(...
        -(cfg.interp_xsize-1)/2, (cfg.interp_xsize-1)/2, cfg.interp_xsize);
    y_coords = linspace(...
        -(cfg.interp_ysize-1)/2, (cfg.interp_ysize-1)/2, cfg.interp_ysize);
    
    [grid_x_coords, grid_y_coords] = meshgrid(x_coords, y_coords);
    
    % Fit the sampled cross-correlation map to a paraboloid using
    % the linear regression routine.
    A = [ones(numel(imX_center_log), 1), ...
        grid_x_coords(:), grid_y_coords(:), ...
        grid_x_coords(:).^2, grid_y_coords(:).^2, ...
        grid_x_coords(:) .* grid_y_coords(:)];
    
    p = lscov(A, imX_center_log(:));
    
    % Calculate the center point of the paraboloid
    max_x = (p(3)*p(6) - 2*p(2)*p(5))/(4*p(4)*p(5)-p(6)^2);
    max_y = (p(2)*p(6) - 2*p(3)*p(4))/(4*p(4)*p(5)-p(6)^2);
    
    % Adjust the integral displacement by interpolated value
    displ(1) = intg_disp_x + max_x;
    displ(2) = intg_disp_y + max_y;
    weight = imX_center_max;
    
    % If the displacement exceeds maximum, then trim it back
    if abs(displ(1)) > cfg.max_disp || abs(displ(2)) > cfg.max_disp
        displ(1) = min(max(displ(1), -cfg.max_disp), cfg.max_disp);
        displ(2) = min(max(displ(2), -cfg.max_disp), cfg.max_disp);
        weight = 0;
        return;
    end
end

function im = detrendImage(im)
% This function removes any gradients from the passed image.
    
    b_x = polyfit(1:size(im,2), mean(im,1), 1);
    b_y = polyfit(1:size(im,1), mean(im,2)', 1);
    
    [x,y] = meshgrid(1:size(im,2), 1:size(im,1));
  
    im = im - b_x(1)*x - b_y(1)*y;

    im = im - mean(im(:));
end
