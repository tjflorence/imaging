function stimulus = nia_createSparseNoise(width, height, params)
%NIA_CREATESPARSENOISE Create sparse noise stimulus
%   stimulus = nia_createSparseNoise(width, height, params) creates a movie
%   with a sparse noise stimulus. In this stimulus, a small, fixed number
%   of grid elements are increased or decreased in intensity while the
%   remainder of the grid maintains a constant luminosity. The position of
%   the first point is drawn at random without replacement until all frames
%   have been generated or until all grid points have been used.  Once all
%   grid points have been used, the random sampling is reset.  The other
%   active points within each frame are drawn at random for each frame with
%   no guarantees outside of not overlapping with each other. All intensity
%   values fall in the range between zero and one.
%
%   The function accepts the following arguments:
%
%       width - Width of the output images
%       height - Height of output images
%       params - Structure of parameters
%
%   The params structure must have the following parameters:
%
%       num_samples - Number of sampled configurations
%       num_active_points - Number of active points per frame
%       background - Background intensity
%       force_binary - True if active points should have fixed
%           values rather being allowed to have an arbitrary
%           value between zero and one.
%       binary_first - If force binary is true, the first active
%           point will be forced to this intensity value.
%       binary_high - If force_binary is true, half of the non-first
%           active points will be forced to this intensity value.
%       binary_low - If force_binary is true, half of the  non-first
%           active points will be forced to this intensity value.
%           Note that binary_low may be greater than or equal to
%           binary_high if desired.
%
%   The function returns:
%
%       stimulus - A 3D array of output images, with frames
%           in teh last dimension.
%
%   Example:
%
%       params.num_samples = 100;
%       params.num_active_points = 2;
%       params.background = 0.5;
%       params.force_binary = true;
%       params.binary_first = 1;
%       params.binary_high = 1;
%       params.binary_low = 1;
%       stimulus = nia_createSparseNoise(20, 20, params);

% Check inputs
if ~nia_isScalarInteger(width) || width <= 1
    error 'The argument ''width'' must be a scalar integer greater than one';
end

if ~nia_isScalarInteger(height) || height <= 1
    error 'The argument ''height'' must be a scalar integer greater than one';
end

if ~isstruct(params) || length(params) ~= 1
    error 'The argument ''params'' must be a structure of parameters';
end

params_allowed = { ...
    'num_samples', ...
    'num_active_points', ...
    'background', ...
    'force_binary', ...
    'binary_first', ...
    'binary_high', ...
    'binary_low'
    };
params_required = params_allowed;

[params_ok, params_msg] = nia_hasValidFieldNames(...
    params, params_allowed, params_required);
if ~params_ok
    error(params_msg, 'params');
end

if ~nia_isScalarInteger(params.num_samples) || params.num_samples < 1
    error 'The argument ''params.num_samples'' has an invalid value';
end

if ~nia_isScalarInteger(params.num_active_points) || params.num_active_points < 1
    error 'The argument ''params.active_points'' has an invalid value';
end

if ~isfloat(params.background) || length(params.background) ~= 1 || ...
        params.background < 0 || params.background > 1
    error 'The argument ''params.background'' has an invalid value';
end

if ~islogical(params.force_binary) || length(params.force_binary) ~= 1
    error 'The argument ''params.force_binary'' has an invalid value';
end

if ~isfloat(params.binary_first) || length(params.binary_first) ~= 1 || ...
        params.binary_first < 0 || params.binary_first > 1
    error 'The argument ''params.binary_first'' has an invalid value';
end

if ~isfloat(params.binary_high) || length(params.binary_high) ~= 1 || ...
        params.binary_high < 0 || params.binary_high > 1
    error 'The argument ''params.binary_high'' has an invalid value';
end

if ~isfloat(params.binary_low) || length(params.binary_low) ~= 1 || ...
        params.binary_low < 0 || params.binary_low > 1
    error 'The argument ''params.binary_low'' has an invalid value';
end

% Allocate memory for stimulus
stimulus = zeros(height, width, params.num_samples);

% Create list of first points
grid_nelem = width * height;
first_pts = zeros(1, params.num_samples);

for first_idx=1:floor(params.num_samples / grid_nelem)
    first_pts((first_idx-1)*grid_nelem+1:first_idx*grid_nelem) = ...
        randperm(grid_nelem);
end

trailing_num = mod(params.num_samples, grid_nelem);
if trailing_num > 0
    first_pts(params.num_samples-trailing_num+1:end) = ...
        randperm(grid_nelem, trailing_num);
end

% Traverse samples
for sample_idx=1:params.num_samples
    
    % Set background intensity
    frame = params.background * ones(height, width);
    
    % Determine active points index
    active_pts_list = zeros(1, params.num_active_points);
    active_pts_list(1) = first_pts(sample_idx);
    active_pts_count = 1;
    
    while active_pts_count < params.num_active_points
        active_pts_list(active_pts_count+1:end) = ...
            randperm(grid_nelem, params.num_active_points - ...
            active_pts_count);
        
        uniq_vals = unique(active_pts_list);
        active_pts_count = length(uniq_vals);
        active_pts_list(1:active_pts_count) = uniq_vals;
    end
    
    % Place each active point
    vals = rand(1, params.num_active_points);
    if params.force_binary
        mask = vals < 0.5;       
        vals(mask) = params.binary_low;
        vals(~mask) = params.binary_high;
        vals(1) = params.binary_first;
    end
    
    frame(active_pts_list) = vals;
    
    % Push frame out to 3D matrix
    stimulus(:,:,sample_idx) = frame;
end

end

