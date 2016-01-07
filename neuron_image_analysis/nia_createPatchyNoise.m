function [stimulus, density] = nia_createPatchyNoise(width, height, ...
    frames_per_epoch, params)
%NIA_CREATEPATCHYNOISE Create a dense noise stimulus
%   [stimulus, density] = nia_createPatchyNoise(width, height,
%   frames_per_epoch, params) creates a movie with a patchy white noise
%   stimulus. The patchy white noise stimulus uses patches produced by the
%   Voronoi diagram of a set of randomly positioned points. Frames are
%   generated in "epochs" in which all frames use the same Voronoi diagram,
%   although the intensity of patches are randomly assigned between frames.
%   All intensity values fall in the range between zero and one.
%
%   The function accepts the following arguments:
%
%       width - Width of the output images
%       height - Height of output images
%       frames_per_epoch - Number of frames in each epoch
%       params - Structure of parameters
%
%   The params structure must have the following parameters:
%
%       space_density - Array of point densities with one
%           element per epoch
%
%   The function returns:
%
%       stimulus - A 3D array of output images, with frames
%           in the last dimensions
%       density - A row vector of the point density used in
%           each epoch, given in the order in which each 
%           epoch appears in the output stimulus movie
%
%   Example:
%
%       params.space_density = [0.1*ones(1,10), 0.01*ones(1,10)];
%       stimulus = nia_createPatchyNoise(20, 20, 2, params);

large_number = 100;

% Check input arguments
if ~nia_isScalarInteger(width) || width <= 1
    error 'The argument ''width'' must be a scalar integer greater than one';
end

if ~nia_isScalarInteger(height) || height <= 1
    error 'The argument ''height'' must be a scalar integer greater than one';
end

if ~nia_isScalarInteger(frames_per_epoch) || frames_per_epoch <= 0
    error 'The argument ''duration'' must be a scalar integer greater than zero';
end

if ~isstruct(params) || length(params) ~= 1
    error 'The argument ''params'' must be a structure of parameters';
end

params_allowed = {'space_density'};
params_required = params_allowed;
[params_ok, params_msg] = nia_hasValidFieldNames(...
    params, params_allowed, params_required);
if ~params_ok
    error(params_msg, 'params');
end

if ~isfloat(params.space_density) || ~isreal(params.space_density) || ...
        ~isvector(params.space_density) || isempty(params.space_density)
    error 'The argument ''params.space_density'' has an invalid type';
end

% Allocate memory for stimulus
stimulus = zeros(height, width, ...
    frames_per_epoch * length(params.space_density));

% Randomize the epochs
space_density_idx = randperm(length(params.space_density));

% Iterate through each of the epochs
for epoch_idx=1:length(params.space_density)
    
    % Find points for Voronoi diagram
    cur_space_density = params.space_density(...
        space_density_idx(epoch_idx));
    
    points_num = ceil(width * height * cur_space_density);
    points_x = (width-1)*rand(points_num,1) + 1;
    points_y = (height-1)*rand(points_num,1) + 1;
    
    % Add some dummy points far from actual image region.
    % These points are defined such that they belong to
    % the convex hull of the point set provided that
    % large_number is greater than one. Thus, they are
    % at the center of the unbounded voronoi cells.  Later
    % on, we discard any unbounded voronoi cells, since we
    % know those cannot be legitimate cells.
    points_x_padded = zeros(points_num+4, 1);
    points_y_padded = zeros(points_num+4, 1);
    
    points_x_padded(1:points_num) = points_x;
    points_x_padded(points_num+1) = -large_number*width;
    points_x_padded(points_num+2) = -large_number*width;
    points_x_padded(points_num+3) = +large_number*width;
    points_x_padded(points_num+4) = +large_number*width;
    
    points_y_padded(1:points_num) = points_y;
    points_y_padded(points_num+1) = -large_number*height;
    points_y_padded(points_num+2) = +large_number*height;
    points_y_padded(points_num+3) = +large_number*height;
    points_y_padded(points_num+4) = -large_number*height;
    
    % Construct the Voronoi diagram
    del_tess = delaunayTriangulation(points_x_padded, points_y_padded);
    [vor_vertices, vor_regions] = voronoiDiagram(del_tess);
    
    % Iterate through each frame of epoch
    for frame_idx=1:frames_per_epoch
        frame = zeros(height, width);
        region_vals = rand(length(vor_regions), 1);
        
        % Iterate through each region
        for region_idx=1:length(vor_regions)
            region_x = vor_vertices(vor_regions{region_idx}, 1);
            region_y = vor_vertices(vor_regions{region_idx}, 2);
            
            % Skip over bounding points
            if nnz(isinf(region_x)) > 0
                continue;
            end
            
            region_mask = poly2mask(region_x, region_y, height, width);
            frame(region_mask) = region_vals(region_idx);
            
            global_frame_idx = (epoch_idx-1)*frames_per_epoch + frame_idx;
        end
        
        stimulus(:,:,global_frame_idx) = frame;
    end
end

density = params.space_density(space_density_idx);
    
end
