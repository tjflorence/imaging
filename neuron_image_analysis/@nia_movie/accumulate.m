function [accum, count] = accumulate(obj, rnge, transform)
%ACCUMULATE Accumulate slices from a movie
%   [accum, count] = accumulate(rnge, transform) traverses the passed movie,
%   selects slices that meet the passed criteria, applies a transform to
%   each image, and then sums the result. This function is highly general,
%   and may be used to perform operations such as finding the mean value
%   for each pixel of a movie, or finding the standard deviaion for each
%   pixel of a movie.
%
%   This function takes the following arguments:
%
%       rnge - The range of slices and channels to be used for the
%           accumulate operation. See isValidRange() for more information
%           about the format of this argument.
%
%       transform - This argument specifies the transform that should be
%           performed prior to performing the accumulate operation. If
%           this argument is a real, scalar floating-point value, then the 
%           accumulate operation is performed on an image in which each
%           pixels has been transformed as pixel_value^(transform). If
%           this argument is not a scalar floating-point value, then it
%           must be a function handle.  In which case, this function is
%           invoked as fcn(slice_image, info) for each of the selected
%           slices and the returned value is accumulated. The function is
%           passed the slice image, and a structure with the fields
%           time, pos, and ch describing the time, position, and channel
%           of the slice, respectively. The type and dimensions of the
%           returned value must be identical for each of the selected
%           images.
%
%   This function returns the following values:
%
%       accum - The mean value of the transformed values
%       count - The number of images that were summed.
%
%   Examples:
% 
%       To find the mean value of all slices in channel 1
%       with a time value between 0 and 1:
%
%           rnge.time = [0; 1];
%           rnge.channels = 1;
%           mean_val = mov.accumulate(rnge, 1);
%
%       To find the standard deviation of all slices in channel 2:
%       
%           rnge.channels = 2;
%           mean_val = mov.ccumulate(rnge, 1);
%           sq_val = mov.accumulate(rnge, 2);
%           std_val = sqrt(sq_val) - mean_val.^2;
%
%       To find the mean of the log-transformed slices where the
%       first element the first position vector has a value
%       of two and the channel identifier is equal to one:
%
%           rnge.pos_vec = 1;
%           rnge.pos_lims = [1; 2];
%           rnge.channels = 1;
%           log_val = mov.accumulate(rnge, @(im,info) log(im));

% Check input arguments
[rnge_ok, rnge_msg] = obj.isValidRange(rnge);
if ~rnge_ok
    error(rnge_msg, 'rnge');
end

if ~isfloat(transform) || ~isreal(transform) || ~isscalar(transform)
    if ~isa(transform, 'function_handle')
        error 'The argument ''transform'' has an invalid type';
    end
end

% Normalize the input range
rnge = obj.normalizeRange(rnge);

% Traverse each slice
has_first = false;
count = 0;
accum = [];

for slice_idx=1:length(obj.slices)
    
    % Skip non-matching slices
    if ~obj.inSliceRange(slice_idx, rnge)
        continue;
    end
    
    % Traverse each channel in slice
    for ch_idx=1:length(obj.slices(slice_idx).channels)
 
        % Skip non-matching channels
        if ~obj.inChannelRange(slice_idx, ch_idx, rnge)
            continue;
        end
        
        % Everything matched. Let's transform that image.
        im = obj.slices(slice_idx).channels(ch_idx).image;
        
        if isa(transform, 'function_handle')
            info.time = obj.slices(slice_idx).time;
            info.pos = obj.slices(slice_idx).pos;
            info.ch = obj.slices(slice_idx).channels(ch_idx).ch;
            
            trans_im = transform(im, info);
            if ~isnumeric(trans_im)
                error 'The transformed value must be numeric';
            end
            
        elseif transform == 1.0
            trans_im = im;
        elseif transform == 2.0
            trans_im = im .* im;
        else
            trans_im = im.^transform;
        end
        
        count = count + 1;
        
        % Check if the returned value has a valid type, if so
        % then add it up
        if has_first
            if ~strcmp(class(accum), class(trans_im)) || ...
               ndims(accum) ~= ndims(trans_im) || ...
               nnz(size(accum) - size(trans_im)) > 0
                error 'The transformed value has an incompatible type';
            end
            
            accum = accum + trans_im;
        else
            accum = trans_im;
        end
        
        has_first = true;
    end
end

if count > 0
    accum = accum ./ count;
end

end

