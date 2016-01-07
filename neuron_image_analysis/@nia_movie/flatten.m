function [flat, info] = flatten(obj, rnge, transform)
%FLATTEN Flatten movie into an N-dimensional array
%   [flat, info] = flatten(rnge, transform) traverses the passed movie,
%   selects slices and channels in the passed range, applies a transform to
%   each image, and concatenates the result. This function is highly
%   general, and may be used to perform operations such as converting a
%   movie to a 3D array, or finding the mean value with an ROI for each
%   time point.
%
%   This function takes the following arguments:
%
%       rnge - The range of slices and channels to be used for the flatten
%           operation. See isValidRange() for more information about the
%           format of this argument.
%
%       transform - This argument specifies the transform that should be
%           performed prior to performing the concatenate operation. If
%           this argument is empty, then each image is simply concatenated.
%           Otherwise, this argument must be a function handle, and is
%           invoked as fcn(slice_image, info) for each of the selected
%           slices and the returned value is concatenated. The function
%           is passed the slice image, and a structure with the fields
%           time, pos, and ch describing the time, position, and channel
%           of the slice, respectively. The type and dimensions of the
%           returned value must be identical for each of the selected
%           images. The returned value may be an empty array, in which
%           case all returned values must be empty, and the result is 
%           any empty array of the same type.
%
%   This function returns the following values:
%
%       flat - The concatenated transformed values.
%
%       info - The information associated with each slice, in
%           the same format used for the transform function above.
%
%   Examples:
%
%       To create a 3D matrix from a movie of 2D images:
%
%           rnge.channels = 1;
%           flat = mov.flatten(rnge, []);
%
%       The above may be equivalently implemented as:
%
%           range.channels = 1;
%           flat = mov.flatten(rnge, @(im, info) im);
%
%       To find the time response of a value within an ROI:
%
%           flat = mov.flatten(rnge, @(im, info)
%               mean(im(roi_mask));


% Check input arguments
[rnge_ok, rnge_msg] = obj.isValidRange(rnge);
if ~rnge_ok
    error(rnge_msg, 'rnge');
end

if ~isempty(transform)
    if ~isa(transform, 'function_handle')
        error 'The argument ''transform'' has an invalid type';
    end
end

% Normalize the input range
rnge = obj.normalizeRange(rnge);

% Traverse each slice to get a full slice count
count = 0;

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
        
        % Everything matched. Increment count
        count = count + 1;
    end
end

% Handle trivial case
if count == 0
    flat = [];
    info = [];
    return;
end

% Initialize output to be empty
has_first = false;
flat = [];

% Pre-allocate memory for info
info = [];
info(count).time = [];
info(count).pos = [];
info(count).ch = [];

% Traverse each slice
cur_idx = 1;
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
        
        % Everything matched, transform image
        im = obj.slices(slice_idx).channels(ch_idx).image;
        info_elem.time = obj.slices(slice_idx).time;
        info_elem.pos = obj.slices(slice_idx).pos;
        info_elem.ch = obj.slices(slice_idx).channels(ch_idx).ch;
        
        if ~isempty(transform)
            trans_im = transform(im, info_elem);
        else
            trans_im = im;
        end
        
        % Check if the returned value has a valid type, if so
        % then add it to the list
        if has_first
            if isempty(flat)
                if ~isempty(trans_im)
                    error 'The transformed value has an incompatible type';
                end
            else
                if ~strcmp(class(flat), class(trans_im)) || ...
                    ndims(flat) ~= ndims(trans_im) + 1
                    error 'The transformed value has an incompatible type';
                end
            
                flat_dims = size(flat);
                trans_im_dims = size(trans_im);
            
                if nnz(flat_dims(1:end-1) - trans_im_dims(1:end-1)) > 0
                    error 'The transformed value has an incompatible type';
                end
                
                % The following assignment must work with trans_im
                % of arbitrary data. this should be implemented with
                % subsasgn but that function was implemented stupidly
                % and performance is absolutely pathetic (R2013a).
                % instead we must access the data using 1-dimensional
                % projection as follows
                first_idx = 1 + (cur_idx-1) * numel(trans_im);
                last_idx = cur_idx * numel(trans_im);
                flat(first_idx:last_idx) = trans_im(:);
            end
        else
            if ~isnumeric(trans_im)
                error 'The transformed value must be a numeric type';
            end
            
            if ~isempty(trans_im)
                % Preallocate the output array
                flat = zeros([size(trans_im), count], 'like', trans_im);

                % See similar assignment in has_first case for explanation
                % of assignment logic
                first_idx = 1 + (cur_idx-1) * numel(trans_im);
                last_idx = cur_idx * numel(trans_im);
                flat(first_idx:last_idx) = trans_im(:);
            else
                % Create an empty array with appropriate type
                flat = zeros(0, 'like', trans_im);
            end
            
            has_first = true;
        end
        
        info(cur_idx) = info_elem; %#ok<AGROW>
        
        cur_idx = cur_idx + 1;
    end
end

end

