function [udata, stop] = analyze(obj, rnge, transform, udata)
%ANALYZE Modify user data based on each frame in range
%   udata = analyze(rnge, transform, udata) traverses the passed movie,
%   selects slices and channels in the passed range, applies a transform
%   to each that accepts a user-defined variable, and overwrites the
%   user-defined variable with the result. Since the transform and the
%   user-defined variable are arbitrary, this function can be used to
%   perform arbitrarily complicated analyses on the movie.
%
%   This function takes the following arguments:
%
%       rnge - The range of slices and channels to be used for the analysis
%           operation. See isValidRange() for more information about the
%           format of this argument.
%
%       transform - This argument specifies the transform that should be
%           performed to produce the modified user data. This argument
%           must be a function handle and is invoked as fcn(image, info,
%           udata) for each of selected images. This function must return a
%           modified user-data 'udata' that is used for subsequent images,
%           and a 'stop' value that is non-zero if analysis should be
%           discontinued. This function is passed the image data, a
%           structure of information, and the current user data. The
%           information structure contains the fields 'time', 'pos', and
%           'ch' that describe the time, position, and channel of the
%           image, respectively.
%
%       udata - This argument is passed to the transform given above. It
%           has arbitrary type and value.
%
%   This function returns the following values:
%
%       udata - The modified user data.
%
%       stop - True if analysis was stopped prematurely, false
%           otherwise.


% Check input arguments
[rnge_ok, rnge_msg] = obj.isValidRange(rnge);
if ~rnge_ok
    error(rnge_msg, 'rnge');
end

if ~isa(transform, 'function_handle')
    error 'The argument ''transform'' has an invalid type';
end

% Normalize the input range
rnge = obj.normalizeRange(rnge);

for slice_idx=1:length(obj.slices)
    
    % Skip non-matching slices
    if ~obj.inSliceRange(slice_idx, rnge)
        continue;
    end
    
    info.time = obj.slices(slice_idx).time;
    info.pos = obj.slices(slice_idx).pos;
    info.ch = 1;
    
    % Traverse each channel in slice
    for ch_idx=1:length(obj.slices(slice_idx).channels)
 
        % Skip non-matching channels
        if ~obj.inChannelRange(slice_idx, ch_idx, rnge)
            continue;
        end
        
        % Everything matched apply the transform
        info.ch = obj.slices(slice_idx).channels(ch_idx).ch;
        im = obj.slices(slice_idx).channels(ch_idx).image;
        
        [udata, stop] = transform(im, info, udata);
        if stop
            break;
        end
    end
end

end

