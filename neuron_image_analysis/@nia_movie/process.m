function process(obj, rnge, slice_transform, ch_transform)
%PROCESS Transform metadata for movie
%   process(range, slice_transform, ch_transform) traverses the movie
%   modifying the metadata for each slice and channel. Note that processing
%   is performed in-place. While this is generally more efficient, in the
%   event of an error (if the transforms return invalid results) this
%   function must destroy the movie to avoid creating a movie with an
%   invalid state. As such when interactively working with data sets it may
%   be wise to first copy the movie using the copy() function and then
%   invoke this function on the copy. That way, if the transform fails, you
%   still have a valid copy of the movie. Of course, this copy should be
%   omitted when this function is called from scripts or functions.
%   Finally, note that this function is able to modify the number of
%   elements in the position vector, unlike functions that operate on one
%   slice at a time.
%
%   This function accepts the following arguments:
%
%       rnge - The range of slices and channels to be used for the
%           transform operations. See isValidRange() for more information
%           about the format of this argument. Note that the slice
%           transform is performed when the first matching channel is
%           encountered. If no channel matches, then the slice is not
%           transformed.
%
%       slice_transform - Function handle that accepts a structure
%           containing information on the current slice, and returns a
%           structure with modified information. The input structure has
%           the fields:
%
%               slice - The index of the slice in the movie.
%               time - The time of the slice in the movie.
%               pos - The position of the slice in the movie.
%
%           The function must return a structure with the fields:
%
%               time - The time of the slice in the modified movie.
%               pos - The position of the slice in the modified movie.
%
%           The argument slice_transform may also be an empty array,
%           in which case no transform is performed.
%
%       ch_transform - Function handle that accepts a structure
%           containing information on the current channel, and returns
%           a structure with modified information. The input structure
%           has the same fields as slice_transform, and the fields:
%
%               ch - The channel of the image in the movie.
%
%           The function must return a structure with the fields:
%
%               ch - The channel of the image in the modified movie.
%
%           The argument ch_transform may also be an empty array,
%           in which case no transform is performed.
%

% Check input arguments
[rnge_ok, rnge_msg] = obj.isValidRange(rnge);
if ~rnge_ok
    error(rnge_msg, 'rnge');
end

% Test the transforms
if ~(isempty(slice_transform) || isa(slice_transform, 'function_handle'))
    error 'The argument ''slice_transform'' has an invalid type';
end

if ~(isempty(ch_transform) || isa(ch_transform, 'function_handle'))
    error 'The argument ''ch_transform'' has an invalid type';
end

% Normalize the input range
rnge = obj.normalizeRange(rnge);

% Traverse each slice in movie
for slice_idx=1:length(obj.slices)
    
    % Skip non-matching slices
    if ~obj.inSliceRange(slice_idx, rnge)
        continue;
    end
    
    in_slice_info = [];
    
    % Traverse each channel in slice
    for ch_idx=1:length(obj.slices(slice_idx).channels)
 
        % Skip non-matching channels
        if ~obj.inChannelRange(slice_idx, ch_idx, rnge)
            continue;
        end
        
        % Everything matched. Append the slice to the list. If this is
        % the first channel to match, then we need to transform the slice.
        if isempty(in_slice_info)
            
            in_slice_info.slice = slice_idx;
            in_slice_info.time = obj.slices(slice_idx).time;
            in_slice_info.pos = obj.slices(slice_idx).pos;
            
            if ~isempty(slice_transform)
                try
                    out_slice_info = slice_transform(in_slice_info);
                catch err
                    nukeMovie(obj);
                    rethrow(err);
                end
                
                if ~isValidSliceInfo(out_slice_info)
                    nukeMovie(obj);
                    error 'The slice transform returned an invalid object';
                end
                
                obj.slices(slice_idx).time = out_slice_info.time;
                obj.slices(slice_idx).pos = out_slice_info.pos;
            end
        end
        
        if ~isempty(ch_transform)
            in_ch_info = in_slice_info;
            in_ch_info.ch = obj.slices(slice_idx).channels(ch_idx).ch;
            
            try
                out_ch_info = ch_transform(in_ch_info);
            catch err
                nukeMovie(obj);
                rethrow(err);
            end
            
            if ~isValidChannelInfo(out_ch_info)
                nukeMovie(obj); 
                error 'The channel transform returned an invalid object';
            end
            
            obj.slices(slice_idx).channels(ch_idx).ch = out_ch_info.ch;
        end
    end
end

% Go back through the movie to make sure to update the channel list
% and check the position vector length
ch_list = [];
new_pos_lens = [];

for slice_idx=1:length(obj.slices)
    slice_pos_lens = cellfun(@length, obj.slices(slice_idx).pos);
    
    % Check the length of the position vector
    if isempty(new_pos_lens)
        new_pos_lens = slice_pos_lens;
    elseif any(new_pos_lens ~= slice_pos_lens)
        nukeMovie(obj);
        error 'The transform returned an invalid position vector';
    end
    
    % Update the channel list
    for ch_idx=1:length(obj.slices(slice_idx).channels)
        new_ch = obj.slices(slice_idx).channels(ch_idx).ch;
        ch_list = unique([ch_list, new_ch]);
    end
end


% Update remaining fields
obj.ch_list = ch_list;
obj.pos_lens = new_pos_lens;

% Remove existing cache information
if ~isempty(slice_transform)
    obj.pos_lu = {};
    obj.pos_ranges = {};
end

if ~isempty(ch_transform)
    obj.hist_info = [];
end

end

function out = isValidSliceInfo(info)
% This function returns true if the slice info is valid and false
% otherwise.

% Check if info is a struct
if ~isstruct(info) || length(info) ~= 1
    out = false;
    return;
end

% Check if info has valid time field
if ~isfield(info, 'time') || ~isfloat(info.time) || ...
    ~isreal(info.time) || ~isscalar(info.time)
    out = false;
    return;
end

% Check if info has pos field with valid type
if isfield(info, 'pos')
    % We intentionally do not test the lengths of the position
    % vectors, since those are allowed to change, and will get
    % checked at the very end for consistency
    
    if ~iscell(info.pos) || ~nia_isRow(info.pos)
        out = false;
        return;
    end
    
    for vec_idx=1:length(info.pos)
        cur_pos = info.pos{vec_idx};
        
        if ~isfloat(cur_pos) || ~isreal(cur_pos) || ~nia_isRow(cur_pos)
            out = false;
            return;
        end
        
        if nnz(abs(cur_pos - round(cur_pos))) > 0 || nnz(cur_pos < 0) > 0
            out = false;
            return;
        end
    end
end

out = true;

end

function out = isValidChannelInfo(info)
% This function returns true if the channel info is valid and false
% otherwise.

% Check if info is a struct
if ~isstruct(info) || length(info) ~= 1
    out = false;
    return;
end

% Check if info has valid ch field
if ~isfield(info, 'ch') || ~isfloat(info.ch) || ...
    ~isreal(info.ch) || ~isscalar(info.ch) || info.ch < 0
    out = false;
    return;
end

out = true;

end

function nukeMovie(obj)
% This function empties the contents of the movie. This function
% is called when the user provided invalid input and the movie must
% be destroyed to avoid creating a movie with an inconsistent state.

obj.slices = [];
obj.ch_list = [];
obj.pos_lens = [];
obj.pos_lu = {};
obj.pos_ranges = {};
obj.hist_info = [];
end
