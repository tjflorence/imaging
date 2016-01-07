function append(obj, mov2, rnge, transform)
%APPEND Append the passed movie to the end of current movie
%   append(mov2, rnge, transform) appends the passed movie to the
%   end of the current movie. In order to produce more semantically
%   meaningful concatenations, this function can also modify the time
%   and position values of the appended movie prior to appending it to
%   the current movie. This function accepts the following arguments:
%
%       mov2 - Movie to append to end of current movie. The length of the
%           position vector in mov2 must be compatible with the length of
%           the position vector in the current movie. Specifically, after
%           the position vector of mov2 is transformed, both movies must
%           have position vectors of identical length. There is an an
%           important exception to this rule: if the current movie is empty
%           then the position vector in mov2 may have any length.
%
%       rnge - The range of slices and channels to be select from mov2 for
%           append operation. See isValidRange() for more information about
%           the format of this argument.
%
%       transform - Each of the slices of the appended movie is 
%           transformed, which allows the concatentation to produce
%           more semantically meaningful movies. If this argument is
%           empty, than no transform is applied. If this argument is
%           non-empty then it must be a structure with one or more
%           of the following fields:
%
%           time_offset - The time stamp for each slice of the appended
%               movie is incremented by this value. If this field is an
%               empty array, then the time is unaltered.
%            
%           pos_val - A [3x1] array that specifies how the position
%               vector of the appended movie should be altered prior to
%               concatentation. The first element specifies the index of
%               the position vector to alter, the second element specifies
%               the index of the element within the vector to alter, and
%               the third element specifies the new value. If the index is
%               one greater than the current length of the position vector
%               of mov2, then the position vector of mov2 is increased in
%               size. Otherwise, the index must be less than the number of
%               position dimensions in mov2.
%


% Check input arguments
if ~isa(mov2, 'nia_movie')
    error 'The argument ''mov2'' must be a nia_movie object';
end

if obj == mov2
    error 'The argument ''mov2'' cannot be the same as target movie';
end

[rnge_ok, rnge_msg] = mov2.isValidRange(rnge);
if ~rnge_ok
    error(rnge_msg, 'rnge');
end

% Check that the transform has a valid form
checkTransform(obj, mov2, transform);

% Check that position dimensions are compatible, we only need to
% check for the case that the transform is empty, since the function
% checkTransform() will have already checked the non-empty case
if isempty(transform)
    if ~isempty(obj.slices)
        if length(obj.pos_lens) ~= length(mov2.pos_lens) || ...
                any(obj.pos_lens ~= mov2.pos_lens)
            error 'The argument mov2 has an incomparticle position vector';
        end
    end
end

% Normalize the input range
rnge = mov2.normalizeRange(rnge);

% Traverse each slice to get a full slice count
append_count = 0;
ch_count_per_slice = zeros(1, length(mov2.slices));

for mov2_slice_idx=1:length(mov2.slices)
    ch_count = 0;
    
    % Skip non-matching slices
    if ~mov2.inSliceRange(mov2_slice_idx, rnge)
        continue;
    end
    
    % Traverse each channel in slice
    for mov2_ch_idx=1:length(mov2.slices(mov2_slice_idx).channels)
 
        % Skip non-matching channels
        if ~mov2.inChannelRange(mov2_slice_idx, mov2_ch_idx, rnge)
            continue;
        end
        
        % Everything matched. Increment counts
        ch_count = ch_count + 1;
    end
    
    append_count = append_count + (ch_count > 0);
    ch_count_per_slice(mov2_slice_idx) = ch_count;
end

% Check for trivial case
if append_count == 0
    return;
end

% Allocate memory for new mov1 slices array
total_count = length(obj.slices) + append_count;

obj.slices(total_count).time = [];
obj.slices(total_count).channels  = [];
obj.slices(total_count).pos = [];

% Retrieve transform info
if ~isempty(transform) && isfield(transform, 'time_offset') ...
        && ~isempty(transform.time_offset)
    time_offset = transform.time_offset;
else
    time_offset = 0;
end

if ~isempty(transform) && isfield(transform, 'pos_val') ...
        && ~isempty(transform.pos_val)
    pos_val = transform.pos_val;
else
    pos_val = [];
end

% Traverse each slice in mov2
obj_slice_idx = total_count - append_count + 1;
ch_list = obj.ch_list;
new_pos = [];

for mov2_slice_idx=1:length(mov2.slices)
    
    % Skip non-matching slices
    if ~mov2.inSliceRange(mov2_slice_idx, rnge)
        continue;
    end
    
    obj_ch_idx = 1;
    
    % Traverse each channel in slice
    for mov2_ch_idx=1:length(mov2.slices(mov2_slice_idx).channels)
 
        % Skip non-matching channels
        if ~mov2.inChannelRange(mov2_slice_idx, mov2_ch_idx, rnge)
            continue;
        end
        
        % Everything matched. Append the slice to the list.  If this is
        % the first channel to match, then we need to prep the slice.
        if obj_ch_idx == 1
            ch_count = ch_count_per_slice(mov2_slice_idx);
            
            obj.slices(obj_slice_idx).time = ...
                mov2.slices(mov2_slice_idx).time + time_offset;
            
            % Preallocate the channels array
            obj.slices(obj_slice_idx).channels(ch_count).ch = [];
            obj.slices(obj_slice_idx).channels(ch_count).image = [];
            
            new_pos = mov2.slices(mov2_slice_idx).pos;
            
            if ~isempty(pos_val)
                new_pos{pos_val(1)}(pos_val(2)) = pos_val(3);
            end
            
            obj.slices(obj_slice_idx).pos = new_pos;
        end
        
        % Copy image to new structure
        new_ch = mov2.slices(mov2_slice_idx).channels(mov2_ch_idx).ch;
        ch_list = unique([ch_list, new_ch]);
        
        obj.slices(obj_slice_idx).channels(obj_ch_idx).ch = new_ch;
        obj.slices(obj_slice_idx).channels(obj_ch_idx).image = ...
            mov2.slices(mov2_slice_idx).channels(mov2_ch_idx).image;

        obj_ch_idx = obj_ch_idx + 1;
    end
    
    obj_slice_idx = obj_slice_idx + 1;
end

% Update remaining mov1 fields
obj.ch_list = ch_list;
obj.pos_lens = mov2.pos_lens;
    
if ~isempty(pos_val)
    obj.pos_lens(pos_val(1)) = length(new_pos{pos_val(1)});
end

% Remove outdated fields
obj.pos_lu = {};
obj.pos_ranges = {};
obj.ch_list = [];
obj.hist_info = [];

end

function checkTransform(obj, mov2, transform)
% This function checks if the transform has a valid format. If there
% is a problem, it reports the error, otherwise it takes no action.

% Check type of the transform variable
if ~isempty(transform)
    if ~isstruct(transform) || length(transform) ~= 1
        error('The argument ''transform'' has an invalid type');
    end
    
    transform_names = {'time_offset', 'pos_val'};
    [transform_ok, transform_msg] = nia_hasValidFieldNames(...
        transform, transform_names, {});
    if ~transform_ok
        error(transform_msg, 'transform');
    end
    
    if isfield(transform, 'time_offset') && ~isempty(transform.time_offset)
        if ~isfloat(transform.time_offset) || ~isreal(transform.time_offset) ...
                || ~isscalar(transform.time_offset)
            error 'The argument ''transform.time_offset'' has an invalid type';
        end
    end
    
    if isfield(transform, 'pos_val') && ~isempty(transform.pos_val)
        if ~isfloat(transform.pos_val) || ~isreal(transform.pos_val) ...
                || ~iscolumn(transform.pos_val) || length(transform.pos_val) ~= 3
            error 'The argument ''transform.pos_val'' has an invalid type';
        end
        
        if nnz(abs(transform.pos_val - round(transform.pos_val))) > 0
            error 'The argument ''transform.pos_val'' has an invalid value';
        end
        
        if ~isempty(obj.slices)
            vec_idx = transform.pos_val(1);
            
            if length(obj.pos_lens) ~= length(mov2.pos_lens)
                error 'The two movies have incompatible position vectors';
            end
            
            if vec_idx < 1 || vec_idx > length(obj.pos_lens)
                error 'The argument ''transform.pos_val(1)'' has an invalid value';
            end
            
            vec_mask = (1:length(obj.pos_lens)) ~= vec_idx;
            
            if any(obj.pos_lens(vec_mask) ~= mov2.pos_lens(vec_mask))
                error 'The two movies have incomparticle position vectors';
            end
            
            if transform.pos_val(2) == mov2.pos_lens(vec_idx)+1
                if obj.pos_lens(vec_idx) ~= mov2.pos_lens(vec_idx)+1
                    error 'The argument ''transform.pos_val(2)'' has an invalid value';
                end
            else
                if obj.pos_lens(vec_idx) ~= mov2.pos_lens(vec_idx)
                    error 'The argument ''transform.pos_val(2)'' has an invalid value'
                end
            end
        else
            vec_idx = transform.pos_val(1);
            
            if transform.pos_val(2) < 1 || transform.pos_val(2) > mov2.pos_lens(vec_idx)+1
                error 'The argument ''transform.pos_val(2)'' has an invalid value';
            end
        end
        
        if transform.pos_val(3) < 0
            error 'The argument ''transform.pos_val(3)'' has an invalid value';
        end
    end
end

end
