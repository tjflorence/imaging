function slice_idx = appendSlice(obj, slice_info)
%APPENDSLICE Append a slice to the movie
%   slice_idx = appendSlice(slice_info) adds a slice to the given movie and
%   returns the slice index for the new entry. The passed slice_info a must
%   be a structure with the fields 'time' and 'pos'. The 'time' field
%   contains the timestamp, and the 'pos' field contains the slice
%   position. If an error occurs, then the movie is unmodified.

% Check input arguments
if ~isstruct(slice_info) || length(slice_info) ~= 1
    error 'The argument ''slice_info'' has an invalid type';
end

slice_info_names = {'time', 'pos'};
[slice_info_ok, slice_info_msg] = nia_hasValidFieldNames(...
    slice_info, slice_info_names, slice_info_names);

if ~slice_info_ok
    error(slice_info_msg, 'slice_info');
end

if ~isfloat(slice_info.time) || ~isreal(slice_info.time) || ...
        ~isscalar(slice_info.time)
    error 'The argument ''slice_info.time'' has an invalid value';
end

if ~obj.isEmpty()
    was_empty = false;
    
    if ~iscell(slice_info.pos) || ~nia_isRow(slice_info.pos) || ...
            length(slice_info.pos) ~= length(obj.pos_lens)
        error 'The argument ''slice_info.pos'' has an invalid value';
    end

    for vec_idx=1:length(slice_info.pos)
        cur_pos = slice_info.pos{vec_idx};

        if ~isfloat(cur_pos) || ~isreal(cur_pos) || ~nia_isRow(cur_pos) || ...
                length(cur_pos) ~= obj.pos_lens(vec_idx)
            error 'The argument ''slice_info.pos'' has an invalid value';
        end

        if nnz(abs(cur_pos - round(cur_pos))) > 0 || nnz(cur_pos < 0) > 0
            error 'The argument ''slice_info.pos'' must have integer values';
        end
    end
else
    was_empty = true;
    
    if ~iscell(slice_info.pos) || ~nia_isRow(slice_info.pos)
        error 'The argument ''slice_info.pos'' has an invalid value';
    end

    for vec_idx=1:length(slice_info.pos)
        cur_pos = slice_info.pos{vec_idx};

        if ~isfloat(cur_pos) || ~isreal(cur_pos) || ~nia_isRow(cur_pos)
            error 'The argument ''slice_info.pos'' has an invalid value';
        end

        if nnz(abs(cur_pos - round(cur_pos))) > 0 || nnz(cur_pos < 0) > 0
            error 'The argument ''slice_info.pos'' must have integer values';
        end
    end
end

% Append to end of slices array
slice_idx = length(obj.slices) + 1;

obj.slices(slice_idx).time = slice_info.time;
obj.slices(slice_idx).channels = [];
obj.slices(slice_idx).pos = slice_info.pos;

% Initialize positions if necessary
if was_empty
    obj.pos_lens = cellfun(@length, slice_info.pos);
end

% Invalidate caches
obj.pos_lu = {};
obj.pos_ranges = {};
obj.ch_list = [];
obj.hist_info = [];

end