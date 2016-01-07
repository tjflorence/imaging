function setSliceInfo(obj, slice_idx, slice_info, skip_caches)
%SETSLICEINFO Retrieves metadata for a single slice
%   setSliceInfo(slice_idx, slice_info, skip_update) modifies the metadata
%   for a single slice of the movie. The input slice_info must be a
%   structure with the fields 'time' and 'pos'. The 'time' field is a
%   scalar double that contains the timestamp and the 'pos' field is a cell
%   array of row vectors describing the slice positions. Any of these
%   fields may be omitted, in which case that element of the slice info is
%   unmodified. If skip_caches is true, then the cache information is left
%   unaltered, and subsequent calls to getSliceIndex() will return the
%   slice index as if the call to setSliceInfo() had not occurred. The
%   function invalidatePosCache() should ultimately be called to force an
%   update of the position cache information. If skip_caches is false, or
%   omitted then subsequent calls to getSliceIndex() will use the current
%   information. Note that this forces the cache to be recomputed and can
%   incur a very substantial performance penalty if setSliceInfo() and
%   getSliceIndex() are used together in a small loop. If an error occurs,
%   then the slice information is unmodified.

% Assign default arguments
if nargin < 4
    skip_caches = false;
end

% Check input arguments
if ~nia_isScalarInteger(slice_idx) || ...
        slice_idx < 1 || slice_idx > length(obj.slices)
    error 'The argument ''slice_idx'' has an invalid value';
end

if ~isstruct(slice_info) || length(slice_info) ~= 1
    error 'The argument ''slice_info'' has an invalid type';
end

slice_info_names = {'time', 'pos'};
[slice_info_ok, slice_info_msg] = nia_hasValidFieldNames(...
    slice_info, slice_info_names, {});

if ~slice_info_ok
    error(slice_info_msg, 'slice_info');
end

if isfield(slice_info, 'time')
    if ~isfloat(slice_info.time) || ~isreal(slice_info.time) || ...
            ~isscalar(slice_info.time)
        error 'The argument ''slice_info.time'' has an invalid value';
    end
end

if isfield(slice_info, 'pos')
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
end

if ~islogical(skip_caches) || ~isscalar(skip_caches)
    error 'The argument ''force_update'' must be a scalar logical';
end

% Retrieve slice info
if isfield(slice_info, 'time')
   obj.slices(slice_idx).time = slice_info.time;
end

if isfield(slice_info, 'pos')
    obj.slices(slice_idx).pos = slice_info.pos;
    
    if ~skip_caches
        obj.pos_lu = {};
        obj.pos_ranges = {};
    end
end

end