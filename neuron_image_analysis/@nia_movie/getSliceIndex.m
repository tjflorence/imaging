function slice_idx = getSliceIndex(obj, pos_vec, pos)
%GETSLICEINDEX Retrieve first slice index for a position
%   slice_idx = getSliceIndex(pos_vec, pos) retrieves the first slice in
%   the movie with a given position value. This function identifies the
%   first slice in the 'pos_vec' position vector with the position value
%   equal to 'pos'. If that position does not appear in the movie, then
%   this function returns zero. This information is cached, so that the
%   first call to this function may require some time to calculate the
%   position lookup table, but subsequent calls will return the stored
%   indices quickly. The cache can be manually loaded using the function
%   loadPosCache(). Some modifications to the movie may cause the position
%   lookup to become out-of-date, at which point the cached lookup table
%   will be recomputed.

% Check input arguments
if ~nia_isScalarInteger(pos_vec) || ...
        pos_vec < 1 || pos_vec > length(obj.pos_lens)
    error 'The argument ''pos_vec'' has an invalid value';
end

if ~isfloat(pos) || ~isreal(pos) || ~nia_isRow(pos) || ...
        length(pos) ~= obj.pos_lens(pos_vec)
    error 'The argument ''pos'' has an invalid type';
end

% Check that position is an integer
if nnz(abs(pos-round(pos))) > 0 || nnz(pos < 0) > 0
    error 'The argument ''pos'' must have integer values';
end

% Compute the lookup table
if isempty(obj.pos_lu)
    obj.loadPosCache();
end

% Check that position is within range
if nnz(pos < obj.pos_ranges{pos_vec}(1,:) | ...
       pos > obj.pos_ranges{pos_vec}(2,:)) > 0
    slice_idx = 0;
    return;
end

% Use hash table to get the correct slice
code = nia_getHashCodeUInt64(uint64(pos));
bin_idx = mod(code, length(obj.pos_lu{pos_vec})) + 1;
bin = obj.pos_lu{pos_vec}{bin_idx};

if isempty(bin)
    slice_idx = 0;
    return;
end

match_mask = all(bsxfun(@eq, bin(:, 2:end), pos), 2);
match_idx = find(match_mask, 1, 'first');

if isempty(match_idx)
    slice_idx = 0;
    return;
end

slice_idx = bin(match_idx, 1);

end