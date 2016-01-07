function loadPosCache(obj)
%LOADPOSCACHE Prepare the position lookup information for movie
%   loadPosCache() scans through the current movie contents and prepares
%   the position lookup cache. If the position lookup cache has been
%   previously loaded, then this function returns immediately.

% Check if we already have position information
if ~isempty(obj.pos_lu) && ~isempty(obj.pos_ranges)
    return;
end

obj.pos_lu = {};
obj.pos_ranges = {};

% Handle trivial case first
if isempty(obj.slices)
    return;
end

% Traverse each position set
nvecs = length(obj.pos_lens);
obj.pos_lu = cell(1, nvecs);
obj.pos_ranges = cell(1, nvecs);

for vec_idx=1:nvecs
    [obj.pos_lu{vec_idx}, obj.pos_ranges{vec_idx}] = ...
        populatePositionSet(obj, vec_idx);
end

end

function [pos_lu_entry, pos_ranges_entry] = populatePositionSet(obj, pos_vec)
% This function constructs the pos_lu and pos_ranges entrys
% for a single position set.

cfg_hash_len = 4096;

% Handle trivial case
if obj.pos_lens(pos_vec) == 0
    pos_lu_entry = [];
    pos_ranges_entry = [];
    return;
end

% Find the position ranges
pos_ranges_entry = zeros(2, obj.pos_lens(pos_vec));
pos_ranges_entry(1,:) = Inf;
pos_ranges_entry(2,:) = -Inf;

for slice_idx=1:length(obj.slices)
    cur_pos = obj.slices(slice_idx).pos{pos_vec};
    
    pos_ranges_entry(1,:) = min([pos_ranges_entry(1,:); cur_pos], [], 1);
    pos_ranges_entry(2,:) = max([pos_ranges_entry(2,:); cur_pos], [], 1);
end

% Allocate the position lookup table. The position lookup table is
% a cell array of matrices. Each of those matrices is a [MxN+1] matrix
% where M is the number of entries in the bin and N is the length of the
% position vector. The first element of each row is the corresponding
% slice index.
pos_lu_entry = cell(1, cfg_hash_len);

% Populate the position lookup
for slice_idx=1:length(obj.slices)
    cur_pos = obj.slices(slice_idx).pos{pos_vec};
    code = nia_getHashCodeUInt64(uint64(cur_pos));
    bin_idx = mod(code, cfg_hash_len) + 1;
    
    bin = pos_lu_entry{bin_idx};
    
    if isempty(bin)
        % Empty bin, use as first element
        pos_lu_entry{bin_idx} = [slice_idx, cur_pos];
        continue;       
    end
    
    if nnz(all(bsxfun(@eq, bin(:, 2:end), cur_pos), 2)) > 0
        % Already have this position in lookup, skip
        continue;
    end
    
    % Append the entry to the end of the list
    bin(end+1,:) = [slice_idx, cur_pos]; %#ok<AGROW>
    
    pos_lu_entry{bin_idx} = bin;
end

end

