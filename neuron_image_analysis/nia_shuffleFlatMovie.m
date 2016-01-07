function [mov,order] = nia_shuffleFlatMovie(mov)
%NIA_SHUFFLEFLATMOVIE Shuffle the order of frames in a flat movie
%   [mov, order] = nia_shuffleFlatMovie(mov) shuffles the order of the
%   frames in the passed flat movie.  The movie must be passed a 3D array,
%   with frames in the last dimension.
%
%   This function has the following outputs:
%
%       mov - Shuffled movie array
%       order - Array of indices.  The first element holds
%           the index of the first frame in the shuffled
%           movie as it occured in the original movie, and
%           so forth.

% Check arguments
if ~isfloat(mov) || ~isreal(mov) || ...
        ndims(mov) ~= 3 || isempty(mov)
    error 'The argument ''mov'' has an invalid type';
end

% Create the new world order
nwo = randperm(size(mov,3));

% Perform the shuffle (we do this using
% swaps to avoid allocating another 3D
% matrix)
needs_swap = true(1, size(mov,3));

while true
    % Find the first frame needing re-positioning
    target_idx = find(needs_swap, 1, 'first');
    if isempty(target_idx)
        % If empty, then we are all done
        break;
    end
    
    % Find source for this frame
    source_idx = nwo(target_idx);
    if target_idx == source_idx
        % If target and source are same, then skip
        needs_swap(target_idx) = false;
        continue;
    end
    
    % Overwrite target with source
    swap_frame = mov(:,:,target_idx);
    mov(:,:,target_idx) = mov(:,:,source_idx);
    needs_swap(target_idx) = false;
    
    % Use overwritten frame as the new target
    % and perform cascading swap
    while nnz(needs_swap) ~= 0
        target_idx = find(nwo == target_idx);
        if ~needs_swap(target_idx)
            break;
        end
        
        tmp = mov(:,:,target_idx);
        mov(:,:,target_idx) = swap_frame;
        needs_swap(target_idx) = false;
        swap_frame = tmp;
    end
end

order = nwo;

end

