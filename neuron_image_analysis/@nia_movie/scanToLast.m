function out = scanToLast(obj, pos_vec, home_pos, scan_extents)
%SCANTOLAST Find the last position with a valid slice in scan
%   out = scanToLast(pos_vec, home_pos, scan_extents) scans through
%   the positions to find the last valid position in the scan.
%
%       pos_vec - A scalar integer that contains the index of the position 
%           vector to use for the scan.
%
%       home_pos - A [1xN] array where N is the number of position
%           dimensions, that specifies the position from which to begin the
%           scan.
%
%       scan_extents - A [2x1] array that specifies how to scan through the
%           position data. The first element specifies the index into the
%           position vector to scan, and the second element specifies the
%           value to scan to (the first value is given in the home_pos
%           argument and must be less than or equal to the value given
%           here). If the second element is NaN, then the scan proceeds
%           through all valid position values.
%
%   This function has the following outputs:
%
%       out - The value of the scanned position element that is
%           immediately after the last valid position. If no valid elements
%           are found, this is equal to the first scanned position value in
%           the scan. If all positions in the scan are valid, then this is
%           equal to one plus the last scanned position value in the scan.

% This function exists solely to document the mex-function
% that does all of the actual work, and to provide a very
% small amount of input processing.

% Check input arguments
if ~nia_isScalarInteger(pos_vec) || ...
        pos_vec < 1 || pos_vec > length(obj.pos_lens)
    error 'The argument ''pos_vec'' has an invalid value';
end

if ~isfloat(home_pos) || ~isreal(home_pos) || ...
        ~nia_isRow(home_pos) || length(home_pos) ~= obj.pos_lens(pos_vec)
    error 'The argument ''home_pos'' has an invalid value';
end

if nnz(abs(home_pos - round(home_pos))) > 0 || nnz(home_pos < 0) > 0
    error 'The argument ''home_pos'' has an invalid value';
end

% Note that the following specifically is specifically designed to
% allow NaN values for scan_extents(2)
if ~isfloat(scan_extents) || ~isreal(scan_extents) || ...
        ~iscolumn(scan_extents) || length(scan_extents) ~= 2 || ...
        nnz(abs(scan_extents - round(scan_extents)) > 0) > 0 || ...
        isnan(scan_extents(1))
    error 'The argument ''scan_extents'' has an invalid value';
end

if scan_extents(1) < 0 || scan_extents(1) > obj.pos_lens(pos_vec)
    error 'The argument ''scan_extents'' has an out of range value';
end

% Load the cache if required
if isempty(obj.pos_lu)
    obj.loadPosCache();
end

% Load default value for scan extents if requested
if isnan(scan_extents(2))
    scan_extents(2) = obj.pos_ranges{pos_vec}(2, scan_extents(1));
end

% Convert the positions to unsigned integers
home_pos = uint64(home_pos);
scan_extents = uint64(scan_extents);

out = nia_scanToLastNiaMovieImpl(obj.slices, obj.pos_lu{pos_vec}, ...
    home_pos, scan_extents);

out = double(out);

end
