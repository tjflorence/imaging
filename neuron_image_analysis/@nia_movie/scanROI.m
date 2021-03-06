function out = scanROI(obj, pos_vec, home_pos, scan_extents, ...
    channel_id, mask)
%SCANROI Find mean value within roi across frames
%   out = scanROI(pos_vec, home_pos, scan_extents, channel_id, mask) finds
%   the mean value of the pixels within the passed ROI for each frame
%   within the passed movie. For each position vector value this function
%   will use the first frame in the movie with a matching position. This
%   function accepts the following arguments:
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
%           up to the last valid position value.
%
%       channel_id - A scalar integer specifying the channel ID to use for the
%           ROI.
%
%       mask - A N-dimensional logical array the same size as the image
%           data with true values for pixels that should be included in the
%           ROI. For efficiency, all image data must have the dimensions.
%           If image data is encountered with a different dimension, then
%           that slice is given a value of NaN.
%
%   This function has the following outputs:
%
%       out - A [3xN] vector where N is equal to the number of processed
%           slices. The first row contains the slice index from which the
%           ROI was calculated, and the second row contains the time
%           associated with the slice, and the third row contains the
%           calculated ROI value.

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

if ~nia_isScalarInteger(channel_id) || channel_id < 0
    error 'The argument ''channel'' has an out of range value';
end

% Load the cache if required
if isempty(obj.pos_lu)
    obj.loadPosCache();
end

% Load default value for scan extents if requested
if isnan(scan_extents(2))
    after_pos_val = obj.scanToLast(pos_vec, home_pos, ...
        [scan_extents(1); NaN]);
    
    if after_pos_val == home_pos(scan_extents(1))
        out = [];
        return;
    end
    
    scan_extents(2) = after_pos_val - 1;
end

% Convert the positions to unsigned integers
home_pos = uint64(home_pos);
scan_extents = uint64(scan_extents);
channel_id = uint64(channel_id);

out = nia_scanROINiaMovieImpl(obj.slices, obj.pos_lu{pos_vec}, ...
    home_pos, scan_extents, channel_id, mask);
    
end
