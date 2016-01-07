function pos_ranges = getPositionRanges(obj)
%GETPOSITIONRANGES Retrieve valid ranges for position vector
%   pos_ranges = getPositionRanges() retrieves the ranges for the values in
%   the position vector of the movie. The range is returned as a [1xM] cell
%   array of [2xN] arrays where M is the number of position vectors in the
%   movie and N is the number of dimensions each position vector. Each
%   column contains the range for a dimension of the position vector. The
%   first row contains the minimum value for the position and the second
%   row contains the maximum value, both given inclusively. This
%   information is cached, so that the first call to this function may
%   require some time to calculate the position ranges, but subsequent
%   calls will return the stored indices quickly. This cache can be
%   manually loaded using the function loadPosCache(). Some modifications
%   to the movie may cause the position ranges to become out-of-date, at
%   which point the cached data will be recomputed.

if isempty(obj.pos_ranges)
    obj.loadPosCache();
end

pos_ranges = obj.pos_ranges;
end