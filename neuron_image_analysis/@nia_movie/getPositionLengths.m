function pos_lens = getPositionLengths(obj)
%GETPOSITIONLENGTHS Get the number of dimensions in position vector
%   pos_lens = getPositionLengths() retrieves the number of dimensions in
%   each position vector as a [1xM] array, were M is the number of position
%   vectors. This information is always available, so the cost of invoking
%   this function is always negligible.

pos_lens = obj.pos_lens;
end