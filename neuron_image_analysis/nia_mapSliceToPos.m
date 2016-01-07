function out_info = nia_mapSliceToPos(in_info)
%NIA_MAPSLICETOPOS Transforms the input slice info
%   out_info = nia_mapSliceToPos(in_info) is a function intended for use as
%   a slice transform for the function nia_movie.process().  When passed as
%   a slice transform this function will add an element to the end of
%   the position vector list and give it a value that is equal to the slice
%   index number.
%
%   Example:
%
%       mov.process([], [], @nia_mapSliceToPos);

out_info.time = in_info.time;
out_info.pos = [in_info.pos, in_info.slice];

end

