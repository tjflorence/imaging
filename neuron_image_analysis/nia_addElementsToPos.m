function out_info = nia_addElementsToPos(num_elem, in_info)
%NIA_addElementsToPos Transforms the input slice info
%   out_info = nia_addElementsToPos(in_info) is a function intended for use
%   as a slice transform for the function nia_movie.process().  When passed
%   as a slice transform this function will append a new position vector to
%   the list of position vectors and initialize it with num_elem zeros.
%
%   Example:
%
%       mov.process([], [], @(x) nia_addElementsToPos(1, x) );

out_info.time = in_info.time;
out_info.pos = [in_info.pos, zeros(1, num_elem)];

end

