function state = isEmpty(obj)
%ISEMPTY Test if movie is empty
%   state = isEmpty() returns true if the current movie is empty,
%   which is to say that it has no slices.

state = isempty(obj.slices);

end

