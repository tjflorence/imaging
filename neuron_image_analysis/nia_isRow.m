function out = nia_isRow(input)
%NIA_ISROW Returns true if input is a row vector
%   out = nia_isRow(input) returns true if the input is a row vector. It is
%   very similar to the builtin isrow() function except that the builtin
%   function returns false if the vector is empty. This function will
%   return true for an empty vector, which is almost certainly the correct
%   generalization, and makes it much easier to use with code for which an
%   empty vector is valid.

out = isempty(input) || (ismatrix(input) && size(input,1) == 1);

end

