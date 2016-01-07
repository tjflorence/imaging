function out = nia_isColumn(input)
%NIA_ISCOLUMN Returns true if input is a column vector
%   out = nia_isColumn(input) returns true if the input is a column vector.
%   It is very similar to the builtin iscolumn() function except that the
%   builtin function returns false if the vector is empty. This function
%   will return true for an empty vector, which is almost certainly the
%   correct generalization, and makes it much easier to use with code for
%   which an empty vector is valid.

out = isempty(input) || (ismatrix(input) && size(input,2) == 1);

end

