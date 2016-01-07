function out = nia_isString(input)
%NIA_ISSTRING Return true if input is a string.

out = ischar(input) && ismatrix(input) && size(input,1) < 2;

end