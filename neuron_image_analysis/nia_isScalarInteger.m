function out = nia_isScalarInteger(x)
%NIA_ISSCALARINTEGER Test if input is a scalar integer.
%   out = nia_isScalarInteger(x) returns true if the input 'x' is a
%   scalar floating point type with an integer value.

out = isfloat(x) && isreal(x) && isscalar(x) && ...
    abs(x-round(x)) == 0;
end