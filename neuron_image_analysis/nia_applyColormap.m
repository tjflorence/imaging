function im = nia_applyColormap(input, zrange, cmap, nan_color)
%NIA_APPLYCOLORMAP Apply colormap to float-point image
%   im = nia_applyColormap(input, zrange, cmap, nan_color) transforms the
%   passed input image into a color image using the passed colormap.
%
%   This function accepts the following arguments:
%
%       input - Input black and white image (2D matrix)
%
%       zrange - Intensity range for colormap.  Must be a [1x2] matrix,
%           with minimum value as the first element and maximum value as
%           the second element.
%
%       cmap - Colormap to use.  Must be an [nx2] matrix with values
%           compatible with the colormap() function.
%
%       nan_color - Color to use for NaN values. Must be a [1x3] matrix,
%           with elements specifying the red, blue, and green elements on
%           a scale from zero to one.
%

% Check arguments
if ~isfloat(input) || ~isreal(input) || ~ismatrix(input)
    error 'The argument ''input'' must be a floating-point matrix';
end

if ~isfloat(zrange) || ~isreal(zrange) || ~ismatrix(zrange) || ...
        size(zrange, 1) ~= 1 || size(zrange, 2) ~= 2
    error 'The argument ''zrange'' must be a [1x2] matrix';
end

if ~isfloat(cmap) || ~isreal(cmap) || ~ismatrix(cmap) || isempty(cmap) || ...
        size(cmap, 2) ~= 3
    error 'The argument ''cmap'' must be an [nx3] matrix';
end

if min(cmap(:)) < 0 || max(cmap(:)) > 1.0
    error 'The argument ''cmap'' must have elements between zero and one';
end

if ~isfloat(nan_color) || ~isreal(nan_color) || ~ismatrix(nan_color) || ...
        size(nan_color, 1) ~= 1 || size(nan_color, 2) ~= 3
    error 'The argument ''nan_color'' must be a [1x3] matrix';
end

if min(nan_color(:)) < 0 || max(nan_color(:)) > 1.0
    error 'The argument ''nan_color'' must have elements between zero and one';
end

nan_mask = isnan(input(:));

% Interpolate each pixel intensity value in cmap
if abs(zrange(1) - zrange(2)) > 10*eps
    
    % Threshold high and low values
    input(input < min(zrange)) = min(zrange);
    input(input > max(zrange)) = max(zrange);

    % Interpolate remaining values
    pixels = interp1(linspace(zrange(1), zrange(2), size(cmap,1)), ...
        cmap, input(:));
else
    T = input(:);
    pixels = zeros(length(input(:)), 3);
    T = T < zrange(1);
    
    for idx=1:3
        pixels(T, idx) = cmap(1,idx);
        pixels(~T, :) = cmap(end,idx);
    end
    
    clear T;
end

% Apply nan_color to NaN pixels
pixels(nan_mask, :) = repmat(nan_color, nnz(nan_mask), 1);

% Reshape the pixels into n image
im = reshape(pixels, size(input,1), size(input,2), 3);

end
