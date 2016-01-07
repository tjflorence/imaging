function out = nia_isColor(input)
%NIA_ISCOLOR Return true if input is a valid ColorSpec color.

if isfloat(input)
    % Check as triplet specfication
    if ~isreal(input) || ~isvector(input) || length(input) ~= 3 || ...
            min(input) < 0 || max(input) > 1
        out = false;
        return
    end
    
    out = true;
    return;
else
    % Check as named color
    if ~ischar(input) || ~ismatrix(input) || size(input,1) > 1
        out = false;
        return;
    end
    
    known_colors = {'blue', 'red', 'green', 'yellow', 'magenta', ...
        'cyan', 'white', 'black'};
    is_known = nnz(cellfun(@(x) strcmp(x, input), known_colors)) ~= 0;
    
    out = is_known;
    return;
end

end