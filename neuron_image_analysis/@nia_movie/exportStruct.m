function sdata = exportStruct(obj)
%EXPORTSTRUCT Export movie as structure
%   sdata = exportStruct() creates a structure that contains all of the
%   movie data. This structure may be subsequently imported using the
%   importStruct() function.

% You're not the boss of me, matlab
s = warning('off', 'MATLAB:structOnObject');

% Convert the structure
sdata = struct(obj);

% Return to previous warning state
warning(s);

end

