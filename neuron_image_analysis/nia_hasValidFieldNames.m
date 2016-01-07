function [out, msg] = nia_hasValidFieldNames(input, allowed, required)
%NIA_HASVALIDFIELDNAMES Return true if structure has valid form.
%   out = nia_hasValidFieldNames(input, allowed, required) scans
%   the argument 'input' to determine if it contains only fields that
%   are in the 'allowed' field names cell array, and has all of the field
%   names in the 'required' field names cell array. If it does then this
%   function returns true, otherwise this function returns false. If the
%   structure is invalid, then msg is set to a string that can be used
%   to generate a specific error message by calling the function
%   error(msg, INPUT_NAME). If the structure is valid, then msg is an
%   empty array.

% Check input arguments
if ~isstruct(input)
    error 'The argument ''input'' must be a struct';
end

if ~iscell(allowed) || ~isvector(allowed)
    error 'The argument ''allowed'' must be a cell array';
end

for idx=1:length(allowed)
    if ~nia_isString(allowed{idx})
        error 'The argument ''allowed'' must be a cell array of strings';
    end
end

if ~iscell(required) || ~(isempty(required) || isvector(required))
    error 'The argument ''required'' must be a cell array';
end

for idx=1:length(required)
    if ~nia_isString(required{idx})
        error 'The argument ''required'' must be a cell array of strings';
    end
end

% Create allowed structure
allowed_struct = cell2struct(cell(1, length(allowed)), allowed, 2);

% Check if all fields in input are allowed
names = fieldnames(input);
for idx=1:length(names)
    if ~isfield(allowed_struct, names{idx})
        out = false;
        msg = ['The argument ''%s'' contains an invalid field ''', names{idx}, ''''];
        return;
    end
end

% Check if all fields in required are in input
for idx=1:length(required)
    if ~isfield(input, required{idx})
        out = false;
        msg = ['The argument ''%s'' is missing the field ''', required{idx}, ''''];
        return;
    end
end

out = true;
msg = [];

end
