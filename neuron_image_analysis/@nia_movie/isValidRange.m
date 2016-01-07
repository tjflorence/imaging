function [out, msg] = isValidRange(obj, rnge)
%ISVALIDRANGE Check if the passed variable is a valid range
%   out = isValidRange(rnge) checks if the passed argument is a valid
%   range for the given movie. It returns true if it is a valid
%   range and false otherwise. If the input is not a valid range, then
%   the ouput 'msg' can be used to generate a meaningful error using
%   the call error(msg, INPUT_NAME). If the input is valid, then the
%   output 'msg' is an empty array.
%
%       rnge - The range to test. The range may be an empty array,
%           indicating that all slices and channels are accepted. If the
%           range is not empty, it should be structure with one or more of
%           the following fields:
%
%           time - A [2x1] array that specifies the times included in the
%               range. The first element specifies the start time and the
%               second element the end time, both are given inclusively. A
%               NaN value indicates that the relevant end of the range is
%               unbounded.  If an empty array is passed, then no time
%               selection is applied.
%
%           pos_vec - A scalar integer that specifies the position vector
%               that should be used when applying the position selector. If
%               this argument is omitted or set to the empty array, then
%               this value default to one.
%
%           pos_lims - A [2xN] array that specifies the positions included
%               in the range, where N is the length of the movie position
%               vector. The first row of the 'pos' field specifies the
%               minimum position value in the range, and the second row of
%               the 'pos' field specifies the maximum value in the range,
%               both are given inclusively. A NaN value indicates that the
%               relevant end of the range is unbounded.  If an empty array
%               is passed, then no position selection is applied.
%
%           channels - A [1xM] array that specifies the channels included
%               in the range. The array should contain the identifier of
%               every selected channel. If an empty array is passed, then
%               no channels selection is applied.


% Check if input is a valid structure
if isempty(rnge)
    out = true;
    msg = [];
    return;
end

if ~isstruct(rnge) || length(rnge) ~= 1
    out = false;
    msg = 'The argument ''%s'' must be a scalar structure';
    return;
end

rnge_names = {'time', 'pos_vec', 'pos_lims', 'channels'};
[rnge_ok, rnge_msg] = nia_hasValidFieldNames(rnge, rnge_names, {});

if ~rnge_ok
    out = rnge_ok;
    msg = rnge_msg;
    return;
end

% Check the time field
if isfield(rnge, 'time') && ~isempty(rnge.time);
    if ~isfloat(rnge.time) || ~isreal(rnge.time) || ...
            ~iscolumn(rnge.time) || size(rnge.time, 1) ~= 2
        out = false;
        msg = 'The argument ''%s.time'' has an invalid type';
        return;
    end
    
    % Note behavior if either value is NaN
    if rnge.time(1) > rnge.time(2)
        out = false;
        msg = 'The argument ''%s.time'' has an invalid value';
        return;
    end
end

% Check the pos_vec field
if isfield(rnge, 'pos_vec') && ~isempty(rnge.pos_vec)
    if isempty(obj.pos_lens)
        out = false;
        msg = 'The argument ''%s.pos_vec'' has an invalid value';
        return;
    end
    
    if ~nia_isScalarInteger(rnge.pos_vec)
        out = false;
        msg = 'The argument ''%s.pos_vec'' has an invalid value';
        return;
    end
    
    if rnge.pos_vec < 1 || rnge.pos_vec > length(obj.pos_lens)
        out = false;
        msg = 'The argument ''%s.pos_vec'' has an invalid value';
        return;
    end
    
    effective_pos_vec = rnge.pos_vec;
else
    if ~isempty(obj.pos_lens)
        effective_pos_vec = 1;
    else
        effective_pos_vec = [];
    end
end

% Check the pos_lims field
if isfield(rnge, 'pos_lims') && ~isempty(rnge.pos_lims)
    if isempty(obj.pos_lens)
        out = false;
        msg = 'The argument ''%s.pos_lims'' has an invalid value';
        return;
    end
    
    if ~isfloat(rnge.pos_lims) || ~isreal(rnge.pos_lims) || ...
            ~ismatrix(rnge.pos_lims) || size(rnge.pos_lims, 1) ~= 2 || ...
            size(rnge.pos_lims, 2) ~= obj.pos_lens(effective_pos_vec)
        out = false;
        msg = 'The argument ''%s.pos_lims'' has an invalid type';
        return;
    end
    
    % Note behavior if either value is NaN
    if nnz(rnge.pos_lims(:) < 0) > 0 || ...
            nnz(rnge.pos_lims(1,:) > rnge.pos_lims(2,:)) > 0
        
        out = false;
        msg = 'The argument ''%s.pos_lims'' has an invalid value';
        return;
    end
end

% Check the channels field
if isfield(rnge, 'channels') && ~isempty(rnge.channels)
    if ~isfloat(rnge.channels) || ~isreal(rnge.channels) || ...
            ~isrow(rnge.channels)
        out = false;
        msg = 'The argument ''%s.channels'' has an invalid type';
        return;
    end
end

out = true;
msg = [];

end
