function rnge = normalizeRange(obj, rnge)
%NORMALIZERANGE Normalize the passed range
%   rnge = normalizeRange(rnge) normalizes the passed range by adding any
%   fields that have been omitted or left empty with NaN ranges.
%
%   This function accepts the following arguments:
%
%       rnge - The range of slices to to use for the test. See the 
%           function isValidRange() for a description of the format
%           of the range.
%
%   Note that for efficiency this function performs no error checking
%   on the input arguments.  The input mus pass isValidRange().

if isempty(rnge)
    rnge.time = [];
    rnge.pos_vec = [];
    rnge.pos_lims = [];
    rnge.channels = [];
end

if ~isfield(rnge, 'time') || isempty(rnge.time)
    rnge.time = [NaN; NaN];
end

if ~isfield(rnge, 'pos_vec') || isempty(rnge.pos_vec)
    if ~isempty(obj.pos_lens)
        rnge.pos_vec = 1;
    else
        rnge.pos_vec = [];
    end
end

if ~isfield(rnge, 'pos_lims') || isempty(rnge.pos_lims)
    if ~isempty(obj.pos_lens)
        rnge.pos_lims = zeros(2, obj.pos_lens(rnge.pos_vec)) * NaN;
    else
        rnge.pos_lims = [];
    end
end

if ~isfield(rnge, 'channels') || isempty(rnge.channels)
    obj.loadChCache();
    rnge.channels = obj.ch_list;
end
        
end
