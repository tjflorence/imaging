function importStruct(obj, sdata)
%IMPORTSTRUCT Import structure as movie
%   importStruct(sdata) loads the structure data into the movie class. This
%   structure should have the same format as that returned by the function
%   exportStruct(). Note that for efficiency this function performs no
%   error checking on the input structure. Invalid inputs are highly likely
%   to produce incorrect results.

obj.slices = sdata.slices;
obj.ch_list = sdata.ch_list;
obj.pos_lens = sdata.pos_lens;
obj.pos_lu = sdata.pos_lu;
obj.pos_ranges = sdata.pos_ranges;
obj.hist_info = sdata.hist_info;

% imported structure is current, so we leave
% obj.api_ver at coded value

obj.safe_mode = sdata.safe_mode;

end

