function out_info = nia_mapChannel(old_val, new_val, in_info)
%NIA_MAPCHANNEL Transforms the input channel info
%   out_info = nia_mapChannel(old_val, new_val, in_info) is a function
%   intended for use as a channel transform for the function
%   nia_movie.process().  When passed as a channel transform this function
%   relabels all channels with the old_val value to the new_val value.
%
%   Example:
%
%       mov.process([], [], @(info) nia_map_channel(2, 1, info));

if in_info.ch == old_val
    out_info.ch = new_val;
else
    out_info.ch = in_info.ch;
end

end

