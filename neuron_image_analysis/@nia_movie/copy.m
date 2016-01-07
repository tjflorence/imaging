function out = copy(obj)
%COPY Create a copy of the movie
%   out = copy() creates a deep copy of the passed movie. Since nia_movie
%   is a handle object, when it is assigned to multiple variables all
%   variables act on the same underlying data. To create a completely
%   independent copy of the data, you must call this function instead. Note
%   that because MATLAB uses copy on write semantics, the actual data copy
%   may be delayed until the data is actually modified.

out = nia_movie;
out.slices = obj.slices;
out.ch_list = obj.ch_list;
out.pos_lens = obj.pos_lens;
out.pos_lu = obj.pos_lu;
out.pos_ranges = obj.pos_ranges;
out.hist_info = obj.hist_info;
end

