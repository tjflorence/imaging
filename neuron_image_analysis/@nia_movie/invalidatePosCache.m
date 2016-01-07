function invalidatePosCache(obj)
%INVALIDATEPOSCACHE Mark current position as requiring an update
%   invalidatePosCache() marks the current position as out of date. The
%   position ranges and lookup cache will be updated the next time that
%   they are accessed.

obj.pos_lu = {};
obj.pos_ranges = {};

end
