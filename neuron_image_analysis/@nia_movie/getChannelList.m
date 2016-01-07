function ch_list = getChannelList(obj)
%GETCHANNELLIST Get the list of channel identifiers
%   ch_list = getChannelList() retrieves the list of channel identifiers
%   that are referenced by the movie. The position list is returned as a
%   row vector where each element contains a channel identifier number.
%   This information is cached, so that the first call to this function may
%   require some time to calculate the channel list, but subsequent calls
%   will return the stored list quickly. The cache can be manually loaded
%   using the function loadChInfo(). Some modifications to the movie may
%   cause the list to become out-of-date, at which point the cached list
%   will be recomputed.

if isempty(obj.ch_list)
    obj.loadChCache();
end

ch_list = obj.ch_list;
end