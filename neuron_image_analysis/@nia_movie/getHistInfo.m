function hist_info = getHistInfo(obj)
%GETHISTINFO Retrieve histogram information
%   hist_info = getHistInfo() retrieves a structure array containing
%   information on this histograms for every channel of the movie. This
%   information is cached, so that the first call to this function may
%   require some time to calculate the histogram, but subsequent calls will
%   return the stored histogram quickly. The cache can be manually loaded
%   using the function loadHistInfo(). Some modifications to the movie may
%   cause the histogram to become out-of-date, at which point the cached
%   histogram will be recomputed.
%
%   The returned value is a structure array that specifies the histogram
%   information for each channel of the movie. It has a number of elements
%   equal to the number of channels each with the fields 'min', 'max', and
%   'dist'. The 'min' field is scalar storing the minimum pixel intensity
%   value, the 'max' field is a scalar storing the maximum pixel intensity
%   value, and 'dist' field is a [2xN] array specifying the histogram count
%   for a series of evenly spaced bins. Each column of 'dist' field
%   describes a bin, with the center position of the bin in the first row,
%   and the count in the second row.

if isempty(obj.hist_info)
    obj.loadHistCache();
end

hist_info = obj.hist_info;
end