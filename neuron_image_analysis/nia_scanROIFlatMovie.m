function out = nia_scanROIFlatMovie(mov, mask, chan)
%NIA_SCANROIFLATMOVIE Find mean value within roi across frames
%   out = nia_scanROIFlatMovie(mov, mask) finds the mean value of the
%   pixels within the passed ROI for each frame within the passed movie.
%   The movie must be passed as a 4D array in which the first two
%   dimensions corresponds to the width and height, the third dimension
%   corresponds to frames, and the fourth dimension corresponds to
%   channels.  Trailing singleton dimensions may be removed, so if the
%   movie maybe represented as a three dimensional matrix if it contains
%   only a single channel, and as a two dimensional matrix if it contains
%   onaly a single frame. The mask must be a two dimensional logical
%   array with the same dimensions as the first two dimensions of movie.
%   The channel must a be a scalar integer index and is optional. The
%   output is a column vector with a length equal to the length of
%   the third dimension of movie. 

% This function exists solely to document the mex-function
% that does all of the actual work, and to provide a very
% small amount of input processing.

    if nargin < 3
        if ~nia_isScalarInteger(chan) || chan < 1 || chan > size(mov,4)
            error 'The argument ''chan'' must be a valid index';
        end
    else
        chan = 1;
    end
        
    chan_idx = uint64(chan);

    out = nia_scanROIFlatMovieImpl(mov, mask, chan_idx);
end
