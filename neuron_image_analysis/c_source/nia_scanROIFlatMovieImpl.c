/* File: nia_scanROIFlatMovieImpl.c 
 *
 * This function is roughly equivalent to the command
 *    out = squeeze(sum(sum(bsxfun(@times, movie, mask), 1), 2)) / sum(mask(:));
 *
 * However, as of Matlab 2013a, this function runs 5-20X faster
 * depending on inputs and uses half as much memory.
 * 
 *
 * This function can be compiled with one of the following commands:
 *
 * 1) on 64-bit linux with gcc (tested R2013a)
 * mex -v -largeArrayDims CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" nia_scanROIFlatMovieImpl.c
 *
 * 2) on 64-bit linux with icc (tested R2013a)
 * mex -v -largeArrayDims CC=icc CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" LDOPTIMFLAGS="-static-intel" nia_scanROIFlatMovieImpl.c
 *
 * 3) on 64-bit mac with gcc (tested R2013a), after downloading Xcode 4.6.3, mathworks has done an absolutely pathetic job of implementing this (as usual), so getting this to build is incredibly twitchy
 * mex -v -largeArrayDims CC="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" CFLAGS="-std=c99 -DNIA_WORKAROUND_MW_CHAR16_BUG -fno-common -arch x86_64 -isysroot /Applications/Xcode_4.6.3.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk -mmacosx-version-min=10.7  -fexceptions" COPTIMFLAGS="-O3 -fstrict-aliasing" LD="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" nia_scanROIFlatMovieImpl.c
 */

#include <stdint.h>

/* With Matlab 2014a on the mac platform there is a bug that
 * causes compilation to fail because mathworks' mex.h assumes
 * that chart16_t is defined without checking. We define it here
 * so that compilation proceeds.
 */
#ifdef NIA_WORKAROUND_MW_CHAR16_BUG
typedef uint16_t char16_t;
#endif

#include "mex.h"

#define NIA_RESTRICT restrict


/* Include source code for nia_getHashCodeUint64 function.
 * Including a static function is a somewhat lame way of
 * implementing a library, but we need only need a single
 * function, and MATLAB has a really awkward build process,
 * so this is more or less rational.
 */
#include "nia_getHashCodeUInt64_src.c"


/**
 * The following function finds the mean value of the
 * elements of the matrix indicated by the non-zero
 * elements of the mask. The matrix must be a two
 * dimensional array in column major format, with
 * nrows rows and ncols columns. The mask must be
 * a one-dimensional array, with nrows elements.
 * The mean value is stored in the passed output array.
 * This array must be a one-dimensional array with
 * ncols elements.
 *
 * Note that this function uses the Matlab column
 * major convention rather than the C row major
 * convention.
 *
 * \param nrows Number of rows in matrix.
 * \param ncols Number of columns in matrix.
 * \param matrix Matrix of values.
 * \param mask Mask for matrix.
 * \param output Output values.
 * \returns Zero on success, one on error.
 */
static
int
nia_scanROIFlatMovieImpl(
  size_t nrows,
  size_t ncols,
  double const* NIA_RESTRICT matrix,
  mxLogical* mask,
  double* NIA_RESTRICT output)
{
   size_t* nz_elems_vals;
   size_t nz_elems_len;
   size_t nz_elems_alloc;

   size_t idx;

   /* make a list of non-zero elements in mask */
   nz_elems_vals = NULL;
   nz_elems_len = 0;
   nz_elems_alloc = 0;

   for (idx = 0; idx < nrows; idx++) {
      if (mask[idx]) {
         /* resize the list if necessary */
         if (nz_elems_len == nz_elems_alloc) {
            size_t* tmp;

            if (nz_elems_alloc == 0) {
               nz_elems_alloc = 4;
            } else {
               nz_elems_alloc *= 2;
            }

            tmp = mxRealloc(nz_elems_vals,
              nz_elems_alloc * sizeof(size_t));
            if (tmp == NULL) {
               mxFree(nz_elems_vals);
               return 1;
            }

            nz_elems_vals = tmp;
         }

         nz_elems_vals[nz_elems_len] = idx;
         nz_elems_len += 1;
      }
   }

   /* scan through each column */
   for (idx = 0; idx < ncols; idx++) {
      double const* NIA_RESTRICT row_vals;
      double accum;
      size_t idx2;

      row_vals = matrix + nrows * idx;
      accum = 0.0;

      for (idx2 = 0; idx2 < nz_elems_len; idx2++) {
         accum += row_vals[nz_elems_vals[idx2]];
      }

      accum /= nz_elems_len;

      output[idx] = accum;
   }

   if (nz_elems_vals != NULL) {
      mxFree(nz_elems_vals);
   }

   return 0;
}


/**
 * The following function is the gateway from the matlab
 * interface.
 *
 * \params nlhs Number of output arguments
 * \params plhs Array of output arguments
 * \param nrhs Number of input arguments
 * \params prhs Array of input arguments
 */
void mexFunction(
   int nlhs, mxArray *plhs[], 
   int nrhs, const mxArray *prhs[])
{
   const mwSize* movie_dims;
   mwSize movie_dim0;
   mwSize movie_dim1;
   mwSize movie_dim2;
   mwSize movie_dim3;

   const mwSize* mask_dims;
   double* movie_vals;
   uint64_t* chan_idx_vals;
   double* output_vals;
   mxLogical* mask_vals;
   int status = 0;

   /* check input arguments */
   if (nlhs != 1) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovie:nlhs",
        "Must have one output argument");
   }

   if (nrhs != 3) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovie:nrhs",
        "Must have three input arguments");
   }

   if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || 
         mxGetNumberOfDimensions(prhs[0]) < 2 ||
         mxGetNumberOfDimensions(prhs[0]) > 4) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovie:movie",
        "The argument 'movie' must be a 4D matrix of floating-point values");
   }

   if (!mxIsLogical(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovie:mask",
        "The argument 'mask' must be a 2D matrix of logical values");
   }

   if (mxGetClassID(prhs[2]) != mxUINT64_CLASS || mxGetNumberOfDimensions(prhs[2]) != 2 ||
       mxGetNumberOfElements(prhs[2]) != 1) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovie:chan_idx",
        "The argument 'chan_idx' must be a 2D matrix of uint64 values");
   }

   movie_dims = mxGetDimensions(prhs[0]);
   mask_dims = mxGetDimensions(prhs[1]);

   /* expand the movie dimensions to allow for collapsed
    * trailing singleton dimensions.
    */
   movie_dim0 = movie_dims[0];
   movie_dim1 = movie_dims[1];

   if (mxGetNumberOfDimensions(prhs[0]) > 2) {
      movie_dim2 = movie_dims[2];
   } else {
      movie_dim2 = 1;
   }

   if (mxGetNumberOfDimensions(prhs[0]) > 3) {
      movie_dim3 = movie_dims[3];
   } else {
      movie_dim3 = 1;
   }

   /* check that mask has correct size */
   if (movie_dim0 != mask_dims[0] || 
       movie_dim1 != mask_dims[1]) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovie:mask",
        "The dimensions of 'mask' must match dimensions of 'movie'");
   }

   movie_vals = mxGetPr(prhs[0]);
   mask_vals = mxGetLogicals(prhs[1]);
   chan_idx_vals = mxGetData(prhs[2]); 

   /* check that channel index is valid */
   if (*chan_idx_vals < 1 || *chan_idx_vals - 1 > movie_dim3) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovie:chan_idx",
        "The channel exceeds the movie dimensions");
   }

   /* movie the movie pointer to the appropriate channel */
   movie_vals += (*chan_idx_vals - 1) * movie_dim0 * movie_dim1 * movie_dim2;
 
   /* allocate the output array */
   plhs[0] = mxCreateDoubleMatrix(movie_dim2, 1, mxREAL);
   if (plhs[0] == NULL) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovieImpl:internal",
        "An internal error has occurred");
   }

   output_vals = mxGetPr(plhs[0]);

   /* evaluate the ROIs. note that the function we call
    * assumes a 2D input matrix rather than a 3D input
    * matrix. however, all of the data is in the right
    * position, so ignoring this amounts to an implicit
    * reshape of the matrix.
    */
   status = nia_scanROIFlatMovieImpl(movie_dim0 * movie_dim1,
     movie_dim2, movie_vals, mask_vals, output_vals);
   if (status) {
      mexErrMsgIdAndTxt("nia_scanROIFlatMovieImpl:internal",
        "An internal error has occurred");
   }
}


