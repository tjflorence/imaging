/* File: nia_scanToLastNiaMovieImpl.c 
 *
 * This function is used to implement nia_movie/scanToLast
 *
 * This function can be compiled with one of the following commands:
 *
 * 1) on 64-bit linux with gcc (tested R2013a)
 * mex -v -largeArrayDims CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" nia_scanToLastNiaMovieImpl.c
 *
 * 2) on 64-bit linux with icc (tested R2013a)
 * mex -v -largeArrayDims CC=icc CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" LDOPTIMFLAGS="-static-intel" nia_scanToLastNiaMovieImpl.c
 *
 * 3) on 64-bit mac with gcc (tested R2013a), after downloading Xcode 4.6.3, mathworks has done an absolutely pathetic job of implementing this (as usual), so getting this to build is incredibly twitchy
 * mex -v -largeArrayDims CC="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" CFLAGS="-std=c99 -DNIA_WORKAROUND_MW_CHAR16_BUG -fno-common -arch x86_64 -isysroot /Applications/Xcode_4.6.3.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk -mmacosx-version-min=10.7  -fexceptions" COPTIMFLAGS="-O3 -fstrict-aliasing" LD="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" nia_scanToLastNiaMovieImpl.c
 */

#include <stdint.h>
#include <math.h>

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


/** The following function scans through the movie to find
 * the last position in the given range that has a valid
 * entry. The format and type of the passed variables is given
 * by the nia_movie object.  Specifically, slices must be a row
 * structure array with the fields 'time', 'channels', and
 * 'pos'. The argument pos_lu must be a hash stored as a
 * cell array where each element contains a [Mx(N+1)] array
 * of doubles, where M is the number slices with the given
 * hash code and N is the length of the position vector.
 * The home_pos is the starting position vector and must be
 * a row vector of uint64 values. The argument scan_extents
 * must be a column vector with two elements, where the
 * first is the index of the position vector to iterate and
 * the second is the last position value. The error checking in
 * the mexFunction that invokes this function checks only these
 * features of the variables. Evaluation of nested variables
 * should confirm that they have appropriate types and values.
 * The result is stored in the argument out, it is assumed to
 * be uninitialized when this function is invoked.
 *
 * \param slices Slices matrix.
 * \param pos_lu Position lookup.
 * \param home_pos Home position.
 * \param scan_extents Scan extents.
 * \param out Output value.
 * \returns Zero on success, one on error.
 */
int
nia_scanToLastNiaMovieImpl(
  const mxArray* slices,
  const mxArray* pos_lu,
  const mxArray* home_pos,
  const mxArray* scan_extents,
  mxArray** out)
{
   uint64_t* home_pos_vals;
   size_t home_pos_len;

   uint64_t* scan_extents_vals;

   size_t scan_idx;
   uint64_t scan_final_val;
   uint64_t scan_cur_val;
   uint64_t* scan_pos = NULL;
   size_t scan_init_idx;

   int status = 0;

   /* initialize the output */
   *out = NULL;

   /* retrieve the home position value */
   home_pos_vals = (uint64_t*)mxGetData(home_pos);
   home_pos_len = mxGetNumberOfElements(home_pos);

   if (home_pos_len == 0) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:home_pos",
        "The argument 'home_pos' has an invalid value");
   }

   /* retrieve the scan extents values */
   scan_extents_vals = (uint64_t*)mxGetData(scan_extents);

   /* check if scan extents are valid */
   if (scan_extents_vals[0] < 1 || scan_extents_vals[0] > home_pos_len) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:scan_extents",
        "The argument 'scan_extents(1)' has an invalid value");
   }

   scan_idx = scan_extents_vals[0] - 1;

   if (scan_extents_vals[1] < home_pos_vals[scan_idx]) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:scan_extents",
        "The argument 'scan_extents(2)' has an invalid value");
   }

   scan_final_val = scan_extents_vals[1];

   /* prepare the iterated position */
   scan_pos = mxMalloc(home_pos_len * sizeof(uint64_t));
   if (scan_pos == NULL) {
      status = 1;
      goto finished;
   }

   for (scan_init_idx = 0; scan_init_idx < home_pos_len; scan_init_idx++) {
      scan_pos[scan_init_idx] = home_pos_vals[scan_init_idx];
   }

   /* traverse positions backwards looking for the first valid element.
    * we do the loop check within the for loop since scan_cur_val is an
    * unsigned integer and we do not want to wrap it if
    * home_pos_vals[scan_idx] is equal to zero. in this logic the iterator
    * starts the loop set to the index of the element that was found to
    * be missing in the previous iteration.  when this loop terminates
    * scan_cur_val should contain the index immediately after the last
    * valid position in the scan, the start position if no valid
    * positions are found, or one plus the last position if all positions
    * are valid.
    */
   scan_cur_val = scan_final_val + 1;
   for (;;) {
      uint64_t hash;
      const mxArray* bin;
      double* bin_vals;
      size_t bin_nrows;
      size_t bin_ncols;
      size_t bin_row_idx;

      /* check if scan terminated without finding a single
       * valid position.
       */
      if (scan_cur_val == home_pos_vals[scan_idx]) {
         break;
      }

      scan_cur_val--;
 
      /* create the current position vector */
      scan_pos[scan_idx] = scan_cur_val;

      /* get the hash code for the current position */
      hash = nia_getHashCodeUInt64(home_pos_len, scan_pos);

      /* get bin and check type */
      bin = mxGetCell(pos_lu, hash % mxGetNumberOfElements(pos_lu));
     
      if (bin == NULL) {
         continue;
      }
 
      if (!mxIsDouble(bin) || mxIsComplex(bin) ||
          mxGetNumberOfDimensions(bin) != 2) {
         mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:pos_lu",
           "The argument 'pos_lu' has an invalid value");
      }

      /* if bin is empty, no slice was found */
      if (mxGetNumberOfElements(bin) == 0) {
         continue;
      }

      /* if bin is not empty, check that it has correct dimensions */
      if (mxGetN(bin) != home_pos_len+1) {
         mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:pos_lu",
           "The argument 'pos_lu' has an invalid value");
      }

      /* scan through bin looking for matching row */
      bin_nrows = mxGetM(bin);
      bin_ncols = mxGetN(bin);
      bin_vals = mxGetPr(bin);

      for (bin_row_idx = 0; bin_row_idx < bin_nrows; bin_row_idx++) {
         size_t bin_col_idx;

         for (bin_col_idx = 1; bin_col_idx < bin_ncols; bin_col_idx++) {
            double bin_flt;
            uint64_t bin_int;

            bin_flt = bin_vals[bin_row_idx + bin_col_idx * bin_nrows];
            bin_int = floor(bin_flt + 0.5);

            if (bin_int != scan_pos[bin_col_idx-1]) {
               break;
            }
         }

         /* if we made it through the above for-loop, then it is
          * a match
          */
         if (bin_col_idx == bin_ncols) {
            break;
         }
      }

      /* if we made through the above for-loop, then there was
       * no match
       */
      if (bin_row_idx == bin_nrows) {
         continue;
      }

      /* we did get a match, stopping search */
      scan_cur_val += 1;
      break;
   }

   *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
   *(uint64_t*)mxGetData(*out) = scan_cur_val;

finished:
   if (scan_pos != NULL) {
      mxFree(scan_pos);
   }

   return status;
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
   const mxArray* slices;
   const mxArray* pos_lu;
   const mxArray* home_pos;
   const mxArray* scan_extents;
   mxArray* out;

   int status = 0;

   /* check input arguments */
   if (nlhs != 1) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:nlhs",
        "Must have one output argument");
   }

   if (nrhs != 4) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:nrhs",
        "Must have seven input arguments");
   }

   slices = prhs[0];
   pos_lu = prhs[1];
   home_pos = prhs[2];
   scan_extents = prhs[3];

   if (!mxIsStruct(slices) || 
       mxGetNumberOfDimensions(slices) != 2 ||
       mxGetM(slices) != 1 ||
       mxGetNumberOfFields(slices) != 3 ||
       mxGetFieldNumber(slices, "time") == -1 ||
       mxGetFieldNumber(slices, "channels") == -1 ||
       mxGetFieldNumber(slices, "pos") == -1) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:slices",
        "The argument 'slices' has an invalid type");
   }

   if (!mxIsCell(pos_lu) ||
       mxGetNumberOfDimensions(pos_lu) != 2 ||
       mxGetM(pos_lu) != 1 ||
       mxGetN(pos_lu) < 1) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:pos_lu",
        "The argument 'pos_lu' has an invalid type");
   }

   if (mxGetClassID(home_pos) != mxUINT64_CLASS ||
       mxGetNumberOfDimensions(home_pos) != 2 ||
       mxGetM(home_pos) != 1) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:mask",
        "The argument 'home_pos' has an invalid type");
   }

   if (mxGetClassID(scan_extents) != mxUINT64_CLASS ||
       mxGetNumberOfDimensions(scan_extents) != 2 ||
       mxGetM(scan_extents) != 2 ||
       mxGetN(scan_extents) != 1) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:scan_extents",
        "The argument 'scan_extents' has an invalid type");
   }

   /* scan the movie */
   status = nia_scanToLastNiaMovieImpl(slices, pos_lu,
     home_pos, scan_extents, &out);
   if (status) {
      mexErrMsgIdAndTxt("nia_scanToLastNiaMovieImpl:internal",
        "An internal error has occurred");
      return;
   }

   plhs[0] = out;
}


