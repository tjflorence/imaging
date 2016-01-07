/* File: nia_scanROINiaMovieImpl.c 
 *
 * This function is used to implement nia_movie/scanROI
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
 * mex -v -largeArrayDims CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" nia_scanROINiaMovieImpl.c
 *
 * 2) on 64-bit linux with icc (tested R2013a)
 * mex -v -largeArrayDims CC=icc CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" LDOPTIMFLAGS="-static-intel" nia_scanROINiaMovieImpl.c
 *
 * 3) on 64-bit mac with gcc (tested R2013a), after downloading Xcode 4.6.3, mathworks has done an absolutely pathetic job of implementing this (as usual), so getting this to build is incredibly twitchy
 * mex -v -largeArrayDims CC="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" CFLAGS="-std=c99 -DNIA_WORKAROUND_MW_CHAR16_BUG -fno-common -arch x86_64 -isysroot /Applications/Xcode_4.6.3.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk -mmacosx-version-min=10.7  -fexceptions" COPTIMFLAGS="-O3 -fstrict-aliasing" LD="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" nia_scanROINiaMovieImpl.c
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


/**
 * The following function scans through the passed matrix of
 * logical values to create an array of non-zero values.
 * 
 * \param mask Logicals mask.
 * \param mask_len Number of elements in mask.
 * \param nz_elems_vals_ptr Output non-zero elements indices
 * \param nz_elems_len_ptr Length of output non-zero elements array
 * \returns Zero on success, an error code on failure.
 */
static
int nia_getNZElems(
  mxLogical* mask,
  size_t mask_len,
  size_t** nz_elems_vals_ptr,
  size_t* nz_elems_len_ptr)
{
   size_t* nz_elems_vals;
   size_t nz_elems_len;
   size_t nz_elems_alloc;

   size_t idx;

   /* make a list of non-zero elements in mask */
   nz_elems_vals = NULL;
   nz_elems_len = 0;
   nz_elems_alloc = 0;

   for (idx = 0; idx < mask_len; idx++) {
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

   /* trim the list to a minimal size */
   if (nz_elems_len > 0) {
      size_t* tmp;

      nz_elems_alloc = nz_elems_len;
      tmp = mxRealloc(nz_elems_vals,
        nz_elems_alloc * sizeof(size_t));
      if (tmp == NULL) {
         mxFree(nz_elems_vals);
         return 1;
      }
   }

   *nz_elems_vals_ptr = nz_elems_vals;
   *nz_elems_len_ptr = nz_elems_len;
   return 0;
}


/**
 * The following function finds the mean value of the
 * elements of the movie indicated by the non-zero
 * elements of the mask. The format and type of the
 * passed variables is given by the nia_movie object.
 * Specifically, slices must be a row structure array
 * with the fields 'time', 'channels', and 'pos'.
 * The argument pos_lu must be a hash stored as a cell
 * array where each element contains a [Mx(N+1)] array of
 * doubles, where M is the number slices with the given
 * hash code and N is the length of the position vector.
 *  The home_pos is the starting position vector
 * and must be a row vector of uint64 values. The argument
 * scan_extents must be a column vector with two elements,
 * where the first is the index of the position vector to
 * iterate and the second is the last position value. The
 * argument channel is a scalar uint64 value that holds the
 * channel identifier to integrate. The mask is a
 * N-dimensional array of logicals that must be the same
 * size as the image data at the home * position. The error
 * checking in the mexFunction that invokes this function
 * checks only these features of the variables. Evaluation
 * of nested variables should confirm that they have
 * appropriate types and values. The result is stored in
 * the argument out, it is assumed to be uninitialized when
 * this function is invoked.
 *
 * \param slices Slices matrix.
 * \param pos_lu Position lookup.
 * \param home_pos Home position.
 * \param scan_extents Scan extents.
 * \param channel_id Evaluated channel.
 * \param mask Image mask.
 * \returns Zero on success, one on error.
 */
int
nia_scanROINiaMovieImpl(
  const mxArray* slices,
  const mxArray* pos_lu,
  const mxArray* home_pos,
  const mxArray* scan_extents,
  const mxArray* channel_id,
  const mxArray* mask,
  mxArray** out)
{
   size_t* nz_elems_vals;
   size_t nz_elems_len;

   uint64_t* home_pos_vals;
   size_t home_pos_len;

   uint64_t* scan_extents_vals;

   size_t scan_idx;
   uint64_t scan_final_val;
   uint64_t scan_cur_val;
   uint64_t* scan_pos = NULL;
   size_t scan_init_idx;

   uint64_t desired_ch_id;

   size_t mask_ndims;
   const mwSize* mask_dims;

   double* out_vals;
   size_t out_idx;

   int status = 0;

   /* create a list of non-zero elements in the mask */
   status = nia_getNZElems(mxGetLogicals(mask),
     mxGetNumberOfElements(mask),
     &nz_elems_vals, &nz_elems_len);
   if (status) {
      return status;
   }

   /* retrieve the home position value */
   home_pos_vals = (uint64_t*)mxGetData(home_pos);
   home_pos_len = mxGetNumberOfElements(home_pos);

   if (home_pos_len == 0) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:home_pos",
        "The argument 'home_pos' has an invalid value");
   }

   /* retrieve the scan extents values */
   scan_extents_vals = (uint64_t*)mxGetData(scan_extents);

   /* check if scan extents are valid */
   if (scan_extents_vals[0] < 1 || scan_extents_vals[0] > home_pos_len) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:scan_extents",
        "The argument 'scan_extents(1)' has an invalid value");
   }

   scan_idx = scan_extents_vals[0] - 1;

   if (scan_extents_vals[1] < home_pos_vals[scan_idx]) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:scan_extents",
        "The argument 'scan_extents(2)' has an invalid value");
   }

   scan_final_val = scan_extents_vals[1];

   /* retrieve the mask dimensions */
   mask_ndims = mxGetNumberOfDimensions(mask);
   mask_dims = mxGetDimensions(mask);

   /* retrieve the channel identifier */
   desired_ch_id = *(uint64_t*)mxGetData(channel_id);

   /* prepare the iterated position */
   scan_pos = mxMalloc(home_pos_len * sizeof(uint64_t));
   if (scan_pos == NULL) {
      status = 1;
      goto finished;
   }

   for (scan_init_idx = 0; scan_init_idx < home_pos_len; scan_init_idx++) {
      scan_pos[scan_init_idx] = home_pos_vals[scan_init_idx];
   }

   /* allocate the output array */
   *out = mxCreateDoubleMatrix(3, scan_final_val - home_pos_vals[scan_idx] + 1,
     mxREAL);
   if (*out == NULL) {
      status = 1;
      goto finished;
   }

   out_vals = (double*)mxGetPr(*out);
   out_idx = 0;

   /* begin traversing the positions */
   for (scan_cur_val = home_pos_vals[scan_idx];
        scan_cur_val <= scan_final_val;
        scan_cur_val++, out_idx++) {

      uint64_t hash;
      const mxArray* bin;
      double* bin_vals;
      size_t bin_nrows;
      size_t bin_ncols;
      size_t bin_row_idx;

      double slice_idx_flt;
      size_t slice_idx;
      const mxArray* time;
      double time_flt;
      const mxArray* channels;
      int ch_field_number;
      int image_field_number;
      size_t num_chs;
      size_t ch_idx;
      const mwSize* image_dims;
      double* image_vals;
      size_t pos_idx;
      double accum;
      const mxArray* image;
      size_t mask_idx;
 
      /* create the current position vector */
      scan_pos[scan_idx] = scan_cur_val;

      /* get the hash code for the current position */
      hash = nia_getHashCodeUInt64(home_pos_len, scan_pos);

      /* get bin and check type */
      bin = mxGetCell(pos_lu, hash % mxGetNumberOfElements(pos_lu));

      if (bin == NULL) {
         out_vals[3*out_idx] = NAN;
         out_vals[3*out_idx+1] = NAN;
         out_vals[3*out_idx+2] = NAN;
         continue;
      }
      
      if (!mxIsDouble(bin) || mxIsComplex(bin) ||
          mxGetNumberOfDimensions(bin) != 2) {
         mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:pos_lu",
           "The argument 'pos_lu' has an invalid value");
      }

      /* if bin is empty, no slice was found */
      if (mxGetNumberOfElements(bin) == 0) {
         out_vals[3*out_idx] = NAN;
         out_vals[3*out_idx+1] = NAN;
         out_vals[3*out_idx+2] = NAN;
         continue;
      }

      /* if bin is not empty, check that it has correct dimensions */
      if (mxGetN(bin) != home_pos_len+1) {
         mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:pos_lu",
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
         out_vals[3*out_idx] = NAN;
         out_vals[3*out_idx+1] = NAN;
         out_vals[3*out_idx+2] = NAN;
         continue;
      }

      /* we did get a match, so take the slice index from first column */
      slice_idx_flt = bin_vals[bin_row_idx];

      if (slice_idx_flt < 0) {
         mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:pos_lu",
           "The argument 'pos_lu' has an invalid value");
      }

      slice_idx = floor(slice_idx_flt + 0.5);

      if (slice_idx < 1 || slice_idx > mxGetNumberOfElements(slices)) {
         mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:pos_lu",
           "The argument 'pos_lu' has an invalid value");
      }

      /* use zero-based indexing */
      slice_idx -= 1;

      /* retrieve the structure information */
      time = mxGetField(slices, slice_idx, "time");
      if (time != NULL) {
         if (!mxIsDouble(time) || 
             mxGetNumberOfDimensions(time) != 2 ||
             mxGetNumberOfElements(time) != 1) {

            mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:slices",
              "The argument 'slices' has an invalid value");
         }

         time_flt = *mxGetPr(time);
      } else {
         time_flt = NAN;
      }

      channels = mxGetField(slices, slice_idx, "channels");
      if (channels == NULL) {
         out_vals[3*out_idx] = NAN;
         out_vals[3*out_idx+1] = NAN;
         out_vals[3*out_idx+2] = NAN;
         continue;
      }

      if (!mxIsStruct(channels) || 
          mxGetNumberOfDimensions(channels) != 2 ||
          mxGetM(channels) != 1 ||
          mxGetNumberOfFields(channels) != 2 ||
	  (ch_field_number = mxGetFieldNumber(channels, "ch")) == -1 ||
	  (image_field_number = mxGetFieldNumber(channels, "image")) == -1) {
         mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:slices",
           "The argument 'slices' has an invalid type");
      }

      /* scan for a matching channel */
      num_chs = mxGetNumberOfElements(channels); 
      for (ch_idx = 0; ch_idx < num_chs; ch_idx++) {
         const mxArray* ch_entry;
         double ch_val;
         uint64_t ch_id;

         ch_entry = mxGetFieldByNumber(channels, ch_idx, ch_field_number);

         /* we skimp on the error checking just a little bit, and only
          * check what we need to assure that the code doesn't crash.
          * this is done for performance, and because errors in these
          * entries would be internal errors elsewhere.
          */
         if (ch_entry == NULL || !mxIsDouble(ch_entry) ||
             mxGetNumberOfElements(ch_entry) != 1) {
            mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:slices",
              "The argument 'slices' has an invalid type");
         }

         ch_val = *mxGetPr(ch_entry);
         ch_id = floor(ch_val + 0.5);

         if (desired_ch_id == ch_id) {
            break;
         }
      }

      /* if we made it to the end of the above for-loop, then we
       * didn't find a matching channel.
       */
      if (ch_idx == num_chs) {
         out_vals[3*out_idx] = NAN;
         out_vals[3*out_idx+1] = NAN;
         out_vals[3*out_idx+2] = NAN;
         continue;
      }

      /* retrieve the associated image data */
      image = mxGetFieldByNumber(channels, ch_idx, image_field_number);

      /* skip image data with incorrect size */
      if (image == NULL ||
          !mxIsDouble(image) ||
          mxGetNumberOfDimensions(image) != mask_ndims) {
         out_vals[3*out_idx] = NAN;
         out_vals[3*out_idx+1] = NAN;
         out_vals[3*out_idx+2] = NAN;
         continue;
      }

      /* check that image has correct dimensions */
      image_dims = mxGetDimensions(image);

      for (pos_idx = 0; pos_idx < mask_ndims; pos_idx++) {
         if (mask_dims[pos_idx] != image_dims[pos_idx]) {
            break;
         }
      }

      /* if we didn't make it to the end of the above for
       * loop then the dimensions didn't match, skip frame
       */
      if (pos_idx != mask_ndims) {
         out_vals[3*out_idx] = NAN;
         out_vals[3*out_idx+1] = NAN;
         out_vals[3*out_idx+2] = NAN;
         continue;
      }

      /* okay, we have a frame, let's find the mean value 
       * for the non-zero elements of the mask
       */
      accum = 0.0;
      image_vals = mxGetPr(image);

      for (mask_idx = 0; mask_idx < nz_elems_len; mask_idx++) {
         accum += image_vals[nz_elems_vals[mask_idx]];
      }

      accum /= nz_elems_len;
      out_vals[3*out_idx] = slice_idx_flt;
      out_vals[3*out_idx+1] = time_flt;
      out_vals[3*out_idx+2] = accum;
   }

finished:
   if (nz_elems_vals != NULL) {
      mxFree(nz_elems_vals);
   }

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
   const mxArray* mask;
   const mxArray* home_pos;
   const mxArray* scan_extents;
   const mxArray* channel_id;
   mxArray* out;

   int status = 0;

   /* check input arguments */
   if (nlhs != 1) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:nlhs",
        "Must have one output argument");
   }

   if (nrhs != 6) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:nrhs",
        "Must have six input arguments");
   }

   slices = prhs[0];
   pos_lu = prhs[1];
   home_pos = prhs[2];
   scan_extents = prhs[3];
   channel_id = prhs[4];
   mask = prhs[5];

   if (!mxIsStruct(slices) || 
       mxGetNumberOfDimensions(slices) != 2 ||
       mxGetM(slices) != 1 ||
       mxGetNumberOfFields(slices) != 3 ||
       mxGetFieldNumber(slices, "time") == -1 ||
       mxGetFieldNumber(slices, "channels") == -1 ||
       mxGetFieldNumber(slices, "pos") == -1) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:slices",
        "The argument 'slices' has an invalid type");
   }

   if (!mxIsCell(pos_lu) ||
       mxGetNumberOfDimensions(pos_lu) != 2 ||
       mxGetM(pos_lu) != 1 ||
       mxGetN(pos_lu) < 1) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:pos_lu",
        "The argument 'pos_lu' has an invalid type");
   }

   if (mxGetClassID(home_pos) != mxUINT64_CLASS ||
       mxGetNumberOfDimensions(home_pos) != 2 ||
       mxGetM(home_pos) != 1) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:mask",
        "The argument 'home_pos' has an invalid type");
   }

   if (mxGetClassID(scan_extents) != mxUINT64_CLASS ||
       mxGetNumberOfDimensions(scan_extents) != 2 ||
       mxGetM(scan_extents) != 2 ||
       mxGetN(scan_extents) != 1) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:scan_extents",
        "The argument 'scan_extents' has an invalid type");
   }

   if (mxGetClassID(channel_id) != mxUINT64_CLASS ||
       mxGetNumberOfDimensions(channel_id) != 2 ||
       mxGetM(channel_id) != 1 ||
       mxGetN(channel_id) != 1) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:channel_id",
        "The argument 'channel_id' has an invalid type");
   }

   if (!mxIsLogical(mask)) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:mask",
        "The argument 'mask' has an invalid type");
   }

   /* evaluate the ROIs  */
   status = nia_scanROINiaMovieImpl(slices, pos_lu,
     home_pos, scan_extents, channel_id, mask, &out);
   if (status) {
      mexErrMsgIdAndTxt("nia_scanROINiaMovieImpl:internal",
        "An internal error has occurred");
      return;
   }

   plhs[0] = out;
}


