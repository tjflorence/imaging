/* File: nia_getHashCodeUInt64.c 
 *
 * This function implements the FNV-1a hash function on an array of doubles.
 *
 * This function can be compiled with one of the following commands:
 *
 * 1) on 64-bit linux with gcc (tested R2013a)
 * mex -v -largeArrayDims CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" nia_getHashCodeUInt64.c
 *
 * 2) on 64-bit linux with icc (tested R2013a)
 * mex -v -largeArrayDims CC=icc CFLAGS="-std=c99 -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3 -fstrict-aliasing" LDOPTIMFLAGS="-static-intel" nia_getHashCodeUInt64.c
 *
 * 3) on 64-bit mac with gcc (tested R2013a), after downloading Xcode 4.6.3, mathworks has done an absolutely pathetic job of implementing this (as usual), so getting this to build is incredibly twitchy
 * mex -v -largeArrayDims CC="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" CFLAGS="-std=c99 -DNIA_WORKAROUND_MW_CHAR16_BUG -fno-common -arch x86_64 -isysroot /Applications/Xcode_4.6.3.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk -mmacosx-version-min=10.7  -fexceptions" COPTIMFLAGS="-O3 -fstrict-aliasing" LD="/Applications/Xcode_4.6.3.app/Contents/Developer/usr/bin/xcrun -sdk macosx10.7 clang" nia_getHashCodeUInt64.c
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
   uint64_t *input_vals;
   size_t input_len;

   const mwSize output_dims[1] = {1};
   uint64_t *output_vals;

   /* check input arguments */
   if (nlhs != 1) {
      mexErrMsgIdAndTxt("nia_getHashCodeUInt64:nlhs",
        "Must have one output argument");
   }

   if (nrhs != 1) {
      mexErrMsgIdAndTxt("nia_getHashCodeUInt64:nrhs",
        "Must have one input arguments");
   }

   if (mxGetClassID(prhs[0]) != mxUINT64_CLASS) {
      mexErrMsgIdAndTxt("nia_getHashCodeUInt64:input",
        "The argument 'input' must be a matrix of uint64 values");
   }

   input_vals = (uint64_t*)mxGetData(prhs[0]);
   input_len = mxGetNumberOfElements(prhs[0]);

   /* allocate the output array */
   plhs[0] = mxCreateNumericArray(1, output_dims, mxUINT64_CLASS, mxREAL);
   if (plhs[0] == NULL) {
      mexErrMsgIdAndTxt("nia_getHashCodeUInt64:internal",
        "An internal error has occurred");
      return;
   }

   output_vals = mxGetData(plhs[0]);

   /* compute the hash */
   *output_vals = nia_getHashCodeUInt64(input_len, input_vals);
}


