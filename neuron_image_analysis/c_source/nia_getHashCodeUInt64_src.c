
/**
 * The following function computes the FNV-1a hash.
 *
 * \param input_len Number of elements in input array.
 * \param input_vals Input array values.
 * \output Output hash code.
 */
static
uint64_t
nia_getHashCodeUInt64(
  size_t input_len,
  uint64_t* input_vals)
{
   const uint64_t prime = 1099511628211ull;
   const uint64_t offset = 14695981039346656037ull;

   size_t i;
   uint64_t code;

   code = offset;
   for (i = 0; i < input_len; i++) {
      int j;
      uint64_t elem;

      elem = input_vals[i];

      for (j = 0; j < 64; j += 4) {
         uint64_t octet;

         octet = (elem >> j) & 0xf;
         code ^= octet;
         code *= prime;
      }
   }

   return code;
}
