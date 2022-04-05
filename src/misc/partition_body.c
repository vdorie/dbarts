#if defined(__USE_AVX2__)
#  define __USE_SIMD__ 1

#  define _mm256_cmpge_epu16(a, b) _mm256_cmpeq_epi16(_mm256_max_epu16(a, b), a)
#  define _mm256_cmple_epu16(a, b) _mm256_cmpge_epu16(b, a)
#  define _mm256_cmpgt_epu16(a, b) _mm256_xor_si256(_mm256_cmple_epu16(a, b), _mm256_set1_epi16(-1))
#  define _mm256_cmplt_epu16(a, b) _mm256_cmpgt_epu16(b, a)

#  define set1(_X_) _mm256_set1_epi16((misc_xint_t) (_X_))

#  define loadLHComp(_X_) \
    (values = _mm256_set_epi16(getDataAt(_X_ + 15), \
                               getDataAt(_X_ + 14), \
                               getDataAt(_X_ + 13), \
                               getDataAt(_X_ + 12), \
                               getDataAt(_X_ + 11), \
                               getDataAt(_X_ + 10), \
                               getDataAt(_X_ +  9), \
                               getDataAt(_X_ +  8), \
                               getDataAt(_X_ +  7), \
                               getDataAt(_X_ +  6), \
                               getDataAt(_X_ +  5), \
                               getDataAt(_X_ +  4), \
                               getDataAt(_X_ +  3), \
                               getDataAt(_X_ +  2), \
                               getDataAt(_X_ +  1), \
                               getDataAt(_X_     )), \
     _mm256_cmpgt_epu16(values, cut_vec))

#  define loadRHComp(_X_) \
    (values = _mm256_set_epi16(getDataAt(_X_ - 15), \
                               getDataAt(_X_ - 14), \
                               getDataAt(_X_ - 13), \
                               getDataAt(_X_ - 12), \
                               getDataAt(_X_ - 11), \
                               getDataAt(_X_ - 10), \
                               getDataAt(_X_ -  9), \
                               getDataAt(_X_ -  8), \
                               getDataAt(_X_ -  7), \
                               getDataAt(_X_ -  6), \
                               getDataAt(_X_ -  5), \
                               getDataAt(_X_ -  4), \
                               getDataAt(_X_ -  3), \
                               getDataAt(_X_ -  2), \
                               getDataAt(_X_ -  1), \
                               getDataAt(_X_     )), \
     _mm256_cmple_epu16(values, cut_vec))

#  define movemask _mm256_movemask_epi8
#  define cmp_width 16
#  define cmp_t __m256i

#elif defined(__USE_SSE2__) || defined(__USE_SSE4_1__)
#  define __USE_SIMD__ 1

#  ifdef __USE_SSE4_1__
#    define _mm_cmpge_epu16(a, b) _mm_cmpeq_epi16(_mm_max_epu16(a, b), a)
#    define _mm_cmple_epu16(a, b) _mm_cmpge_epu16(b, a)
#    define _mm_cmpgt_epu16(a, b) _mm_xor_si128(_mm_cmple_epu16(a, b), _mm_set1_epi16(-1))
#    define _mm_cmplt_epu16(a, b) _mm_cmpgt_epu16(b, a)
#  else
#    define _mm_cmpgt_epu16(a, b) _mm_cmpgt_epi16(_mm_add_epi16(a, _mm_set1_epi16((uint16_t) 0x8000u)), \
                                                  _mm_add_epi16(b, _mm_set1_epi16((uint16_t) 0x8000u)))
#    define _mm_cmplt_epu16(a, b) _mm_cmpgt_epu16(b, a)
#    define _mm_cmpge_epu16(a, b) _mm_xor_si128(_mm_cmplt_epu16(a, b), _mm_set1_epi16(-1))
#    define _mm_cmple_epu16(a, b) _mm_cmpge_epu16(b, a)
#  endif

#  define set1(_X_) _mm_set1_epi16((misc_xint_t) (_X_))

#  define loadLHComp(_X_) \
    (values = _mm_set_epi16(getDataAt(_X_ + 7), \
                            getDataAt(_X_ + 6), \
                            getDataAt(_X_ + 5), \
                            getDataAt(_X_ + 4), \
                            getDataAt(_X_ + 3), \
                            getDataAt(_X_ + 2), \
                            getDataAt(_X_ + 1), \
                            getDataAt(_X_    )), \
     _mm_cmpgt_epu16(values, cut_vec))

#  define loadRHComp(_X_) \
    (values = _mm_set_epi16(getDataAt(_X_ - 7), \
                            getDataAt(_X_ - 6), \
                            getDataAt(_X_ - 5), \
                            getDataAt(_X_ - 4), \
                            getDataAt(_X_ - 3), \
                            getDataAt(_X_ - 2), \
                            getDataAt(_X_ - 1), \
                            getDataAt(_X_    )), \
     _mm_cmple_epu16(values, cut_vec))

#  define movemask _mm_movemask_epi8
#  define cmp_width 8
#  define cmp_t __m128i

#elif defined(__USE_NEON__)
#  define __USE_SIMD__ 1

//#  if PARTITION_RANGE == 1
/* #    define loadLHComp(_X_) \
       (values = vcombine_u16(vcreate_u16(*((uint64_t*) (x + _X_))), vcreate_u16(*(((uint64_t*) (x + _X_ + 4))))), \
          vcgtq_u16(values, cut_vec))
#    define loadRHComp(_X_) \
       (values = vrev64q_u16(vcombine_u16(vcreate_u16(*((uint64_t*) (x + _X_ - 3))), vcreate_u16(*(((uint64_t*) (x + _X_ - 7)))))), \
          vcleq_u16(values, cut_vec)) */
/* #    define loadLHComp(_X_) \
       (values = (((uintptr_t) (x + _X_)) % (8 * sizeof(double))) == 0 ? \
         vld1q_u16(x + _X_) : \
         vcombine_u16(vcreate_u16(*((uint64_t*) (x + _X_))), vcreate_u16(*(((uint64_t*) (x + _X_ + 4))))), \
         vcgtq_u16(values, cut_vec))
#    define loadRHComp(_X_) \
       (values = (((uintptr_t) (x + _X_ - 7)) % (8 * sizeof(double))) == 0 ? \
         vld1q_u16(x + _X_ - 7) : \
         vcombine_u16(vcreate_u16(*((uint64_t*) (x + _X_ - 7))), vcreate_u16(*((uint64_t*) (x + _X_ - 3)))), \
         vcleq_u16(values, cut_vec))  */
//#  else
#    define vset_u16(_X7_, _X6_, _X5_, _X4_, _X3_, _X2_, _X1_, _X0_) \
       vcombine_u16( \
         vcreate_u16(((uint64_t) _X0_) + (((uint64_t) _X1_) << 16) + (((uint64_t) _X2_) << 32) + (((uint64_t) _X3_) << 48)), \
         vcreate_u16(((uint64_t) _X4_) + (((uint64_t) _X5_) << 16) + (((uint64_t) _X6_) << 32) + (((uint64_t) _X7_) << 48)))
#    define loadLHComp(_X_) \
       (values = vset_u16(getDataAt(_X_ + 7), \
                          getDataAt(_X_ + 6), \
                          getDataAt(_X_ + 5), \
                          getDataAt(_X_ + 4), \
                          getDataAt(_X_ + 3), \
                          getDataAt(_X_ + 2), \
                          getDataAt(_X_ + 1), \
                          getDataAt(_X_    )), \
        vcgtq_u16(values, cut_vec))

#    define loadRHComp(_X_) \
       (values = vset_u16(getDataAt(_X_ - 7), \
                          getDataAt(_X_ - 6), \
                          getDataAt(_X_ - 5), \
                          getDataAt(_X_ - 4), \
                          getDataAt(_X_ - 3), \
                          getDataAt(_X_ - 2), \
                          getDataAt(_X_ - 1), \
                          getDataAt(_X_    )), \
        vcleq_u16(values, cut_vec))
//#  endif

#  define set1(_X_) vdupq_n_u16((misc_xint_t) (_X_))

#  define movemask(_X_) vmovemask_u8(vreinterpretq_u8_u16(_X_))
#  define cmp_width 8
#  define cmp_t uint16x8_t

#endif

#if PARTITION_RANGE == 1

#  define getDataAt(_I_) x[_I_]

  size_t lengthOfLeft;
  
  for (size_t i = 0; i < length; ++i) indices[i] = i;
  
  size_t lh = 0, rh = length - 1;

#  ifdef __USE_SIMD__

  if (lh + 2 * cmp_width < rh) {
    
    cmp_t lh_comp, rh_comp, values, cut_vec;
    cut_vec = set1(cut);
    uint_least8_t lh_sub = 0, rh_sub = 0;
    unsigned int lh_mask = 0, rh_mask = 0;
    
    lh_comp = loadLHComp(lh);
    lh_mask = (unsigned int) movemask(lh_comp);
    rh_comp = loadRHComp(rh);
    rh_mask = (unsigned int) movemask(rh_comp);
    
    while (true) {
      while (lh_mask == 0 && lh + 2 * cmp_width < rh) {
        lh += cmp_width;
        lh_comp = loadLHComp(lh);
        lh_mask = (unsigned int) movemask(lh_comp);
        lh_sub = 0;
      }
      while (rh_mask == 0 && lh + 2 * cmp_width < rh) {
        rh -= cmp_width;
        rh_comp = loadRHComp(rh);
        rh_mask = (unsigned int) movemask(rh_comp);
        rh_sub = 0;
      }
      if (lh + 2 * cmp_width >= rh) {
        lh += lh_sub;
        rh -= rh_sub;
        break;
      }
      
      do {
        unsigned int zeros = (unsigned int) countTrailingZeros(lh_mask);
        lh_mask >>= zeros;
        lh_sub += (uint_least8_t) zeros / 2;
        
        zeros = (unsigned int) countTrailingZeros(rh_mask);
        rh_mask >>= zeros;
        rh_sub += (uint_least8_t) zeros / 2;
        
        indices[rh - rh_sub] = lh + lh_sub;
        indices[lh + lh_sub] = rh - rh_sub;
                  
        lh_mask >>= 2;
        rh_mask >>= 2;
        ++lh_sub;
        ++rh_sub;
        
      } while (lh_mask != 0 && rh_mask != 0);
    }
  }

#  endif // __USE_SIMD__
  
  while (true) {
    while (x[lh] <= cut && lh < rh) ++lh;
    while (x[rh]  > cut && lh < rh) --rh;
    
    if (lh >= rh) break;
    
    indices[rh] = lh;
    indices[lh] = rh;
    
    ++lh;
    --rh;
  }
  
  lengthOfLeft = x[indices[lh]] <= cut ? lh + 1 : lh;
  
  return lengthOfLeft;

#else // PARTITION_RANGE == 1 above, 0 below

#  define getDataAt(_I_) x[indices[_I_]]
  
  if (length == 0) return 0;
  
  size_t lengthOfLeft;
  
  size_t lh = 0, rh = length - 1;
  
#  ifdef __USE_SIMD__
  if (lh + 2 * cmp_width < rh) {
    
    cmp_t lh_comp, rh_comp, values, cut_vec;
    cut_vec = set1(cut);
    uint_least8_t lh_sub = 0, rh_sub = 0;
    unsigned int lh_mask = 0, rh_mask = 0;
    
    lh_comp = loadLHComp(lh);
    lh_mask = (unsigned int) movemask(lh_comp);
    rh_comp = loadRHComp(rh);
    rh_mask = (unsigned int) movemask(rh_comp);
    
    while (true) {
      while (lh_mask == 0 && lh + 2 * cmp_width < rh) {
      lh += cmp_width;
        lh_comp = loadLHComp(lh);
        lh_mask = (unsigned int) movemask(lh_comp);
        lh_sub = 0;
      }
      while (rh_mask == 0 && lh + 2 * cmp_width < rh) {
        rh -= cmp_width;
        rh_comp = loadRHComp(rh);
        rh_mask = (unsigned int) movemask(rh_comp);
        rh_sub = 0;
      }
      if (lh + 2 * cmp_width >= rh) {
        lh += lh_sub;
        rh -= rh_sub;
        break;
      }
      
      do {
        unsigned int zeros = (unsigned int) countTrailingZeros(lh_mask);
        lh_mask >>= zeros;
        lh_sub += (uint_least8_t) zeros / 2;
        
        zeros = (unsigned int) countTrailingZeros(rh_mask);
        rh_mask >>= zeros;
        rh_sub += (uint_least8_t) zeros / 2;
        
        size_t temp = indices[rh - rh_sub];
        indices[rh - rh_sub] = indices[lh + lh_sub];
        indices[lh + lh_sub] = temp;
        
        lh_mask >>= 2;
        rh_mask >>= 2;
        ++lh_sub;
        ++rh_sub;
      } while (lh_mask != 0 && rh_mask != 0);
    }
  }

#  endif // __USE_SIMD__

  while (true) {
    while (x[indices[lh]] <= cut && lh < rh) ++lh;
    while (x[indices[rh]]  > cut && lh < rh) --rh;
    
    
    if (lh >= rh) break;
    
    size_t temp = indices[rh];
    indices[rh] = indices[lh];
    indices[lh] = temp;
    
    ++lh;
    --rh;
  }
  
  lengthOfLeft = x[indices[lh]] <= cut ? lh + 1 : lh;
  
  return lengthOfLeft;

#endif // PARTITION_RANGE

#undef cmp_type
#undef cmp_width
#undef movemask
#undef getDataAt
#undef loadRHComp
#undef loadLHComp
#undef _mm_cmplt_epu16
#undef _mm_cmpgt_epu16
#undef _mm_cmple_epu16
#undef _mm_cmpge_epu16
#undef _mm256_cmplt_epu16
#undef _mm256_cmpgt_epu16
#undef _mm256_cmple_epu16
#undef _mm256_cmpge_epu16
#undef __USE_SIMD__

