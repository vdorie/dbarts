#if defined(__USE_AVX2__)
#  define __USE_SIMD__ 1

#  define _mm256_cmpge_epu16(a, b) _mm256_cmpeq_epi16(_mm256_max_epu16(a, b), a)
#  define _mm256_cmple_epu16(a, b) _mm256_cmpge_epu16(b, a)
#  define _mm256_cmpgt_epu16(a, b) _mm256_xor_si256(_mm256_cmple_epu16(a, b), _mm256_set1_epi16(-1))
#  define _mm256_cmplt_epu16(a, b) _mm256_cmpgt_epu16(b, a)

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
                               getDataAt(_X_    )), \
     _mm256_cmpgt_epu16(values, _mm256_set1_epi16(cut)))

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
                               getDataAt(_X_    )), \
     _mm256_cmple_epu16(values, _mm256_set1_epi16(cut)))

#  define movemask _mm256_movemask_epi8
#  define cmp_width 16
#  define cmp_type __m256i

#elif defined(__USE_SSE2__) || defined(__USE_SSE4_1__)
#  define __USE_SIMD__ 1

#  ifdef __USE_SSE4_1__
#    define _mm_cmpge_epu16(a, b) _mm_cmpeq_epi16(_mm_max_epu16(a, b), a)
#  else
#    define _mm_cmpge_epu16(a, b) ~_mm_cmplt_epi16(_mm_add_epi16(a, _mm_set1_epi16((uint16_t) 0x8000u)), \
                                                   _mm_add_epi16(b, _mm_set1_epi16((uint16_t) 0x8000u)))
#  endif

#  define _mm_cmple_epu16(a, b) _mm_cmpge_epu16(b, a)
#  define _mm_cmpgt_epu16(a, b) _mm_xor_si128(_mm_cmple_epu16(a, b), _mm_set1_epi16(-1))
#  define _mm_cmplt_epu16(a, b) _mm_cmpgt_epu16(b, a)

#  define loadLHComp(_X_) \
    (values = _mm_set_epi16(getDataAt(_X_ + 7), \
                            getDataAt(_X_ + 6), \
                            getDataAt(_X_ + 5), \
                            getDataAt(_X_ + 4), \
                            getDataAt(_X_ + 3), \
                            getDataAt(_X_ + 2), \
                            getDataAt(_X_ + 1), \
                            getDataAt(_X_    )), \
     _mm_cmpgt_epu16(values, _mm_set1_epi16((misc_xint_t) cut)))

#  define loadRHComp(_X_) \
    (values = _mm_set_epi16(getDataAt(_X_ - 7), \
                            getDataAt(_X_ - 6), \
                            getDataAt(_X_ - 5), \
                            getDataAt(_X_ - 4), \
                            getDataAt(_X_ - 3), \
                            getDataAt(_X_ - 2), \
                            getDataAt(_X_ - 1), \
                            getDataAt(_X_    )), \
     _mm_cmple_epu16(values, _mm_set1_epi16((misc_xint_t) cut)))

#  define movemask _mm_movemask_epi8
#  define cmp_width 8
#  define cmp_type __m128i

#endif

#if PARTITION_RANGE == 1

#  define getDataAt(_I_) x[_I_]


  size_t lengthOfLeft;
  
  for (size_t i = 0; i < length; ++i) indices[i] = i;
  
  size_t lh = 0, rh = length - 1;

#  ifdef __USE_SIMD__

  if (lh + 2 * cmp_width < rh) {
    
    cmp_type lh_comp, rh_comp, values;
    uint8_t lh_sub = 0, rh_sub = 0;
    uint16_t lh_mask = 0, rh_mask = 0;
    
    lh_comp = loadLHComp(lh);
    lh_mask = movemask(lh_comp);
    rh_comp = loadRHComp(rh);
    rh_mask = movemask(rh_comp);
    
    while (true) {
      while (lh_mask == 0 && lh + 2 * cmp_width < rh) {
        lh += cmp_width;
        lh_comp = loadLHComp(lh);
        lh_mask = movemask(lh_comp);
        lh_sub = 0;
      }
      while (rh_mask == 0 && lh + 2 * cmp_width < rh) {
        rh -= cmp_width;
        rh_comp = loadRHComp(rh);
        rh_mask = movemask(rh_comp);
        rh_sub = 0;
      }
      if (lh + 2 * cmp_width >= rh) {
        lh += lh_sub;
        rh -= rh_sub;
        break;
      }
      
      do {
        int zeros = countTrailingZeros(lh_mask);
        lh_mask >>= zeros;
        lh_sub += zeros / 2;
        
        zeros = countTrailingZeros(rh_mask);
        rh_mask >>= zeros;
        rh_sub += zeros / 2;
        
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
    
    cmp_type lh_comp, rh_comp, values;
    uint8_t lh_sub = 0, rh_sub = 0;
    uint16_t lh_mask = 0, rh_mask = 0;
    
    lh_comp = loadLHComp(lh);
    lh_mask = movemask(lh_comp);
    rh_comp = loadRHComp(rh);
    rh_mask = movemask(rh_comp);
    
    while (true) {
      while (lh_mask == 0 && lh + 2 * cmp_width < rh) {
      lh += cmp_width;
        lh_comp = loadLHComp(lh);
        lh_mask = movemask(lh_comp);
        lh_sub = 0;
      }
      while (rh_mask == 0 && lh + 2 * cmp_width < rh) {
        rh -= cmp_width;
        rh_comp = loadRHComp(rh);
        rh_mask = movemask(rh_comp);
        rh_sub = 0;
      }
      if (lh + 2 * cmp_width >= rh) {
        lh += lh_sub;
        rh -= rh_sub;
        break;
      }
      
      do {
        int zeros = countTrailingZeros(lh_mask);
        lh_mask >>= zeros;
        lh_sub += zeros / 2;
        
        zeros = countTrailingZeros(rh_mask);
        rh_mask >>= zeros;
        rh_sub += zeros / 2;
        
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

