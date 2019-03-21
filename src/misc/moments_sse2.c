#include "config.h"

#include <math.h>
#include <stdint.h>

#include <misc/intrinsic.h>

double misc_computeUnrolledMean_sse2(const double* x, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = offset == 0 ? 0 : sizeof(double) - offset;
  
  if (prefix > length) prefix = length;
  
  double result = 0.0;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    result += x[i];
  
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  if (suffix > prefix) {
    __m128d result_vec = _mm_setzero_pd();
    
    for ( ; i < suffix; i += 16) {
      result_vec = _mm_add_pd(result_vec,
        _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_load_pd(x + i     ), _mm_load_pd(x + i +  2)),
                              _mm_add_pd(_mm_load_pd(x + i +  4), _mm_load_pd(x + i +  6))),
                   _mm_add_pd(_mm_add_pd(_mm_load_pd(x + i +  8), _mm_load_pd(x + i + 10)),
                              _mm_add_pd(_mm_load_pd(x + i + 12), _mm_load_pd(x + i + 14)))));
    }
    
    double result_arr[2];
    _mm_storeu_pd(result_arr, result_vec);
    
    result += result_arr[0] + result_arr[1];
  }
  
  for ( ; i < length; ++i)
    result += x[i];
  
  return result / (double) length;
}

double misc_computeIndexedUnrolledMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t i = 0;
  size_t lengthMod12 = length % 12;
  
  double result = 0.0;
  if (lengthMod12 != 0) {
    for ( ; i < lengthMod12; ++i) result += x[indices[i]];
    if (length < 12) return result / (double) length;
  }
    
  __m128d result_vec = _mm_setzero_pd();
  
  for ( ; i < length; i += 12) {
    result_vec = _mm_add_pd(result_vec,
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_set_pd(x[indices[i     ]], x[indices[i +  1]]),
                                       _mm_set_pd(x[indices[i +  2]], x[indices[i +  3]])),
                            _mm_add_pd(_mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]),
                                       _mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]))),
                 _mm_add_pd(           _mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]),
                                       _mm_set_pd(x[indices[i + 10]], x[indices[i + 11]]))));
  }
  
  double result_arr[2];
  _mm_storeu_pd(result_arr, result_vec);
  
  result += result_arr[0] + result_arr[1];
  
  return result / (double) length;
}

double misc_computeOnlineUnrolledMean_sse2(const double* x, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = 2 + offset == 0 ? 0 : sizeof(double) - offset;
  
  if (prefix > length) prefix = length;
  
  double result = x[0];
  size_t i = 1;
  for ( ; i < prefix; ++i)
    result += (x[i] - result) / (double) (i + 1);
  
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  if (suffix > prefix) {
    double sum_arr[2];
    
    for ( ; i < suffix; i += 16) {
      __m128d sum_vec =
        _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_load_pd(x + i     ), _mm_load_pd(x + i +  2)),
                              _mm_add_pd(_mm_load_pd(x + i +  4), _mm_load_pd(x + i +  6))),
                   _mm_add_pd(_mm_add_pd(_mm_load_pd(x + i +  8), _mm_load_pd(x + i + 10)),
                              _mm_add_pd(_mm_load_pd(x + i + 12), _mm_load_pd(x + i + 14))));
      
      _mm_storeu_pd(sum_arr, sum_vec);
      
      result += ((sum_arr[0] - 8.0 * result) + (sum_arr[1] - 8.0 * result)) / (double) (i + 16);
    }
  }
  
  for ( ; i < length; ++i)
    result += (x[i] - result) / (double) (i + 1);
  
  return result;
}

double misc_computeIndexedOnlineUnrolledMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t i = 1;
  size_t lengthMod12 = (length - 1) % 12;
  
  double result = x[indices[0]];
  if (lengthMod12++ != 0) {
    for ( ; i < lengthMod12; ++i) result += (x[indices[i]] - result) / (double) (i + 1);
    if (length < 12) return result;
  }
  
  double sum_arr[2];
  for ( ; i < length; i += 12) {
    __m128d sum_vec =
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_set_pd(x[indices[i     ]], x[indices[i +  1]]),
                                       _mm_set_pd(x[indices[i +  2]], x[indices[i +  3]])),
                            _mm_add_pd(_mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]),
                                       _mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]))),
                 _mm_add_pd(           _mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]),
                                       _mm_set_pd(x[indices[i + 10]], x[indices[i + 11]])));
    
    _mm_storeu_pd(sum_arr, sum_vec);
    
    result += ((sum_arr[0] - 6.0 * result) + (sum_arr[1] - 6.0 * result)) / (double) (i + 12);
  }
    
  return result;
}

double misc_computeUnrolledVarianceForKnownMean_sse2(const double* x, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = offset == 0 ? 0 : sizeof(double) - offset;
  
  if (prefix > length) prefix = length;
  
  double result = 0.0;
  size_t i = 0;
  for ( ; i < prefix; ++i)
    result += (x[i] - mean) * (x[i] - mean);
  
  size_t suffix = prefix + 12 * ((length - prefix) / 12);
  
  if (suffix > prefix) {
    __m128d m = _mm_set1_pd(mean);
    __m128d result_vec = _mm_setzero_pd();

    for ( ; i < suffix; i += 12) {
      __m128d a = _mm_load_pd(x + i     ), b = _mm_load_pd(x + i +  2),
              c = _mm_load_pd(x + i +  4), d = _mm_load_pd(x + i +  6),
              e = _mm_load_pd(x + i +  8), f = _mm_load_pd(x + i + 10);
      
      a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
      c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
      e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
      
      result_vec = _mm_add_pd(result_vec,
        _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, a),
                                         _mm_mul_pd(b, b)),
                              _mm_add_pd(_mm_mul_pd(c, c),
                                         _mm_mul_pd(d, d))),
                   _mm_add_pd(           _mm_mul_pd(e, e),
                                         _mm_mul_pd(f, f))));
    }
    
    double result_arr[2];
    _mm_storeu_pd(result_arr, result_vec);
    
    result += result_arr[0] + result_arr[1];
  }
  
  for ( ; i < length; ++i)
    result += (x[i] - mean) * (x[i] - mean);
    
  return result / (double) (length - 1);
}

double misc_computeIndexedUnrolledVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 0;
  size_t lengthMod12 = length % 12;
  
  double result = 0.0;
  if (lengthMod12 != 0) {
    for ( ; i < lengthMod12; ++i) result += (x[indices[i]] - mean) * (x[indices[i]] - mean);
    if (length < 12) return result / (double) (length - 1);
  }
  
  __m128d m = _mm_set1_pd(mean);
  __m128d result_vec = _mm_setzero_pd();
  for ( ; i < length; i += 12) {
    __m128d a = _mm_set_pd(x[indices[i     ]], x[indices[i +  1]]),
            b = _mm_set_pd(x[indices[i +  2]], x[indices[i +  3]]),
            c = _mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]),
            d = _mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]),
            e = _mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]),
            f = _mm_set_pd(x[indices[i + 10]], x[indices[i + 11]]);
    
    a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
    c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
    e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
    
    result_vec = _mm_add_pd(result_vec,
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, a),
                                       _mm_mul_pd(b, b)),
                            _mm_add_pd(_mm_mul_pd(c, c),
                                       _mm_mul_pd(d, d))),
                 _mm_add_pd(           _mm_mul_pd(e, e),
                                       _mm_mul_pd(f, f))));
  }
  
  double result_arr[2];
  _mm_storeu_pd(result_arr, result_vec);
  
  result += result_arr[0] + result_arr[1];
  
  return result / (double) (length - 1);
}

double misc_computeOnlineUnrolledVarianceForKnownMean_sse2(const double* x, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = 2 + offset == 0 ? 0 : sizeof(double) - offset;
  
  if (prefix > length) prefix = length;
  
  double result = (x[0] - mean) * (x[0] - mean) + (x[1] - mean) * (x[1] - mean);
  size_t i = 2;
  for ( ; i < prefix; ++i)
    result += ((x[i] - mean) * (x[i] - mean) - result) / (double) i;
  
  size_t suffix = prefix + 12 * ((length - prefix) / 12);
  
  if (suffix > prefix) {
    __m128d m = _mm_set1_pd(mean);
    
    for ( ; i < suffix; i += 12) {
      __m128d a = _mm_load_pd(x + i     ), b = _mm_load_pd(x + i +  2),
              c = _mm_load_pd(x + i +  4), d = _mm_load_pd(x + i +  6),
              e = _mm_load_pd(x + i +  8), f = _mm_load_pd(x + i + 10);
      
      a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
      c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
      e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
      
      __m128d sum_vec = 
        _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, a),
                                         _mm_mul_pd(b, b)),
                              _mm_add_pd(_mm_mul_pd(c, c),
                                         _mm_mul_pd(d, d))),
                   _mm_add_pd(           _mm_mul_pd(e, e),
                                         _mm_mul_pd(f, f)));
      double sum_arr[2];
      _mm_storeu_pd(sum_arr, sum_vec);
      
      result += ((sum_arr[0] - 6.0 * result) + (sum_arr[1] - 6.0 * result)) / (double) (i + 11);
    }
  }
  
  for ( ; i < length; ++i)
    result += ((x[i] - mean) * (x[i] - mean) - result) / (double) i;
    
  return result;
}

double misc_computeIndexedOnlineUnrolledVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 2;
  size_t lengthMod12 = (length - 2) % 12;
  
  double result = (x[indices[0]] - mean) * (x[indices[0]] - mean) + (x[indices[1]] - mean) * (x[indices[1]] - mean);
  if (lengthMod12 != 0) {
    for ( ; i < lengthMod12 + 2; ++i) result += ((x[indices[i]] - mean) * (x[indices[i]] - mean) - result) / (double) i;
    if (length < 12 + 2) return result;
  }
  
  __m128d m = _mm_set1_pd(mean);
  
  for ( ; i < length; i += 12) {
    __m128d a = _mm_set_pd(x[indices[i     ]], x[indices[i +  1]]),
            b = _mm_set_pd(x[indices[i +  2]], x[indices[i +  3]]),
            c = _mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]),
            d = _mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]),
            e = _mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]),
            f = _mm_set_pd(x[indices[i + 10]], x[indices[i + 11]]);
    
    a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
    c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
    e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
    
    __m128d sum_vec =
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, a),
                                       _mm_mul_pd(b, b)),
                            _mm_add_pd(_mm_mul_pd(c, c),
                                       _mm_mul_pd(d, d))),
                 _mm_add_pd(           _mm_mul_pd(e, e),
                                       _mm_mul_pd(f, f)));
    
    double sum_arr[2];
    _mm_storeu_pd(sum_arr, sum_vec);
    
    result += ((sum_arr[0] - 6.0 * result) + (sum_arr[1] - 6.0 * result)) / (double) (i + 11);
  }
  
  return result;
}

double misc_computeUnrolledWeightedMean_sse2(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t w_offset = ((uintptr_t) w) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = x_offset == 0 ? 0 : sizeof(double) - x_offset;
  
  if (prefix > length) prefix = length;
  
  double result = 0.0, n = 0.0;
  size_t i = 0;
  for ( ; i < prefix; ++i) {
    result += x[i] * w[i];
    n += w[i];
  }
  
  size_t suffix = prefix + 12 * ((length - prefix) / 12);
  
  if (suffix > prefix) {
    __m128d sum_vec = _mm_setzero_pd();
    __m128d n_vec   = _mm_setzero_pd();
    
    if (x_offset == w_offset) {
      for ( ; i < suffix; i += 12) {
        __m128d w_0 = _mm_load_pd(w + i     ), w_1 = _mm_load_pd(w + i +  2);
        __m128d w_2 = _mm_load_pd(w + i +  4), w_3 = _mm_load_pd(w + i +  6);
        __m128d w_4 = _mm_load_pd(w + i +  8), w_5 = _mm_load_pd(w + i + 10);
        
        sum_vec = _mm_add_pd(sum_vec,
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i     ), w_0),
                                           _mm_mul_pd(_mm_load_pd(x + i +  2), w_1)),
                                _mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i +  4), w_2),
                                           _mm_mul_pd(_mm_load_pd(x + i +  6), w_3))),
                     _mm_add_pd(           _mm_mul_pd(_mm_load_pd(x + i +  8), w_4),
                                           _mm_mul_pd(_mm_load_pd(x + i + 10), w_5))));
        
        n_vec = _mm_add_pd(n_vec, _mm_add_pd(_mm_add_pd(_mm_add_pd(w_0, w_1), _mm_add_pd(w_2, w_3)), _mm_add_pd(w_4, w_5)));
      }
      
      double arr[2];
      _mm_storeu_pd(arr, sum_vec);
      result += arr[0] + arr[1];
      
      _mm_storeu_pd(arr, n_vec);
      n += arr[0] + arr[1];
    } else {
      for ( ; i < suffix; i += 12) {
        __m128d w_0 = _mm_loadu_pd(w + i     ), w_1 = _mm_loadu_pd(w + i +  2);
        __m128d w_2 = _mm_loadu_pd(w + i +  4), w_3 = _mm_loadu_pd(w + i +  6);
        __m128d w_4 = _mm_loadu_pd(w + i +  8), w_5 = _mm_loadu_pd(w + i + 10);
        
        sum_vec = _mm_add_pd(sum_vec,
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i     ), w_0),
                                           _mm_mul_pd(_mm_load_pd(x + i +  2), w_1)),
                                _mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i +  4), w_2),
                                           _mm_mul_pd(_mm_load_pd(x + i +  6), w_3))),
                     _mm_add_pd(           _mm_mul_pd(_mm_load_pd(x + i +  8), w_4),
                                           _mm_mul_pd(_mm_load_pd(x + i + 10), w_5))));
        
        n_vec = _mm_add_pd(n_vec, _mm_add_pd(_mm_add_pd(_mm_add_pd(w_0, w_1), _mm_add_pd(w_2, w_3)), _mm_add_pd(w_4, w_5)));
      }
      
      double arr[2];
      _mm_storeu_pd(arr, sum_vec);
      result += arr[0] + arr[1];
      
      _mm_storeu_pd(arr, n_vec);
      n += arr[0] + arr[1];
    }
  }
  
  for ( ; i < length; ++i) {
    result += x[i] * w[i];
    n += w[i];
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result / n;
}

double misc_computeIndexedUnrolledWeightedMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t i = 0;
  size_t lengthMod12 = length % 12;
  
  double result = 0.0;
  double n = 0.0;
  if (lengthMod12 != 0) {
    for ( ; i < lengthMod12; ++i) {
      result += x[indices[i]] * w[indices[i]];
      n += w[indices[i]];
    }
    if (length < 12) { if (nPtr != NULL) *nPtr = n; return result / n; }
  }
  
  __m128d sum_vec = _mm_setzero_pd();
  __m128d n_vec   = _mm_setzero_pd();
  for ( ; i < length; i += 12) {
    __m128d w_0 = _mm_set_pd(w[indices[i     ]], w[indices[i +  1]]);
    __m128d w_1 = _mm_set_pd(w[indices[i +  2]], w[indices[i +  3]]);
    __m128d w_2 = _mm_set_pd(w[indices[i +  4]], w[indices[i +  5]]);
    __m128d w_3 = _mm_set_pd(w[indices[i +  6]], w[indices[i +  7]]);
    __m128d w_4 = _mm_set_pd(w[indices[i +  8]], w[indices[i +  9]]);
    __m128d w_5 = _mm_set_pd(w[indices[i + 10]], w[indices[i + 11]]);
    
    sum_vec = _mm_add_pd(sum_vec,
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set_pd(x[indices[i     ]], x[indices[i +  1]]), w_0),
                                       _mm_mul_pd(_mm_set_pd(x[indices[i +  2]], x[indices[i +  3]]), w_1)),
                            _mm_add_pd(_mm_mul_pd(_mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]), w_2),
                                       _mm_mul_pd(_mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]), w_3))),
                 _mm_add_pd(           _mm_mul_pd(_mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]), w_4),
                                       _mm_mul_pd(_mm_set_pd(x[indices[i + 10]], x[indices[i + 11]]), w_5))));
        
    n_vec = _mm_add_pd(n_vec, _mm_add_pd(_mm_add_pd(_mm_add_pd(w_0, w_1), _mm_add_pd(w_2, w_3)), _mm_add_pd(w_4, w_5)));
  }
  
  double arr[2];
  _mm_storeu_pd(arr, sum_vec);
  result += arr[0] + arr[1];
  
  _mm_storeu_pd(arr, n_vec);
  n += arr[0] + arr[1];
  
  if (nPtr != NULL) *nPtr = n;
  return result / n;
}

double misc_computeOnlineUnrolledWeightedMean_sse2(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t w_offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = 2 + x_offset == 0 ? 0 : sizeof(double) - x_offset;
  
  if (prefix > length) prefix = length;
  
  size_t i = 1;
  
  double n = w[0];
  double result = x[0];
  for ( ; i < prefix; ++i) {
    n += w[i];
    result += (x[i] - result) * (w[i] / n);
  }
  
  size_t suffix = prefix + 12 * ((length - prefix) / 12);
  
  if (suffix > prefix) {
    
    double arr[2];
    
    if (x_offset == w_offset) {
      for ( ; i < suffix; i += 12) {
        __m128d w_0 = _mm_load_pd(w + i     ), w_1 = _mm_load_pd(w + i +  2);
        __m128d w_2 = _mm_load_pd(w + i +  4), w_3 = _mm_load_pd(w + i +  6);
        __m128d w_4 = _mm_load_pd(w + i +  8), w_5 = _mm_load_pd(w + i + 10);
        
        __m128d n_vec = _mm_add_pd(_mm_add_pd(_mm_add_pd(w_0, w_1), _mm_add_pd(w_2, w_3)), _mm_add_pd(w_4, w_5));
        __m128d sum_vec =
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i     ), w_0),
                                           _mm_mul_pd(_mm_load_pd(x + i +  2), w_1)),
                                _mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i +  4), w_2),
                                           _mm_mul_pd(_mm_load_pd(x + i +  6), w_3))),
                     _mm_add_pd(           _mm_mul_pd(_mm_load_pd(x + i +  8), w_4),
                                           _mm_mul_pd(_mm_load_pd(x + i + 10), w_5)));
        _mm_storeu_pd(arr, n_vec);
        double delta_n = arr[0] + arr[1];
        
        _mm_storeu_pd(arr, sum_vec);
         
        n += delta_n;
        result += (arr[0] + arr[1] - delta_n * result) / n;
      }
    } else {
      for ( ; i < suffix; i += 12) {
        __m128d w_0 = _mm_loadu_pd(w + i     ), w_1 = _mm_loadu_pd(w + i +  2);
        __m128d w_2 = _mm_loadu_pd(w + i +  4), w_3 = _mm_loadu_pd(w + i +  6);
        __m128d w_4 = _mm_loadu_pd(w + i +  8), w_5 = _mm_loadu_pd(w + i + 10);
        
        __m128d n_vec = _mm_add_pd(_mm_add_pd(_mm_add_pd(w_0, w_1), _mm_add_pd(w_2, w_3)), _mm_add_pd(w_4, w_5));
        __m128d sum_vec =
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i     ), w_0),
                                           _mm_mul_pd(_mm_load_pd(x + i +  2), w_1)),
                                _mm_add_pd(_mm_mul_pd(_mm_load_pd(x + i +  4), w_2),
                                           _mm_mul_pd(_mm_load_pd(x + i +  6), w_3))),
                     _mm_add_pd(           _mm_mul_pd(_mm_load_pd(x + i +  8), w_4),
                                           _mm_mul_pd(_mm_load_pd(x + i + 10), w_5)));
        _mm_storeu_pd(arr, n_vec);
        double delta_n = arr[0] + arr[1];
        
        _mm_storeu_pd(arr, sum_vec);
         
        n += delta_n;
        result += (arr[0] + arr[1] - delta_n * result) / n;
      }
    }
  }
  for ( ; i < length; ++i) {
    n += w[i];
    result += (x[i] - result) * (w[i] / n);
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result;
}

double misc_computeIndexedOnlineUnrolledWeightedMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t lengthMod12 = (length - 1) % 12;
  
  double n = w[indices[0]];
  double result = x[indices[0]];
  
  size_t i = 1;
  if (lengthMod12++ != 0) {
    for ( ; i < lengthMod12; ++i) {
      n += w[indices[i]];
      result += (x[indices[i]] - result) * (w[indices[i]] / n);
    }
    if (length < 12 + 1) {
      if (nPtr != NULL) *nPtr = n;
      return result;
    }
  }
  
  for ( ; i < length; i += 12) {
    __m128d w_0 = _mm_set_pd(w[indices[i     ]], w[indices[i +  1]]);
    __m128d w_1 = _mm_set_pd(w[indices[i +  2]], w[indices[i +  3]]);
    __m128d w_2 = _mm_set_pd(w[indices[i +  4]], w[indices[i +  5]]);
    __m128d w_3 = _mm_set_pd(w[indices[i +  6]], w[indices[i +  7]]);
    __m128d w_4 = _mm_set_pd(w[indices[i +  8]], w[indices[i +  9]]);
    __m128d w_5 = _mm_set_pd(w[indices[i + 10]], w[indices[i + 11]]);
    
    __m128d sum_vec = 
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_set_pd(x[indices[i     ]], x[indices[i +  1]]), w_0),
                                       _mm_mul_pd(_mm_set_pd(x[indices[i +  2]], x[indices[i +  3]]), w_1)),
                            _mm_add_pd(_mm_mul_pd(_mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]), w_2),
                                       _mm_mul_pd(_mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]), w_3))),
                 _mm_add_pd(           _mm_mul_pd(_mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]), w_4),
                                       _mm_mul_pd(_mm_set_pd(x[indices[i + 10]], x[indices[i + 11]]), w_5)));
        
    __m128d n_vec = _mm_add_pd(_mm_add_pd(_mm_add_pd(w_0, w_1), _mm_add_pd(w_2, w_3)), _mm_add_pd(w_4, w_5));
    
    double arr[2];
    _mm_storeu_pd(arr, n_vec);
    
    double delta_n = arr[0] + arr[1];
    
    _mm_storeu_pd(arr, sum_vec);
    
    n += delta_n;
    result += (arr[0] + arr[1] - delta_n * result) / n;
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result;
}

double misc_computeUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t w_offset = ((uintptr_t) w) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = x_offset == 0 ? 0 : sizeof(double) - x_offset;
  
  if (prefix > length) prefix = length;
  
  double result = 0.0;
  size_t i = 0;
  for ( ; i < prefix; ++i)
    result += w[i] * (x[i] - mean) * (x[i] - mean);
  
  size_t suffix = prefix + 12 * ((length - prefix) / 12);
  
  if (suffix > prefix) {
    __m128d m = _mm_set1_pd(mean);
    __m128d result_vec = _mm_setzero_pd();
    
    if (x_offset == w_offset) {
      for ( ; i < suffix; i += 12) {
        __m128d a = _mm_load_pd(x + i     ), b = _mm_load_pd(x + i +  2),
                c = _mm_load_pd(x + i +  4), d = _mm_load_pd(x + i +  6),
                e = _mm_load_pd(x + i +  8), f = _mm_load_pd(x + i + 10);
        
        a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
        c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
        e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
        
        a = _mm_mul_pd(a, a); b = _mm_mul_pd(b, b);
        c = _mm_mul_pd(c, c); d = _mm_mul_pd(d, d);
        e = _mm_mul_pd(e, e); f = _mm_mul_pd(f, f);
        
        result_vec = _mm_add_pd(result_vec,
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, _mm_load_pd(w + i     )),
                                           _mm_mul_pd(b, _mm_load_pd(w + i +  2))),
                                _mm_add_pd(_mm_mul_pd(c, _mm_load_pd(w + i +  4)),
                                           _mm_mul_pd(d, _mm_load_pd(w + i +  6)))),
                     _mm_add_pd(           _mm_mul_pd(e, _mm_load_pd(w + i +  8)),
                                           _mm_mul_pd(f, _mm_load_pd(w + i + 10)))));
      }
    } else {
      for ( ; i < suffix; i += 12) {
        __m128d a = _mm_load_pd(x + i     ), b = _mm_load_pd(x + i +  2),
                c = _mm_load_pd(x + i +  4), d = _mm_load_pd(x + i +  6),
                e = _mm_load_pd(x + i +  8), f = _mm_load_pd(x + i + 10);
        
        a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
        c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
        e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
        
        a = _mm_mul_pd(a, a); b = _mm_mul_pd(b, b);
        c = _mm_mul_pd(c, c); d = _mm_mul_pd(d, d);
        e = _mm_mul_pd(e, e); f = _mm_mul_pd(f, f);
        
        result_vec = _mm_add_pd(result_vec,
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, _mm_loadu_pd(w + i     )),
                                           _mm_mul_pd(b, _mm_loadu_pd(w + i +  2))),
                                _mm_add_pd(_mm_mul_pd(c, _mm_loadu_pd(w + i +  4)),
                                           _mm_mul_pd(d, _mm_loadu_pd(w + i +  6)))),
                     _mm_add_pd(           _mm_mul_pd(e, _mm_loadu_pd(w + i +  8)),
                                           _mm_mul_pd(f, _mm_loadu_pd(w + i + 10)))));
      }
    }
    
    double result_arr[2];
    _mm_storeu_pd(result_arr, result_vec);
    
    result += result_arr[0] + result_arr[1];
  }
  
  for ( ; i < length; ++i)
    result += w[i] * (x[i] - mean) * (x[i] - mean);
    
  return result / (double) (length - 1);
}

double misc_computeIndexedUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
    
  size_t i = 0;
  size_t lengthMod12 = length % 12;
  
  double result = 0.0;
  if (lengthMod12 != 0) {
    for ( ; i < lengthMod12; ++i) result += w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean);
    if (length < 12) return result / (double) (length - 1);
  }
  
  __m128d m = _mm_set1_pd(mean);
  __m128d result_vec = _mm_setzero_pd();
  
  for ( ; i < length; i += 12) {
    __m128d a = _mm_set_pd(x[indices[i     ]], x[indices[i +  1]]);
    __m128d b = _mm_set_pd(x[indices[i +  2]], x[indices[i +  3]]);
    __m128d c = _mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]);
    __m128d d = _mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]);
    __m128d e = _mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]);
    __m128d f = _mm_set_pd(x[indices[i + 10]], x[indices[i + 11]]);
    
    a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
    c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
    e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
    
    a = _mm_mul_pd(a, a); b = _mm_mul_pd(b, b);
    c = _mm_mul_pd(c, c); d = _mm_mul_pd(d, d);
    e = _mm_mul_pd(e, e); f = _mm_mul_pd(f, f);
    
    result_vec = _mm_add_pd(result_vec,
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, _mm_set_pd(w[indices[i     ]], w[indices[i +  1]])),
                                       _mm_mul_pd(b, _mm_set_pd(w[indices[i +  2]], w[indices[i +  3]]))),
                            _mm_add_pd(_mm_mul_pd(c, _mm_set_pd(w[indices[i +  4]], w[indices[i +  5]])),
                                       _mm_mul_pd(d, _mm_set_pd(w[indices[i +  6]], w[indices[i +  7]])))),
                 _mm_add_pd(           _mm_mul_pd(e, _mm_set_pd(w[indices[i +  8]], w[indices[i +  9]])),
                                       _mm_mul_pd(f, _mm_set_pd(w[indices[i + 10]], w[indices[i + 11]])))));
  }
  
  double result_arr[2];
  _mm_storeu_pd(result_arr, result_vec);
  
  result += result_arr[0] + result_arr[1];
      
  return result / (double) (length - 1);
}

double misc_computeOnlineUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double)) / sizeof(double);
  size_t w_offset = ((uintptr_t) w) % (2 * sizeof(double)) / sizeof(double);
  size_t prefix = 2 + x_offset == 0 ? 0 : sizeof(double) - x_offset;
  
  if (prefix > length) prefix = length;
  
  double result = w[0] * (x[0] - mean) * (x[0] - mean) + w[1] * (x[1] - mean) * (x[1] - mean);
  size_t i = 2;
  for ( ; i < prefix; ++i)
    result += (w[i] * (x[i] - mean) * (x[i] - mean) - result) / (double) i;
  
  size_t suffix = prefix + 12 * ((length - prefix) / 12);
  
  if (suffix > prefix) {
    __m128d m = _mm_set1_pd(mean);
    double sum_arr[2];
    
    if (x_offset == w_offset) {
      for ( ; i < suffix; i += 12) {
        __m128d a = _mm_load_pd(x + i     ), b = _mm_load_pd(x + i +  2),
                c = _mm_load_pd(x + i +  4), d = _mm_load_pd(x + i +  6),
                e = _mm_load_pd(x + i +  8), f = _mm_load_pd(x + i + 10);
        
        a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
        c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
        e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
        
        a = _mm_mul_pd(a, a); b = _mm_mul_pd(b, b);
        c = _mm_mul_pd(c, c); d = _mm_mul_pd(d, d);
        e = _mm_mul_pd(e, e); f = _mm_mul_pd(f, f);
        
        __m128d sum_vec =
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, _mm_load_pd(w + i     )),
                                           _mm_mul_pd(b, _mm_load_pd(w + i +  2))),
                                _mm_add_pd(_mm_mul_pd(c, _mm_load_pd(w + i +  4)),
                                           _mm_mul_pd(d, _mm_load_pd(w + i +  6)))),
                     _mm_add_pd(           _mm_mul_pd(e, _mm_load_pd(w + i +  8)),
                                           _mm_mul_pd(f, _mm_load_pd(w + i + 10))));
        
        _mm_storeu_pd(sum_arr, sum_vec);
        result += ((sum_arr[0] - 6.0 * result) + (sum_arr[1] - 6.0 * result)) / (double) (i + 11);
      }
    } else {
      for ( ; i < suffix; i += 12) {
        __m128d a = _mm_load_pd(x + i     ), b = _mm_load_pd(x + i +  2),
                c = _mm_load_pd(x + i +  4), d = _mm_load_pd(x + i +  6),
                e = _mm_load_pd(x + i +  8), f = _mm_load_pd(x + i + 10);
        
        a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
        c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
        e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
        
        a = _mm_mul_pd(a, a); b = _mm_mul_pd(b, b);
        c = _mm_mul_pd(c, c); d = _mm_mul_pd(d, d);
        e = _mm_mul_pd(e, e); f = _mm_mul_pd(f, f);
        
        __m128d sum_vec =
          _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, _mm_loadu_pd(w + i     )),
                                           _mm_mul_pd(b, _mm_loadu_pd(w + i +  2))),
                                _mm_add_pd(_mm_mul_pd(c, _mm_loadu_pd(w + i +  4)),
                                           _mm_mul_pd(d, _mm_loadu_pd(w + i +  6)))),
                     _mm_add_pd(           _mm_mul_pd(e, _mm_loadu_pd(w + i +  8)),
                                           _mm_mul_pd(f, _mm_loadu_pd(w + i + 10))));
        
        _mm_storeu_pd(sum_arr, sum_vec);
        result += ((sum_arr[0] - 6.0 * result) + (sum_arr[1] - 6.0 * result)) / (double) (i + 11);
      }
    }
  }
  
  for ( ; i < length; ++i)
    result += (w[i] * (x[i] - mean) * (x[i] - mean) - result) / (double) i;
    
  return result;
}

double misc_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
    
  size_t i = 2;
  size_t lengthMod12 = (length - 2) % 12;
  
  double result = w[indices[0]] * (x[indices[0]] - mean) * (x[indices[0]] - mean) + w[indices[1]] * (x[indices[1]] - mean) * (x[indices[1]] - mean);
  if (lengthMod12 != 0) {
    for ( ; i < lengthMod12 + 2; ++i) result += (w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean) - result) / (double) i;
    if (length < 12 + 2) return result;
  }
  
  __m128d m = _mm_set1_pd(mean);
  double sum_arr[2];
    
  for ( ; i < length; i += 12) {
    __m128d a = _mm_set_pd(x[indices[i     ]], x[indices[i +  1]]);
    __m128d b = _mm_set_pd(x[indices[i +  2]], x[indices[i +  3]]);
    __m128d c = _mm_set_pd(x[indices[i +  4]], x[indices[i +  5]]);
    __m128d d = _mm_set_pd(x[indices[i +  6]], x[indices[i +  7]]);
    __m128d e = _mm_set_pd(x[indices[i +  8]], x[indices[i +  9]]);
    __m128d f = _mm_set_pd(x[indices[i + 10]], x[indices[i + 11]]);
    
    a = _mm_sub_pd(a, m); b = _mm_sub_pd(b, m);
    c = _mm_sub_pd(c, m); d = _mm_sub_pd(d, m);
    e = _mm_sub_pd(e, m); f = _mm_sub_pd(f, m);
    
    a = _mm_mul_pd(a, a); b = _mm_mul_pd(b, b);
    c = _mm_mul_pd(c, c); d = _mm_mul_pd(d, d);
    e = _mm_mul_pd(e, e); f = _mm_mul_pd(f, f);
    
    __m128d sum_vec =
      _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(a, _mm_set_pd(w[indices[i     ]], w[indices[i +  1]])),
                                       _mm_mul_pd(b, _mm_set_pd(w[indices[i +  2]], w[indices[i +  3]]))),
                            _mm_add_pd(_mm_mul_pd(c, _mm_set_pd(w[indices[i +  4]], w[indices[i +  5]])),
                                       _mm_mul_pd(d, _mm_set_pd(w[indices[i +  6]], w[indices[i +  7]])))),
                 _mm_add_pd(           _mm_mul_pd(e, _mm_set_pd(w[indices[i +  8]], w[indices[i +  9]])),
                                       _mm_mul_pd(f, _mm_set_pd(w[indices[i + 10]], w[indices[i + 11]]))));
    
    _mm_storeu_pd(sum_arr, sum_vec);
    result += ((sum_arr[0] - 6.0 * result) + (sum_arr[1] - 6.0 * result)) / (double) (i + 11);
  }
  
  return result;
}

