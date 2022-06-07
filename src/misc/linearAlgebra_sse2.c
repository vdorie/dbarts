#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/intrinsic.h>

/* void misc_addVectors_sse2(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (2 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (2 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 8 * ((length - prefix) / 8);

  for ( ; i < prefix; ++i)
    z[i] = y[i] + x[i];
  
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
      _mm_stream_pd(z + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_load_pd(x + i + 2)));
      _mm_stream_pd(z + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_load_pd(x + i + 4)));
      _mm_stream_pd(z + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_load_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
    }
  } else {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_loadu_pd(y + i    ), _mm_loadu_pd(x + i    )));
      _mm_stream_pd(z + i + 2, _mm_add_pd(_mm_loadu_pd(y + i + 2), _mm_loadu_pd(x + i + 2)));
      _mm_stream_pd(z + i + 4, _mm_add_pd(_mm_loadu_pd(y + i + 4), _mm_loadu_pd(x + i + 4)));
      _mm_stream_pd(z + i + 6, _mm_add_pd(_mm_loadu_pd(y + i + 6), _mm_loadu_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_loadu_pd(y + i    ), _mm_loadu_pd(x + i    )));
    }
  }

  
  for ( ; i < length; ++i)
    z[i] = y[i] + x[i];
 }

void misc_subtractVectors_sse2(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (2 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (2 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 8 * ((length - prefix) / 8);

  
  for ( ; i < prefix; ++i)
    z[i] = y[i] - x[i];
  
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(z + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
      _mm_stream_pd(z + i + 2, _mm_sub_pd(_mm_load_pd(y + i + 2), _mm_load_pd(x + i + 2)));
      _mm_stream_pd(z + i + 4, _mm_sub_pd(_mm_load_pd(y + i + 4), _mm_load_pd(x + i + 4)));
      _mm_stream_pd(z + i + 6, _mm_sub_pd(_mm_load_pd(y + i + 6), _mm_load_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(z + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
    }
  } else {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(z + i    , _mm_sub_pd(_mm_loadu_pd(y + i    ), _mm_loadu_pd(x + i    )));
      _mm_stream_pd(z + i + 2, _mm_sub_pd(_mm_loadu_pd(y + i + 2), _mm_loadu_pd(x + i + 2)));
      _mm_stream_pd(z + i + 4, _mm_sub_pd(_mm_loadu_pd(y + i + 4), _mm_loadu_pd(x + i + 4)));
      _mm_stream_pd(z + i + 6, _mm_sub_pd(_mm_loadu_pd(y + i + 6), _mm_loadu_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(z + i    , _mm_sub_pd(_mm_loadu_pd(y + i    ), _mm_loadu_pd(x + i    )));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] - x[i];
}

void misc_addVectorsWithMultiplier_sse2(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (2 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (2 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 8 * ((length - prefix) / 8);

  for ( ; i < prefix; ++i)
    z[i] = y[i] + alpha * x[i];
  
  __m128d alpha_vec = _mm_load1_pd(&alpha);
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
      _mm_stream_pd(z + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_mul_pd(_mm_load_pd(x + i + 2), alpha_vec)));
      _mm_stream_pd(z + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_mul_pd(_mm_load_pd(x + i + 4), alpha_vec)));
      _mm_stream_pd(z + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_mul_pd(_mm_load_pd(x + i + 6), alpha_vec)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
    }
  } else {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_loadu_pd(y + i    ), _mm_mul_pd(_mm_loadu_pd(x + i    ), alpha_vec)));
      _mm_stream_pd(z + i + 2, _mm_add_pd(_mm_loadu_pd(y + i + 2), _mm_mul_pd(_mm_loadu_pd(x + i + 2), alpha_vec)));
      _mm_stream_pd(z + i + 4, _mm_add_pd(_mm_loadu_pd(y + i + 4), _mm_mul_pd(_mm_loadu_pd(x + i + 4), alpha_vec)));
      _mm_stream_pd(z + i + 6, _mm_add_pd(_mm_loadu_pd(y + i + 6), _mm_mul_pd(_mm_loadu_pd(x + i + 6), alpha_vec)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(z + i    , _mm_add_pd(_mm_loadu_pd(y + i    ), _mm_mul_pd(_mm_loadu_pd(x + i    ), alpha_vec)));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + alpha * x[i];
} */

void misc_addVectorsInPlace_sse2(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) y[i] += x[i];
  
  for ( ; i < length; i += 4) {
    y[i    ] += x[i    ];
    y[i + 1] += x[i + 1];
    y[i + 2] += x[i + 2];
    y[i + 3] += x[i + 3];
  }
  /*
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (2 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;

  size_t i = 0;
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  for ( ; i < prefix; ++i)
    y[i] += x[i];
  
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
      _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_load_pd(x + i + 2)));
      _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_load_pd(x + i + 4)));
      _mm_stream_pd(y + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_load_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
    }
  } else {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_loadu_pd(x + i    )));
      _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_loadu_pd(x + i + 2)));
      _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_loadu_pd(x + i + 4)));
      _mm_stream_pd(y + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_loadu_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_loadu_pd(x + i    )));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += x[i]; */
}

void misc_subtractVectorsInPlace_sse2(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) y[i] -= x[i];
  
  for ( ; i < length; i += 4) {
    y[i    ] -= x[i    ];
    y[i + 1] -= x[i + 1];
    y[i + 2] -= x[i + 2];
    y[i + 3] -= x[i + 3];
  }
  /*
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (2 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;

  size_t i = 0;
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  for ( ; i < prefix; ++i)
    y[i] -= x[i];
  
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(y + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
      _mm_stream_pd(y + i + 2, _mm_sub_pd(_mm_load_pd(y + i + 2), _mm_load_pd(x + i + 2)));
      _mm_stream_pd(y + i + 4, _mm_sub_pd(_mm_load_pd(y + i + 4), _mm_load_pd(x + i + 4)));
      _mm_stream_pd(y + i + 6, _mm_sub_pd(_mm_load_pd(y + i + 6), _mm_load_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(y + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
    }
  } else {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(y + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_loadu_pd(x + i    )));
      _mm_stream_pd(y + i + 2, _mm_sub_pd(_mm_load_pd(y + i + 2), _mm_loadu_pd(x + i + 2)));
      _mm_stream_pd(y + i + 4, _mm_sub_pd(_mm_load_pd(y + i + 4), _mm_loadu_pd(x + i + 4)));
      _mm_stream_pd(y + i + 6, _mm_sub_pd(_mm_load_pd(y + i + 6), _mm_loadu_pd(x + i + 6)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(y + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_loadu_pd(x + i    )));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] -= x[i]; */
}

void misc_addVectorsInPlaceWithMultiplier_sse2(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) y[i] += alpha * x[i];
  
  for ( ; i < length; i += 4) {
    y[i    ] += alpha * x[i    ];
    y[i + 1] += alpha * x[i + 1];
    y[i + 2] += alpha * x[i + 2];
    y[i + 3] += alpha * x[i + 3];
  }
  
  /*
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (2 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;

  size_t i = 0;
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  for ( ; i < prefix; ++i)
    y[i] += alpha * x[i];
  
  __m128d alpha_vec = _mm_load1_pd(&alpha);
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
      _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_mul_pd(_mm_load_pd(x + i + 2), alpha_vec)));
      _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_mul_pd(_mm_load_pd(x + i + 4), alpha_vec)));
      _mm_stream_pd(y + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_mul_pd(_mm_load_pd(x + i + 6), alpha_vec)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
    }
  } else {
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_loadu_pd(x + i    ), alpha_vec)));
      _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_mul_pd(_mm_loadu_pd(x + i + 2), alpha_vec)));
      _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_mul_pd(_mm_loadu_pd(x + i + 4), alpha_vec)));
      _mm_stream_pd(y + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_mul_pd(_mm_loadu_pd(x + i + 6), alpha_vec)));
    }
    
    suffix = prefix + 2 * ((length - prefix) / 2);
    
    for ( ; i < suffix; i += 2) {
      _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_loadu_pd(x + i    ), alpha_vec)));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += alpha * x[i]; */
}

/* void misc_addAlignedVectorsInPlace_sse2(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t suffix = 8 * (length / 8);
  
  for ( ; i < suffix; i += 8) {
    _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
    _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_load_pd(x + i + 2)));
    _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_load_pd(x + i + 4)));
    _mm_stream_pd(y + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_load_pd(x + i + 6)));
  }
  
  suffix = 2 * (length / 2);
  
  for ( ; i < suffix; i += 2) {
    _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
  }
    
  for ( ; i < length; ++i)
    y[i] += x[i];
}

void misc_subtractAlignedVectorsInPlace_sse2(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;

  size_t i = 0;
  size_t suffix = 8 * (length / 8);
  
  for ( ; i < suffix; i += 8) {
    _mm_stream_pd(y + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
    _mm_stream_pd(y + i + 2, _mm_sub_pd(_mm_load_pd(y + i + 2), _mm_load_pd(x + i + 2)));
    _mm_stream_pd(y + i + 4, _mm_sub_pd(_mm_load_pd(y + i + 4), _mm_load_pd(x + i + 4)));
    _mm_stream_pd(y + i + 6, _mm_sub_pd(_mm_load_pd(y + i + 6), _mm_load_pd(x + i + 6)));
  }
  
  suffix = 2 * (length / 2);
  
  for ( ; i < suffix; i += 2) {
    _mm_stream_pd(y + i    , _mm_sub_pd(_mm_load_pd(y + i    ), _mm_load_pd(x + i    )));
  }
    
  for ( ; i < length; ++i)
    y[i] -= x[i];
}

void misc_addAlignedVectorsInPlaceWithMultiplier_sse2(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t suffix = 8 * (length / 8);
  
  __m128d alpha_vec = _mm_load1_pd(&alpha);
  for ( ; i < suffix; i += 8) {
    _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
    _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_mul_pd(_mm_load_pd(x + i + 2), alpha_vec)));
    _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_mul_pd(_mm_load_pd(x + i + 4), alpha_vec)));
    _mm_stream_pd(y + i + 6, _mm_add_pd(_mm_load_pd(y + i + 6), _mm_mul_pd(_mm_load_pd(x + i + 6), alpha_vec)));
  }
  
  suffix = 2 * (length / 2);
  
  for ( ; i < suffix; i += 2) {
    _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
  }
    
  for ( ; i < length; ++i)
    y[i] += x[i];
} */

void misc_addScalarToVectorInPlace_sse2(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) x[i] += alpha;
  
  for ( ; i < length; i += 4) {
    x[i    ] += alpha;
    x[i + 1] += alpha;
    x[i + 2] += alpha;
    x[i + 3] += alpha;
  }
  /*
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (2 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] += alpha;
  
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  __m128d alpha_vec = _mm_load1_pd(&alpha);
  for ( ; i < suffix; i += 8) {
    _mm_stream_pd(x + i    , _mm_add_pd(_mm_load_pd(x + i    ), alpha_vec));
    _mm_stream_pd(x + i + 2, _mm_add_pd(_mm_load_pd(x + i + 2), alpha_vec));
    _mm_stream_pd(x + i + 4, _mm_add_pd(_mm_load_pd(x + i + 4), alpha_vec));
    _mm_stream_pd(x + i + 6, _mm_add_pd(_mm_load_pd(x + i + 6), alpha_vec));
  }
  
  suffix = prefix + 2 * ((length - prefix) / 2);
  
  for ( ; i < suffix; i += 2) {
    _mm_stream_pd(x + i    , _mm_add_pd(_mm_load_pd(x + i    ), alpha_vec));
  }
  
  for ( ; i < length; ++i)
    x[i] += alpha;
  */
}

void misc_setVectorToConstant_sse2(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) x[i] = alpha;
  
  for ( ; i < length; i += 4) {
    x[i    ] = alpha;
    x[i + 1] = alpha;
    x[i + 2] = alpha;
    x[i + 3] = alpha;
  }
  /* 
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (2 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] = alpha;
  
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  __m128d alpha_vec = _mm_load1_pd(&alpha);
  for ( ; i < suffix; i += 8) {
    _mm_stream_pd(x + i    , alpha_vec);
    _mm_stream_pd(x + i + 2, alpha_vec);
    _mm_stream_pd(x + i + 4, alpha_vec);
    _mm_stream_pd(x + i + 6, alpha_vec);
  }

  suffix = prefix + 2 * ((length - prefix) / 2);
  
  for ( ; i < suffix; i += 2) {
    _mm_stream_pd(x + i    , alpha_vec);
  }
  
  for ( ; i < length; ++i)
    x[i] = alpha;
  */
}


static inline void transposeMatrixBlock(const double* restrict x, size_t ldx, double* restrict y, size_t ldy)
{
  // x:  0  4  8 12
  //     1  5  9 13
  //     2  6 10 14
  //     3  7 11 15
  // y:  0  1  2  3
  //     4  5  6  7
  //     8  9 10 11
  //    12 13 14 15
  
  __m128d temp0, temp1;

  temp0 = _mm_loadu_pd(x              ); // 0, 1
  temp1 = _mm_loadu_pd(x +         ldx); // 4, 5
  _mm_storeu_pd(y              , _mm_shuffle_pd(temp0, temp1, 0x0)); // 0, 4
  _mm_storeu_pd(y +         ldy, _mm_shuffle_pd(temp0, temp1, 0x3)); // 1 5
  
  temp0 = _mm_loadu_pd(x + 2          ); // 2, 3
  temp1 = _mm_loadu_pd(x + 2 +     ldx); // 6, 7
  _mm_storeu_pd(y +     2 * ldy, _mm_shuffle_pd(temp0, temp1, 0x0));
  _mm_storeu_pd(y +     3 * ldy, _mm_shuffle_pd(temp0, temp1, 0x3));
  
  temp0 = _mm_loadu_pd(x +     2 * ldx); // 8, 9
  temp1 = _mm_loadu_pd(x +     3 * ldx); // 12, 13
  _mm_storeu_pd(y + 2          , _mm_shuffle_pd(temp0, temp1, 0x0));
  _mm_storeu_pd(y + 2 +     ldy, _mm_shuffle_pd(temp0, temp1, 0x3));
  
  temp0 = _mm_loadu_pd(x + 2 + 2 * ldx); // 10, 11
  temp1 = _mm_loadu_pd(x + 2 + 3 * ldx); // 14, 15
  _mm_storeu_pd(y + 2 + 2 * ldy, _mm_shuffle_pd(temp0, temp1, 0x0));
  _mm_storeu_pd(y + 2 + 3 * ldy, _mm_shuffle_pd(temp0, temp1, 0x3));
}

void misc_transposeMatrix_sse2(const double* restrict x, size_t numRows, size_t numCols, double* restrict y)
{
  if (numRows == 0 || numCols == 0) return;
  
  // We can't really ensure that loads/stores occur on 16 byte boundaries, since
  // any time there is an odd number of rows that completely screws up
  // block transposing. For fun, we start x on a boundary and transpose
  // the first row explicitly, if necessary.
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = x_offset == 0 ? 0 : (2 * sizeof(double) - x_offset) / sizeof(double);
  prefix = prefix > numRows ? numRows : prefix;
  
  size_t row = 0;
  
  for ( ; row < prefix; ++row) {
    for (size_t col = 0; col < numCols; ++col) {
      y[col + row * numCols] = x[row + col * numRows];
    }
  }
  
  size_t suffix = prefix + 4 * ((numRows - prefix) / 4);
  
  if (suffix > prefix) {
    for ( ; row < suffix; row += 4) {
      size_t col = 0, colEnd = 4 * (numCols / 4);
      for ( ; col < colEnd; col += 4)
        transposeMatrixBlock(x + row + col * numRows, numRows, y + col + row * numCols, numCols);
      
      for (size_t rowInBlock = row; rowInBlock < row + 4; ++rowInBlock) {
        size_t colInBlock = col;
        for ( ; colInBlock < numCols; ++colInBlock) {
          y[colInBlock + rowInBlock * numCols] = x[rowInBlock + colInBlock * numRows];
        }
      }
    }
  }
  
  for ( ; row < numRows; ++row) {
    for (size_t col = 0; col < numCols; ++col) {
      y[col + row * numCols] = x[row + col * numRows];
    }
  }
}
