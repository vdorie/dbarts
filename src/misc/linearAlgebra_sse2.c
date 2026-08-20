#include "config.h"
#include <misc/linearAlgebra.h>

#include <misc/intrinsic.h>

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
  
}

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
}
