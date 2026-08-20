#include "config.h"
#include <misc/linearAlgebra.h>

#include <misc/alloca.h>
#include <misc/intrinsic.h>

void misc_addVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) y[i] += x[i];
  
  for ( ; i < length; i += 8) {
    y[i    ] += x[i    ];
    y[i + 1] += x[i + 1];
    y[i + 2] += x[i + 2];
    y[i + 3] += x[i + 3];
    y[i + 4] += x[i + 4];
    y[i + 5] += x[i + 5];
    y[i + 6] += x[i + 6];
    y[i + 7] += x[i + 7];
  }
}
  

void misc_subtractVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) y[i] -= x[i];
  
  for ( ; i < length; i += 8) {
    y[i    ] -= x[i    ];
    y[i + 1] -= x[i + 1];
    y[i + 2] -= x[i + 2];
    y[i + 3] -= x[i + 3];
    y[i + 4] -= x[i + 4];
    y[i + 5] -= x[i + 5];
    y[i + 6] -= x[i + 6];
    y[i + 7] -= x[i + 7];
  }
}
  

void misc_addVectorsInPlaceWithMultiplier_avx(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) y[i] += alpha * x[i];
  
  for ( ; i < length; i += 8) {
    y[i    ] += alpha * x[i    ];
    y[i + 1] += alpha * x[i + 1];
    y[i + 2] += alpha * x[i + 2];
    y[i + 3] += alpha * x[i + 3];
    y[i + 4] += alpha * x[i + 4];
    y[i + 5] += alpha * x[i + 5];
    y[i + 6] += alpha * x[i + 6];
    y[i + 7] += alpha * x[i + 7];
  }
}
  

void misc_addScalarToVectorInPlace_avx(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) x[i] += alpha;
  
  for ( ; i < length; i += 8) {
    x[i    ] += alpha;
    x[i + 1] += alpha;
    x[i + 2] += alpha;
    x[i + 3] += alpha;
    x[i + 4] += alpha;
    x[i + 5] += alpha;
    x[i + 6] += alpha;
    x[i + 7] += alpha;
  }
  
}


void misc_setVectorToConstant_avx(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) x[i] = alpha;
  
  for ( ; i < length; i += 8) {
    x[i    ] = alpha;
    x[i + 1] = alpha;
    x[i + 2] = alpha;
    x[i + 3] = alpha;
    x[i + 4] = alpha;
    x[i + 5] = alpha;
    x[i + 6] = alpha;
    x[i + 7] = alpha;
  }
}
