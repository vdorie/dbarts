#include <external/stats.h>
#include <external/thread.h>

#include <math.h>

#include <Rmath.h>


size_t ext_drawFromDiscreteDistribution(const double* probabilities, size_t length)
{
  if (length == 0) return EXT_DISCRETE_DRAW_FAILURE;
  
  double u = unif_rand();
  
  size_t result = 0;
  double sum = probabilities[0];
  
  while (sum < u && result < length - 1) {
    sum += probabilities[++result];
  }
  
  if (result == length - 1 && sum < u) return EXT_DISCRETE_DRAW_FAILURE;
  return result;
}

int64_t ext_simulateIntegerUniformInRange(int64_t min_inclusive, int64_t max_exclusive)
{
  double range = fabs((double) (max_exclusive - min_inclusive));
  int64_t actualMin = (min_inclusive < max_exclusive ? min_inclusive : max_exclusive);
  
  double u = unif_rand();
  
  return actualMin + (int64_t) (u * range);
}

uint64_t ext_simulateUnsignedIntegerUniformInRange(uint64_t min_inclusive, uint64_t max_exclusive)
{
  uint64_t actualMin, actualMax;
  if (min_inclusive < max_exclusive) {
    actualMin = min_inclusive;
    actualMax = max_exclusive;
  } else {
    actualMin = max_exclusive;
    actualMax = min_inclusive;
  }
  
  double range = (double) (actualMax - actualMin);
  
  double u = unif_rand();
  
  return actualMin + (uint64_t) (u * range);
}

/*
double ext_computeAndSumSquaresOfResiduals(const double* x, size_t length, double x_hat)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) result += (x[i] - x_hat) * (x[i] - x_hat);
    if (length < 5) return result;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    result += (x[i] - x_hat) * (x[i] - x_hat) +
              (x[i + 1] - x_hat) * (x[i + 1] - x_hat) +
              (x[i + 2] - x_hat) * (x[i + 2] - x_hat) +
              (x[i + 3] - x_hat) * (x[i + 3] - x_hat) +
              (x[i + 4] - x_hat) * (x[i + 4] - x_hat);
  }
  return result;
} */
