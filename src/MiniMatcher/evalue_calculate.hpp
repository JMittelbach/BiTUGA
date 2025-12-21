#ifndef EVALUE_CALCULATE_HPP
#define EVALUE_CALCULATE_HPP
#include <cstdlib>
#include <cassert>
#include <cmath>

static inline double log_binomial_coefficent(size_t l, size_t k)
{
  assert(k <= l) ;
  return std::lgamma(l+1) - std::lgamma(k+1) - std::lgamma(l-k+1);
}

static inline double evalue_calculate(size_t n, size_t m, size_t l, size_t k)
{
  assert(k <= l);
  static constexpr const double p = 0.25;
  return std::exp(std::log(n) +
                  std::log(m) +
                  log_binomial_coefficent(l,k) +
                  (l-k) * std::log(p) + (k+2) * std::log(1-p));
}
#endif
