#ifndef SEQUENCE_BIAS_GET_HPP
#define SEQUENCE_BIAS_GET_HPP
#include <cstddef>
#include <array>
#include <cassert>
#include "sequences/literate_multiseq.hpp"

static inline double sequence_bias_get(const GttlMultiseq &multiseq)
{
  static constexpr const double bias_factor[] = {.690, .690, .690, .690, .780,
                                                 .850, .900, .933, .966, 1.000};

  constexpr const int undefined_rank = 4;
  static constexpr const char nucleotides_upper_lower[] = "Aa|Cc|Gg|TtUu";
  LiterateMultiseq<nucleotides_upper_lower,undefined_rank>
    literate_multiseq(multiseq);
  const std::array<size_t,undefined_rank+1> &rank_dist
    = literate_multiseq.rank_dist_get();

  const size_t at_count = rank_dist[0] + rank_dist[2];
  const size_t gc_count = rank_dist[1] + rank_dist[3];
  const double ratio = static_cast<double>(std::min(at_count, gc_count))/
                       (at_count + gc_count);
  const size_t bias_index
    = static_cast<int>(std::max(0.0, (ratio + 0.025) * 20.0 - 1.0));
  assert(bias_index < sizeof bias_factor/sizeof bias_factor[0]);
  return bias_factor[bias_index];
}
#endif
