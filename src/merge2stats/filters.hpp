#pragma once
#include <cstdint>
#include <cmath>
#include "complexity.hpp"
#include "kmc_api/kmer_api.h"

struct PrevThresholds {
  uint64_t need_all_counts = 0;
  uint64_t allow_all_counts = 0;
  uint64_t need_min_A = 0;
  uint64_t allow_max_A = 0;
  uint64_t need_min_B = 0;
  uint64_t allow_max_B = 0;
  uint64_t need_all_indiv = 0;
  uint64_t allow_all_indiv = 0;
  uint64_t need_min_A_ind = 0;
  uint64_t allow_max_A_ind = 0;
  uint64_t need_min_B_ind = 0;
  uint64_t allow_max_B_ind = 0;
};

inline PrevThresholds make_thresholds(
  int nA, int nB,
  double prev_min, double prev_max,
  double group_min, double group_max)
{
  PrevThresholds t{};
  const int N = nA + nB;
  auto to_cnt = [](uint64_t indiv)->uint64_t { return 2ull * indiv; };

  t.need_all_indiv  = (uint64_t)std::ceil (prev_min  * (double)N);
  t.allow_all_indiv = (uint64_t)std::floor(prev_max  * (double)N);

  t.need_min_A_ind  = (uint64_t)std::ceil (group_min * (double)nA);
  t.allow_max_A_ind = (uint64_t)std::floor(group_max * (double)nA);
  t.need_min_B_ind  = (uint64_t)std::ceil (group_min * (double)nB);
  t.allow_max_B_ind = (uint64_t)std::floor(group_max * (double)nB);

  t.need_all_counts  = to_cnt(t.need_all_indiv);
  t.allow_all_counts = to_cnt(t.allow_all_indiv);
  t.need_min_A  = to_cnt(t.need_min_A_ind);
  t.allow_max_A = to_cnt(t.allow_max_A_ind);
  t.need_min_B  = to_cnt(t.need_min_B_ind);
  t.allow_max_B = to_cnt(t.allow_max_B_ind);
  return t;
}

struct OverallPrevFilter {
  uint64_t need_all_counts{};
  uint64_t allow_all_counts{};
  bool pass(uint64_t cA, uint64_t cB) const {
    const uint64_t tot = cA + cB;
    return tot >= need_all_counts && tot <= allow_all_counts;
  }
};

struct SymmetricContrastFilter {
  uint64_t need_min_A{}, allow_max_B{}, need_min_B{}, allow_max_A{};
  bool enabled = false;
  bool pass(uint64_t cA, uint64_t cB) const {
    if (!enabled) return true;
    const bool A_high_B_low = (cA >= need_min_A) && (cB <= allow_max_B);
    const bool B_high_A_low = (cB >= need_min_B) && (cA <= allow_max_A);
    return A_high_B_low || B_high_A_low;
  }
};

struct ComplexityFilter {
  bool check_homo = true;
  bool check_di   = true;
  bool check_tri  = true;
  double entropy_min = 1.5;
  double homo_max_frac = 0.33;
}
;
