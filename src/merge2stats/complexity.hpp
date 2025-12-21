#pragma once
#include <cstddef>
#include <cmath>
#include <algorithm>

inline size_t longest_homopolymer_run(const char* s, size_t n) {
  if (n == 0) return 0;
  size_t best = 1, cur = 1;
  for (size_t i = 1; i < n; ++i) {
    if (s[i] == s[i - 1]) {
      cur++;
      if (cur > best) best = cur;
    } else {
      cur = 1;
    }
  }
  return best;
}

inline size_t homopolymer_allow_run(size_t k, double frac_max, size_t hard_cap = 15) {
  if (frac_max <= 0.0) return 0;
  size_t dynamic_cap = static_cast<size_t>(std::floor(frac_max * static_cast<double>(k)));
  if (dynamic_cap < 1) dynamic_cap = 1;
  if (dynamic_cap > hard_cap) dynamic_cap = hard_cap;
  return dynamic_cap;
}

inline bool homopolymer_exceeds(const char* s, size_t n, double frac_max, size_t hard_cap = 15) {
  if (n == 0 || frac_max <= 0.0) return false;
  const size_t limit = homopolymer_allow_run(n, frac_max, hard_cap);
  const size_t run = longest_homopolymer_run(s, n);
  return run > limit;
}

inline bool is_periodic2(const char* s, size_t n) {
  if (n < 2) return false;
  for (size_t i = 2; i < n; ++i) if (s[i] != s[i % 2]) return false;
  return true;
}
inline bool is_periodic3(const char* s, size_t n) {
  if (n < 3) return false;
  for (size_t i = 3; i < n; ++i) if (s[i] != s[i % 3]) return false;
  return true;
}
inline double shannon_entropy_bits(const char* s, size_t n) {
  int cnt[4]={0,0,0,0}; size_t valid=0;
  auto idx=[](char c)->int{
    switch(c){case 'A':case 'a': return 0; case 'C':case 'c': return 1;
              case 'G':case 'g': return 2; case 'T':case 't': return 3;
              default: return -1;}
  };
  for (size_t i=0;i<n;++i) { int k=idx(s[i]); if (k>=0) { cnt[k]++; valid++; } }
  if (!valid) return 0.0;
  double H=0.0;
  for(int i=0;i<4;i++) if (cnt[i]) { double p = double(cnt[i])/valid; H -= p*std::log2(p); }
  return H;
}
