#pragma once
#include <cstdint>
#include <stdexcept>
#include "kmc_index.hpp"

template<class Processor>
inline void merge_two_indices(KMCIndex& A, KMCIndex& B, Processor& proc) {
  if (A.k != B.k) throw std::runtime_error("k mismatch between indices");
  CKmerAPI kA(A.k), kB(B.k);
  uint64_t cA=0, cB=0;
  bool okA = kmc_next(A.fh, kA, cA);
  bool okB = kmc_next(B.fh, kB, cB);
  uint64_t iters=0, common=0, onlyA=0, onlyB=0;

  while (okA && okB) {
    ++iters;
    if (kA == kB) {
      proc.on_pair(kA, cA, cB);
      ++common;
      okA = kmc_next(A.fh, kA, cA);
      okB = kmc_next(B.fh, kB, cB);
    } else if (kA < kB) {
      proc.on_pair(kA, cA, 0);
      ++onlyA;
      okA = kmc_next(A.fh, kA, cA);
    } else {
      proc.on_pair(kB, 0, cB);
      ++onlyB;
      okB = kmc_next(B.fh, kB, cB);
    }
  }
  while (okA) { ++iters; proc.on_pair(kA, cA, 0); ++onlyA; okA = kmc_next(A.fh, kA, cA); }
  while (okB) { ++iters; proc.on_pair(kB, 0, cB); ++onlyB; okB = kmc_next(B.fh, kB, cB); }
  proc.on_finish(iters, common, onlyA, onlyB);
}
