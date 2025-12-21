#pragma once
#include <string>
#include <stdexcept>
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"

struct KMCIndex {
  CKMCFile fh;
  CKMCFileInfo info{};
  uint32_t k = 0;

  explicit KMCIndex(const std::string& basename) {
    if (!fh.OpenForListing(basename)) {
      throw std::runtime_error("cannot open KMC index: " + basename);
    }
    if (!fh.Info(info)) {
      throw std::runtime_error("cannot read KMC Info: " + basename);
    }
    k = info.kmer_length;
    if (k == 0) {
      throw std::runtime_error("k=0 in KMC index: " + basename);
    }
  }
  ~KMCIndex() { fh.Close(); }
};

inline bool kmc_next(CKMCFile& fh, CKmerAPI& kmer, uint64_t& count) {
  return fh.ReadNextKmer(kmer, count);
}
