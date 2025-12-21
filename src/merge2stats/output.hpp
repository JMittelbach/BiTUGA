#pragma once
#include <cstdio>
#include <string>
#include <vector>
#include <stdexcept>
#include <inttypes.h>

struct Histogram {
  std::vector<uint64_t> bins;
  explicit Histogram(int n_individuals) : bins(n_individuals+1, 0) {}
  void add(int prev_indiv) {
    if (prev_indiv < 0) prev_indiv = 0;
    if ((size_t)prev_indiv >= bins.size()) prev_indiv = (int)bins.size()-1;
    bins[(size_t)prev_indiv]++;
  }
  static void write_table(FILE* fp, const char* header, const std::vector<uint64_t>& v){
    std::fprintf(fp, "# %s\nprevalence_individuals\tn_kmers\n", header);
    for (size_t i=0;i<v.size();++i)
      std::fprintf(fp, "%zu\t%" PRIu64 "\n", i, (uint64_t)v[i]);
  }
};

struct FASTAWriter {
  FILE* fp = nullptr;
  std::string tag0;
  std::string tag1;
  FASTAWriter(const std::string& path, const std::string& _tag0, const std::string& _tag1)
    : tag0(_tag0), tag1(_tag1)
  {
    fp = std::fopen(path.c_str(), "w");
    if (!fp) throw std::runtime_error("cannot open FASTA for writing: " + path);
  }
  ~FASTAWriter(){ if (fp) std::fclose(fp); }

  void write_record_buf(uint64_t id, uint32_t prev0, uint32_t prev1, const char* seq) {
    std::fprintf(fp, ">%" PRIu64 "%s%u%s%u\n", (uint64_t)id, tag0.c_str(), prev0, tag1.c_str(), prev1);
    std::fputs(seq, fp);
    std::fputc('\n', fp);
  }
};

inline FILE* open_log(const std::string& path) {
  FILE* fp = std::fopen(path.c_str(), "w");
  if (!fp) throw std::runtime_error("cannot open log file: " + path);
  return fp;
}
