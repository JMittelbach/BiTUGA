#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <stdexcept>
#include <vector>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <sys/resource.h>
#include <sys/time.h>
#include <cctype>
#include <inttypes.h>
#include <memory>

#include "kmc_index.hpp"
#include "merge_engine.hpp"
#include "filters.hpp"
#include "output.hpp"
#include "kmc_api/kmer_api.h"

struct Config {
  std::string results_dir = "results";
  std::string logging_dir = "logging";
  std::string base0 = "";
  std::string base1 = "";
  std::string trait0 = "female";
  std::string trait1 = "male";
  int n0 = -1, n1 = -1;
  double prev_min = 0.25;
  double prev_max = 0.75;
  bool prev_min_set = false;
  bool prev_max_set = false;
  std::string adaptive_prev_mode = "auto"; // auto | off
  double group_min = 0.00;
  double group_max = 1.00;
  bool perfect_only = false;
  bool perfect_extra = false;
  bool homopolymer_filter = true;
  double homopolymer_frac = 0.33;
  double entropy_min = 1.5;
  std::string out_prefix = "";
};

static std::string rel_path(const std::string& path, const std::string& base) {
  if (base.empty()) return path;
  std::string b = base;
  if (!b.empty() && b.back() == '/') b.pop_back();
  if (path.rfind(b, 0) == 0) {
    std::string rest = path.substr(b.size());
    if (rest.empty()) return ".";
    if (rest[0] != '/') rest = "/" + rest;
    return "." + rest;
  }
  return path;
}

static void usage(const char* prog) {
  std::fprintf(stderr,
    "Usage:\n"
    "  %s --n0 INT --n1 INT [options]\n\n"
    "Options:\n"
    "  --output DIR                base directory; results=DIR/results, logging=DIR/logging\n"
    "  --results-dir DIR           default: results\n"
    "  --logging-dir DIR           default: logging\n"
    "  --index0 BASENAME           default: <results>/pooled_kmc_files/pooled_<trait0>\n"
    "  --index1 BASENAME           default: <results>/pooled_kmc_files/pooled_<trait1>\n"
    "  --trait0 NAME               default: female\n"
    "  --trait1 NAME               default: male\n"
    "  --n0 INT                    #individuals for trait_0 (required)\n"
    "  --n1 INT                    #individuals for trait_1 (required)\n"
    "  --prev-min FLOAT            default: 0.25 (global min prevalence)\n"
    "  --prev-max FLOAT            default: 0.75 (global max prevalence)\n"
    "  --group-min FLOAT           default: 0.00 (hi-group min fraction when using contrast)\n"
    "  --group-max FLOAT           default: 1.00 (lo-group max fraction when using contrast)\n"
    "  --perfect-only              emit only kmers where one group is 100%% and the other 0%%\n"
    "  --perfect-extra             also emit those max-contrast kmers to a separate FASTA\n"
    "  --shannon-entropy-kmer FLOAT default: 1.5 (<=0 disables entropy check)\n"
    "  --entropy-min FLOAT         alias for --shannon-entropy-kmer\n"
    "  --homopolymer-max-frac FLOAT default: 0.33 (max allowed consecutive base run, capped at 15)\n"
    "  --no-homopolymer-filter     disable homopolymer filter entirely\n"
    "  --adaptive-prev-thresholds auto|off default: auto (auto adjusts defaults for asymmetric groups)\n"
    "  --disable-adaptive-defaults deprecated alias for --adaptive-prev-thresholds off\n"
    "  --force-static-defaults     deprecated alias for --adaptive-prev-thresholds off\n"
    "  --out-prefix PATH           default: <results>/candidate_kmers\n",
    prog);
}

static std::string sanitize(const std::string& s) {
  std::string t; t.reserve(s.size());
  for (unsigned char c: s) if (std::isalnum(c)) t.push_back(c);
  return t;
}

static std::pair<std::string,std::string> unique_prefixes(std::string a, std::string b) {
  auto S = [](const std::string& x){ std::string y; y.reserve(x.size()); for (unsigned char c: x) if (std::isalnum(c)) y.push_back(c); return y; };
  a = S(a); b = S(b);
  if (a.empty() && b.empty()) return {"0","1"};
  if (a.empty()) a = "0";
  if (b.empty()) b = "1";
  std::string al=a, bl=b;
  for (auto& c: al) c = (char)std::tolower((unsigned char)c);
  for (auto& c: bl) c = (char)std::tolower((unsigned char)c);
  size_t i=1, L=std::max(al.size(), bl.size());
  while (i<=L) {
    std::string pa = al.substr(0,i);
    std::string pb = bl.substr(0,i);
    if (pa != pb) return {a.substr(0,std::min(i,a.size())), b.substr(0,std::min(i,b.size()))};
    i++;
  }
  return {a + "0", b + "1"};
}

static Config parse_args(int argc, char** argv) {
  Config cfg;
  bool results_dir_set=false, logging_dir_set=false, base0_set=false, base1_set=false, out_prefix_set=false;
  std::string output_dir;

  for (int i=1;i<argc;i++) {
    auto need = [&](const char* name)->const char*{
      if (i+1>=argc) { std::fprintf(stderr,"missing value after %s\n",name); usage(argv[0]); std::exit(1); }
      return argv[++i];
    };
    if      (!std::strcmp(argv[i],"--output"))        output_dir = need("--output");
    else if (!std::strcmp(argv[i],"--results-dir")) { cfg.results_dir = need("--results-dir"); results_dir_set=true; }
    else if (!std::strcmp(argv[i],"--logging-dir")) { cfg.logging_dir = need("--logging-dir"); logging_dir_set=true; }
    else if (!std::strcmp(argv[i],"--index0"))      { cfg.base0 = need("--index0"); base0_set=true; }
    else if (!std::strcmp(argv[i],"--index1"))      { cfg.base1 = need("--index1"); base1_set=true; }
    else if (!std::strcmp(argv[i],"--trait0"))      { cfg.trait0 = need("--trait0"); }
    else if (!std::strcmp(argv[i],"--trait1"))      { cfg.trait1 = need("--trait1"); }
    else if (!std::strcmp(argv[i],"--n0"))          cfg.n0 = std::atoi(need("--n0"));
    else if (!std::strcmp(argv[i],"--n1"))          cfg.n1 = std::atoi(need("--n1"));
    else if (!std::strcmp(argv[i],"--prev-min"))  { cfg.prev_min = std::atof(need("--prev-min")); cfg.prev_min_set = true; }
    else if (!std::strcmp(argv[i],"--prev-max"))  { cfg.prev_max = std::atof(need("--prev-max")); cfg.prev_max_set = true; }
    else if (!std::strcmp(argv[i],"--group-min"))   cfg.group_min = std::atof(need("--group-min"));
    else if (!std::strcmp(argv[i],"--group-max"))   cfg.group_max = std::atof(need("--group-max"));
    else if (!std::strcmp(argv[i],"--perfect-only"))  cfg.perfect_only = true;
    else if (!std::strcmp(argv[i],"--perfect-extra")) cfg.perfect_extra = true;
    else if (!std::strcmp(argv[i],"--shannon-entropy-kmer")) { cfg.entropy_min = std::atof(need("--shannon-entropy-kmer")); }
    else if (!std::strcmp(argv[i],"--entropy-min")) { cfg.entropy_min = std::atof(need("--entropy-min")); }
    else if (!std::strcmp(argv[i],"--homopolymer-max-frac")) cfg.homopolymer_frac = std::atof(need("--homopolymer-max-frac"));
    else if (!std::strcmp(argv[i],"--no-homopolymer-filter")) cfg.homopolymer_filter = false;
    else if (!std::strcmp(argv[i],"--adaptive-prev-thresholds")) {
      std::string mode = need("--adaptive-prev-thresholds");
      std::transform(mode.begin(), mode.end(), mode.begin(), [](unsigned char c){ return (char)std::tolower(c); });
      if (mode == "off" || mode == "none" || mode == "disabled") mode = "off";
      else mode = "auto";
      cfg.adaptive_prev_mode = mode;
    }
    else if (!std::strcmp(argv[i],"--disable-adaptive-defaults") || !std::strcmp(argv[i],"--force-static-defaults")) {
      cfg.adaptive_prev_mode = "off";
    }
    else if (!std::strcmp(argv[i],"--out-prefix"))  { cfg.out_prefix = need("--out-prefix"); out_prefix_set=true; }
    else {
      std::fprintf(stderr,"unknown option: %s\n", argv[i]);
      usage(argv[0]); std::exit(1);
    }
  }

  if (cfg.n0 < 0 || cfg.n1 < 0) {
    std::fprintf(stderr,"ERROR: please provide --n0 <N_trait_0> and --n1 <N_trait_1>\n");
    usage(argv[0]); std::exit(1);
  }

  if (cfg.perfect_only && cfg.perfect_extra) {
    std::fprintf(stderr,"ERROR: please choose either --perfect-only or --perfect-extra, not both\n");
    usage(argv[0]); std::exit(1);
  }

  if (!output_dir.empty()) {
    const bool has_slash = !output_dir.empty() && output_dir.back()=='/';
    if (!results_dir_set)  cfg.results_dir = output_dir + (has_slash ? "results" : "/results");
    if (!logging_dir_set)  cfg.logging_dir = output_dir + (has_slash ? "logging" : "/logging");
  }

  if (!out_prefix_set)
    cfg.out_prefix = cfg.results_dir + "/candidate_kmers";

  if (!base0_set) {
    std::string t0 = sanitize(cfg.trait0);
    if (t0.empty()) t0 = "trait0";
    cfg.base0 = cfg.results_dir + "/pooled_kmc_files/pooled_" + t0;
  }
  if (!base1_set) {
    std::string t1 = sanitize(cfg.trait1);
    if (t1.empty()) t1 = "trait1";
    cfg.base1 = cfg.results_dir + "/pooled_kmc_files/pooled_" + t1;
  }

  return cfg;
}

struct Processor {
  KMCIndex& A; KMCIndex& B;
  const int n0, n1;
  Histogram H0, H1;
  Histogram Hdiff;
  OverallPrevFilter f_all;
  SymmetricContrastFilter f_contrast;
  ComplexityFilter f_comp;
  FASTAWriter* main_fasta;
  FILE* logfp;
  FASTAWriter* perfect_fasta;
  bool perfect_only_mode;

  std::string trait0_name;
  std::string trait1_name;
  const size_t k;
  std::vector<char> kbuf;

  uint64_t emitted = 0, id = 0;
  uint64_t perfect_emitted = 0, perfect_id = 0;
  uint64_t pre_bias_0 = 0, pre_bias_1 = 0, pre_ties = 0;
  uint64_t pre_all0_01 = 0, pre_all1_00 = 0;
  uint64_t post_bias_0 = 0, post_bias_1 = 0, post_ties = 0;
  uint64_t post_all0_01 = 0, post_all1_00 = 0;
  uint64_t lc_tested = 0, lc_pass = 0;
  uint64_t lc_reject_homo = 0, lc_reject_di = 0, lc_reject_tri = 0, lc_reject_entropy = 0;
  uint64_t iter_count = 0;
  uint64_t common_kmers = 0;
  uint64_t only0_kmers = 0;
  uint64_t only1_kmers = 0;

  Processor(KMCIndex& _A, KMCIndex& _B,
            int _n0, int _n1,
            const PrevThresholds& thr,
            const ComplexityFilter& comp,
            FASTAWriter* main_writer,
            FILE* _logfp,
            FASTAWriter* perfect_writer,
            bool perfect_only,
            std::string t0_name,
            std::string t1_name)
    : A(_A), B(_B), n0(_n0), n1(_n1),
      H0(_n0), H1(_n1),
      Hdiff(std::max(_n0,_n1)),
      main_fasta(main_writer),
      logfp(_logfp),
      perfect_fasta(perfect_writer),
      perfect_only_mode(perfect_only),
      trait0_name(std::move(t0_name)),
      trait1_name(std::move(t1_name)),
      k(_A.k),
      kbuf(_A.k + 1, '\0')
  {
    f_all.need_all_counts  = thr.need_all_counts;
    f_all.allow_all_counts = thr.allow_all_counts;
    f_contrast.need_min_A = thr.need_min_A;
    f_contrast.allow_max_B = thr.allow_max_B;
    f_contrast.need_min_B = thr.need_min_B;
    f_contrast.allow_max_A = thr.allow_max_A;
    f_comp = comp;
  }

  inline void on_pair(const CKmerAPI& kmer, uint64_t c0, uint64_t c1) {
    const int p0 = static_cast<int>(c0 / 2u);
    const int p1 = static_cast<int>(c1 / 2u);
    H0.add(p0); H1.add(p1);
    const int diff = std::abs(p0 - p1);
    Hdiff.add(diff);
    if (p0 > p1) pre_bias_0++;
    else if (p1 > p0) pre_bias_1++;
    else pre_ties++;
    if (p0 == n0 && p1 == 0) pre_all0_01++;
    if (p1 == n1 && p0 == 0) pre_all1_00++;

    if (!perfect_only_mode) {
      if (!f_all.pass(c0, c1)) return;
      if (!f_contrast.pass(c0, c1)) return;
    }

    kmer.to_string(kbuf.data());

    lc_tested++;
    bool fail = false;
    if (f_comp.check_homo) {
      if (homopolymer_exceeds(kbuf.data(), k, f_comp.homo_max_frac)) {
        lc_reject_homo++;
        fail = true;
      }
    }
    if (f_comp.check_di && is_periodic2(kbuf.data(), k)) { lc_reject_di++; fail = true; }
    if (f_comp.check_tri && is_periodic3(kbuf.data(), k)) { lc_reject_tri++; fail = true; }
    if (f_comp.entropy_min > 0.0 &&
        shannon_entropy_bits(kbuf.data(), k) < f_comp.entropy_min)    { lc_reject_entropy++; fail = true; }

    if (fail) return;
    lc_pass++;

    const bool perfect = (p0 == n0 && p1 == 0) || (p1 == n1 && p0 == 0);
    if (perfect_only_mode && !perfect) return;

    if (p0 > p1) post_bias_0++;
    else if (p1 > p0) post_bias_1++;
    else post_ties++;
    if (p0 == n0 && p1 == 0) post_all0_01++;
    if (p1 == n1 && p0 == 0) post_all1_00++;

    ++id;
    if (main_fasta) {
      main_fasta->write_record_buf(id, static_cast<uint32_t>(p0), static_cast<uint32_t>(p1), kbuf.data());
    }
    ++emitted;
    if (perfect && perfect_fasta) {
      perfect_fasta->write_record_buf(++perfect_id, static_cast<uint32_t>(p0), static_cast<uint32_t>(p1), kbuf.data());
      ++perfect_emitted;
    }
  }

  inline void on_finish(uint64_t iters, uint64_t common, uint64_t only0, uint64_t only1) {
    iter_count = iters;
    common_kmers = common;
    only0_kmers = only0;
    only1_kmers = only1;
    const uint64_t union_total = common + only0 + only1;

    std::fprintf(logfp, "\n=== kmer_overlap ===\n");
    std::fprintf(logfp, "k\t%u\n", A.k);
    std::fprintf(logfp, "unique_kmers_%s\t%" PRIu64 "\n", trait0_name.c_str(), (uint64_t)A.info.total_kmers);
    std::fprintf(logfp, "unique_kmers_%s\t%" PRIu64 "\n", trait1_name.c_str(), (uint64_t)B.info.total_kmers);
    std::fprintf(logfp, "total_unique_kmers\t%" PRIu64 "\n", (uint64_t)iters);
    std::fprintf(logfp, "overlap_common\t%" PRIu64 "\n", (uint64_t)common);
    std::fprintf(logfp, "only_%s\t%" PRIu64 "\n", trait0_name.c_str(), (uint64_t)only0);
    std::fprintf(logfp, "only_%s\t%" PRIu64 "\n", trait1_name.c_str(), (uint64_t)only1);
    std::fprintf(logfp, "union_total\t%" PRIu64 "\n", (uint64_t)union_total);

    std::fprintf(logfp, "\n=== prevalence_summary ===\n");
    std::fprintf(logfp, "pre_bias_%s\t%" PRIu64 "\n", trait0_name.c_str(), (uint64_t)pre_bias_0);
    std::fprintf(logfp, "pre_bias_%s\t%" PRIu64 "\n", trait1_name.c_str(), (uint64_t)pre_bias_1);
    std::fprintf(logfp, "pre_ties\t%" PRIu64 "\n",        (uint64_t)pre_ties);
    std::fprintf(logfp, "pre_all%s_0%s\t%" PRIu64 "\n",   trait0_name.c_str(), trait1_name.c_str(), (uint64_t)pre_all0_01);
    std::fprintf(logfp, "pre_all%s_0%s\t%" PRIu64 "\n",   trait1_name.c_str(), trait0_name.c_str(), (uint64_t)pre_all1_00);
    std::fprintf(logfp, "post_bias_%s\t%" PRIu64 "\n", trait0_name.c_str(), (uint64_t)post_bias_0);
    std::fprintf(logfp, "post_bias_%s\t%" PRIu64 "\n", trait1_name.c_str(), (uint64_t)post_bias_1);
    std::fprintf(logfp, "post_ties\t%" PRIu64 "\n",        (uint64_t)post_ties);
    std::fprintf(logfp, "post_all%s_0%s\t%" PRIu64 "\n",   trait0_name.c_str(), trait1_name.c_str(), (uint64_t)post_all0_01);
    std::fprintf(logfp, "post_all%s_0%s\t%" PRIu64 "\n",   trait1_name.c_str(), trait0_name.c_str(), (uint64_t)post_all1_00);

    std::fprintf(logfp, "\n=== lc_summary_on_numeric_pass ===\n");
    std::fprintf(logfp, "lc_tested\t%" PRIu64 "\n",         (uint64_t)lc_tested);
    std::fprintf(logfp, "lc_pass\t%" PRIu64 "\n",           (uint64_t)lc_pass);
    std::fprintf(logfp, "lc_reject_homopolymer\t%" PRIu64 "\n",    (uint64_t)lc_reject_homo);
    std::fprintf(logfp, "lc_reject_dinucleotide_repeat\t%" PRIu64 "\n",      (uint64_t)lc_reject_di);
    std::fprintf(logfp, "lc_reject_trinucleotide_repeat\t%" PRIu64 "\n",     (uint64_t)lc_reject_tri);
    std::fprintf(logfp, "lc_reject_entropy\t%" PRIu64 "\n", (uint64_t)lc_reject_entropy);

    std::fprintf(logfp, "\n=== output_counts ===\n");
    std::fprintf(logfp, "emitted_kmers\t%" PRIu64 "\n", (uint64_t)emitted);
  }
};

int main(int argc, char** argv) {
  try {
    const auto t0 = std::chrono::steady_clock::now();

    const Config cfg = parse_args(argc, argv);

    const auto tags = unique_prefixes(cfg.trait0, cfg.trait1);
    const std::string tag0 = tags.first;
    const std::string tag1 = tags.second;

    const std::string fasta_path = cfg.out_prefix + ".fasta";
    const std::string stats_path = cfg.results_dir + "/statistics/group_prevalence_histogram.txt";
    const std::string log_path   = cfg.logging_dir + "/candidate_kmers.log";
    const std::string perfect_fasta_path = cfg.results_dir + "/kmers_perfect_contrast.fasta";
    const std::string disp_fasta = rel_path(fasta_path, cfg.results_dir);
    const std::string disp_stats = rel_path(stats_path, cfg.results_dir);
    const std::string disp_perfect = rel_path(perfect_fasta_path, cfg.results_dir);

    FILE* logfp = open_log(log_path);
    std::fprintf(logfp, "=== Pass II — candidate k-mer selection (merge2stats) ===\n");

    KMCIndex A(cfg.base0), B(cfg.base1);
    if (A.k != B.k) throw std::runtime_error("k mismatch between indices");

    const int total_samples = cfg.n0 + cfg.n1;
    const int minority = std::min(cfg.n0, cfg.n1);
    const double minority_ratio = (total_samples > 0) ? (double)minority / (double)total_samples : 0.0;
    const double adaptive_min = minority_ratio * 0.5;
    const double adaptive_max = 1.0 - adaptive_min;

    const double base_prev_min = cfg.prev_min_set ? cfg.prev_min : 0.25;
    const double base_prev_max = cfg.prev_max_set ? cfg.prev_max : 0.75;

    double eff_prev_min = base_prev_min;
    double eff_prev_max = base_prev_max;
    if (cfg.adaptive_prev_mode == "auto" && !cfg.prev_min_set && !cfg.prev_max_set && cfg.n0 != cfg.n1) {
      eff_prev_min = adaptive_min;
      eff_prev_max = adaptive_max;
    }

    PrevThresholds thr = make_thresholds(cfg.n0, cfg.n1, eff_prev_min, eff_prev_max, cfg.group_min, cfg.group_max);
    ComplexityFilter comp;
    comp.check_homo  = cfg.homopolymer_filter;
    comp.check_di    = true;  // always on
    comp.check_tri   = true;  // always on
    comp.entropy_min = cfg.entropy_min;
    comp.homo_max_frac = cfg.homopolymer_frac;

    const int N = cfg.n0 + cfg.n1;

    const double eff_total_min_frac = N > 0 ? double(thr.need_all_indiv)  / double(N)  : 0.0;
    const double eff_total_max_frac = N > 0 ? double(thr.allow_all_indiv) / double(N)  : 0.0;
    const double eff_g0_min_frac    = cfg.n0 > 0 ? double(thr.need_min_A_ind)  / double(cfg.n0) : 0.0;
    const double eff_g0_max_frac    = cfg.n0 > 0 ? double(thr.allow_max_A_ind) / double(cfg.n0) : 0.0;
    const double eff_g1_min_frac    = cfg.n1 > 0 ? double(thr.need_min_B_ind)  / double(cfg.n1) : 0.0;
    const double eff_g1_max_frac    = cfg.n1 > 0 ? double(thr.allow_max_B_ind) / double(cfg.n1) : 0.0;

    const bool contrast_enabled = (cfg.group_min > 0.0) && (cfg.group_max < 1.0);
    const char* perfect_mode_str = cfg.perfect_only ? "only" : (cfg.perfect_extra ? "extra" : "off");
    const char* global_filters_applied = cfg.perfect_only ? "no (perfect-only)" : "yes";

    std::unique_ptr<FASTAWriter> main_writer;
    if (!cfg.perfect_only)
      main_writer.reset(new FASTAWriter(fasta_path, tag0, tag1));
    std::unique_ptr<FASTAWriter> perfect_writer;
    if (cfg.perfect_only || cfg.perfect_extra) {
      perfect_writer.reset(new FASTAWriter(perfect_fasta_path, tag0, tag1));
    }

    Processor proc(A, B, cfg.n0, cfg.n1, thr, comp, main_writer.get(), logfp, perfect_writer.get(), cfg.perfect_only, cfg.trait0, cfg.trait1);
    proc.f_contrast.enabled = contrast_enabled;
    merge_two_indices(proc.A, proc.B, proc);

    const uint64_t perfect_bias_trait0 = proc.post_all0_01;
    const uint64_t perfect_bias_trait1 = proc.post_all1_00;
    const uint64_t perfect_total = perfect_bias_trait0 + perfect_bias_trait1;
    const uint64_t union_total = proc.common_kmers + proc.only0_kmers + proc.only1_kmers;

    FILE* statsfp = std::fopen(stats_path.c_str(), "w");
    if (!statsfp) throw std::runtime_error("cannot open stats file: " + stats_path);
    Histogram::write_table(statsfp, cfg.trait0.c_str(), proc.H0.bins);
    Histogram::write_table(statsfp, cfg.trait1.c_str(), proc.H1.bins);
    std::fprintf(statsfp, "# diff_abs(|%s - %s|)\nprevalence_diff\t n_kmers\n", cfg.trait0.c_str(), cfg.trait1.c_str());
    for (size_t i=0;i<proc.Hdiff.bins.size();++i) {
      std::fprintf(statsfp, "%zu\t%" PRIu64 "\n", i, (uint64_t)proc.Hdiff.bins[i]);
    }
    std::fclose(statsfp);

    std::fprintf(logfp, "\n=== perfect_contrast_results ===\n");
    std::fprintf(logfp, "perfect_mode\t%s\n", perfect_mode_str);
    std::fprintf(logfp, "perfect_detected\t%s\n", perfect_total ? "yes" : "no");
    std::fprintf(logfp, "perfect_total\t%" PRIu64 "\n", perfect_total);
    std::fprintf(logfp, "perfect_bias_%s\t%" PRIu64 "\n", cfg.trait0.c_str(), perfect_bias_trait0);
    std::fprintf(logfp, "perfect_bias_%s\t%" PRIu64 "\n", cfg.trait1.c_str(), perfect_bias_trait1);
    if (cfg.perfect_only || cfg.perfect_extra) {
      std::fprintf(logfp, "perfect_emitted\t%" PRIu64 "\n", (uint64_t)proc.perfect_emitted);
    }

    struct rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    double ut = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1e6;
    double st = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1e6;
    double cpu_total = ut + st;
    double rss_mb =
    #ifdef __APPLE__
      ru.ru_maxrss / (1024.0 * 1024.0);
    #else
      ru.ru_maxrss / 1024.0;
    #endif
    const auto t1 = std::chrono::steady_clock::now();
    double wall = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0;

    std::fprintf(logfp, "\n=== summary ===\n");
    std::fprintf(logfp, "kmer_union_total\t%" PRIu64 "\n", union_total);
    std::fprintf(logfp, "kmer_emitted\t%" PRIu64 "\n", (uint64_t)proc.emitted);
    std::fprintf(logfp, "perfect_total\t%" PRIu64 "\n", perfect_total);
    std::fprintf(logfp, "perfect_bias_%s\t%" PRIu64 "\n", cfg.trait0.c_str(), perfect_bias_trait0);
    std::fprintf(logfp, "perfect_bias_%s\t%" PRIu64 "\n", cfg.trait1.c_str(), perfect_bias_trait1);
    std::fprintf(logfp, "perf_wall_sec\t%.3f\n", wall);
    std::fprintf(logfp, "perf_max_rss_MB\t%.3f\n", rss_mb);

    std::fprintf(logfp, "\n=== perf ===\n");
    std::fprintf(logfp, "cpu_user_sec\t%.6f\n", ut);
    std::fprintf(logfp, "cpu_sys_sec\t%.6f\n", st);
    std::fprintf(logfp, "cpu_total_sec\t%.6f\n", cpu_total);
    std::fprintf(logfp, "wall_sec\t%.6f\n", wall);
    std::fprintf(logfp, "max_rss_MB\t%.3f\n", rss_mb);
    std::fprintf(logfp, "end_time_epoch\t%ld\n", std::time(nullptr));

    std::fprintf(logfp, "\n=== config ===\n");
    std::fprintf(logfp, "traits\t%s=%d\t%s=%d\n", cfg.trait0.c_str(), cfg.n0, cfg.trait1.c_str(), cfg.n1);

    std::fprintf(logfp, "\n=== global_prevalence ===\n");
    std::fprintf(logfp, "prev_min\t%.6f\n",    eff_prev_min);
    std::fprintf(logfp, "prev_max\t%.6f\n",    eff_prev_max);
    std::fprintf(logfp, "adaptive_prev_min\t%.6f\n", adaptive_min);
    std::fprintf(logfp, "adaptive_prev_max\t%.6f\n", adaptive_max);
    std::fprintf(logfp, "adaptive_prev_mode\t%s\n", cfg.adaptive_prev_mode.c_str());
    std::fprintf(logfp, "minority_ratio\t%.6f\n", minority_ratio);
    std::fprintf(logfp, "need_min_indiv\t%" PRIu64 "\n", thr.need_all_indiv);
    std::fprintf(logfp, "allow_max_indiv\t%" PRIu64 "\n", thr.allow_all_indiv);
    std::fprintf(logfp, "effective_min_frac\t%.6f\n", eff_total_min_frac);
    std::fprintf(logfp, "effective_max_frac\t%.6f\n", eff_total_max_frac);
    std::fprintf(logfp, "global_filters_applied\t%s\n", global_filters_applied);

    std::fprintf(logfp, "\n=== group_contrast ===\n");
    std::fprintf(logfp, "hi_group_min_frac\t%.6f\n", cfg.group_min);
    std::fprintf(logfp, "lo_group_max_frac\t%.6f\n", cfg.group_max);
    std::fprintf(logfp, "%s_hi_need_indiv\t%" PRIu64 "\n", cfg.trait0.c_str(), thr.need_min_A_ind);
    std::fprintf(logfp, "%s_hi_need_indiv\t%" PRIu64 "\n", cfg.trait1.c_str(), thr.need_min_B_ind);
    std::fprintf(logfp, "%s_lo_allow_indiv_when_%s_hi\t%" PRIu64 "\n",
                 cfg.trait1.c_str(), cfg.trait0.c_str(), thr.allow_max_B_ind);
    std::fprintf(logfp, "%s_lo_allow_indiv_when_%s_hi\t%" PRIu64 "\n",
                 cfg.trait0.c_str(), cfg.trait1.c_str(), thr.allow_max_A_ind);
    std::fprintf(logfp, "%s_hi_need_frac\t%.6f\n", cfg.trait0.c_str(), eff_g0_min_frac);
    std::fprintf(logfp, "%s_hi_need_frac\t%.6f\n", cfg.trait1.c_str(), eff_g1_min_frac);
    std::fprintf(logfp, "%s_lo_allow_frac_when_%s_hi\t%.6f\n",
                 cfg.trait1.c_str(), cfg.trait0.c_str(), eff_g1_max_frac);
    std::fprintf(logfp, "%s_lo_allow_frac_when_%s_hi\t%.6f\n",
                 cfg.trait0.c_str(), cfg.trait1.c_str(), eff_g0_max_frac);
    std::fprintf(logfp, "group_contrast_enabled\t%s\n", contrast_enabled ? "yes" : "no");
    std::fprintf(logfp, "group_contrast_applied\t%s\n",
                 (cfg.perfect_only ? "no (perfect-only)" : (contrast_enabled ? "yes" : "no")));

    const size_t homo_limit = homopolymer_allow_run(proc.k, comp.homo_max_frac, 15);

    std::fprintf(logfp, "\n=== filters ===\n");
    std::fprintf(logfp, "homopolymer\t%s\tmax_frac=%.3f\tcap=15\tmax_run_bases=%zu\n", comp.check_homo ? "on" : "off", comp.homo_max_frac, homo_limit);
    std::fprintf(logfp, "dinucleotide_repeat\t%s\n", comp.check_di ? "on" : "off");
    std::fprintf(logfp, "trinucleotide_repeat\t%s\n", comp.check_tri ? "on" : "off");
    std::fprintf(logfp, "entropy_min_bits\t%.6f\n", comp.entropy_min);

    std::fclose(logfp);

    {
      const std::string main_log_path = cfg.logging_dir + "/main.log";
      if (FILE* mfp = std::fopen(main_log_path.c_str(), "a")) {
        const std::time_t now = std::time(nullptr);
        char iso[32];
        std::strftime(iso, sizeof(iso), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
        char ts_local[32];
        if (std::tm* lt = std::localtime(&now)) {
          std::strftime(ts_local, sizeof(ts_local), "%Y-%m-%d %H:%M:%S", lt);
        } else {
          std::strncpy(ts_local, "1970-01-01 00:00:00", sizeof(ts_local));
          ts_local[sizeof(ts_local) - 1] = '\0';
        }

        auto log_line = [&](const char* fmt, auto... args) {
          std::fprintf(mfp, "[%s] [INFO] ", ts_local);
          std::fprintf(mfp, fmt, args...);
          std::fputc('\n', mfp);
        };

        log_line("[merge2stats] %s", iso);
        log_line("groups: %s (n=%d) vs %s (n=%d); k=%u",
                 cfg.trait0.c_str(), cfg.n0, cfg.trait1.c_str(), cfg.n1, A.k);
        log_line("global: min_indiv=%" PRIu64 " (%.4f), max_indiv=%" PRIu64 " (%.4f); applied=%s; minority_ratio=%.4f; adaptive_min=%.4f; adaptive_max=%.4f; adaptive_mode=%s",
                 thr.need_all_indiv, eff_total_min_frac, thr.allow_all_indiv, eff_total_max_frac,
                 global_filters_applied,
                 minority_ratio, adaptive_min, adaptive_max, cfg.adaptive_prev_mode.c_str());
        if (contrast_enabled) {
          if (cfg.perfect_only) {
            log_line("contrast: enabled=yes but bypassed due to perfect-only");
          } else {
            log_line("contrast: enabled=yes; %s >= %" PRIu64 " (%.4f) && %s <= %" PRIu64 " (%.4f), OR %s >= %" PRIu64 " (%.4f) && %s <= %" PRIu64 " (%.4f)",
                     cfg.trait0.c_str(), thr.need_min_A_ind, eff_g0_min_frac,
                     cfg.trait1.c_str(), thr.allow_max_B_ind, eff_g1_max_frac,
                     cfg.trait1.c_str(), thr.need_min_B_ind, eff_g1_min_frac,
                     cfg.trait0.c_str(), thr.allow_max_A_ind, eff_g0_max_frac);
          }
        } else {
        log_line("contrast: enabled=no");
        }
        log_line("lc: homo=%s(max_frac=%.3f,cap=15,max_run_bases=%zu), di=%s, tri=%s, entropy_min=%.3f",
                 comp.check_homo ? "on" : "off", comp.homo_max_frac,
                 homo_limit,
                 comp.check_di ? "on" : "off",
                 comp.check_tri ? "on" : "off",
                 comp.entropy_min);
        log_line("input_overlap: %s=%" PRIu64 ", %s=%" PRIu64 ", common=%" PRIu64 ", only_%s=%" PRIu64 ", only_%s=%" PRIu64 ", union=%" PRIu64,
                 cfg.trait0.c_str(), (uint64_t)A.info.total_kmers,
                 cfg.trait1.c_str(), (uint64_t)B.info.total_kmers,
                 (uint64_t)proc.common_kmers,
                 cfg.trait0.c_str(), (uint64_t)proc.only0_kmers,
                 cfg.trait1.c_str(), (uint64_t)proc.only1_kmers,
                 (uint64_t)union_total);
        log_line("lc_results: tested=%" PRIu64 ", pass=%" PRIu64 ", reject_homo=%" PRIu64 ", reject_di=%" PRIu64 ", reject_tri=%" PRIu64 ", reject_entropy=%" PRIu64,
                 (uint64_t)proc.lc_tested, (uint64_t)proc.lc_pass,
                 (uint64_t)proc.lc_reject_homo, (uint64_t)proc.lc_reject_di, (uint64_t)proc.lc_reject_tri, (uint64_t)proc.lc_reject_entropy);
        log_line("results: emitted=%" PRIu64 "; bias: %s=%" PRIu64 ", %s=%" PRIu64 ", ties=%" PRIu64,
                 (uint64_t)proc.emitted,
                 cfg.trait0.c_str(), (uint64_t)proc.post_bias_0,
                 cfg.trait1.c_str(), (uint64_t)proc.post_bias_1,
                 (uint64_t)proc.post_ties);
        if (cfg.perfect_only || cfg.perfect_extra) {
          log_line("perfect_contrast: mode=%s; detected=%s; total=%" PRIu64 "; %s_only=%" PRIu64 "; %s_only=%" PRIu64,
                   perfect_mode_str,
                   perfect_total ? "yes" : "no",
                   perfect_total,
                   cfg.trait0.c_str(), perfect_bias_trait0,
                   cfg.trait1.c_str(), perfect_bias_trait1);
          log_line("perfect_contrast_fasta: path=%s; emitted=%" PRIu64,
                   perfect_fasta_path.c_str(),
                   (uint64_t)proc.perfect_emitted);
        }
        if (cfg.perfect_only) {
          log_line("files: fasta=<disabled via perfect-only>; stats=%s", disp_stats.c_str());
        } else {
          log_line("files: fasta=%s; stats=%s", disp_fasta.c_str(), disp_stats.c_str());
        }
        log_line("perf: wall=%.3fs; max_rss=%.1fMB", wall, rss_mb);
        std::fputc('\n', mfp);
        std::fclose(mfp);
      }
    }

    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
