#include <algorithm>
#include <cstdint>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

struct Options {
  std::string unitigs_fa;
  std::string matches_tsv;
  std::string present_out;

  int carrier_threshold = 3;

  int short_max = 37;
  int medium_max = 92;

  double short_frac = 1.0;
  double medium_frac = 0.98;
  double long_frac = 0.95;

  int anchor_bp = 12;

  int short_min_bp = 37;
  int medium_min_bp = 92;
  int long_min_bp = 250;

  int bridge_min_span = 180;

  int short_min_segment = 31;
  int medium_min_segment = 50;
  int long_min_segment = 80;

  int min_bridging_pairs = 2;

  bool full_gapless_long = false;
};

static void die(const std::string &msg) {
  std::cerr << "ERROR: " << msg << "\n";
  std::exit(1);
}

static std::string trim(const std::string &s) {
  std::size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) {
    ++i;
  }
  if (i == s.size()) return std::string();
  std::size_t j = s.size();
  while (j > i &&
         std::isspace(static_cast<unsigned char>(s[j - 1]))) {
    --j;
  }
  return s.substr(i, j - i);
}

static std::vector<std::string> split_tab(const std::string &line) {
  std::vector<std::string> res;
  std::string field;
  for (char c : line) {
    if (c == '\t') {
      res.push_back(field);
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  res.push_back(field);
  return res;
}

static void parse_args(int argc, char **argv, Options &opt) {
  if (argc < 4) {
    die("Usage: postprocess_minimatcher unitigs_fa matches_tsv present_out "
        "[--carrier-threshold INT] [--short-max INT] [--medium-max INT] "
        "[--short-frac FLOAT] [--medium-frac FLOAT] [--long-frac FLOAT] "
        "[--anchor-bp INT] [--short-min-bp INT] [--medium-min-bp INT] "
        "[--long-min-bp INT] [--bridge-min-span INT] "
        "[--short-min-segment INT] [--medium-min-segment INT] "
        "[--long-min-segment INT] [--min-bridging-pairs INT] "
        "[--full-gapless-long]");
  }

  opt.unitigs_fa = argv[1];
  opt.matches_tsv = argv[2];
  opt.present_out = argv[3];

  int i = 4;
  while (i < argc) {
    std::string a = argv[i];
    if (a == "--carrier-threshold" && i + 1 < argc) {
      opt.carrier_threshold = std::stoi(argv[++i]);
    } else if (a == "--short-max" && i + 1 < argc) {
      opt.short_max = std::stoi(argv[++i]);
    } else if (a == "--medium-max" && i + 1 < argc) {
      opt.medium_max = std::stoi(argv[++i]);
    } else if (a == "--short-frac" && i + 1 < argc) {
      opt.short_frac = std::stod(argv[++i]);
    } else if (a == "--medium-frac" && i + 1 < argc) {
      opt.medium_frac = std::stod(argv[++i]);
    } else if (a == "--long-frac" && i + 1 < argc) {
      opt.long_frac = std::stod(argv[++i]);
    } else if (a == "--anchor-bp" && i + 1 < argc) {
      opt.anchor_bp = std::stoi(argv[++i]);
    } else if (a == "--short-min-bp" && i + 1 < argc) {
      opt.short_min_bp = std::stoi(argv[++i]);
    } else if (a == "--medium-min-bp" && i + 1 < argc) {
      opt.medium_min_bp = std::stoi(argv[++i]);
    } else if (a == "--long-min-bp" && i + 1 < argc) {
      opt.long_min_bp = std::stoi(argv[++i]);
    } else if (a == "--bridge-min-span" && i + 1 < argc) {
      opt.bridge_min_span = std::stoi(argv[++i]);
    } else if (a == "--short-min-segment" && i + 1 < argc) {
      opt.short_min_segment = std::stoi(argv[++i]);
    } else if (a == "--medium-min-segment" && i + 1 < argc) {
      opt.medium_min_segment = std::stoi(argv[++i]);
    } else if (a == "--long-min-segment" && i + 1 < argc) {
      opt.long_min_segment = std::stoi(argv[++i]);
    } else if (a == "--min-bridging-pairs" && i + 1 < argc) {
      opt.min_bridging_pairs = std::stoi(argv[++i]);
    } else if (a == "--full-gapless-long") {
      opt.full_gapless_long = true;
    } else {
      die("Unknown or incomplete option: " + a);
    }
    ++i;
  }
}

static std::unordered_map<std::string, int>
load_unitig_lengths(const std::string &fa_path) {
  std::unordered_map<std::string, int> lens;
  std::ifstream in(fa_path);
  if (!in) {
    die("Cannot open unitigs FASTA: " + fa_path);
  }
  std::string line;
  std::string current_id;
  int current_len = 0;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (!current_id.empty()) {
        lens[current_id] = current_len;
      }
      std::string hdr = line.substr(1);
      std::string id;
      for (char c : hdr) {
        if (std::isspace(static_cast<unsigned char>(c))) break;
        id.push_back(c);
      }
      current_id = id;
      current_len = 0;
    } else {
      for (char c : line) {
        if (!std::isspace(static_cast<unsigned char>(c))) {
          ++current_len;
        }
      }
    }
  }
  if (!current_id.empty()) {
    lens[current_id] = current_len;
  }
  return lens;
}

struct CarrierInfo {
  std::unordered_set<int> readinsts;
  int min_start = std::numeric_limits<int>::max();
  int max_end = std::numeric_limits<int>::min();

  void update(int start, int end, int readinst) {
    if (readinst >= 0) {
      readinsts.insert(readinst);
    }
    if (start < min_start) {
      min_start = start;
    }
    if (end > max_end) {
      max_end = end;
    }
  }
};

struct UnitigSupport {
  std::unordered_set<std::string> carriers;
  std::vector<std::pair<int, int>> intervals;
  int max_read_len = 0;
  std::unordered_map<std::string, CarrierInfo> carrier_infos;
};

static std::vector<std::pair<int, int>>
merge_intervals(const std::vector<std::pair<int, int>> &intervals,
                int unitig_len) {
  if (intervals.empty()) {
    return {};
  }
  std::vector<std::pair<int, int>> clipped;
  clipped.reserve(intervals.size());
  for (auto &iv : intervals) {
    int s = iv.first;
    int e = iv.second;
    if (s < 0) s = 0;
    if (s > unitig_len) s = unitig_len;
    if (e < s) e = s;
    if (e > unitig_len) e = unitig_len;
    if (s == e) continue;
    clipped.emplace_back(s, e);
  }
  if (clipped.empty()) {
    return {};
  }
  std::sort(clipped.begin(), clipped.end(),
            [](const std::pair<int, int> &a,
               const std::pair<int, int> &b) {
              if (a.first < b.first) return true;
              if (a.first > b.first) return false;
              return a.second < b.second;
            });
  std::vector<std::pair<int, int>> merged;
  merged.reserve(clipped.size());
  merged.push_back(clipped[0]);
  for (std::size_t i = 1; i < clipped.size(); ++i) {
    int s = clipped[i].first;
    int e = clipped[i].second;
    auto &last = merged.back();
    if (s <= last.second) {
      if (e > last.second) {
        last.second = e;
      }
    } else {
      merged.emplace_back(s, e);
    }
  }
  return merged;
}

static int total_covered_bases(
    const std::vector<std::pair<int, int>> &merged) {
  long long total = 0;
  for (auto &iv : merged) {
    total += static_cast<long long>(iv.second - iv.first);
  }
  if (total > std::numeric_limits<int>::max()) {
    return std::numeric_limits<int>::max();
  }
  return static_cast<int>(total);
}

static double coverage_fraction(
    const std::vector<std::pair<int, int>> &merged, int unitig_len) {
  if (unitig_len <= 0 || merged.empty()) return 0.0;
  int cov = total_covered_bases(merged);
  return static_cast<double>(cov) /
         static_cast<double>(unitig_len);
}

static bool interval_overlaps_window(
    const std::vector<std::pair<int, int>> &intervals,
    int w_start, int w_end) {
  for (auto &iv : intervals) {
    int s = iv.first;
    int e = iv.second;
    if (s >= w_end) {
      break;
    }
    if (e > w_start && s < w_end) {
      return true;
    }
  }
  return false;
}

static bool has_anchor_support(
    const std::vector<std::pair<int, int>> &merged,
    int unitig_len,
    int anchor_bp) {
  if (merged.empty() || unitig_len <= 0) return false;
  int half = (unitig_len > 1) ? (unitig_len / 2) : 1;
  int span = anchor_bp;
  if (span < 1) span = 1;
  if (span > half) span = half;
  int start0 = 0;
  int end0 = span;
  int end1 = unitig_len;
  int start1 = unitig_len - span;
  if (start1 < 0) start1 = 0;
  bool ok_start = interval_overlaps_window(merged, start0, end0);
  bool ok_end = interval_overlaps_window(merged, start1, end1);
  return ok_start && ok_end;
}

static int max_gap_size(const std::vector<std::pair<int, int>> &merged,
                        int unitig_len) {
  if (merged.empty()) {
    return unitig_len;
  }
  int prev_end = 0;
  int max_gap = merged[0].first; 
  for (auto &iv : merged) {
    int s = iv.first;
    int e = iv.second;
    if (s > prev_end) {
      int gap = s - prev_end;
      if (gap > max_gap) {
        max_gap = gap;
      }
    }
    prev_end = e;
  }
  int tail_gap = unitig_len - prev_end;
  if (tail_gap > max_gap) {
    max_gap = tail_gap;
  }
  if (max_gap < 0) max_gap = 0;
  return max_gap;
}

static int allowed_gap(int unitig_len, double frac_threshold) {
  double slack_d =
      static_cast<double>(unitig_len) *
      std::max(0.0, 1.0 - frac_threshold);
  int slack = static_cast<int>(std::round(slack_d));
  if (slack < 0) slack = 0;
  return slack;
}

static double fraction_threshold(int unitig_len,
                                 const Options &opt) {
  if (unitig_len <= opt.short_max) {
    return opt.short_frac;
  }
  if (unitig_len <= opt.medium_max) {
    return opt.medium_frac;
  }
  return opt.long_frac;
}

static int min_bp_threshold(int unitig_len, const Options &opt) {
  if (unitig_len <= opt.short_max) {
    return std::min(unitig_len, opt.short_min_bp);
  }
  if (unitig_len <= opt.medium_max) {
    return std::min(unitig_len, opt.medium_min_bp);
  }
  return std::min(unitig_len, opt.long_min_bp);
}

static int min_segment_threshold(int unitig_len,
                                 const Options &opt) {
  if (unitig_len <= opt.short_max) {
    return std::min(unitig_len, opt.short_min_segment);
  }
  if (unitig_len <= opt.medium_max) {
    return std::min(unitig_len, opt.medium_min_segment);
  }
  return std::min(unitig_len, opt.long_min_segment);
}

static int parse_readinst(const std::string &field) {
  if (field.empty()) return -1;
  std::string digits;
  digits.reserve(field.size());
  for (char c : field) {
    if (std::isdigit(static_cast<unsigned char>(c))) {
      digits.push_back(c);
    } else {
      break;
    }
  }
  if (digits.empty()) return -1;
  try {
    return std::stoi(digits);
  } catch (...) {
    return -1;
  }
}

static int count_bridging_pairs(const UnitigSupport &data,
                                int unitig_len,
                                int min_span) {
  if (data.carrier_infos.empty()) return 0;
  int target_span = min_span;
  if (target_span < 0) target_span = 0;
  if (target_span > unitig_len) target_span = unitig_len;

  int count = 0;
  const int INT_MAXI = std::numeric_limits<int>::max();
  const int INT_MINI = std::numeric_limits<int>::min();
  for (const auto &kv : data.carrier_infos) {
    const CarrierInfo &info = kv.second;
    if (info.readinsts.size() < 2) {
      continue;
    }
    if (info.min_start == INT_MAXI ||
        info.max_end == INT_MINI) {
      continue;
    }
    int span = info.max_end - info.min_start;
    if (span >= target_span) {
      ++count;
    }
  }
  return count;
}

static bool has_paired_information(const UnitigSupport &data) {
  for (const auto &kv : data.carrier_infos) {
    const CarrierInfo &info = kv.second;
    if (info.readinsts.size() >= 2) {
      return true;
    }
  }
  return false;
}

static void read_matches_stream(std::istream &in,
    std::unordered_map<std::string, UnitigSupport> &supports) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    auto parts = split_tab(line);
    if (parts.size() < 3) continue;

    std::string uid = trim(parts[0]);
    if (uid.empty()) continue;

    int start = 0;
    int length = 0;
    try {
      start = std::stoi(parts[1]);
      length = std::stoi(parts[2]);
    } catch (...) {
      continue;
    }
    if (length <= 0) continue;

    std::string carrier_id;
    if (parts.size() > 3) {
      carrier_id = trim(parts[3]);
    }
    std::string readinst_field;
    if (parts.size() > 4) {
      readinst_field = trim(parts[4]);
    }

    int read_match_len = length;
    if (parts.size() > 6) {
      try {
        read_match_len = std::stoi(parts[6]);
      } catch (...) {
        read_match_len = length;
      }
    }
    int read_seq_len = read_match_len;
    if (parts.size() > 7) {
      try {
        read_seq_len = std::stoi(parts[7]);
      } catch (...) {
        read_seq_len = read_match_len;
      }
    }

    UnitigSupport &sup = supports[uid];
    int end = start + length;

    if (!carrier_id.empty()) {
      sup.carriers.insert(carrier_id);
      CarrierInfo &ci = sup.carrier_infos[carrier_id];
      int ri = parse_readinst(readinst_field);
      ci.update(start, end, ri);
    }

    sup.intervals.emplace_back(start, end);

    int max_len_candidate = read_match_len;
    if (read_seq_len > max_len_candidate) {
      max_len_candidate = read_seq_len;
    }
    if (max_len_candidate > sup.max_read_len) {
      sup.max_read_len = max_len_candidate;
    }
  }
}

struct KeySortEntry {
  std::string key;
  bool is_int = false;
  long long val = 0;
};

static std::vector<KeySortEntry> make_sorted_keys(
    const std::unordered_map<std::string, UnitigSupport> &m) {
  std::vector<KeySortEntry> keys;
  keys.reserve(m.size());
  for (const auto &kv : m) {
    const std::string &k = kv.first;
    KeySortEntry e;
    e.key = k;
    bool ok = true;
    if (k.empty()) {
      ok = false;
    } else {
      std::size_t pos = 0;
      if (k[0] == '+' || k[0] == '-') {
        pos = 1;
      }
      for (; pos < k.size(); ++pos) {
        if (!std::isdigit(static_cast<unsigned char>(k[pos]))) {
          ok = false;
          break;
        }
      }
    }
    if (ok) {
      try {
        e.val = std::stoll(k);
        e.is_int = true;
      } catch (...) {
        e.is_int = false;
        e.val = 0;
      }
    }
    keys.push_back(e);
  }
  std::sort(keys.begin(), keys.end(),
            [](const KeySortEntry &a, const KeySortEntry &b) {
              if (a.is_int && b.is_int) {
                if (a.val < b.val) return true;
                if (a.val > b.val) return false;
                return a.key < b.key;
              }
              if (a.is_int != b.is_int) {
                return a.is_int;
              }
              return a.key < b.key;
            });
  return keys;
}

int main(int argc, char **argv) {
  Options opt;
  parse_args(argc, argv, opt);

  auto unitig_lengths = load_unitig_lengths(opt.unitigs_fa);

  std::unordered_map<std::string, UnitigSupport> supports;
  supports.reserve(unitig_lengths.size());

  if (opt.matches_tsv == "-") {
    read_matches_stream(std::cin, supports);
  } else {
    std::ifstream in(opt.matches_tsv);
    if (!in) {
      die("Cannot open matches TSV: " + opt.matches_tsv);
    }
    read_matches_stream(in, supports);
  }

  std::vector<KeySortEntry> keys = make_sorted_keys(supports);
  std::vector<std::string> passing;
  passing.reserve(keys.size());

  for (const auto &e : keys) {
    const std::string &uid = e.key;
    auto it_len = unitig_lengths.find(uid);
    if (it_len == unitig_lengths.end()) {
      continue;
    }
    int unitig_len = it_len->second;
    if (unitig_len <= 0) continue;

    const UnitigSupport &data = supports[uid];

    if (static_cast<int>(data.carriers.size()) <
        opt.carrier_threshold) {
      continue;
    }

    int min_seg = min_segment_threshold(unitig_len, opt);
    std::vector<std::pair<int, int>> usable;
    usable.reserve(data.intervals.size());
    for (auto &iv : data.intervals) {
      int seg_len = iv.second - iv.first;
      if (seg_len >= min_seg) {
        usable.push_back(iv);
      }
    }
    if (usable.empty()) continue;

    auto merged = merge_intervals(usable, unitig_len);
    if (merged.empty()) continue;

    double frac_thr = fraction_threshold(unitig_len, opt);
    double cov_frac = coverage_fraction(merged, unitig_len);
    if (cov_frac < frac_thr) continue;

    int covered_bp = total_covered_bases(merged);
    int min_bp = min_bp_threshold(unitig_len, opt);
    if (covered_bp < min_bp) continue;

    if (!has_anchor_support(merged, unitig_len, opt.anchor_bp)) {
      continue;
    }

    int gap = max_gap_size(merged, unitig_len);
    int gap_allow = 0;
    if (opt.full_gapless_long && unitig_len > opt.medium_max) {
      gap_allow = 0;
    } else {
      gap_allow = allowed_gap(unitig_len, frac_thr);
    }
    if (gap > gap_allow) continue;

    bool needs_bridge = (unitig_len > opt.medium_max);
    if (needs_bridge) {
      if (has_paired_information(data)) {
        int bridge_pairs = count_bridging_pairs(
            data, unitig_len, opt.bridge_min_span);
        int minbp = std::max(1, opt.min_bridging_pairs);
        if (bridge_pairs < minbp) {
          continue;
        }
      } else {
        if (cov_frac < 1.0) {
          continue;
        }
      }
    }

    passing.push_back(uid);
  }

  std::ofstream out(opt.present_out, std::ios::binary);
  if (!out) {
    die("Cannot open present_out: " + opt.present_out);
  }
  uint32_t count = static_cast<uint32_t>(passing.size());
  out.write(reinterpret_cast<const char *>(&count), sizeof(uint32_t));
  for (const auto &uid : passing) {
    uint32_t len = static_cast<uint32_t>(uid.size());
    out.write(reinterpret_cast<const char *>(&len), sizeof(uint32_t));
    out.write(uid.data(), len);
  }

  return 0;
}
