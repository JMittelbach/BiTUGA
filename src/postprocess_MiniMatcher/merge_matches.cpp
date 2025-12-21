#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Args {
  std::string indir;
  std::string metadata;
  std::string out;
  std::string stats;
  std::string unitig_stats;
};

static void die(const std::string &msg) {
  std::cerr << "ERROR: " << msg << "\n";
  std::exit(1);
}

static bool file_exists(const std::string &path) {
  std::ifstream in(path);
  return in.good();
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

static void parse_args(int argc, char **argv, Args &args) {
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--indir" && i + 1 < argc) {
      args.indir = argv[++i];
    } else if (a == "--metadata" && i + 1 < argc) {
      args.metadata = argv[++i];
    } else if (a == "--out" && i + 1 < argc) {
      args.out = argv[++i];
    } else if (a == "--stats" && i + 1 < argc) {
      args.stats = argv[++i];
    } else if (a == "--unitig-stats" && i + 1 < argc) {
      args.unitig_stats = argv[++i];
    } else {
      die("Unknown or incomplete argument: " + a);
    }
  }
  if (args.indir.empty() || args.metadata.empty() || args.out.empty() ||
      args.stats.empty()) {
    die("Required arguments: --indir --metadata --out --stats [--unitig-stats]");
  }
}

static int read_total_unitigs(const std::string &unitig_stats_path) {
  if (unitig_stats_path.empty() || !file_exists(unitig_stats_path)) {
    return -1;
  }
  std::ifstream in(unitig_stats_path);
  if (!in) {
    return -1;
  }
  std::string line;
  int total = -1;
  int after_minlen = -1;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (line.rfind("n_unitigs_after_minlen_", 0) == 0) {
      std::size_t tabpos = line.find('\t');
      if (tabpos != std::string::npos && tabpos + 1 < line.size()) {
        std::string v = line.substr(tabpos + 1);
        try {
          after_minlen = std::stoi(v);
          break;
        } catch (...) {
        }
      }
    } else if (line.rfind("n_unitigs", 0) == 0) {
      std::size_t tabpos = line.find('\t');
      if (tabpos != std::string::npos && tabpos + 1 < line.size()) {
        std::string v = line.substr(tabpos + 1);
        try {
          total = std::stoi(v);
        } catch (...) {
          total = -1;
        }
      }
    }
  }
  if (after_minlen >= 0) {
    return after_minlen;
  }
  return total;
}

struct Metadata {
  std::unordered_map<std::string, std::string> sample_to_trait;
  std::vector<std::string> traits;
  std::vector<int> trait_sample_counts;
};

static Metadata load_metadata(const std::string &meta_path) {
  Metadata md;
  std::ifstream in(meta_path);
  if (!in) {
    die("Cannot open metadata: " + meta_path);
  }
  std::string header;
  if (!std::getline(in, header)) {
    die("Empty metadata file: " + meta_path);
  }
  std::vector<std::string> cols = split_tab(header);
  int idx_sample = -1;
  int idx_trait = -1;
  for (std::size_t i = 0; i < cols.size(); ++i) {
    if (cols[i] == "sample_id") {
      idx_sample = static_cast<int>(i);
    } else if (cols[i] == "trait") {
      idx_trait = static_cast<int>(i);
    }
  }
  if (idx_sample < 0 || idx_trait < 0) {
    die("Metadata must contain columns 'sample_id' and 'trait'");
  }
  std::map<std::string, std::set<std::string>> trait_to_samples_set;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::vector<std::string> parts = split_tab(line);
    if (static_cast<int>(parts.size()) <= std::max(idx_sample, idx_trait)) {
      continue;
    }
    std::string sid = parts[idx_sample];
    std::string tr = parts[idx_trait];
    if (sid.empty() || tr.empty()) {
      continue;
    }
    md.sample_to_trait[sid] = tr;
    trait_to_samples_set[tr].insert(sid);
  }
  for (auto &kv : trait_to_samples_set) {
    md.traits.push_back(kv.first);
  }
  std::sort(md.traits.begin(), md.traits.end());
  md.trait_sample_counts.resize(md.traits.size(), 0);
  for (std::size_t i = 0; i < md.traits.size(); ++i) {
    const std::string &t = md.traits[i];
    md.trait_sample_counts[i] = static_cast<int>(trait_to_samples_set[t].size());
  }
  return md;
}

struct Presence {
  std::unordered_map<std::string, std::vector<int>> counts;
};

static std::vector<std::string> read_present_bin(const std::string &path) {
  std::vector<std::string> ids;
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return ids;
  }
  uint32_t count = 0;
  in.read(reinterpret_cast<char *>(&count), sizeof(uint32_t));
  if (!in) {
    return ids;
  }
  ids.reserve(count);
  for (uint32_t i = 0; i < count; ++i) {
    uint32_t len = 0;
    in.read(reinterpret_cast<char *>(&len), sizeof(uint32_t));
    if (!in) {
      ids.clear();
      return ids;
    }
    std::string uid(len, '\0');
    in.read(uid.data(), len);
    if (!in) {
      ids.clear();
      return ids;
    }
    ids.push_back(std::move(uid));
  }
  return ids;
}

static void build_presence(const std::string &indir, const Metadata &md,
                           Presence &presence) {
  std::unordered_map<std::string, int> trait_index;
  for (std::size_t i = 0; i < md.traits.size(); ++i) {
    trait_index[md.traits[i]] = static_cast<int>(i);
  }

  for (const auto &kv : md.sample_to_trait) {
    const std::string &sid = kv.first;
    const std::string &trait = kv.second;
    auto it_idx = trait_index.find(trait);
    if (it_idx == trait_index.end()) {
      continue;
    }
    int t_idx = it_idx->second;

    std::string present_file = indir;
    if (!present_file.empty() && present_file.back() != '/') {
      present_file.push_back('/');
    }
    present_file += sid;
    present_file += ".present_unitigs.bin";

    if (!file_exists(present_file)) {
      continue;
    }

    std::vector<std::string> ids = read_present_bin(present_file);
    if (ids.empty()) {
      continue;
    }

    std::unordered_set<std::string> seen;
    for (const auto &uid : ids) {
      if (uid.empty()) {
        continue;
      }
      if (!seen.insert(uid).second) {
        continue;
      }
      auto it = presence.counts.find(uid);
      if (it == presence.counts.end()) {
        std::vector<int> vec(md.traits.size(), 0);
        vec[t_idx] = 1;
        presence.counts.emplace(uid, std::move(vec));
      } else {
        it->second[t_idx] += 1;
      }
    }
  }
}

struct KeySortEntry {
  std::string key;
  bool is_int;
  long long val;
};

static std::vector<KeySortEntry> make_sorted_keys(
    const std::unordered_map<std::string, std::vector<int>> &m) {
  std::vector<KeySortEntry> keys;
  keys.reserve(m.size());
  for (const auto &kv : m) {
    const std::string &k = kv.first;
    KeySortEntry e;
    e.key = k;
    e.is_int = false;
    e.val = 0;
    bool ok = true;
    if (!k.empty()) {
      std::size_t pos = 0;
      if (k[0] == '+' || k[0] == '-') {
        pos = 1;
      }
      for (; pos < k.size(); ++pos) {
        if (k[pos] < '0' || k[pos] > '9') {
          ok = false;
          break;
        }
      }
    } else {
      ok = false;
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
  Args args;
  parse_args(argc, argv, args);

  Metadata md = load_metadata(args.metadata);
  if (md.traits.empty()) {
    die("No traits found in metadata");
  }

  Presence presence;
  build_presence(args.indir, md, presence);

  int total_unitigs = read_total_unitigs(args.unitig_stats);
  int present_unitigs = static_cast<int>(presence.counts.size());
  int absent_unitigs = -1;
  if (total_unitigs >= 0) {
    absent_unitigs = total_unitigs - present_unitigs;
  }

  {
    std::ofstream out(args.out);
    if (!out) {
      die("Cannot open output: " + args.out);
    }
    out << "unitig";
    for (const auto &t : md.traits) {
      out << "\t" << t;
    }
    out << "\n";

    std::vector<KeySortEntry> keys = make_sorted_keys(presence.counts);
    for (const auto &e : keys) {
      const std::string &u = e.key;
      const std::vector<int> &vec = presence.counts[u];
      out << u;
      for (int c : vec) {
        out << "\t" << c;
      }
      out << "\n";
    }
  }

  std::map<int, long long> total_presence_hist;
  std::map<int, long long> abs_diff_hist;
  std::map<std::string, long long> trait_bias_specific_only;
  std::map<std::string, long long> trait_bias_biased;
  std::map<std::string, long long> full_prevalence;
  std::map<std::string, long long> trait_specific_only;
  std::map<std::string, long long> full_prevalence_absent_others;

  for (const auto &t : md.traits) {
    trait_bias_specific_only[t] = 0;
    trait_bias_biased[t] = 0;
    full_prevalence[t] = 0;
    trait_specific_only[t] = 0;
    full_prevalence_absent_others[t] = 0;
  }

  for (const auto &kv : presence.counts) {
    const std::vector<int> &vec = kv.second;
    int total_presence = 0;
    for (int c : vec) {
      total_presence += c;
    }
    if (total_presence > 0) {
      total_presence_hist[total_presence] += 1;
    }

    int n_traits = static_cast<int>(md.traits.size());
    if (n_traits == 2) {
      int c1 = vec[0];
      int c2 = vec[1];
      int d = c1 - c2;
      if (d < 0) d = -d;
      abs_diff_hist[d] += 1;
    } else {
      int max_c = 0;
      std::vector<int> vals = vec;
      for (int c : vals) {
        if (c > max_c) max_c = c;
      }
      for (int i = 0; i < n_traits; ++i) {
        for (int j = i + 1; j < n_traits; ++j) {
          if (vals[i] != vals[j]) {
            int d = vals[i] - vals[j];
            if (d < 0) d = -d;
            abs_diff_hist[d] += 1;
          }
        }
      }
      (void)max_c;
    }

    for (int ti = 0; ti < n_traits; ++ti) {
      const std::string &t = md.traits[ti];
      int c_t = vec[ti];
      std::vector<int> other_counts;
      other_counts.reserve(n_traits - 1);
      for (int tj = 0; tj < n_traits; ++tj) {
        if (tj == ti) continue;
        other_counts.push_back(vec[tj]);
      }
      bool all_others_zero = true;
      for (int oc : other_counts) {
        if (oc != 0) {
          all_others_zero = false;
          break;
        }
      }
      if (c_t > 0 && all_others_zero) {
        trait_bias_specific_only[t] += 1;
        trait_specific_only[t] += 1;
      }

      bool strictly_larger = true;
      for (int tj = 0; tj < n_traits; ++tj) {
        if (tj == ti) continue;
        if (!(c_t > vec[tj])) {
          strictly_larger = false;
          break;
        }
      }
      if (c_t > 0 && strictly_larger) {
        trait_bias_biased[t] += 1;
      }
    }

    for (int ti = 0; ti < n_traits; ++ti) {
      const std::string &t = md.traits[ti];
      int c_t = vec[ti];
      int n_t = md.trait_sample_counts[ti];
      if (n_t > 0 && c_t == n_t) {
        full_prevalence[t] += 1;
        bool others_zero = true;
        for (int tj = 0; tj < n_traits; ++tj) {
          if (tj == ti) continue;
          if (vec[tj] != 0) {
            others_zero = false;
            break;
          }
        }
        if (others_zero) {
          full_prevalence_absent_others[t] += 1;
        }
      }
    }
  }

  std::ofstream sf(args.stats);
  if (!sf) {
    die("Cannot open stats output: " + args.stats);
  }

  if (total_unitigs >= 0) {
    sf << "total_unitigs\t" << total_unitigs << "\n";
  }
  sf << "present_unitigs\t" << present_unitigs << "\n";
  if (absent_unitigs >= 0) {
    sf << "absent_unitigs\t" << absent_unitigs << "\n";
  }
  sf << "\n";

  sf << "presence_total_hist\n";
  sf << "total_presence\tn_unitigs\n";
  for (const auto &kv2 : total_presence_hist) {
    sf << kv2.first << "\t" << kv2.second << "\n";
  }
  sf << "\n";

  sf << "abs_diff_hist\n";
  sf << "abs_diff\tn_unitigs\n";
  for (const auto &kv2 : abs_diff_hist) {
    sf << kv2.first << "\t" << kv2.second << "\n";
  }
  sf << "\n";

  for (const auto &t : md.traits) {
    long long spec_only = trait_bias_specific_only[t];
    long long biased = trait_bias_biased[t];
    sf << "TRAIT_BIAS\t" << t << "\tspecific_only=" << spec_only
       << "\tbiased_gt_other=" << biased << "\n";
  }
  sf << "\n";

  for (const auto &t : md.traits) {
    sf << "FULL_PREVALENCE\t" << t << "\t" << full_prevalence[t] << "\n";
  }
  sf << "\n";

  for (const auto &t : md.traits) {
    sf << "TRAIT_SPECIFIC_ONLY\t" << t << "\t" << trait_specific_only[t]
       << "\n";
  }
  sf << "\n";

  for (const auto &t : md.traits) {
    sf << "FULL_PREVALENCE_ABSENT_OTHERS\t" << t << "\t"
       << full_prevalence_absent_others[t] << "\n";
  }

  return 0;
}
