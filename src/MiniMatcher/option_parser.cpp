#include <string>
#include <vector>
#include "utilities/cxxopts.hpp"
#include "threading/threads_output_files.hpp"
#include "option_parser.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << '\n';
  std::cerr << "CXX_VERSION\t" << CXX_VERSION << '\n';
  std::cerr << "CXX_FLAGS\t" << CXX_FLAGS << '\n';
}

// Option parser based on https://github.com/jarro2783/cxxopts

void Options::parse(int argc, char **argv)
{
  ReverseComplementFor default_reverse_complement_for{};

  std::string msg_reverse_complement_for
    = std::string("specify from which sequences the reverse complement is "
                  "taken to compute reverse complement matches; if not "
                  "specified, then only matches on forward strand are "
                  "computed; possible values are ")
      + default_reverse_complement_for.string_values_joined(", ");

  std::string thread_out_prefix_help;

  const char *this_help_line = ThreadsOutputFiles::help_line;
  int spaces = 0;
  for (size_t idx = 0; this_help_line[idx] != '\0'; idx++)
  {
    const char cc = this_help_line[idx];
    if (isspace(cc))
    {
      spaces++;
    } else
    {
      if (spaces > 0)
      {
        thread_out_prefix_help += ' ';
        spaces = 0;
      }
      thread_out_prefix_help += cc;
    }
  }

  cxxopts::Options options(argv[0],
                           gi_joiner_code
                             ? ("run join methods on minimizers generated "
                                "from sequences of database and queryfile")
                             : ("compute maximal matches and alignments "
                                "based on seeds from minimizers"));

  options.set_width(80);
  options.custom_help(std::string("[options] reference-file [query-file1 ..]"));
  options.set_tab_expansion();
  options.add_options()
     ("k,kmer_length", "k-mer length",
      cxxopts::value<int>(kmer_size)->default_value(default_kmer_size))

     ("b,hash_bits", "number of bits used for hashing, cannot be combined with "
                     "--invint",
      cxxopts::value<int>(hash_bits)->default_value(default_hash_bits))

     ("t,threads", "number of threads",
      cxxopts::value<int>(number_of_threads)->
                          default_value(default_number_of_threads))

     ("o,threads_out_prefix",
      thread_out_prefix_help,
      cxxopts::value<std::string>(threads_out_prefix)->
                          default_value(default_threads_out_prefix))

     ("s,thread_pool_max_size",
      "number of work packages stored in the thread pool; if not specified, "
      "the number will equal the number of threads",
      cxxopts::value<int>(thread_pool_max_size)
                          ->default_value(default_thread_pool_max_size))

     ("q,query_split_size", "number of query sequences (from a single file) "
                            "which are processed per thread",
      cxxopts::value<int>(query_split_size)->
                          default_value(default_query_split_size))

     ("max_replicate_ref", "maximum number of replicates of hash value on "
                           "reference sequence to produce a seed; if not set, "
                           "the number of replicates is unbounded",
      cxxopts::value<size_t>(max_replicate_ref)->
                             default_value(default_max_replicate_ref))

     ("max_replicate_qry", "maximum number of replicates of hash value on "
                           "reference sequence to produce a seedl if not set, "
                           "the number of replicates is unbounded",
      cxxopts::value<size_t>(max_replicate_qry)->
                             default_value(default_max_replicate_qry))

     ("index-cache", "path to cache file for reference index (load if present, "
                     "save if missing)",
      cxxopts::value<std::string>(index_cache_path)->
                          default_value(default_index_cache_path))

     ("build-index-only", "build the reference index cache and exit",
      cxxopts::value<bool>(build_index_only)->default_value("false"))

     ("a,all_vs_all", "compare each of the input files against each other, "
                      "requires to specify at least three input files",
      cxxopts::value<bool>(all_vs_all_option)->default_value("false"))

     ("r,read_pairs", "query sequences are provided as pairs of fastq-files, "
                      "i.e. the number of query files must ba a positive "
                      "multiple of 2 and all pairs of fastq-files must "
                      "contain the same number of sequences",
      cxxopts::value<bool>(read_pairs_option)->default_value("false"))

     ("g,group_by_query", "primary order for sorting the output is by query; "
                          "secondary order is by subject",
      cxxopts::value<bool>(group_by_query_option)->default_value("false"))

     ("d,display", ("specify what to display and in which order; possible "
                    "values are ") +
                   this->display_options.display_help_string(),
      cxxopts::value<std::string>(display_options_spec)->default_value(""))

     ("c,reverse_complement",
      msg_reverse_complement_for,
      cxxopts::value<std::string>(reverse_complement_for_from_argv)->
                                  default_value(default_reverse_complement_for
                                                   .str().c_str()))

     ("i,invint",
        gi_joiner_code
          ? "use invertible integer encoding "
          : "use invertible integer encoding and no sampling of k-mers",
      cxxopts::value<bool>(invint_option)->default_value("false"));

  if (gi_joiner_code)
  {
    options.add_options()
    ("w,window_size", "size of the window from which minimizers are drawn",
      cxxopts::value<int>(window_size)->default_value(default_window_size));
  } else
  {
    options.add_options()
     ("x,store_seeds", "store all seeds, before extending them to "
                       "maximal exact matches",
      cxxopts::value<bool>(store_seeds_option)->default_value("false"))

     ("l,minimum_mem_length", "minimum length of maximal exact match (MEM)",
      cxxopts::value<int>(minimum_mem_length)
                        ->default_value(default_minimum_mem_length))

     ("m,minimum_als_length", "minimum length of aligned sequences",
      cxxopts::value<int>(minimum_als_length)
                         ->default_value(default_minimum_als_length))

     ("e,maximum_error_percentage", "maximum percentage of errors of "
                                    "aligned sequences",
      cxxopts::value<double>(maximum_error_percentage)
                            ->default_value(default_maximum_error_percentage))

     ("f,final_polishing", "polish coordinates of aligned substrings obtained "
                           "from extending seeds",
      cxxopts::value<bool>(final_polishing_option)->default_value("false"))

     ("p,polishing_error_percentage", "percentage value for polishing errors, "
                                      "if not specified, then "
                                      "maximum_error_percentage is used",
      cxxopts::value<double>(polishing_error_percentage)
                            ->default_value(default_polishing_error_percentage))

     ("cop", "use constant co-prime distances for sampling seeds",
      cxxopts::value<bool>(cop_option)->default_value("false"));
   }

   options.add_options()
     ("v,verbose", "be verbose, i.e. show 'Options:'-line",
      cxxopts::value<bool>(verbose_option)->default_value("false"))

     ("h,help", "print usage");
  try
  {
    auto result = options.parse(argc, argv);

    if (result.count("help") > 0)
    {
      help_option = true;
      usage(options);
    }
    const std::vector<std::string>& unmatched_args = result.unmatched();
    if (unmatched_args.size() < 1)
    {
      throw std::invalid_argument("missing positional reference file argument");
    }
    for (size_t idx = 0; idx < unmatched_args.size(); idx++)
    {
      this->input_files.push_back(unmatched_args[idx]);
    }
    try
    {
      this->reverse_complement_for.set(reverse_complement_for_from_argv);
    }
    catch (const std::exception &err)
    {
      throw std::invalid_argument(
              std::string("option -c/--reverse_complement: ") + err.what());
    }
    if (this->kmer_size == -1)
    {
      if (gi_joiner_code)
      {
        this->kmer_size = 22;
      } else
      {
        constexpr const int min_len2_kmer_length[] = {
          1 /* 0 */,
          1 /* 1 */,
          1 /* 2 */,
          1 /* 3 */,
          2 /* 4 */,
          3 /* 5 */,
          4 /* 6 */,
          5 /* 7 */,
          6 /* 8 */,
          7 /* 9 */,
          8 /* 10 */,
          9 /* 11 */,
          10 /* 12 */,
          11 /* 13 */,
          12 /* 14 */,
          13 /* 15 */,
          14 /* 16 */,
          14 /* 17 */,
          14 /* 18 */,
          14 /* 19 */,
          16 /* 20 */,
          16 /* 21 */,
          18 /* 22 */,
          19 /* 23 */,
          19 /* 24 */,
          19 /* 25 */
        };
        if (this->minimum_mem_length <
            static_cast<int>(sizeof min_len2_kmer_length/
                             sizeof min_len2_kmer_length[0]))
        {
          this->kmer_size = min_len2_kmer_length[this->minimum_mem_length];
          if (this->reverse_complement_for_get().is(for_reference) &&
              this->kmer_size % 2 == 1)
          {
            this->kmer_size--;
          }
        } else
        {
          this->kmer_size = size_t(22);
        }
      }
    }
    if (gi_joiner_code)
    {
      if (this->kmer_size < 2)
      {
        throw std::invalid_argument("kmer length must be >= 2");
      }
      if (this->window_size <= 0)
      {
        throw std::invalid_argument("window size must be > 0");
      }
    } else
    {
      if (this->kmer_size < 1 || this->kmer_size > this->minimum_mem_length)
      {
        throw std::invalid_argument("kmer length must be > 0 and not larger "
                                    "than minimal length of MEM");
      }
    }
    if (this->hash_bits != -1)
    {
      if (this->invint_option)
      {
        throw std::invalid_argument("if option --invint is used then "
                                    "option -b,--hash_bits cannot be used");
      } else
      {
        if (this->hash_bits < 1)
        {
          throw std::invalid_argument("number of hash bits must be > 0");
        }
      }
    }
    if (this->number_of_threads < 1)
    {
      throw std::invalid_argument("number of threads must be > 0");
    }
    if (this->thread_pool_max_size < 0)
    {
      throw std::invalid_argument("thread_pool_max_size must be > 0");
    }
    if (this->thread_pool_max_size > 0 and this->number_of_threads == 1)
    {
      throw std::invalid_argument("option -s/--thread_pool_max_size can only "
                                  "be used if the number of threads is at "
                                  "least 2");
    }
    if (this->all_vs_all_option_is_set() &&
        this->index_cache_path.size() > 0)
    {
      throw std::invalid_argument("option --index-cache cannot be combined "
                                  "with option -a/--all_vs_all");
    }
    if (this->build_index_only && this->index_cache_path.size() == 0)
    {
      throw std::invalid_argument("option --build-index-only requires option "
                                  "--index_cache to be set");
    }
    if (this->build_index_only && this->all_vs_all_option_is_set())
    {
      throw std::invalid_argument("option --build-index-only cannot be "
                                  "combined with option -a/--all_vs_all");
    }
    this->display_options.set_flags(this->display_options_spec);
    if (this->number_of_threads == 1 &&
        this->threads_out_prefix != default_threads_out_prefix)
    {
      throw std::invalid_argument("option -o/--threads_out_prefix can only be "
                                  "used if the number of threads is at least "
                                  "2");
    }
    if (this->all_vs_all_option_is_set() && this->read_pairs_option_is_set())
    {
      throw std::invalid_argument("options -a/--all_vs_all and -r/--read_pairs "
                                  "exclude each other");
    }
    if (this->number_of_threads > 1)
    {
      if (this->reverse_complement_for_get().is(for_reference_canonical))
      {
        throw std::invalid_argument("option '-c/--reverse_complement "
                                    "for_reference_canonical' can only be "
                                    "used with one thread");
      }
      if (this->display_options.minimizers_display())
      {
        throw std::invalid_argument("display flag minimizers_display "
                                    "can only be used with one thread");
      }
      if (this->display_options.chain_elements_display())
      {
        throw std::invalid_argument("display flag chain_elements "
                                    "can only be used with one thread");
      }
      if (this->display_options.chain_closed_display())
      {
        throw std::invalid_argument("display flag chain_closed "
                                    "can only be used with one thread");
      }
      if (this->input_files.size() == 1)
      {
        throw std::invalid_argument("option -t/--threads "
                                    "can only be used if multiple input files "
                                    "are specified");
      } else
      {
        if (this->all_vs_all_option_is_set())
        {
          throw std::invalid_argument("option -t/--threads "
                                      "can not be used with option "
                                      "-a/--all_vs_all");
        }
      }
    }
    if (not gi_joiner_code and this->polishing_error_percentage == -1.0)
    {
      this->polishing_error_percentage = this->maximum_error_percentage;
    }
    if (this->thread_pool_max_size == 0 and this->number_of_threads > 1)
    {
      this->thread_pool_max_size = this->number_of_threads;
    }
  }
  catch (const cxxopts::exceptions::exception &err)
  {
    usage(options);
    throw std::invalid_argument(err.what());
  }
}

size_t Options::minimum_als_length_get(void) const noexcept
{
  return static_cast<size_t>(this->minimum_als_length);
}

size_t Options::minimum_mem_length_get(void) const noexcept
{
  return static_cast<size_t>(this->minimum_mem_length);
}

double Options::maximum_error_percentage_get(void) const noexcept
{
  return this->maximum_error_percentage;
}

double Options::polishing_error_percentage_get(void) const noexcept
{
  return this->polishing_error_percentage;
}

size_t Options::kmer_size_get(void) const noexcept
{
  return static_cast<size_t>(this->kmer_size);
}

size_t Options::number_of_threads_get(void) const noexcept
{
  return static_cast<size_t>(this->number_of_threads);
}

size_t Options::thread_pool_max_size_get(void) const noexcept
{
  return static_cast<size_t>(this->thread_pool_max_size);
}

size_t Options::query_split_size_get(void) const noexcept
{
  return static_cast<size_t>(this->query_split_size);
}

size_t Options::max_replicate_ref_get(void) const noexcept
{
  return static_cast<size_t>(this->max_replicate_ref);
}

size_t Options::max_replicate_qry_get(void) const noexcept
{
  return static_cast<size_t>(this->max_replicate_qry);
}

int Options::hash_bits_get(void) const noexcept
{
  return this->hash_bits;
}

size_t Options::window_size_get(void) const noexcept
{
  return static_cast<size_t>(this->window_size);
}

const std::vector<std::string> &Options::input_files_get(void) const noexcept
{
  return this->input_files;
}

bool Options::help_option_is_set(void) const noexcept
{
  return this->help_option;
}

bool Options::verbose_option_is_set(void) const noexcept
{
  return this->verbose_option;
}

bool Options::all_vs_all_option_is_set(void) const noexcept
{
  return this->all_vs_all_option;
}

bool Options::store_seeds_option_is_set(void) const noexcept
{
  return this->store_seeds_option;
}

bool Options::cop_option_is_set(void) const noexcept
{
  return this->cop_option;
}

bool Options::invint_option_is_set(void) const noexcept
{
  return this->invint_option;
}

bool Options::read_pairs_option_is_set(void) const noexcept
{
  return this->read_pairs_option;
}

bool Options::final_polishing_option_is_set(void) const noexcept
{
  return this->final_polishing_option;
}

bool Options::group_by_query_option_is_set(void) const noexcept
{
  return this->group_by_query_option;
}

const std::string &Options::threads_out_prefix_get(void) const noexcept
{
  return this->threads_out_prefix;
}

bool Options::build_index_only_is_set(void) const noexcept
{
  return this->build_index_only;
}

const std::string &Options::index_cache_path_get(void) const noexcept
{
  return this->index_cache_path;
}

ReverseComplementFor Options::reverse_complement_for_get(void) const noexcept
{
  return this->reverse_complement_for;
}

const DisplayOptions &Options::display_options_get(void) const noexcept
{
  return this->display_options;
}
