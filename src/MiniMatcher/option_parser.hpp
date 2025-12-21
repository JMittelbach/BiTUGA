#ifndef OPTION_PARSER_HPP
#define OPTION_PARSER_HPP
#include <cstddef>
#include <vector>
#include "display_options.hpp"
#include "reverse_complement_for.hpp"

class Options
{
  private:
  std::vector<std::string> input_files;
  const std::string default_minimum_mem_length = "19";
  const std::string default_minimum_als_length = "100";
  const std::string default_maximum_error_percentage = "15.0";
  const std::string default_polishing_error_percentage = "-1";
  const std::string default_kmer_size = "-1";
  const std::string default_hash_bits = "-1";
  const std::string default_window_size = "1";
  const std::string default_number_of_threads = "1";
  const std::string default_query_split_size = "0";
  const std::string default_thread_pool_max_size = "0";
  const std::string default_threads_out_prefix = "";
  const std::string default_index_cache_path = "";
  const std::string default_max_replicate_ref = "0";
  const std::string default_max_replicate_qry = "0";

  std::string display_options_spec;
  int minimum_mem_length;
  int minimum_als_length;
  int kmer_size;
  int hash_bits;
  int window_size;
  int number_of_threads;
  int thread_pool_max_size;
  int query_split_size;
  size_t max_replicate_ref;
  size_t max_replicate_qry;
  double maximum_error_percentage;
  double polishing_error_percentage;
  bool gi_joiner_code = false;
  bool verbose_option = false;
  bool help_option = false;
  bool read_pairs_option = false;
  bool final_polishing_option = false;
  bool group_by_query_option = false;
  bool all_vs_all_option = false;
  bool store_seeds_option = false;
  bool cop_option = false;
  bool invint_option = false;
  bool build_index_only = false;
  std::string index_cache_path;
  std::string reverse_complement_for_from_argv;
  std::string threads_out_prefix;
  ReverseComplementFor reverse_complement_for{};
  DisplayOptions display_options{};

  public:
  Options(bool _gi_joiner_code)
   : gi_joiner_code(_gi_joiner_code)
  { }

  void parse(int argc, char **argv);
  size_t minimum_mem_length_get(void) const noexcept;
  size_t minimum_als_length_get(void) const noexcept;
  double maximum_error_percentage_get(void) const noexcept;
  double polishing_error_percentage_get(void) const noexcept;
  size_t kmer_size_get(void) const noexcept;
  size_t number_of_threads_get(void) const noexcept;
  size_t thread_pool_max_size_get(void) const noexcept;
  size_t query_split_size_get(void) const noexcept;
  size_t max_replicate_ref_get(void) const noexcept;
  size_t max_replicate_qry_get(void) const noexcept;
  int hash_bits_get(void) const noexcept;
  size_t window_size_get(void) const noexcept;
  const std::vector<std::string> &input_files_get(void) const noexcept;
  bool help_option_is_set(void) const noexcept;
  bool verbose_option_is_set(void) const noexcept;
  bool all_vs_all_option_is_set(void) const noexcept;
  bool read_pairs_option_is_set(void) const noexcept;
  bool final_polishing_option_is_set(void) const noexcept;
  bool group_by_query_option_is_set(void) const noexcept;
  bool store_seeds_option_is_set(void) const noexcept;
  bool cop_option_is_set(void) const noexcept;
  bool invint_option_is_set(void) const noexcept;
  const std::string &threads_out_prefix_get(void) const noexcept;
  bool build_index_only_is_set(void) const noexcept;
  const std::string &index_cache_path_get(void) const noexcept;
  ReverseComplementFor reverse_complement_for_get(void) const noexcept;
  const DisplayOptions &display_options_get(void) const noexcept;
};
#endif
