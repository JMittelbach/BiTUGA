/*
  Copyright (c) 2021-2025 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2021-2025 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#include <algorithm>                                   // for min
#include <array>                                       // for array
#include <bit>                                         // for bit_width
#include <cassert>                                     // for assert
#include <cinttypes>                                   // for uint8_t, UINT8...
#include <cmath>                                       // for ceil
#include <climits>                                     // for UCHAR_MAX, CHA...
#include <cstdio>                                      // for size_t, FILE
#include <cstdlib>                                     // for EXIT_FAILURE
#include <exception>                                   // for exception
#include <format>                                      // for format
#include <filesystem>                                  // for exists
#include <iostream>                                    // for cerr
#include <limits>                                      // for numeric_limits
#include <memory>                                      // for unique_ptr
#include <stdexcept>                                   // for runtime_error
#include <string>                                      // for basic_string
#include <type_traits>                                 // for conditional
#include <utility>                                     // for pair, get
#include <vector>                                      // for vector
#include "display_options.hpp"                         // for DisplayOptions
#include "enum_position_pair_matches.hpp"              // for SelfMatchMode
#include "guide_params.hpp"                            // for GuideParams
#include "match_class.hpp"                             // for Match
#include "option_parser.hpp"                           // for Options
#include "reverse_complement_for.hpp"                  // for EnumReverseCom...
#include "sequences/gttl_fasta_generator.hpp"          // for GttlFastAGener...
#include "sequences/gttl_fastq_generator.hpp"          // for GttlFastQGener...
#include "sequences/gttl_multiseq.hpp"                 // for GttlMultiseq
#include "sequences/guess_if_protein_seq.hpp"          // for guess_if_prote...
#include "sequences/hashed_qgrams.hpp"                 // for HashedQgramsGe...
#include "sequences/multiseq_factory.hpp"              // for GttlMultiseqFa...
#include "sequences/multiseq_generator.hpp"            // for GttlMultiseqGe...
#include "sequences/qgrams_hash_invint.hpp"            // for InvertibleInte...
#include "sequences/qgrams_hash_nthash.hpp"            // for QgramNtHashFwd...
#include "sequences/sorted_match_list.hpp"             // for SortedMatchList
#include "threading/thread_pool_var.hpp"               // for GttlThreadPoolVar
#include "threading/threads_output_files.hpp"          // for ThreadsOutputF...
#include "threading/thread_specific_index.hpp"         // for ThreadSpecificIn..
#include "threading/thread_pool_unknown_tasks.hpp"     // for ThreadPoolUnknow..
#include "utilities/argv2string.hpp"                   // for gttl_argv2string
#include "utilities/argv_concat.hpp"                   // for ArgvConcat
#include "utilities/constexpr_for.hpp"                 // for constexpr_for
#include "utilities/duplicated_filenames.hpp"          // for gttl_duplicate...
#include "utilities/gttl_file_open.hpp"                // for GttlFpType
#include "utilities/has_fasta_or_fastq_extension.hpp"  // for gttl_likely_fa...
#include "utilities/mathsupport.hpp"                   // for mega_bytes
#include "utilities/runtime_class.hpp"                 // for RunTimeClass
#ifndef GI_JOINER_CODE
#include "sequences/sorted_match_list.hpp"
#include "match_class.hpp"
#define IGNORE_CHAINING
#ifndef IGNORE_CHAINING
#include "segment_chaining.hpp"
#endif
#endif

#ifdef GI_JOINER_CODE
#define INSTANTIATE_NT_HASH_STORE_SEEDS_USE_DIAGONAL_ENCODER(FUNC,...)\
        if (options.invint_option_is_set())\
        {\
          constexpr const bool with_nt_hash_local = false;\
          constexpr const bool store_seeds_local = false;\
          constexpr const bool use_diagonal_encoder_local = false;\
          FUNC<with_nt_hash_local,store_seeds_local,\
               use_diagonal_encoder_local>(__VA_ARGS__);\
        } else\
        {\
          constexpr const bool with_nt_hash_local = true;\
          constexpr const bool store_seeds_local = false;\
          constexpr const bool use_diagonal_encoder_local = false;\
          FUNC<with_nt_hash_local,store_seeds_local,\
               use_diagonal_encoder_local>(__VA_ARGS__);\
        }
#else
#define INSTANTIATE_NT_HASH_STORE_SEEDS_USE_DIAGONAL_ENCODER(FUNC,...)\
        if (options.invint_option_is_set())\
        {\
          constexpr const bool with_nt_hash_local = false;\
          if (options.store_seeds_option_is_set())\
          {\
            constexpr const bool store_seeds_local = true;\
            if (!guide_params.perform_qgram_sampling())\
            {\
              constexpr const bool use_diagonal_encoder_local = true;\
              FUNC<with_nt_hash_local,store_seeds_local,\
                   use_diagonal_encoder_local>(__VA_ARGS__);\
            } else\
            {\
              constexpr const bool use_diagonal_encoder_local = false;\
              FUNC<with_nt_hash_local,store_seeds_local,\
                   use_diagonal_encoder_local>(__VA_ARGS__);\
            }\
          } else\
          {\
            constexpr const bool store_seeds_local = false;\
            constexpr const bool use_diagonal_encoder_local = false;\
            FUNC<with_nt_hash_local,store_seeds_local,\
                 use_diagonal_encoder_local>(__VA_ARGS__);\
          }\
        } else\
        {\
          constexpr const bool with_nt_hash_local = true;\
          constexpr const bool use_diagonal_encoder_local = false;\
          if (options.store_seeds_option_is_set())\
          {\
            constexpr const bool store_seeds_local = true;\
            FUNC<with_nt_hash_local,store_seeds_local,\
                 use_diagonal_encoder_local>(__VA_ARGS__);\
          } else\
          {\
            constexpr const bool store_seeds_local = false;\
            FUNC<with_nt_hash_local,store_seeds_local,\
                 use_diagonal_encoder_local>(__VA_ARGS__);\
          }\
        }
#endif

static void check_hash_bits_log(std::vector<std::string> *log_vector,
                                const char *tag,
                                int hash_bits,
                                int sequences_bits,
                                bool with_nt_hash)
{
  constexpr const int max_hashed_qgram_bits = 72;
  const int hashed_qgram_bits = hash_bits + sequences_bits;
  if (hashed_qgram_bits > max_hashed_qgram_bits) /* check_err.py */
  {
    if (with_nt_hash)
    {
      throw std::runtime_error(
              std::format(": number of hash bits is too large. Maximum value "
                          "for option -b,--hash_bits for this data set is {}",
                          max_hashed_qgram_bits - sequences_bits));

    } else
    {
      throw std::runtime_error(
              std::format(": kmer size is too large. Maximum value for "
                          "option -k,--kmer_size for this data set is {}",
                          (max_hashed_qgram_bits - sequences_bits)/2));
    }
  }
  if (log_vector != nullptr)
  {
    log_vector->push_back("hash_bits\t" + std::to_string(hash_bits));
    log_vector->push_back("ref_sequences_bits\t" +
                          std::to_string(sequences_bits));
    log_vector->push_back(std::string(tag) + "_hashed_qgram_bits\t"
                          + std::to_string(hashed_qgram_bits));
    log_vector->push_back(std::string(tag)
                          + "_hashed_qgram_bytes\t"
                          + std::to_string((hashed_qgram_bits + 7)/8));
  }
}

static void log_vector_show(FILE *out_fp,
                            const std::vector<std::string> *log_vector)
{
  assert(log_vector != nullptr);
  for (size_t idx = 0; idx < log_vector->size(); idx++)
  {
    fprintf(out_fp,"# %s\n",log_vector->at(idx).c_str());
  }
}

static void append_cache_log(std::vector<std::string> *log_vector,
                             const std::string &tag,
                             const std::string &cache_path,
                             size_t num_hashed_qgrams,
                             size_t count_all_qgrams,
                             size_t sizeof_unit)
{
  if (log_vector == nullptr)
  {
    return;
  }
  log_vector->push_back(tag + "\t" + cache_path);
  log_vector->push_back("number of hashed kmers\t" +
                        std::to_string(num_hashed_qgrams));
  if (count_all_qgrams > 0)
  {
    const double density
      = static_cast<double>(num_hashed_qgrams) /
        static_cast<double>(count_all_qgrams);
    log_vector->push_back(std::format("hashed kmers density\t{:.2f}", density));
  }
  const double space_in_mega_bytes
    = mega_bytes(num_hashed_qgrams * sizeof_unit);
  log_vector->push_back(std::format("SPACE\thashed kmers (MB)\t{:.0f}",
                                    space_in_mega_bytes));
}

static void check_hash_bits(FILE *out_fp,
                            const char *tag,
                            int hash_bits,
                            int sequences_bits,
                            bool with_nt_hash)
{
  if (out_fp != nullptr)
  {
    std::vector<std::string> log_vector;
    check_hash_bits_log(&log_vector,tag,hash_bits,sequences_bits,with_nt_hash);
    log_vector_show(out_fp,&log_vector);
  } else
  {
    check_hash_bits_log(nullptr,tag,hash_bits,sequences_bits,with_nt_hash);
  }
}

static void check_sum_sequences_bits(int ref_sequence_bits,
                                     int query_sequence_bits)
{
  if (ref_sequence_bits + query_sequence_bits > 64)
  {
    throw std::runtime_error(
      std::format(": need {} bits for sequences of the reference "
                  "and {} bits for sequences of the query; the "
                  "sum {} exceeds the limit of 64; please split "
                  "the reference or the query sequences",
                  ref_sequence_bits,
                  query_sequence_bits,
                  ref_sequence_bits + query_sequence_bits));
  }
}

static bool check_false_sequence_type(const char *progname,
                                      const char *filename)
{
  bool has_err = false, is_protein = false;
  try
  {
    is_protein = guess_if_protein_file(filename);
  }
  catch (const std::exception &err) /* check_err.py */
  {
    std::cerr << progname << ": file \"" << filename << "\""
              << err.what() << '\n';
    has_err = true;
  }
  if (!has_err && is_protein) /* check_err.py */
  {
    std::cerr << progname << ": file \"" << filename << "\""
              << ": can only handle DNA sequences" << '\n';
    has_err = true;
  }
  return has_err;
}

#ifndef GI_JOINER_CODE
template<class MatchTable>
static void output_encoded_matches(FILE *out_fp,
                                   const DisplayOptions &display_options,
                                   const MatchTable &match_table,
                                   const GttlMultiseq &ref_multiseq,
                                   const GttlMultiseq &query_multiseq)
{
  const size_t query_seqnum_offset
    = query_multiseq.sequence_number_offset_get();
  for (size_t idx = 0; idx < match_table.size(); idx++)
  {
    Match this_match(match_table.ref_endpos_get(idx),
                     match_table.length_get(idx),
                     match_table.query_endpos_get(idx),
                     match_table.length_get(idx),
                     0);
    constexpr const bool fast_output = false;
    this_match.display<fast_output>(out_fp,
                                    display_options,
                                    ref_multiseq,
                                    match_table.ref_seqnum_get(idx),
                                    query_multiseq,
                                    match_table.query_seqnum_get(idx),
                                    query_seqnum_offset,
                                    nullptr);
  }
}
#endif

template<SelfMatchMode self_match,
         class RefHashedQgramsClass,
         class QueryHashedQgramsClass,
         bool store_seeds,
         bool use_diagonal_encoder,
         int sizeof_unit_pospair,
         int ref_idx,
         int sizeof_unit_match>
static void extend_and_process_matches_unit_match(
                                     FILE *out_fp,
#ifdef GI_JOINER_CODE
                                     [[maybe_unused]]
#endif
                                     const GuideParams &guide_params,
                                     const EnumPositionPairMatches
                                       <self_match,
                                        RefHashedQgramsClass,
                                        QueryHashedQgramsClass,
                                        store_seeds,
                                        use_diagonal_encoder,
                                        sizeof_unit_pospair>
                                     &enum_pos_pair_enumerator,
#ifdef GI_JOINER_CODE
                                     [[maybe_unused]]
#endif
                                     const GttlMultiseq &ref_multiseq,
#ifdef GI_JOINER_CODE
                                     [[maybe_unused]]
#endif
                                     const GttlMultiseq &query_multiseq)
{
#ifdef GI_JOINER_CODE
  if (guide_params.display_options_get().only_count_seeds_display())
  {
    size_t count_seeds = 0;
    for ([[maybe_unused]] auto const &pp : enum_pos_pair_enumerator)
    {
      count_seeds++;
    }
    if (guide_params.verbose_option_get())
    {
      fprintf(out_fp,"# number of seeds:\t%zu\n",count_seeds);
    }
  } else
  {
    for (auto const &pp : enum_pos_pair_enumerator)
    {
      printf("%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
             pp.seqnum0,
             pp.startpos0,
             pp.seqnum1,
             pp.startpos1);
    }
  }
#else
  constexpr const bool seed_output = false;
  constexpr const bool from_same_sequence = false;
  const DisplayOptions &display_options = guide_params.display_options_get();
  using ThisEnumPosPairMatches = EnumPositionPairMatches
                                       <self_match,
                                        RefHashedQgramsClass,
                                        QueryHashedQgramsClass,
                                        store_seeds,
                                        use_diagonal_encoder,
                                        sizeof_unit_pospair>;
  using ThisSortedMatchList
    = SortedMatchList<GttlMultiseq,
                      ThisEnumPosPairMatches,
                      self_match != None,
                      from_same_sequence,
                      sizeof_unit_match,
                      seed_output,
                      ref_idx>;
  ThisSortedMatchList sorted_match_list(guide_params.qgram_length_get(),
                                        guide_params.minimum_mem_length_get(),
                                        enum_pos_pair_enumerator,
                                        ref_multiseq,
                                        query_multiseq);
  if (guide_params.verbose_option_get())
  {
    const double space
      = mega_bytes(sorted_match_list.number_of_all_matches_get() *
                   sizeof_unit_match);
    fprintf(out_fp,"# SPACE\tencoded_match_list (MB):\t%.0f\n",space);
  }
  if (display_options.chain())
  {
#ifndef IGNORE_CHAINING
    segment_chaining<ThisSortedMatchList,Match>
                    (guide_params,
                     ref_multiseq,
                     query_multiseq,
                     self_match == Regular,
                     sorted_match_list,
                     out_fp);
#endif
  } else
  {
    output_encoded_matches<ThisSortedMatchList>
                          (out_fp,
                           display_options,
                           sorted_match_list,
                           ref_multiseq,
                           query_multiseq);
  }
  if (guide_params.verbose_option_get())
  {
    sorted_match_list.statistics(out_fp);
  }
#endif
}

#define EXTEND_AND_PROCESS_MATCHES_UNIT_MATCH(USE_DIAGONAL_ENCODER,\
                                              SIZEOF_UNIT_MATCH)\
        extend_and_process_matches_unit_match<self_match,\
                                              RefHashedQgramsClass,\
                                              QueryHashedQgramsClass,\
                                              store_seeds,\
                                              USE_DIAGONAL_ENCODER,\
                                              sizeof_unit_pospair,\
                                              ref_idx,\
                                              SIZEOF_UNIT_MATCH>\
                                             (out_fp,\
                                              guide_params,\
                                              enum_pos_pair_enumerator,\
                                              ref_multiseq,\
                                              query_multiseq)

template<SelfMatchMode self_match,
         class RefHashedQgramsClass,
         class QueryHashedQgramsClass,
         bool store_seeds,
         bool use_diagonal_encoder,
         int sizeof_unit_pospair,
         int ref_idx>
static void extend_and_process_matches(FILE *out_fp,
                                       const GuideParams &guide_params,
                                       const EnumPositionPairMatches
                                         <self_match,
                                          RefHashedQgramsClass,
                                          QueryHashedQgramsClass,
                                          store_seeds,
                                          use_diagonal_encoder,
                                          sizeof_unit_pospair>
                                         &enum_pos_pair_enumerator,
                                       const GttlMultiseq &ref_multiseq,
                                       const GttlMultiseq &query_multiseq)
{
  constexpr const int uint64_bits
    = static_cast<int>(sizeof(uint64_t) * CHAR_BIT);
  const int bits_for_sequences
    = enum_pos_pair_enumerator.sequences_bits_sum_get();
  if constexpr (use_diagonal_encoder)
  {
    const size_t maximum_match_length
     = enum_pos_pair_enumerator.maximum_match_length_get();
    assert(maximum_match_length >= guide_params.qgram_length_get());
    const int bits_for_match_length
      = static_cast<int>(std::bit_width(maximum_match_length -
                                        guide_params.qgram_length_get()));
    if (bits_for_sequences + bits_for_match_length <= uint64_bits)
    {
      EXTEND_AND_PROCESS_MATCHES_UNIT_MATCH(true,8);
    } else
    {
      if (bits_for_sequences + bits_for_match_length <= uint64_bits + CHAR_BIT)
      {
        EXTEND_AND_PROCESS_MATCHES_UNIT_MATCH(true,9);
      } else
      {
        EXTEND_AND_PROCESS_MATCHES_UNIT_MATCH(true,10);
      }
    }
  } else
  {
    bool use9bytes_version = false;
    if (bits_for_sequences <= uint64_bits)
    {
      try
      {
        EXTEND_AND_PROCESS_MATCHES_UNIT_MATCH(false,8);
      }
      catch (const std::overflow_error &msg)
      {
        if (guide_params.verbose_option_get())
        {
          fprintf(out_fp,"# %s\n",msg.what());
        }
        use9bytes_version = true;
      }
    } else
    {
      use9bytes_version = true;
    }
    bool use10bytes_version = false;
    if (use9bytes_version)
    {
      try
      {
        EXTEND_AND_PROCESS_MATCHES_UNIT_MATCH(false,9);
      }
      catch (const std::overflow_error &msg)
      {
        if (guide_params.verbose_option_get())
        {
          fprintf(out_fp,"# %s\n",msg.what());
        }
        use10bytes_version = true;
      }
    }
    if (use10bytes_version)
    {
      EXTEND_AND_PROCESS_MATCHES_UNIT_MATCH(false,10);
    }
  }
}

static void apply_wilcard_replacement(std::array<uint8_t,UINT8_MAX+1>
                                        &wildcard_replacement_arr,
                                      char *char_seq,size_t len)
{
  uint8_t *byte_seq = reinterpret_cast<uint8_t *>(char_seq);
  for (size_t idx = 0; idx < len; idx++)
  {
    byte_seq[idx] = wildcard_replacement_arr[static_cast<int>(char_seq[idx])];
  }
}

static constexpr const char nucleotides_upper_lower[] = "Aa|Cc|Gg|TtUu";

static void wildcard_replacement(GttlMultiseq *multiseq,uint8_t undef)
{
  std::array<uint8_t,UINT8_MAX+1> wildcard_replacement_arr;
  wildcard_replacement_arr.fill(undef);
  constexpr const size_t spec_size = sizeof nucleotides_upper_lower - 1;
  for (size_t idx = 0; idx < spec_size; idx++)
  {
    if (nucleotides_upper_lower[idx] != '|')
    {
      const uint8_t r = static_cast<uint8_t>(nucleotides_upper_lower[idx]);
      assert(r != undef);
      wildcard_replacement_arr[r] = r;
    }
  }
  multiseq->transformer<std::array<uint8_t,UINT8_MAX+1>,
                        apply_wilcard_replacement>
                       (wildcard_replacement_arr);
}

template<bool rc_opt>
using NtHashHashedQgramsGeneric
  = typename std::conditional<rc_opt,
                              QgramNtHashIterator4,
                              QgramNtHashFwdIterator4>::type;

template<bool rc_opt>
using InvertibleIntegerHashedQgramsGeneric
  = typename std::conditional<rc_opt,
                              InvertibleIntegercode2Iterator4,
                              InvertibleIntegercodeIterator4>::type;

template<int sizeof_unit,bool rc_opt,bool with_nt_hash>
using HashedQgrams
  = HashedQgramsGeneric<sizeof_unit,
                        typename std::conditional
                          <with_nt_hash,
                           NtHashHashedQgramsGeneric<rc_opt>,
                           InvertibleIntegerHashedQgramsGeneric<rc_opt>>::type>;

/* As the minimizer concept is based on Nthash and this only allows
   to process DNA sequences, it is safe to assume that we can use value >= 4
   as passing characters. So we can use 0,1,2,3 for the nucleotides. */
static constexpr const uint8_t smallest_undef_char = 4;

template<SelfMatchMode self_match,
         class RefHashedQgramsClass,
         class QueryHashedQgramsClass,
         bool store_seeds,
         bool use_diagonal_encoder,
         int sizeof_unit_pospair,
         int ref_idx>
static void hashed_qgrams_seeded_reference_vs_query_matches_generic(
                   FILE *out_fp,
                   const RefHashedQgramsClass &rhqg,
                   const QueryHashedQgramsClass &qhqg,
                   const GuideParams &guide_params)
{
  RunTimeClass rt_merge_extension{};
  if (guide_params.verbose_option_get())
  {
    fprintf(out_fp,"# enumerate position pairs matches for %zu hashed qgrams "
                   "from the reference multiseq versus %zu hashed qgrams from "
                   "the query\n",rhqg.size(),qhqg.size());
  }
  EnumPositionPairMatches<self_match,
                          RefHashedQgramsClass,
                          QueryHashedQgramsClass,
                          store_seeds,
                          use_diagonal_encoder,
                          sizeof_unit_pospair>
    enum_pos_pair_enumerator(rhqg,qhqg,out_fp);
  if constexpr (store_seeds)
  {
    if (guide_params.verbose_option_get())
    {
      const double space = mega_bytes(enum_pos_pair_enumerator.size_in_RAM());
      fprintf(out_fp,"# SPACE\tseeds.ref_vs_query (MB):\t%.0f\n",space);
      fprintf(out_fp,"# %s\n",
              rt_merge_extension.to_string("generation of seeds").c_str());
    }
    rt_merge_extension.reset();
  }
  extend_and_process_matches<self_match,
                             RefHashedQgramsClass,
                             QueryHashedQgramsClass,
                             store_seeds,
                             use_diagonal_encoder,
                             sizeof_unit_pospair,
                             ref_idx>
                            (out_fp,
                             guide_params,
                             enum_pos_pair_enumerator,
                             rhqg.multiseq,
                             qhqg.multiseq);
  if (guide_params.verbose_option_get())
  {
    fprintf(out_fp,"# %s\n",
            rt_merge_extension
              .to_string(store_seeds ? "extension and output of matches"
                                     : ("generation of seeds, extension "
                                        "and output of matches")).c_str());
  }
}

template<class RefHashedQgramsClass,
         class QueryHashedQgramsClass,
         bool store_seeds,
         bool use_diagonal_encoder,
         int sizeof_unit_pospair,
         int ref_idx>
static void hashed_qgrams_seeded_reference_vs_query_matches(
                   SelfMatchMode var_self_match,
                   FILE *out_fp,
                   const RefHashedQgramsClass &rhqg,
                   GttlMultiseq *query_multiseq,
                   const GuideParams &guide_params)
{
  if (var_self_match == Regular)
  {
    constexpr const SelfMatchMode self_match_local = Regular;
    hashed_qgrams_seeded_reference_vs_query_matches_generic
               <self_match_local,
                RefHashedQgramsClass,
                RefHashedQgramsClass,
                store_seeds,
                use_diagonal_encoder,
                sizeof_unit_pospair,
                ref_idx>
               (out_fp,
                rhqg,
                rhqg,
                guide_params);
  } else
  {
    if (var_self_match == WithReverseComplement)
    {
      constexpr const SelfMatchMode self_match_local = WithReverseComplement;
      hashed_qgrams_seeded_reference_vs_query_matches_generic
               <self_match_local,
                RefHashedQgramsClass,
                RefHashedQgramsClass,
                store_seeds,
                use_diagonal_encoder,
                sizeof_unit_pospair,
                ref_idx>
               (out_fp,
                rhqg,
                rhqg,
                guide_params);
    } else
    {
      /* merge join algorithm requires sorted streams with values sorted by
         the hash-value. So we set sort_by_hashvalue to true */
      constexpr const bool sort_by_hashvalue = true;
      std::vector<std::string> log_vector;
      QueryHashedQgramsClass qhqg(*query_multiseq,
                                  guide_params.number_of_threads_get(),
                                  guide_params.qgram_length_get(),
                                  guide_params.query_window_size_get(),
                                  guide_params.hash_bits_get(),
                                  sort_by_hashvalue,
                                  guide_params.at_constant_distance_get(),
                                  guide_params.max_replicate_qry_get(),
                                  guide_params.verbose_option_get()
                                    ? &log_vector : nullptr);
      if (guide_params.verbose_option_get())
      {
        log_vector_show(out_fp,&log_vector);
        log_vector.clear();
      }
      if (guide_params.display_options_get().minimizers_display())
      {
        qhqg.show();
      }
      assert(rhqg.multiseq.padding_char_get() !=
             query_multiseq->padding_char_get());
      /* We only need to perform replacement of wildcards if
         the reference contains wildcards (which have already been replaced)
         and the query contains wildcards */
      if (rhqg.sequence_has_wildcards() &&
          qhqg.sequence_has_wildcards())
      {
        wildcard_replacement(query_multiseq,smallest_undef_char + 1);
      }
      constexpr const SelfMatchMode self_match_local = None;
      hashed_qgrams_seeded_reference_vs_query_matches_generic
                 <self_match_local,
                  RefHashedQgramsClass,
                  QueryHashedQgramsClass,
                  store_seeds,
                  use_diagonal_encoder,
                  sizeof_unit_pospair,
                  ref_idx>
                 (out_fp,
                  rhqg,
                  qhqg,
                  guide_params);
    }
  }
}

static constexpr const bool store_sequence = true;

template<class RefHashedQgramsClass,
         bool store_seeds,
         bool use_diagonal_encoder>
static void process_single_query_multiseq(FILE *out_fp,
                                          const int var_ref_idx,
                                          const RefHashedQgramsClass &rhqg,
                                          GttlMultiseq *query_multiseq,
                                          const GuideParams &guide_params)
{
  constexpr const bool with_nt_hash
    = RefHashedQgramsClass::possible_false_positive_matches;
  /* Both HashedQgramClass need to use the same hashing method */
  check_hash_bits(guide_params.verbose_option_get() ? out_fp : nullptr,
                  "query",
                  guide_params.hash_bits_get(),
                  query_multiseq->sequences_bits_get(),
                  with_nt_hash);
  check_sum_sequences_bits(rhqg.multiseq.sequences_bits_get(),
                           query_multiseq->sequences_bits_get());
  const double space = std::ceil(mega_bytes(query_multiseq->size_in_bytes()));
  if (guide_params.verbose_option_get())
  {
    fprintf(out_fp,"# SPACE\tquery_multiseq (MB)\t%.0f\n",space);
  }
  if (guide_params.display_options_get().q_seqid_display())
  {
    query_multiseq->short_header_cache_create<'|','|'>();
  }
  const int query_sizeof_unit_hashed_qgram
    = (guide_params.hash_bits_get() +
       query_multiseq->sequences_bits_get() <= 64) ? 8 : 9;
  constexpr_for<8,9+1,1>([&](auto ct_query_sizeof_unit_hashed_qgram)
  {
    if (query_sizeof_unit_hashed_qgram == ct_query_sizeof_unit_hashed_qgram)
    {
      assert(var_ref_idx == 0 or var_ref_idx == 1);
      constexpr_for<0,1+1,1>([&](auto ct_ref_idx)
      {
        if (var_ref_idx == ct_ref_idx)
        {
          constexpr const bool query_rc_opt = false;
          using QueryHashedQgrams
            = HashedQgrams<ct_query_sizeof_unit_hashed_qgram,
                           query_rc_opt,
                           with_nt_hash>;
          constexpr const SelfMatchMode self_match_local = None;
          hashed_qgrams_seeded_reference_vs_query_matches<RefHashedQgramsClass,
                                                          QueryHashedQgrams,
                                                          store_seeds,
                                                          use_diagonal_encoder,
                                                          8,
                                                          ct_ref_idx>
                                                         (self_match_local,
                                                          out_fp,
                                                          rhqg,
                                                          query_multiseq,
                                                          guide_params);
        }
      });
    }
  });
}

template<class SequenceGeneratorClass,
         bool fastq_paired,
         class RefHashedQgramsClass,
         bool store_seeds,
         bool use_diagonal_encoder>
static void match_query_against_reference(FILE *out_fp,
                                          const std::vector<std::string>
                                            &query_inputfiles,
                                          int var_ref_idx,
                                          const RefHashedQgramsClass &rhqg,
                                          const GuideParams &guide_params,
                                          bool store_header,
                                          size_t max_num_sequences,
                                          bool padding_char,
                                          uint8_t grant_owner_ship)
{
  auto ms_gen = GttlMultiseqGenerator<SequenceGeneratorClass,fastq_paired>
                                     (query_inputfiles,
                                      store_header,
                                      max_num_sequences,
                                      padding_char,
                                      grant_owner_ship);
  for (GttlMultiseq *query_multiseq : ms_gen)
  {
    if (guide_params.verbose_option_get())
    {
      fprintf(out_fp, "# process query_multiseq with seqnum offset %zu of size "
                      "%.0f (MB)\n",
                      query_multiseq->sequence_number_offset_get(),
                      mega_bytes(query_multiseq->size_in_bytes()));
    }
    process_single_query_multiseq<RefHashedQgramsClass,
                                  store_seeds,
                                  use_diagonal_encoder>
                                 (out_fp,
                                  var_ref_idx,
                                  rhqg,
                                  query_multiseq,
                                  guide_params);
  }
}

template<class SequenceGeneratorClass,
         bool fastq_paired,
         class RefHashedQgramsClass,
         bool store_seeds,
         bool use_diagonal_encoder>
static void threaded_match_query_against_reference(
                                          ThreadsOutputFiles
                                            *threads_output_files,
                                          size_t max_size_of_queue,
                                          const std::vector<std::string>
                                            &query_inputfiles,
                                          int var_ref_idx,
                                          const RefHashedQgramsClass &rhqg,
                                          const GuideParams &guide_params,
                                          bool store_header,
                                          size_t max_num_sequences,
                                          uint8_t padding_char)
{
  const bool grant_owner_ship = false;
  auto ms_gen = GttlMultiseqGenerator<SequenceGeneratorClass,fastq_paired>
                                     (query_inputfiles,
                                      store_header,
                                      max_num_sequences,
                                      padding_char,
                                      grant_owner_ship);
  const size_t num_threads = guide_params.number_of_threads_get();
  ThreadPoolUnknownTasks thread_pool(num_threads);
  ThreadSpecificIndex thread_specific_index(num_threads);
  const std::vector<FILE *> output_filepointers
    = threads_output_files->filepointers_vector_get();
  std::vector<size_t> solved_tasks(num_threads,0);
  for (GttlMultiseq *query_multiseq : ms_gen)
  {
    while (thread_pool.size_of_queue() > max_size_of_queue - 1)
    {
      using namespace std::chrono_literals;
      std::this_thread::sleep_for(40ms); // NOLINT(misc-include-cleaner)
    }
    thread_pool.enqueue([&thread_specific_index,
                         &output_filepointers,
                         query_multiseq,
                         var_ref_idx,
                         rhqg,
                         guide_params,
                         &solved_tasks] {
      const size_t thread_idx = thread_specific_index.get();
      assert(thread_idx < output_filepointers.size());
      solved_tasks[thread_idx] += query_multiseq->sequences_total_length_get();
      FILE *out_fp = output_filepointers[thread_idx];
      if (guide_params.verbose_option_get())
      {
        fprintf(out_fp, "# start processing query_multiseq with seqnum offset "
                        "%zu of size %.0f (MB)\n",
                        query_multiseq->sequence_number_offset_get(),
                        mega_bytes(query_multiseq->size_in_bytes()));
      }
      process_single_query_multiseq<RefHashedQgramsClass,
                                    store_seeds,
                                    use_diagonal_encoder>
                                   (out_fp,
                                    var_ref_idx,
                                    rhqg,
                                    query_multiseq,
                                    guide_params);
      if (guide_params.verbose_option_get())
      {
        fprintf(out_fp, "# process query_multiseq with seqnum offset %zu "
                        "finished\n",
                        query_multiseq->sequence_number_offset_get());
      }
      delete query_multiseq;
    });
  }
  while (thread_pool.size_of_queue() > 0)
  {
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(40ms); // NOLINT(misc-include-cleaner)
  }
  if (guide_params.verbose_option_get())
  {
    size_t sum_solved = 0;
    for (auto s : solved_tasks)
    {
      sum_solved += s;
    }
    for (size_t thread_idx = 0; thread_idx < num_threads; thread_idx++)
    {
      const size_t s = solved_tasks[thread_idx];
      const double s_ratio = static_cast<double>(s)/sum_solved;
      printf("# thread %zu\t%zu\t%.2f%%\n",thread_idx,s, 100.0 * s_ratio);
    }
  }
}

template<class RefHashedQgramsClass,bool store_seeds,bool use_diagonal_encoder>
static void process_single_query_unit_threaded(size_t thread_id,
                                               size_t task_num,
                                               const std::vector<FILE *>
                                                 &output_filepointers,
                                               const RefHashedQgramsClass &rhqg,
                                               const std::vector<std::string>
                                                 *query_files,
                                               const GttlMultiseqFactory
                                                 *multiseq_factory,
                                               const GuideParams &guide_params,
                                               bool read_pairs_option,
                                               int var_ref_idx)
{
  RunTimeClass rt_create_query_multiseq{};
  GttlMultiseq *query_multiseq;

  if (query_files != nullptr)
  {
    const bool store_header
      = guide_params.display_options_get().q_seqid_display();
    assert(task_num <= UINT8_MAX - 1);
    const uint8_t padding_char = UINT8_MAX - 1 - task_num;
    constexpr const bool with_reverse_complement = false;
    query_multiseq
      = read_pairs_option
          ? new GttlMultiseq(query_files->at(2*task_num), /* CONSTRUCTOR */
                             query_files->at(2*task_num+1),
                             store_header,
                             store_sequence,
                             padding_char)
          : new GttlMultiseq(query_files->at(task_num), /* CONSTRUCTOR */
                             store_header,
                             store_sequence,
                             padding_char,
                             with_reverse_complement);
  } else
  {
    assert(multiseq_factory != nullptr);
    query_multiseq = multiseq_factory->at(task_num);
  }
  FILE *out_fp = output_filepointers[thread_id];
  if (guide_params.verbose_option_get())
  {
    fprintf(out_fp,"# thread\t%zu\ttask\t%zu\n",thread_id,task_num);
    fprintf(out_fp,"# %s\n", rt_create_query_multiseq
                             .to_string("create query multiseq").c_str());
  }
  process_single_query_multiseq<RefHashedQgramsClass,
                                store_seeds,
                                use_diagonal_encoder>
                               (out_fp,
                                var_ref_idx,
                                rhqg,
                                query_multiseq,
                                guide_params);
  if (query_files != nullptr)
  {
    delete query_multiseq; /* as this was created in the if case above */
  }
  /* in the else case, the query_multiseq-objects are owned by the
     pointer multiseq_factory_ptr to the MultiseqFactory,
     which takes care of the deletion */
}

template<bool ref_rc_opt,bool with_nt_hash,bool store_seeds,
         bool use_diagonal_encoder>
static void iterate_over_queryfiles_generic(
                                    int var_ref_sizeof_unit_hashed_qgram,
                                    int var_ref_idx,
                                    GttlMultiseq *ref_multiseq,
                                    const GuideParams &guide_params,
                                    bool read_pairs_option,
                                    const std::vector<std::string> &query_files,
                                    std::vector<std::string> *log_vector,
                                    const std::string &threads_out_prefix,
                                    const std::string &index_cache_path,
                                    bool build_index_only)
{
  constexpr const int buf_size = 1 << 14;
  assert(guide_params.hash_bits_get() != -1);
  constexpr_for<8,9+1,1>([&](auto ref_sizeof_unit_hashed_qgram)
  {
    if (ref_sizeof_unit_hashed_qgram == var_ref_sizeof_unit_hashed_qgram)
    {
      using RefHashedQgrams = HashedQgrams<ref_sizeof_unit_hashed_qgram,
                                           ref_rc_opt,
                                           with_nt_hash>;
      std::unique_ptr<RefHashedQgrams> rhqg_ptr;
      const bool have_cache = not index_cache_path.empty();
      if (have_cache && not build_index_only &&
          std::filesystem::exists(index_cache_path))
      {
        try
        {
          rhqg_ptr = std::make_unique<RefHashedQgrams>(
                       RefHashedQgrams::load_cache(
                         index_cache_path,
                         *ref_multiseq,
                         guide_params.qgram_length_get(),
                         guide_params.hash_bits_get(),
                         guide_params.ref_window_size_get(),
                         guide_params.at_constant_distance_get()));
          append_cache_log(log_vector,
                           "index_cache_load",
                           index_cache_path,
                           rhqg_ptr->size(),
                           rhqg_ptr->count_all_qgrams_get(),
                           static_cast<size_t>(RefHashedQgrams
                                                ::sizeof_unit_const));
        }
        catch (const std::exception &err)
        {
          std::cerr << "# warning: failed to load index cache from "
                    << index_cache_path << ": "
                    << err.what() << '\n';
        }
      }
      if (!rhqg_ptr)
      {
        constexpr const bool sort_by_hashvalue = true;
        rhqg_ptr = std::make_unique<RefHashedQgrams>(
                     *ref_multiseq,
                     guide_params.number_of_threads_get(),
                     guide_params.qgram_length_get(),
                     guide_params.ref_window_size_get(),
                     guide_params.hash_bits_get(),
                     sort_by_hashvalue,
                     guide_params.at_constant_distance_get(),
                     guide_params.max_replicate_ref_get(),
                     log_vector);
        if (have_cache)
        {
          if (!rhqg_ptr->save_cache(index_cache_path,
                                    guide_params.ref_window_size_get(),
                                    guide_params.at_constant_distance_get()))
          {
            std::cerr << "# warning: failed to save index cache to "
                      << index_cache_path << '\n';
          } else
          {
            append_cache_log(log_vector,
                             "index_cache_save",
                             index_cache_path,
                             rhqg_ptr->size(),
                             rhqg_ptr->count_all_qgrams_get(),
                             static_cast<size_t>(RefHashedQgrams
                                                   ::sizeof_unit_const));
          }
        }
      }
      RefHashedQgrams &rhqg = *rhqg_ptr;
      if (build_index_only)
      {
        if (guide_params.verbose_option_get() && log_vector != nullptr)
        {
          log_vector_show(stdout,log_vector);
        }
        return;
      }
      if (guide_params.display_options_get().minimizers_display())
      {
        rhqg.show();
      }
      /* We replace the wildcards in the reference by a different character than
         the one used for the query in case this contains wildcards. */
      if (rhqg.sequence_has_wildcards())
      {
        RunTimeClass rt_wc_replacement{};
        wildcard_replacement(ref_multiseq,smallest_undef_char);
        if (log_vector != nullptr)
        {
          log_vector->push_back(rt_wc_replacement
                                  .to_string("replacement of wildcards"));
        }
      }
      const bool store_header
        = guide_params.display_options_get().q_seqid_display();
      const uint8_t padding_char = ref_multiseq->padding_char_get() - 1;
      if (guide_params.number_of_threads_get() == 1)
      {
        const size_t step = read_pairs_option ? 2 : 1;
        constexpr const bool grant_owner_ship = true;
        const size_t max_num_sequences
          = guide_params.query_split_size_get() == 0
              ? std::numeric_limits<size_t>::max()
              : guide_params.query_split_size_get();
        FILE *out_fp = stdout;
        if (log_vector != nullptr)
        {
          log_vector_show(out_fp,log_vector);
          log_vector->clear();
        }
        for (size_t qfidx = 0; qfidx < query_files.size(); qfidx += step)
        {
          if (read_pairs_option)
          {
            for (size_t this_idx = qfidx; this_idx < qfidx + step; this_idx++)
            {
              if (gttl_likely_fasta_format(query_files[this_idx]))
              {
                throw std::format("the fasta file {} cannot be processed with "
                                  "the option --read_pairs",
                                  query_files[this_idx]);
              }
            }
            constexpr const bool fastq_paired = true;
            const std::vector<std::string>
              paired_files{query_files[qfidx], query_files[qfidx+1]};
            match_query_against_reference
                <GttlFastQGenerator<buf_size>,
                 fastq_paired,
                 RefHashedQgrams,
                 store_seeds,
                 use_diagonal_encoder>
                (out_fp,
                 paired_files,
                 var_ref_idx,
                 rhqg,
                 guide_params,
                 store_header,
                 max_num_sequences,
                 padding_char,
                 grant_owner_ship);
          } else
          {
            constexpr const bool fastq_paired = false;
            const std::vector<std::string> query_inputfiles{query_files[qfidx]};
            (gttl_likely_fasta_format(query_files[qfidx])
              ? match_query_against_reference
                  <GttlFastAGenerator<buf_size>,
                   fastq_paired,
                   RefHashedQgrams,
                   store_seeds,
                   use_diagonal_encoder>
              : match_query_against_reference
                  <GttlFastQGenerator<buf_size>,
                   fastq_paired,
                   RefHashedQgrams,
                   store_seeds,
                   use_diagonal_encoder>)
                (out_fp,
                 query_inputfiles,
                 var_ref_idx,
                 rhqg,
                 guide_params,
                 store_header,
                 max_num_sequences,
                 padding_char,
                 grant_owner_ship);
          }
        }
      } else /* guide_params.number_of_threads_get() > 1 */
      {
        ThreadsOutputFiles threads_output_files("nt_mini_matcher",
                                                threads_out_prefix.empty()
                                                  ? nullptr
                                                  : threads_out_prefix.c_str(),
                                                guide_params
                                                  .number_of_threads_get());
        if (guide_params.verbose_option_get())
        {
          if (threads_out_prefix.empty())
          {
            /* then finally the output will be shown on stdout, so
               we output it here */
            log_vector_show(stdout,log_vector);
          } else
          {
            for (size_t thd_idx = 0;
                 thd_idx < guide_params.number_of_threads_get();
                 thd_idx++)
            {
              log_vector_show(threads_output_files.filepointer(thd_idx),
                              log_vector);
            }
          }
          log_vector->clear();
        }
        const size_t units_to_process
          = read_pairs_option ? (query_files.size()/2) : query_files.size();
        if (guide_params.query_split_size_get() > 0 && units_to_process == 1)
        {
          assert (query_files.size() == 1 or query_files.size() == 2);
#define WITH_MULTISEQ_GENERATOR
#ifdef WITH_MULTISEQ_GENERATOR
          const size_t max_size_of_queue
            = guide_params.thread_pool_max_size_get();
          const size_t max_num_sequences = guide_params.query_split_size_get();
          (query_files.size() == 2
            ? threaded_match_query_against_reference
                <GttlFastQGenerator<buf_size>,
                 true,
                 RefHashedQgrams,
                 store_seeds,
                 use_diagonal_encoder>
            : threaded_match_query_against_reference
                <GttlFastAGenerator<buf_size>,
                 false,
                 RefHashedQgrams,
                 store_seeds,
                 use_diagonal_encoder>)
                (&threads_output_files,
                 max_size_of_queue,
                 query_files,
                 var_ref_idx,
                 rhqg,
                 guide_params,
                 store_header,
                 max_num_sequences,
                 padding_char);
#else
          /* if we have only one unit (i.e. on file or a pair of files)
             and we have to split the queries, we create a MultiseqFactory,
             which allows to split the input into multiple queries */
          constexpr const uint8_t query_padding_char = UINT8_MAX-1;
          const bool short_header = store_header;
          GttlMultiseqFactory *multiseq_factory_ptr
            = query_files.size() == 2
                ? new GttlMultiseqFactory (query_files[0],
                                           query_files[1],
                                           0,
                                           0,
                                           guide_params.query_split_size_get(),
                                           query_padding_char,
                                           store_header,
                                           short_header)
                : new GttlMultiseqFactory (query_files[0],
                                           0,
                                           0,
                                           guide_params.query_split_size_get(),
                                           query_padding_char,
                                           store_header,
                                           short_header);
          /* Thread the function given as third argument using the
             Parameters which follow */
          GttlThreadPoolVar(guide_params.number_of_threads_get(),
                            multiseq_factory_ptr->size(),
                            process_single_query_unit_threaded
                              <RefHashedQgrams,
                               store_seeds,
                               use_diagonal_encoder>,
                            threads_output_files.filepointers_vector_get(),
                            rhqg,
                            nullptr, /* =query_files, as these are in factory*/
                            multiseq_factory_ptr,
                            guide_params,
                            read_pairs_option,
                            var_ref_idx);
          delete multiseq_factory_ptr;
#endif
        } else /* guide_params.query_split_size_get() == 0 or
                  units_to_process > 1 */
        {
          if (guide_params.number_of_threads_get() > units_to_process)
          {
            throw std::runtime_error(
                    std::format(": the number of threads ({}) must not be "
                                "larger than the number of {}query files "
                                "({}) to process",
                                guide_params.number_of_threads_get(),
                                read_pairs_option ? "pairs of " : "",
                                units_to_process));
          }
          /* Thread the function given as third argument using the
             Parameters which follow */
          GttlThreadPoolVar(guide_params.number_of_threads_get(),
                            units_to_process,
                            process_single_query_unit_threaded
                               <RefHashedQgrams,
                                store_seeds,
                                use_diagonal_encoder>,
                            threads_output_files.filepointers_vector_get(),
                            rhqg,
                            &query_files,
                            nullptr, /* multiseq_factory_ptr */
                            guide_params,
                            read_pairs_option,
                            var_ref_idx);
        }
      }
    }
  });
}

template<bool with_nt_hash,bool store_seeds,bool use_diagonal_encoder>
static void iterate_over_queryfiles(bool var_ref_rc_opt,
                                    int var_ref_sizeof_unit_hashed_qgram,
                                    int var_ref_idx,
                                    GttlMultiseq *ref_multiseq,
                                    const GuideParams &guide_params,
                                    bool read_pairs_option,
                                    const std::vector<std::string> &query_files,
                                    std::vector<std::string> *log_vector,
                                    const std::string &threads_out_prefix,
                                    const std::string &index_cache_path,
                                    bool build_index_only)
{
  (var_ref_rc_opt
     ? iterate_over_queryfiles_generic<true,
                                       with_nt_hash,
                                       store_seeds,
                                       use_diagonal_encoder>
     : iterate_over_queryfiles_generic<false,
                                       with_nt_hash,
                                       store_seeds,
                                       use_diagonal_encoder>)
                                      (var_ref_sizeof_unit_hashed_qgram,
                                       var_ref_idx,
                                       ref_multiseq,
                                       guide_params,
                                       read_pairs_option,
                                       query_files,
                                       log_vector,
                                       threads_out_prefix,
                                       index_cache_path,
                                       build_index_only);
}

#define ALL_VS_ALL_ENUMPOSITIONPAIRMATCHES(SIZE_UNIT_I,I_IDX,SIZE_UNIT_J,J_IDX)\
        if (guide_params.verbose_option_get()) \
        {\
          fprintf(out_fp,"# compute matches of minimum length %zu from "\
                         "minimizer based seeds for k=%zu and w=%zu between "\
                         "sequences in files \t%s\t%s\n",\
                         guide_params.minimum_mem_length_get(),\
                         guide_params.qgram_length_get(),\
                         guide_params.ref_window_size_get(),\
                         input_files[I_IDX].c_str(),\
                         input_files[J_IDX].c_str());\
        }\
        using HashedQgrams_j = HashedQgrams<SIZE_UNIT_J,\
                                            ref_rc_opt,\
                                            with_nt_hash>;\
        constexpr const bool query_rc_opt = false;\
        using HashedQgrams_i = HashedQgrams<SIZE_UNIT_I,\
                                            query_rc_opt,\
                                            with_nt_hash>;\
        const HashedQgrams_j &j_mi = get##SIZE_UNIT_J(J_IDX);\
        assert(i_mi.multiseq.padding_char_get() != \
               j_mi.multiseq.padding_char_get());\
        RunTimeClass rt_merge_extension{};\
        constexpr const SelfMatchMode self_match_local = None;\
        EnumPositionPairMatches<self_match_local,\
                                HashedQgrams_i,\
                                HashedQgrams_j,\
                                store_seeds,\
                                use_diagonal_encoder,\
                                8> enum_pos_pair_enumerator(i_mi,j_mi,out_fp);\
        if constexpr (store_seeds)\
        {\
          if (guide_params.verbose_option_get()) \
          {\
            const double space \
              = mega_bytes(enum_pos_pair_enumerator.size_in_RAM());\
            fprintf(out_fp,"# SPACE\tseeds.all_vs_all (MB):\t%.0f\n",space);\
            fprintf(out_fp,"# %s\n",rt_merge_extension\
                                    .to_string("generation of seeds").c_str());\
          }\
          rt_merge_extension.reset();\
        }\
        extend_and_process_matches<self_match_local,\
                                   HashedQgrams_i,\
                                   HashedQgrams_j,\
                                   store_seeds,\
                                   use_diagonal_encoder,\
                                   sizeof_unit_pospair,\
                                   ref_idx>\
                                  (out_fp,\
                                   guide_params,\
                                   enum_pos_pair_enumerator,\
                                   i_mi.multiseq,\
                                   j_mi.multiseq);\
        if (guide_params.verbose_option_get()) \
        {\
          fprintf(out_fp,"# %s\n",rt_merge_extension.to_string(\
                                  store_seeds \
                                    ? "extension and output of matches"\
                                    : ("generation of seeds, "\
                                       "extension and output of "\
                                       "matches")).c_str());\
        }

template<bool ref_rc_opt,bool with_nt_hash>
class HashedQgramsTable
{
  private:
  static constexpr const bool at_constant_distance = false;
  using HashedQgrams8 = HashedQgrams<8,ref_rc_opt,with_nt_hash>;
  using HashedQgrams9 = HashedQgrams<9,ref_rc_opt,with_nt_hash>;
  std::vector<GttlMultiseq *> multiseq_vec;
  std::vector<std::pair<bool,size_t>> idx8_or_9;
  std::vector<HashedQgrams8> mi8_vec;
  std::vector<HashedQgrams9> mi9_vec;
  bool is_eight(size_t idx) const noexcept
  {
    return std::get<0>(idx8_or_9[idx]);
  }
  const HashedQgrams8 &get8(size_t idx) const noexcept
  {
    auto result = idx8_or_9[idx];
    assert(std::get<0>(result));
    return mi8_vec[std::get<1>(result)];
  }
  const HashedQgrams9 &get9(size_t idx) const noexcept
  {
    auto result = idx8_or_9[idx];
    assert(!std::get<0>(result));
    return mi9_vec[std::get<1>(result)];
  }
  public:
  HashedQgramsTable(const std::vector<std::string> &input_files,
                    GuideParams *guide_params)
  {
    const bool store_header
      = guide_params->display_options_get().s_seqid_display();
    uint8_t padding_char = UINT8_MAX;
    size_t files_with_wildcards = 0;
    constexpr const bool with_reverse_complement = false;
    int remaining_hash_bits = INT_MAX;

    FILE *out_fp = stdout;
    for (auto &&input_file : input_files)
    {
      if (guide_params->verbose_option_get())
      {
        fprintf(out_fp,"# create multiseq for sequences in file\t%s\n",
                input_file.c_str());
      }
      GttlMultiseq *multiseq_ptr = nullptr;
      try
      {
        multiseq_ptr = new GttlMultiseq(input_file, /* CONSTRUCTOR */
                                        store_header,
                                        store_sequence,
                                        padding_char,
                                        with_reverse_complement);
        if (guide_params->display_options_get().s_seqid_display())
        {
          multiseq_ptr->short_header_cache_create<'|','|'>();
        }
      }
      catch (const std::exception &err)
      {
        for (auto &&ptr : multiseq_vec)
        {
          delete ptr;
        }
        throw err;
      }
      if (guide_params->hash_bits_get() == -1)
      {
        if (multiseq_ptr->sequences_bits_get() >= 64)
        {
          remaining_hash_bits = 0;
        } else
        {
          remaining_hash_bits
            = std::min(remaining_hash_bits,
                       64 - multiseq_ptr->sequences_bits_get());
        }
      }
      multiseq_vec.push_back(multiseq_ptr);
      padding_char--;
    }
    assert(guide_params->hash_bits_get() != -1 or
           remaining_hash_bits != INT_MAX);
    guide_params->hash_bits_set_if_undefined(
      HashedQgrams8::possible_false_positive_matches,remaining_hash_bits);
    for (size_t file_idx = 0; file_idx < input_files.size(); file_idx++)
    {
      GttlMultiseq *multiseq_ptr = multiseq_vec[file_idx];
      if (guide_params->verbose_option_get())
      {
        fprintf(out_fp,"# extract minimimizers from sequences in file\t%s\n",
                input_files[file_idx].c_str());
      }
      check_hash_bits(guide_params->verbose_option_get() ? out_fp : nullptr,
                      "all_vs_all",
                      guide_params->hash_bits_get(),
                      multiseq_ptr->sequences_bits_get(),
                      HashedQgrams8::possible_false_positive_matches);
      constexpr const bool multiseq_rc_opt = false;
      if (guide_params->hash_bits_get() + multiseq_ptr->sequences_bits_get()
            <= 64)
      {
        idx8_or_9.push_back({true,mi8_vec.size()});
        mi8_vec.push_back(HashedQgrams<8,multiseq_rc_opt,with_nt_hash>
                                      (*multiseq_ptr,
                                       guide_params->number_of_threads_get(),
                                       guide_params->qgram_length_get(),
                                       guide_params->ref_window_size_get(),
                                       guide_params->hash_bits_get(),
                                       true,
                                       guide_params->at_constant_distance_get(),
                                       guide_params->max_replicate_ref_get(),
                                       nullptr));
        if (mi8_vec.back().sequence_has_wildcards())
        {
          files_with_wildcards++;
        }
      } else
      {
        constexpr const int max_hashed_qgram_bits = 72;
        if (guide_params->hash_bits_get() + multiseq_ptr->sequences_bits_get()
              <= max_hashed_qgram_bits)
        {
          idx8_or_9.push_back({false,mi9_vec.size()});
          mi9_vec.push_back(HashedQgrams<9,multiseq_rc_opt,with_nt_hash>
                                        (*multiseq_ptr,
                                         guide_params->number_of_threads_get(),
                                         guide_params->qgram_length_get(),
                                         guide_params->ref_window_size_get(),
                                         guide_params->hash_bits_get(),
                                         true,
                                         guide_params
                                           ->at_constant_distance_get(),
                                         guide_params->max_replicate_ref_get(),
                                         nullptr));
          if (mi9_vec.back().sequence_has_wildcards())
          {
            files_with_wildcards++;
          }
        } else
        {
          throw std::runtime_error(
                  std::format("number of hash bits ({}) + number of sequences "
                              "bits ({}) for file {} is too large",
                              guide_params->hash_bits_get(),
                              multiseq_ptr->sequences_bits_get(),
                              input_files[file_idx].c_str()));
        }
      }
    }
    const uint8_t smallest_padding_char = padding_char + 1;
    if (files_with_wildcards >= 2)
    {
      std::array<bool,UINT8_MAX+1> valid_char;
      valid_char.fill(false);
      constexpr const size_t spec_size = sizeof nucleotides_upper_lower - 1;
      for (size_t idx = 0; idx < spec_size; idx++)
      {
        if (nucleotides_upper_lower[idx] != '|')
        {
          valid_char[static_cast<int>(nucleotides_upper_lower[idx])] = true;
        }
      }
      std::vector<uint8_t> invalid_reservoir;
      for (uint8_t cidx = smallest_undef_char; cidx < smallest_padding_char;
           cidx++)
      {
        if (!valid_char[cidx])
        {
          invalid_reservoir.push_back(cidx);
        }
      }
      for (size_t idx = 0; idx < multiseq_vec.size(); idx++)
      {
        if (is_eight(idx))
        {
          if (get8(idx).sequence_has_wildcards())
          {
            const uint8_t current_undef_char = invalid_reservoir.back();
            invalid_reservoir.pop_back();
            wildcard_replacement(multiseq_vec[idx],current_undef_char);
          }
        } else
        {
          if (get9(idx).sequence_has_wildcards())
          {
            assert(invalid_reservoir.size() > 0);
            const uint8_t current_undef_char = invalid_reservoir.back();
            invalid_reservoir.pop_back();
            wildcard_replacement(multiseq_vec[idx],current_undef_char);
          }
        }
      }
    }
  }
  template<int sizeof_unit_pospair,int ref_idx,
           bool store_seeds,bool use_diagonal_encoder>
  void all_vs_all(FILE *out_fp,
                  const std::vector<std::string> &input_files,
                  const GuideParams &guide_params)
  {
    /* This function currently computes the hashed_qgrams
       for each sequence (possibly in multiple threads)
       and performs a comparison for each pair (using only one thread). */
    assert(multiseq_vec.size() > 0);
    for (size_t i = 0; i < multiseq_vec.size() - 1; i++)
    {
      if (is_eight(i))
      {
        const HashedQgrams8 &i_mi = get8(i);
        for (size_t j = i+1; j < multiseq_vec.size(); j++)
        {
          if (is_eight(j))
          {
            ALL_VS_ALL_ENUMPOSITIONPAIRMATCHES(8,i,8,j);
          } else
          {
            ALL_VS_ALL_ENUMPOSITIONPAIRMATCHES(8,i,9,j);
          }
        }
      } else
      {
        const HashedQgrams9 &i_mi = get9(i);
        for (size_t j = i+1; j < multiseq_vec.size(); j++)
        {
          if (is_eight(j))
          {
            ALL_VS_ALL_ENUMPOSITIONPAIRMATCHES(9,i,8,j);
          } else
          {
            ALL_VS_ALL_ENUMPOSITIONPAIRMATCHES(9,i,9,j);
          }
        }
      }
    }
  }
  ~HashedQgramsTable(void)
  {
    for (auto &&ptr : multiseq_vec)
    {
      delete ptr;
    }
  }
};

template<bool rc_opt,
         bool with_nt_hash,
         bool store_seeds,
         bool use_diagonal_encoder,
         int ref_sizeof_unit_hashed_qgram>
static void hashed_qgrams_seeded_self_matches_generic(
                                              FILE *out_fp,
                                              const GuideParams &guide_params,
                                              const GttlMultiseq &ref_multiseq)
{
  constexpr const int sizeof_unit_pospair = 8;
  constexpr const size_t self_match_ref_idx = 0;
  using RefHashedQgrams = HashedQgrams<ref_sizeof_unit_hashed_qgram,
                                       rc_opt,
                                       with_nt_hash>;
  std::vector<std::string> log_vector;
  constexpr const bool sort_by_hashvalue = true;
  RefHashedQgrams rhqg(ref_multiseq,
                       guide_params.number_of_threads_get(),
                       guide_params.qgram_length_get(),
                       guide_params.ref_window_size_get(),
                       guide_params.hash_bits_get(),
                       sort_by_hashvalue,
                       guide_params.at_constant_distance_get(),
                       guide_params.max_replicate_ref_get(),
                       guide_params.verbose_option_get() ? &log_vector
                                                         : nullptr);
  if (guide_params.verbose_option_get())
  {
    log_vector_show(out_fp,&log_vector);
    log_vector.clear();
  }
  if (guide_params.display_options_get().minimizers_display())
  {
    rhqg.show();
  }
  hashed_qgrams_seeded_reference_vs_query_matches_generic
    <rc_opt ? WithReverseComplement : Regular,
     RefHashedQgrams,
     RefHashedQgrams,
     store_seeds,
     use_diagonal_encoder,
     sizeof_unit_pospair,
     self_match_ref_idx>
    (out_fp,rhqg,rhqg,guide_params);
}

template<bool with_nt_hash,
         bool store_seeds,
         bool use_diagonal_encoder>
static void hashed_qgrams_seeded_self_matches(FILE *out_fp,
                                              bool rc_opt,
                                              const GuideParams &guide_params,
                                              const GttlMultiseq &ref_multiseq,
                                              const char *reference_file)

{
  const int bits_for_hashed_sequence
    = guide_params.hash_bits_get() + ref_multiseq.sequences_bits_get();
  if (bits_for_hashed_sequence > 72)
  {
    throw std::runtime_error(
            std::format("the number of hash bits ({}) + the number of bits for "
                        "the sequence coordinates ({}) for sequences in file "
                        "{} is larger than the maximum value of 72",
                        guide_params.hash_bits_get(),
                        ref_multiseq.sequences_bits_get(),
                        reference_file));
  }
  const int bytes_for_hashed_sequence
    = bits_for_hashed_sequence <= 64 ? 8 : 9;
  constexpr_for<8,9+1,1>([&](auto ct_bytes_for_hashed_sequence)
  {
    if (bytes_for_hashed_sequence == ct_bytes_for_hashed_sequence)
    {
      const int rc_opt_int = rc_opt ? 1 : 0;
      constexpr_for<0,1+1,1>([&](auto ct_rc_opt_int)
      {
        if (rc_opt_int == ct_rc_opt_int)
        {
          constexpr const bool ct_rc_opt = ct_rc_opt_int == 1 ? true : false;
          hashed_qgrams_seeded_self_matches_generic
            <ct_rc_opt,
             with_nt_hash,
             store_seeds,
             use_diagonal_encoder,
             ct_bytes_for_hashed_sequence>
            (out_fp,guide_params,ref_multiseq);
        }
      });
    }
  });
}

template<bool with_nt_hash,bool store_seeds,bool use_diagonal_encoder>
static void compare_all_vs_all(FILE *out_fp,
                               const std::vector<std::string> &input_files,
                               GuideParams *guide_params)
{
  constexpr const bool ref_rc_opt = false;
  HashedQgramsTable<ref_rc_opt,with_nt_hash> hqg_tab(input_files,
                                                     guide_params);
  constexpr const int ref_idx = 0;
  hqg_tab.template all_vs_all<8,ref_idx,store_seeds,use_diagonal_encoder>
                             (out_fp,
                              input_files,
                              *guide_params);
}

template <class OptionClass>
static GttlMultiseq *reference_multiseq_get(std::vector<std::string>
                                              *log_vector,
                                            const std::string &reference_file,
                                            const OptionClass &options)
{
  const bool store_header = options.display_options_get().s_seqid_display();
  RunTimeClass rt_create_ref_multiseq{};
  const bool with_reverse_complement
    = options.reverse_complement_for_get().is(for_reference_canonical) ||
      options.reverse_complement_for_get().is(for_reference);
  constexpr const uint8_t padding_char = UINT8_MAX;
  GttlMultiseq *ref_multiseq
    = new GttlMultiseq(reference_file, /* CONSTRUCTOR */
                       store_header,
                       store_sequence,
                       padding_char,
                       with_reverse_complement);
  if (options.display_options_get().s_seqid_display())
  {
    ref_multiseq->short_header_cache_create<'|','|'>();
  }
  if (log_vector != nullptr)
  {
    std::string log_msg = std::string("create reference multiseq");
    if (with_reverse_complement)
    {
      log_msg += std::string(" including reverse complement sequences");
    }
    log_vector->push_back(rt_create_ref_multiseq.to_string(log_msg));
    const double space = std::ceil(mega_bytes(ref_multiseq->size_in_bytes()));
    log_vector->push_back(std::format("SPACE\tref_multiseq (MB)\t{:.0f}",
                                      space));
  }
  return ref_multiseq;
}

int main(int argc, char *argv[])
{
  RunTimeClass rt_total{};
#ifdef GI_JOINER_CODE
  Options options{true};
#else
  Options options{false};
#endif
  static constexpr const bool replace_single_hyphen_options = false;
  ArgvConcat my_argv(argc,argv," ",replace_single_hyphen_options,
                     std::vector<std::string>{"-d","--display"});
  try
  {
    options.parse(static_cast<int>(my_argv.size()),my_argv.concat_get());
  }
  catch (const std::invalid_argument &err) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << err.what() << '\n';
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  if (options.display_options_get().strand_display() &&
      options.display_options_get().q_read_pairs_display())
  {
    std::cerr << argv[0] << (": display option 'strand' and 'q.read_pairs' "
                             "cannot be combined") << '\n';
    return EXIT_FAILURE;
  }
  const std::string option_string = gttl_argv2string("Options:",argc,argv);
  std::vector<std::string> log_vector;
  log_vector.push_back(option_string.c_str());
  log_vector.push_back(options.display_options_get().to_string().c_str());
  const char *progname = argv[0];
  const std::vector<std::string> &input_files = options.input_files_get();
  const std::string &reference_file = input_files[0];
  bool has_err = false;
  for (auto &&input_file : input_files)
  {
    if (options.read_pairs_option_is_set())
    {
      GttlFpType in_fp = gttl_fp_type_open(input_file.c_str(),"rb");

      if (in_fp == nullptr)
      {
        std::cerr << progname << ": cnnot open file: " << input_file << '\n';
        has_err = true;
      }
      gttl_fp_type_close(in_fp);
    } else
    {
      has_err = check_false_sequence_type(progname,input_file.c_str());
    }
    if (has_err)
    {
      break;
    }
  }
  if (!has_err)
  {
    try
    {
      gttl_duplicated_filenames(input_files);
    }
    catch (const std::exception &err) /* check_err.py */
    {
      std::cerr << progname << ": " << err.what() << '\n';
      has_err = true;
    }
  }
  if (!has_err && options.all_vs_all_option_is_set() &&
      (input_files.size() == 2 || input_files.size() == 1))
  {
    std::cerr << progname /* check_err.py */
              << std::format(": option -a/--all_vs_all can only be used "
                             "if at least three input files are specified; "
                             "you likely want to omit option -a to compare "
                             "the given reference {} file against {}",
                             input_files[0].c_str(),
                             input_files.size() == 2 ? input_files[1].c_str()
                                                     : "itself")
              << '\n';
    has_err = true;
  }
  if (!has_err)
  {
    GuideParams guide_params(options.kmer_size_get(),
                             options.hash_bits_get(),
#ifdef GI_JOINER_CODE
                             options.window_size_get(),
                             options.window_size_get(),
                             options.number_of_threads_get(),
                             options.query_split_size_get(),
#else
                             options.minimum_mem_length_get(),
                             options.minimum_als_length_get(),
                             options.number_of_threads_get(),
                             options.thread_pool_max_size_get(),
                             options.query_split_size_get(),
                             options.max_replicate_ref_get(),
                             options.max_replicate_qry_get(),
                             options.maximum_error_percentage_get(),
                             options.polishing_error_percentage_get(),
                             options.cop_option_is_set(),
                             options.final_polishing_option_is_set(),
                             options.group_by_query_option_is_set(),
                             options.verbose_option_is_set(),
#endif
                             options.display_options_get());

    if (input_files.size() == 1)
    {
      FILE *out_fp = stdout;
      if (options.build_index_only_is_set())
      {
        try
        {
          GttlMultiseq *ref_multiseq
            = reference_multiseq_get<Options>(options.verbose_option_is_set()
                                                ? &log_vector : nullptr,
                                              reference_file,
                                              options);
          guide_params.hash_bits_set_if_undefined(
             not options.invint_option_is_set(),
             ref_multiseq->sequences_bits_get());
          check_hash_bits_log(options.verbose_option_is_set()
                                ? &log_vector : nullptr,
                              "ref",
                              guide_params.hash_bits_get(),
                              ref_multiseq->sequences_bits_get(),
                              not options.invint_option_is_set());
          const int var_ref_sizeof_unit_hashed_qgram
            = (guide_params.hash_bits_get() +
               ref_multiseq->sequences_bits_get() <= 64) ? 8 : 9;
          const int var_ref_idx = 0;
          const bool read_pairs_option = false;
          const std::vector<std::string> query_files{};
          INSTANTIATE_NT_HASH_STORE_SEEDS_USE_DIAGONAL_ENCODER(
            iterate_over_queryfiles,
            options.reverse_complement_for_get().is(for_reference_canonical),
            var_ref_sizeof_unit_hashed_qgram,
            var_ref_idx,
            ref_multiseq,
            guide_params,
            read_pairs_option,
            query_files,
            options.verbose_option_is_set() ? &log_vector : nullptr,
            options.threads_out_prefix_get(),
            options.index_cache_path_get(),
            options.build_index_only_is_set());
          delete ref_multiseq;
          if (options.verbose_option_is_set())
          {
            log_vector_show(stdout,&log_vector);
            log_vector.clear();
          }
        }
        catch (const std::exception &err) /* check_err.py */
        {
          std::cerr << progname << ": file \"" << reference_file << "\" "
                    << err.what() << '\n';
          has_err = true;
        }
      } else
      {
        assert(options.number_of_threads_get() == 1);
        try
        {
          GttlMultiseq *ref_multiseq
            = reference_multiseq_get<Options>(options.verbose_option_is_set()
                                                ? &log_vector : nullptr,
                                              reference_file,
                                              options);
          if (options.display_options_get().s_seqid_display())
          {
            ref_multiseq->short_header_cache_create<'|','|'>();
          }
          assert (2 * ref_multiseq->sequences_bits_get() <= 64);
          guide_params.hash_bits_set_if_undefined(
             not options.invint_option_is_set(),
             ref_multiseq->sequences_bits_get());
          check_hash_bits(options.verbose_option_is_set() ? out_fp : nullptr,
                          "ref",
                          guide_params.hash_bits_get(),
                          ref_multiseq->sequences_bits_get(),
                          not options.invint_option_is_set());
          INSTANTIATE_NT_HASH_STORE_SEEDS_USE_DIAGONAL_ENCODER(
            hashed_qgrams_seeded_self_matches,
            out_fp,
            options.reverse_complement_for_get().is(for_reference),
            guide_params,
            *ref_multiseq,
            reference_file.c_str());
        }
        catch (const std::exception &err) /* check_err.py */
        {
          std::cerr << progname << ": file \"" << reference_file << "\" "
                    << err.what() << '\n';
          has_err = true;
        }
      }
      if (!has_err and options.verbose_option_is_set())
      {
        log_vector_show(out_fp,&log_vector);
        log_vector.clear();
      }
    } else
    {
      assert(input_files.size() > 1);
      if (input_files.size() > 2 && options.all_vs_all_option_is_set())
      {
        assert(options.number_of_threads_get() == 1);
        try
        {
          FILE *out_fp = stdout;
          INSTANTIATE_NT_HASH_STORE_SEEDS_USE_DIAGONAL_ENCODER(
            compare_all_vs_all,
            out_fp,
            input_files,
            &guide_params);
        }
        catch (const std::exception &err) /* check_err.py */
        {
          std::cerr << progname << ": " << err.what() << '\n';
          has_err = true;
        }
      } else
      {
        const bool read_pairs_option = options.read_pairs_option_is_set();
        auto query_files = std::vector<std::string>(input_files.begin() + 1,
                                                    input_files.end());
        if (read_pairs_option && query_files.size() % 2 != 0)
        {
          std::cerr << progname
                    << std::string(": if option -r/--read_pairs is used, "
                                   "there must be an even number of query "
                                   "files")
                    << '\n';
          has_err = true;
        } else
        {
          GttlMultiseq *ref_multiseq = nullptr;
          try
          {
            ref_multiseq
              = reference_multiseq_get<Options>(options.verbose_option_is_set()
                                                 ? &log_vector : nullptr,
                                                reference_file,
                                                options);
            guide_params.hash_bits_set_if_undefined(
              not options.invint_option_is_set(),
              ref_multiseq->sequences_bits_get());
            check_hash_bits_log(options.verbose_option_is_set()
                                  ? &log_vector : nullptr,
                                "ref",
                                guide_params.hash_bits_get(),
                                ref_multiseq->sequences_bits_get(),
                                not options.invint_option_is_set());
            const int var_ref_sizeof_unit_hashed_qgram
              = (guide_params.hash_bits_get() +
                 ref_multiseq->sequences_bits_get() <= 64) ? 8 : 9;
            const int var_ref_idx = guide_params.group_by_query_get() ? 1 : 0;

            INSTANTIATE_NT_HASH_STORE_SEEDS_USE_DIAGONAL_ENCODER(
              iterate_over_queryfiles,
              options.reverse_complement_for_get().is(for_reference_canonical),
              var_ref_sizeof_unit_hashed_qgram,
              var_ref_idx,
              ref_multiseq,
              guide_params,
              read_pairs_option,
              query_files,
              options.verbose_option_is_set() ? &log_vector : nullptr,
              options.threads_out_prefix_get(),
              options.index_cache_path_get(),
              options.build_index_only_is_set());
          }
          catch (const std::exception &err) /* check_err.py */
          {
            std::cerr << progname << ": file \"" << reference_file << "\""
                      << err.what() << '\n';
            has_err = true;
          }
          delete ref_multiseq;
        }
      }
    }
  }
  if (!has_err)
  {
    rt_total.show("total");
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
