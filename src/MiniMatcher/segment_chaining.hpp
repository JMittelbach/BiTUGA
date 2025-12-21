#ifndef SEGMENT_CHAINING_HPP
#define SEGMENT_CHAINING_HPP
/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

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
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <algorithm>
#include <cstdio>
#include <array>
#include <cmath>
#include <type_traits>
#include <iomanip>
#include <ios>

#include "utilities/mathsupport.hpp"
#include "sequences/gttl_multiseq.hpp"
#ifndef NDEBUG
#include "sequences/outsenseedist_inplace.hpp"
#endif
#include "sequences/matching_characters.hpp"
#include "sequences/local_chainer.hpp"
#include "sequences/close_gaps.hpp"
#include "sequences/pol_pref_and_suff_ext.hpp"
#include "sequences/non_redundant_matches.hpp"

#include "display_options.hpp"
#include "guide_params.hpp"
#include "segmentation.hpp"
#include "sequence_bias_get.hpp"

template<class SeedTable,class MatchClass>
static void show_chain_elements(const std::vector<uint32_t> &chain_elements,
                                const SeedTable &seed_table,
                                size_t segment_start,
                                const GttlMultiseq &ref_multiseq,
                                const GttlMultiseq &query_multiseq,
                                const DisplayOptions &display_options)
{
  const size_t query_seqnum_offset
    = query_multiseq.sequence_number_offset_get();
  for (size_t idx = chain_elements.size(); idx > 0; idx--)
  {
    const size_t mtidx = segment_start + chain_elements[idx-1];
    assert(mtidx < seed_table.size());
    MatchClass this_match(seed_table.ref_endpos_get(mtidx),
                          seed_table.length_get(mtidx),
                          seed_table.query_endpos_get(mtidx),
                          seed_table.length_get(mtidx),
                          0);
    if (!display_options.swallow_display())
    {
      this_match.display(stdout,
                         display_options,
                         ref_multiseq,
                         seed_table.ref_seqnum_get(mtidx),
                         query_multiseq,
                         seed_table.query_seqnum_get(mtidx),
                         query_seqnum_offset,
                         nullptr);
    }
  }
}

template<class MatchClass>
class SegmentMatches
{
  struct CompleteMatch
  {
    size_t s_seqnum, q_seqnum;
    const MatchClass match;
    size_t primary_startpos;
    double weight;
    private:
    size_t primary_startpos_evaluate(const GttlMultiseq &ref_multiseq,
                                     const GttlMultiseq &query_multiseq)
    {
      const bool forward_strand
        = not ref_multiseq.has_reverse_complement_is_set()
          or (s_seqnum % 2 == 0);
      if (forward_strand)
      {
        return match.q_startpos_get();
      }
      const size_t query_seqlen = query_multiseq.sequence_length_get(q_seqnum);
      return query_seqlen - (match.q_endpos_get() + 1);
    }
    public:
    CompleteMatch(const GttlMultiseq &ref_multiseq,
                  const GttlMultiseq &query_multiseq,
                  size_t _s_seqnum, size_t _q_seqnum,const MatchClass &_match)
      : s_seqnum(_s_seqnum)
      , q_seqnum(_q_seqnum)
      , match(_match)
      , primary_startpos(primary_startpos_evaluate(ref_multiseq,query_multiseq))
      , weight(0.5 * match.aligned_len_get() *
               (100.0 - match.match_error_percentage_get()))
    {}
    size_t primary_startpos_get(void) const noexcept
    {
      return primary_startpos;
    }
    size_t primary_len_get(void) const noexcept
    {
      return match.q_len_get();
    }
    double weight_get(void) const noexcept
    {
      return weight;
    }
    bool superior_weight (const CompleteMatch& other) const noexcept
    {
      return this->weight_get() > other.weight_get();
    }
    bool superior_weight_tie_primary_startpos(const CompleteMatch& other)
      const noexcept
    {
      return this->superior_weight(other) or
             (this->weight_get() == other.weight_get() and
              this->primary_startpos_get() > other.primary_startpos_get());
    }
    std::string to_string(void) const noexcept
    {
      return std::to_string(primary_startpos_get()) + "\t" +
             std::to_string(primary_len_get()) + "\t" +
             std::to_string(weight_get());
    }
  };
  std::vector<CompleteMatch> store;
  const DisplayOptions &display_options;
  const GttlMultiseq &ref_multiseq;
  const GttlMultiseq &query_multiseq;
  const int shift_if_complement_for_reference,
            shift_if_query_readpairs;
  FILE *out_fp;
  size_t real_s_seqnum(size_t seqnum) const noexcept
  {
    return seqnum >> shift_if_complement_for_reference;
  }
  size_t real_q_seqnum(size_t seqnum) const noexcept
  {
    return seqnum >> shift_if_query_readpairs;
  }
  public:
  SegmentMatches(const GuideParams &_guide_params,
                 const GttlMultiseq &_ref_multiseq,
                 const GttlMultiseq &_query_multiseq,
                 FILE *_out_fp)
    : store({})
    , display_options(_guide_params.display_options_get())
    , ref_multiseq(_ref_multiseq)
    , query_multiseq(_query_multiseq)
    , shift_if_complement_for_reference(_ref_multiseq.
                                        has_reverse_complement_is_set() ? 1 : 0)
    , shift_if_query_readpairs(_query_multiseq.has_read_pairs_is_set() ? 1 : 0)
    , out_fp(_out_fp)
  {}
  void display(void) const noexcept
  {
    if (store.size() > 0)
    {
      if (display_options.non_redundant_display())
      {
        printf("# segment %lu %lu\n",store[0].s_seqnum,store[0].q_seqnum);
      }
      const size_t query_seqnum_offset
        = query_multiseq.sequence_number_offset_get();
      NonRedundantMatches<CompleteMatch> non_redundant_matches(store);
      for (size_t idx = 0; idx < store.size(); idx++)
      {
        if (display_options.non_redundant_display())
        {
          printf("%s\t",non_redundant_matches[idx] ? "++" : "--");
        }
        if (display_options.non_redundant_display() or
            non_redundant_matches[idx])
        {
          const CompleteMatch &cm = store[idx];
          if (cm.match.distance_get() > 0)
          {
            const char *ref_seq = ref_multiseq.sequence_ptr_get(cm.s_seqnum);
            const char *query_seq
              = query_multiseq.sequence_ptr_get(cm.q_seqnum);
            Eoplist eoplist = cm.match.match2eoplist(ref_seq,
                                                     cm.s_seqnum,
                                                     query_seq,
                                                     cm.q_seqnum);
            cm.match.display(out_fp,
                             display_options,
                             ref_multiseq,
                             cm.s_seqnum,
                             query_multiseq,
                             cm.q_seqnum,
                             query_seqnum_offset,
                             &eoplist);
          } else
          {
            cm.match.display(out_fp,
                             display_options,
                             ref_multiseq,
                             cm.s_seqnum,
                             query_multiseq,
                             cm.q_seqnum,
                             query_seqnum_offset,
                             nullptr);
          }
        }
      }
    }
  }
  void append(size_t s_seqnum,size_t q_seqnum, const MatchClass &match)
  {
    if (store.size() > 0)
    {
      const CompleteMatch &cm = store[store.size()-1];
      if (real_s_seqnum(cm.s_seqnum) != real_s_seqnum(s_seqnum) or
          real_q_seqnum(cm.q_seqnum) != real_q_seqnum(q_seqnum))
      {
        display();
        store.clear();
      }
    }
    store.emplace_back(CompleteMatch(ref_multiseq,query_multiseq,
                                     s_seqnum,q_seqnum,match));
  }
};

template<bool self_match,class SeedTable,bool (*match_method)(char,char),
         class MatchClass>
static void single_segment_chaining_template(
                const GuideParams &guide_params,
                const AlignmentPolishing &alignment_polishing,
                const GttlMultiseq &ref_multiseq,
                size_t ref_seqnum,
                const char *ref_seq,
                size_t ref_len,
                const GttlMultiseq &query_multiseq,
                size_t query_seqnum,
                const char *query_seq,
                size_t query_len,
                /* reference to table of matches seed_table must provide the
                   methods.  size()
                             ref_endpos_get(i), i is index in seed_table
                             query_endpos_get(i), is index in seed_table
                             length_get(i)
                             ref_seqnum_get(i)
                             query_seqnum_get(i)
                   Matches must be sorted in ascending order
                   of order_endpos_get() */
                const SeedTable &seed_table,
                size_t segment_start, /* index of first match in segment,
                                         refers to seed_table */
                size_t segment_length,/* number of matches in segment beginning
                                         with segment_start */
                SegmentMatches<MatchClass> *segment_matches,
                FILE *out_fp)
{
  //bool had_output_for_segment = false;
  if (seed_table.size() == 0)
  {
    return;
  }
  using GapType = uint32_t;
  using ScoreType = uint32_t;
  using PredecessorType = uint32_t;
  using ThisLocalChainer
    = LocalChainer<SeedTable,GapType,ScoreType,PredecessorType>;
  const DisplayOptions &display_options = guide_params.display_options_get();

  ThisLocalChainer local_chainer(seed_table,segment_start,segment_length,
                                 guide_params.max_previous_get());
  MatchClass extended_match{};
  const size_t query_seqnum_offset
    = query_multiseq.sequence_number_offset_get();
  for (auto &&this_chain : local_chainer) /* desc. score order */
  {
    if (this_chain.score_get() == 0)
    {
      break; /* as all following chains will have score <= 0 */
    }
    if constexpr (self_match)
    {
      if (ref_seqnum == query_seqnum && this_chain.self_overlap())
      {
        continue;
      }
    }
    const size_t distance
      = close_gaps_between_chain_elements<SeedTable,
                                          ThisLocalChainer,
                                          char,
                                          match_method>
                                         (seed_table,
                                          segment_start,
                                          local_chainer,
                                          ref_seq,
                                          query_seq,
                                          this_chain.from_element_get(),
                                          this_chain.size(),
                                          ref_seqnum,
                                          query_seqnum);
    MatchClass inner_match(this_chain.ref_endpos_get(),
                           this_chain.ref_match_length_get(),
                           this_chain.query_endpos_get(),
                           this_chain.query_match_length_get(),
                           distance);
    if (!display_options.swallow_display() &&
        (display_options.chain_elements_display() ||
         display_options.chain_closed_display()))
    {
      std::vector<uint32_t> chain_elements
        = local_chainer.local_chain_get(this_chain.from_element_get(),
                                        this_chain.size());
      assert (this_chain.size() == chain_elements.size() &&
              this_chain.size() > 0);
      if (this_chain.size() > 1 && display_options.chain_elements_display())
      {
        printf("chain\t%lu",this_chain.size());
        for (size_t idx = this_chain.size() - 1; idx >= 1; idx--)
        {
          auto ce = chain_elements[idx-1];
          printf("\t%u\t%u",
                 static_cast<unsigned int>(local_chainer
                                            .ref_gap_length_get(ce)),
                 static_cast<unsigned int>(local_chainer
                                            .query_gap_length_get(ce)));
        }
        printf("\t%u\n",this_chain.score_get());
      }
      if (display_options.chain_closed_display())
      {
        printf("# %s_match\t", this_chain.size() == 1 ? "singleton"
                                                      : "chain");
        inner_match.display(stdout,
                            display_options,
                            ref_multiseq,
                            ref_seqnum,
                            query_multiseq,
                            query_seqnum,
                            query_seqnum_offset,
                            nullptr);
      }
      if (this_chain.size() > 1 && display_options.chain_elements_display())
      {
        show_chain_elements<SeedTable,MatchClass>
                           (chain_elements,
                            seed_table,
                            segment_start,
                            ref_multiseq,
                            query_multiseq,
                            display_options);
      }
    }
    if (polished_prefix_and_suffix_extension<MatchClass>
                                            (&extended_match,
                                             alignment_polishing,
                                             guide_params
                                               .minimum_als_length_get(),
                                             guide_params
                                               .maximum_error_percentage_get(),
                                             ref_seq,
                                             ref_len,
                                             query_seq,
                                             query_len,
                                             ref_seqnum,
                                             query_seqnum,
                                             inner_match))
    {
      if (!display_options.swallow_display())
      {
        if (display_options.chain_elements_display() ||
            display_options.chain_closed_display())
        {
          printf("# ext_match\t");
        }
        if (extended_match.distance_get() > 0)
        {
          Eoplist eoplist = extended_match.match2eoplist(ref_seq,
                                                         ref_seqnum,
                                                         query_seq,
                                                         query_seqnum);
          extended_match.correct_coordinates(eoplist);
          if (guide_params.final_polishing_get())
          {
            segment_matches->append(ref_seqnum,
                                    query_seqnum,
                                    extended_match);
          } else
          {
            /*if (not had_output_for_segment)
            {
              printf("# segment %lu %lu\n",ref_seqnum,query_seqnum);
              had_output_for_segment = true;
            }*/
            extended_match.display(out_fp,
                                   display_options,
                                   ref_multiseq,
                                   ref_seqnum,
                                   query_multiseq,
                                   query_seqnum,
                                   query_seqnum_offset,
                                   &eoplist);
          }
        } else
        {
          if (guide_params.final_polishing_get())
          {
            segment_matches->append(ref_seqnum,
                                    query_seqnum,
                                    extended_match);
          } else
          {
            /*
            if (not had_output_for_segment)
            {
              printf("# segment %lu %lu\n",ref_seqnum,query_seqnum);
              had_output_for_segment = true;
            }
            */
            extended_match.display(out_fp,
                                   display_options,
                                   ref_multiseq,
                                   ref_seqnum,
                                   query_multiseq,
                                   query_seqnum,
                                   query_seqnum_offset,
                                   nullptr);
          }
        }
#ifndef NDEBUG
        static constexpr const size_t d_max = 0;
        const size_t edist
          = fastedist_inplace<char,size_t,false,lcplen_fwd<match_method,true>>
                             (d_max,
                              ref_seq + extended_match.s_startpos_get(),
                              extended_match.s_len_get(),
                              query_seq + extended_match.q_startpos_get(),
                              extended_match.q_len_get(),
                              ref_seqnum,query_seqnum);
        if (edist > extended_match.distance_get())
        {
          std::cerr << "Unexpected edist=" << edist << " > "
                    << extended_match.distance_get() << '\n';
          exit(EXIT_FAILURE);
        }
#endif
      }
    }
  }
}

template<bool self_match,class SeedTable,bool (*match_method)(char,char),
         class MatchClass>
static void segment_chaining_template(const GuideParams &guide_params,
                                      const GttlMultiseq &ref_multiseq,
                                      const GttlMultiseq &query_multiseq,
                                      const SeedTable &seed_table,
                                      FILE *out_fp)
{
  if (seed_table.size() == 0)
  {
    return;
  }
  std::map<size_t,size_t> segment_length_dist{};
  Segmentation<SeedTable> segmentation(seed_table);

  const double matchscore_bias = sequence_bias_get(ref_multiseq);
  constexpr const size_t history_size = 64;
  AlignmentPolishing alignment_polishing(guide_params
                                           .polishing_error_percentage_get(),
                                         history_size,
                                         matchscore_bias);
  SegmentMatches<MatchClass> segment_matches(guide_params,
                                             ref_multiseq,
                                             query_multiseq,
                                             out_fp);
  for (auto &&[segment_start,segment_length] : segmentation)
  {
    segment_length_dist[segment_length]++;
    assert(segment_length > 0);
    const size_t ref_seqnum = seed_table.ref_seqnum_get(segment_start);
    const size_t query_seqnum = seed_table.query_seqnum_get(segment_start);
    const char *ref_seq = ref_multiseq.sequence_ptr_get(ref_seqnum);
    const char *query_seq = query_multiseq.sequence_ptr_get(query_seqnum);
    const size_t ref_len = ref_multiseq.sequence_length_get(ref_seqnum);
    const size_t query_len = query_multiseq.sequence_length_get(query_seqnum);
    single_segment_chaining_template<self_match,SeedTable,match_method,
                                     MatchClass>
                                    (guide_params,
                                     alignment_polishing,
                                     ref_multiseq,
                                     ref_seqnum,
                                     ref_seq,
                                     ref_len,
                                     query_multiseq,
                                     query_seqnum,
                                     query_seq,
                                     query_len,
                                     seed_table,
                                     segment_start,
                                     segment_length,
                                     guide_params.final_polishing_get()
                                       ? &segment_matches
                                       : nullptr,
                                     out_fp);
  }
#ifndef NDEBUG
  size_t sum_segment_lengths = 0;
  for (auto &&element : segment_length_dist)
  {
    sum_segment_lengths += std::get<0>(element) * std::get<1>(element);
  }
  assert(sum_segment_lengths == seed_table.size());
#endif
  segment_matches.display();
}

#define SEGMENT_CHAINING_TEMPLATE(SELF_MATCH,MATCH_METHOD)\
        segment_chaining_template<SELF_MATCH,SeedTable,MATCH_METHOD,MatchClass>\
                                 (guide_params,\
                                  ref_multiseq,\
                                  query_multiseq,\
                                  seed_table,\
                                  out_fp)

template<class SeedTable,class MatchClass>
void segment_chaining(const GuideParams &guide_params,
                      const GttlMultiseq &ref_multiseq,
                      const GttlMultiseq &query_multiseq,
                      bool var_self_match,
                      const SeedTable &seed_table,
                      FILE *out_fp)
{
  assert(var_self_match == true || var_self_match == false);

  if (var_self_match)
  {
    SEGMENT_CHAINING_TEMPLATE(true,matching_characters_wc);
  } else
  {
    SEGMENT_CHAINING_TEMPLATE(false,matching_characters);
  }
}
#endif
