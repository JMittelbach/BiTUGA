#ifndef MATCH_CLASS_HPP
#define MATCH_CLASS_HPP
#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <cmath>
#include "utilities/mathsupport.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/outsenseedist_aln.hpp"
#include "sequences/matching_characters.hpp"
#include "sequences/prefix_suffix_extension.hpp"
#include "utilities/string_values_join.hpp"
#include "display_options.hpp"
#include "evalue_calculate.hpp"

#define DECL_GETTER(ELEM) size_t ELEM##_get(void) const noexcept {return ELEM;}
#define DECL_SETTER(ELEM)\
        void ELEM##_set(size_t _##ELEM) noexcept { ELEM=_##ELEM; }

class Match
{
  size_t s_endpos;
  size_t s_len;
  size_t q_endpos;
  size_t q_len;
  size_t distance;

  DECL_SETTER(s_endpos)
  DECL_SETTER(s_len)
  DECL_SETTER(q_endpos)
  DECL_SETTER(q_len)
  DECL_SETTER(distance)

  public:
  Match(size_t _s_endpos, size_t _s_len,
        size_t _q_endpos, size_t _q_len,
        size_t _distance)
    : s_endpos(_s_endpos)
    , s_len(_s_len)
    , q_endpos(_q_endpos)
    , q_len(_q_len)
    , distance(_distance)
  { }

  DECL_GETTER(s_endpos)
  DECL_GETTER(s_len)
  DECL_GETTER(q_endpos)
  DECL_GETTER(q_len)
  DECL_GETTER(distance)

  size_t s_startpos_get(void) const noexcept
  {
    return s_endpos - s_len + 1;
  }

  size_t q_startpos_get(void) const noexcept
  {
    return q_endpos - q_len + 1;
  }

  size_t aligned_len_get(void) const noexcept
  {
    return s_len_get() + q_len_get();
  }

  double match_error_percentage_get(void) const noexcept
  {
    return error_percentage_get(distance_get(), aligned_len_get());
  }

  double evalue_get(size_t ref_seqlen, size_t query_seqlen) const noexcept
  {
    return evalue_calculate(ref_seqlen, query_seqlen,
                            0.5 * aligned_len_get(),
                            distance);
  }

  Eoplist match2eoplist(const char *ref_seq,
                        size_t s_seqnum,
                        const char *query_seq,
                        size_t q_seqnum) const noexcept
  {
    static constexpr const bool track_eop = true;
    TrackEditoperations track_editoperations{};
    prefix_or_suffix_extension_generic<track_eop,
                                       TrackEditoperations,
                                       FrontValueTrace,
                                       lcplen_fwd<matching_characters_wc,
                                                  true>
                                      >
                                     (&track_editoperations,
                                      ref_seq + s_startpos_get(),
                                      s_len_get(),
                                      query_seq + q_startpos_get(),
                                      q_len_get(),
                                      s_seqnum,
                                      q_seqnum);
    //track_editoperations.show();
    Eoplist eoplist = track_editoperations.traceback_one(s_len_get(),
                                                         q_len_get());
    [[maybe_unused]] bool has_cut_off = eoplist.cut_off_unpolished_tail();
    return eoplist;
  }
  void correct_coordinates(const Eoplist &eoplist)
  {
    distance_set(eoplist.errors_get());
    const size_t this_s_len = eoplist.aligned_len_u_get();
    const size_t s_startpos = s_startpos_get();
    const size_t this_q_len = eoplist.aligned_len_v_get();
    const size_t q_startpos = q_startpos_get();
    s_len_set(this_s_len);
    s_endpos_set(s_startpos + this_s_len - 1);
    q_len_set(this_q_len);
    q_endpos_set(q_startpos + this_q_len - 1);
  }

  template<bool fast_version>
  void display(FILE *out_fp,
               const DisplayOptions &display_options,
               const GttlMultiseq &ref_multiseq,
               size_t s_seqnum,
               const GttlMultiseq &query_multiseq,
               size_t q_seqnum,
               size_t query_seqnum_offset,
               const Eoplist *eoplist_ptr) const noexcept
  {
    bool forward_strand = true;
    if constexpr (fast_version)
    {
      const size_t query_seqlen = query_multiseq.sequence_length_get(q_seqnum);
      size_t this_s_start;
      if (not ref_multiseq.has_reverse_complement_is_set() or
          s_seqnum % 2 == 0)
      {
        this_s_start = s_startpos_get();
      } else
      {
        const size_t ref_seqlen = ref_multiseq.sequence_length_get(s_seqnum);
        this_s_start = ref_seqlen - (s_endpos + 1);
        forward_strand = false;
      }
      const bool halve = query_multiseq.has_reverse_complement_is_set() or
                         display_options.q_read_pairs_display();
      fprintf(out_fp,"%zu\t%zu\t%zu\t%zu\t%zu%c\t%zu\t%zu\t%zu\n",
              ref_multiseq.has_reverse_complement_is_set() ? s_seqnum/2
                                                           : s_seqnum,
              this_s_start,
              s_len,
              query_seqnum_offset + (halve ? q_seqnum/2
                                           : q_seqnum),
              1 + q_seqnum % 2, forward_strand ? '+' : '-',
              forward_strand ? q_startpos_get() : query_seqlen - (q_endpos + 1),
              q_len,
              query_seqlen);
    } else
    {
      assert(query_seqnum_offset
             == query_multiseq.sequence_number_offset_get());
      if (display_options.swallow_display())
      {
        return;
      }
      assert(out_fp != nullptr);
      const size_t ref_seqlen = ref_multiseq.sequence_length_get(s_seqnum);
      const size_t query_seqlen = query_multiseq.sequence_length_get(q_seqnum);
      std::vector<std::string> line;
      if (display_options.s_seqid_display())
      {
        size_t sh_offset, sh_len;
        std::tie(sh_offset, sh_len) = ref_multiseq.short_header_get(s_seqnum);
        const std::string_view header = ref_multiseq.header_get(s_seqnum);
        line.push_back(std::string(header.begin() + sh_offset,
                                   header.begin() + sh_offset + sh_len));
      } else
      {
        line.push_back(std::format("{}",
                                   ref_multiseq.has_reverse_complement_is_set()
                                     ? s_seqnum/2
                                     : s_seqnum));
      }
      if (ref_multiseq.has_reverse_complement_is_set())
      {
        if (s_seqnum % 2 == 0)
        {
          line.push_back(std::format("{}", s_startpos_get()));
        } else
        {
          forward_strand = false;
          assert(ref_seqlen > s_endpos);
          line.push_back(std::format("{}", ref_seqlen - (s_endpos + 1)));
        }
      } else
      {
        line.push_back(std::format("{}", s_startpos_get()));
      }
      line.push_back(std::format("{}", s_len));
      if (display_options.strand_display())
      {
        line.push_back(forward_strand ? "+" : "-");
      }
      if (display_options.q_seqid_display())
      {
        size_t sh_offset, sh_len;
        std::tie(sh_offset, sh_len) = query_multiseq.short_header_get(q_seqnum);
        const std::string_view header = query_multiseq.header_get(q_seqnum);
        line.push_back(std::string(header.begin() + sh_offset,
                                   header.begin() + sh_offset + sh_len));
      } else
      {
        const bool halve = query_multiseq.has_reverse_complement_is_set() or
                           display_options.q_read_pairs_display();
        line.push_back(std::format("{}", query_seqnum_offset +
                                         (halve ? q_seqnum/2
                                                : q_seqnum)));
      }
      if (display_options.q_read_pairs_display())
      {
        line.push_back(std::format("{}{}", 1 + q_seqnum % 2,
                                          forward_strand ? '+' : '-'));
      }
      if (forward_strand)
      {
        line.push_back(std::format("{}", q_startpos_get()));
      } else
      {
#undef POS_ON_ORIG
#ifdef POS_ON_ORIG
        line.push_back(std::format("{}", q_startpos_get()));
#else
        assert(query_seqlen > q_endpos);
        line.push_back(std::format("{}", query_seqlen - (q_endpos + 1)));
#endif
      }
#ifdef WITH_DIAGONAL
      line.push_back(std::format("{}", static_cast<int64_t>(q_startpos_get() -
                                                           s_startpos_get())));
#endif
      line.push_back(std::format("{}", q_len));
      if (display_options.distance_display())
      {
        line.push_back(std::format("{}", distance));
      }
      if (display_options.error_percentage_display())
      {
        line.push_back(std::format("{:.2f}", match_error_percentage_get()));
      }
      if (display_options.evalue_display())
      {
        line.push_back(std::format("{:.1e}",
                                   evalue_get(ref_seqlen, query_seqlen)));
      }
      if (display_options.q_seqlen_display())
      {
        line.push_back(std::format("{}", query_seqlen));
      }
      if (display_options.s_seq_display())
      {
        const char *ref_seq = ref_multiseq.sequence_ptr_get(s_seqnum);
        line.push_back(std::string(ref_seq + s_startpos_get(), s_len));
#ifndef NDEBUG
        const char *query_seq = query_multiseq.sequence_ptr_get(q_seqnum);
        assert(distance > 0 or
               std::string(ref_seq + s_startpos_get(), s_len) ==
               std::string(query_seq + q_startpos_get(), q_len));
#endif
      }
      if (display_options.q_seq_display())
      {
        const char *query_seq = query_multiseq.sequence_ptr_get(q_seqnum);
        line.push_back(std::string(query_seq + q_startpos_get(), q_len));
  #ifndef NDEBUG
        const char *ref_seq = ref_multiseq.sequence_ptr_get(s_seqnum);
        assert(distance > 0 or
               std::string(ref_seq + s_startpos_get(), s_len) ==
               std::string(query_seq + q_startpos_get(), q_len));
  #endif
      }
      if (display_options.cigar_string_display())
      {
        if (distance == 0)
        {
          line.push_back(std::format("{}=", s_len));
        } else
        {
          assert(eoplist_ptr != nullptr and
                 eoplist_ptr->errors_get() <= distance);
          constexpr const bool distinguish_mismatch_match = true;
          if (eoplist_ptr != nullptr)
          {
            const std::string cigar_string
              = eoplist_ptr->cigar_string_get(distinguish_mismatch_match);
            line.push_back(cigar_string);
          }
        }
      }
      fprintf(out_fp, "%s\n", string_values_join("\t", line.begin(), line.end())
                              .c_str());
    }
  }
};
#endif
