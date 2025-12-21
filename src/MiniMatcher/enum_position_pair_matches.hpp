#ifndef ENUM_POSITION_PAIR_MATCHES_HPP
#define ENUM_POSITION_PAIR_MATCHES_HPP
#include <cstdint>
#include <vector>
#include <typeinfo>
#include <bit>
#include <map>
#include "utilities/bytes_unit.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/bitpacker.hpp"
#include "utilities/is_big_endian.hpp"
#include "segmentation.hpp"

enum SelfMatchMode
{
  Regular,
  WithReverseComplement,
  None
};

class DiagonalEncoder
{
  private:
    const size_t ref_max_len;
  public:
    DiagonalEncoder(size_t _ref_max_len)
      : ref_max_len(_ref_max_len)
    { }
    uint64_t encode(size_t ref_pos,uint32_t query_pos) const noexcept
    {
      assert(ref_max_len + static_cast<uint64_t>(query_pos) >= ref_pos);
      return static_cast<uint64_t>(ref_max_len +
                                   static_cast<size_t>(query_pos) - ref_pos);
    }
    uint64_t decode(uint64_t value,uint64_t query_pos) const noexcept
    {
      assert(ref_max_len + query_pos >= value);
      return ref_max_len + query_pos - value;
    }
};

template <class HashedQgramsClass>
static inline size_t match_pos_pair_determine_length(
    const HashedQgramsClass &hashed_qgrams,
    uint64_t first_code,
    size_t start_range)
{
  size_t idx;

  for (idx = start_range + 1;
       idx < hashed_qgrams.size() &&
       first_code == hashed_qgrams.hash_value_get(idx);
       idx++)
    /* Nothing */;
  return idx - start_range;
}

template <typename StoreType,class HashedQgramsClass>
static inline void qgram_store_extract_coords(
    std::vector<StoreType> *dec,
    const HashedQgramsClass &hashed_qgrams,
    uint64_t first_code,
    size_t start_range)
{
  dec->clear();
  for (size_t idx = start_range;
       idx < hashed_qgrams.size() &&
       first_code == hashed_qgrams.hash_value_get(idx);
       idx++)
  {
    assert(hashed_qgrams.sequence_number_get(idx) <= UINT32_MAX &&
           hashed_qgrams.startpos_get(idx) <= UINT32_MAX);
    dec->emplace_back(
      StoreType(static_cast<uint32_t>(hashed_qgrams.sequence_number_get(idx)),
                static_cast<uint32_t>(hashed_qgrams.startpos_get(idx)))
    );
  }
}

template <SelfMatchMode self_match,
          class RefHashedQgramsClass,
          class QueryHashedQgramsClass,
          bool store_seeds,
          bool use_diagonal_encoder,
          int sizeof_unit_pospair>
class EnumPositionPairMatches
{
  using PositionPairEncoding = BytesUnit<sizeof_unit_pospair,4>;
  using PosPairPacker = GttlBitPacker<sizeof_unit_pospair,4>;
  class PositionPairStore
  {
    std::vector<PositionPairEncoding> vec;
    const PosPairPacker &pos_pair_packer;
    const bool all_same_segment_var;
    const bool has_complement_for_reference;
    const bool has_query_read_pairs;
    public:
    PositionPairStore(const PosPairPacker &_pos_pair_packer,
                      const bool _all_same_segment_var,
                      const bool _has_complement_for_reference,
                      const bool _has_query_read_pairs)
      : pos_pair_packer(_pos_pair_packer)
      , all_same_segment_var(_all_same_segment_var)
      , has_complement_for_reference(_has_complement_for_reference)
      , has_query_read_pairs(_has_query_read_pairs)
    { }
    size_t size(void) const noexcept
    {
      return vec.size();
    }
    void emplace_back(PositionPairEncoding elem)
    {
      vec.emplace_back(elem);
    }
    PositionPairEncoding operator [](size_t idx) const noexcept
    {
      assert(idx < vec.size());
      return vec[idx];
    }
    uint8_t *byte_array(void)
    {
      return reinterpret_cast<uint8_t *>(vec.data());
    }
    size_t size_in_RAM(void) const noexcept
    {
      return sizeof(PositionPairEncoding) * this->size();
    }
    bool all_same_segment(void) const noexcept
    {
      return all_same_segment_var;
    }
    bool same_segment(size_t i, size_t j) const noexcept
    {
      const PositionPairEncoding &bu_i = vec[i],
                                 &bu_j = vec[j];
#undef STRICT_BY_SEQUENCE
#ifdef STRICT_BY_SEQUENCE
      return (bu_i.template decode_at<0>(pos_pair_packer) ==
              bu_j.template decode_at<0>(pos_pair_packer)) &&
             (bu_i.template decode_at<1>(pos_pair_packer) ==
              bu_j.template decode_at<1>(pos_pair_packer));
#else
      const bool same_refer
        = has_complement_for_reference
            ? (bu_i.template decode_at<0>(pos_pair_packer)/2 ==
               bu_j.template decode_at<0>(pos_pair_packer)/2)
            : (bu_i.template decode_at<0>(pos_pair_packer) ==
               bu_j.template decode_at<0>(pos_pair_packer));
      if (not same_refer)
      {
        return false;
      }
      const bool same_query
         = has_query_read_pairs
             ? (bu_i.template decode_at<1>(pos_pair_packer)/2 ==
                bu_j.template decode_at<1>(pos_pair_packer)/2)
             : (bu_i.template decode_at<1>(pos_pair_packer) ==
                bu_j.template decode_at<1>(pos_pair_packer));
      return same_query;
#endif
    }
  };
  template<typename T>
  struct PositionPair
  {
    T seqnum0,
      seqnum1,
      startpos1,
      startpos0,
      length;
    PositionPair(void)
    { }
    PositionPair(T _seqnum0,
                 T _seqnum1,
                 T _startpos1,
                 T _startpos0,
                 T _length = 0)
     : seqnum0(_seqnum0)
     , seqnum1(_seqnum1)
     , startpos1(_startpos1)
     , startpos0(_startpos0)
     , length(_length)
    { }
  };
  static constexpr const size_t ref_first_seqnum = 0;
  static constexpr const size_t query_first_seqnum = 0;

  struct DecodedEnvCoords
  {
    uint32_t seqnum, startpos;
    DecodedEnvCoords(uint32_t _seqnum, uint32_t _startpos)
      : seqnum(_seqnum)
      , startpos(_startpos)
    { }
  };

  class IteratorStreamed
  {
    PositionPair<uint64_t> position_pair;
    std::vector<DecodedEnvCoords> query_dec;
    std::vector<PositionPair<uint64_t>> pos_pair_buf;
    size_t pos_pair_buf_max_size;
    const RefHashedQgramsClass &ref_hashed_qgrams;
    const QueryHashedQgramsClass &query_hashed_qgrams;
    size_t ref_start_range, query_start_range, buf_index;
    bool exhausted;
    FILE *out_fp;
#undef WITH_COMBINATION_HISTOGRAM
#ifdef WITH_COMBINATION_HISTOGRAM
    std::map<std::pair<size_t,size_t>,size_t> combination_histogram;
#endif

    public:
    void assign_next_match_to_pp(void) noexcept
    {
      assert(!exhausted && ref_first_seqnum % 2 == 0 &&
                           query_first_seqnum % 2 == 0);
      while (true)
      {
        if (buf_index < pos_pair_buf.size())
        {
          position_pair = pos_pair_buf[buf_index++];
          break;
        }
        if (ref_start_range >= ref_hashed_qgrams.size() ||
            query_start_range >= query_hashed_qgrams.size())
        {
          exhausted = true;
          break;
        }
        const uint64_t ref_first_code
          = ref_hashed_qgrams.hash_value_get(ref_start_range);
        const uint64_t query_first_code
          = query_hashed_qgrams.hash_value_get(query_start_range);

        if (ref_first_code < query_first_code)
        {
          ref_start_range
            += match_pos_pair_determine_length<RefHashedQgramsClass>
                                              (ref_hashed_qgrams,
                                               ref_first_code,
                                               ref_start_range);
        } else
        {
          if (ref_first_code > query_first_code)
          {
            query_start_range
              += match_pos_pair_determine_length<QueryHashedQgramsClass>
                                                (query_hashed_qgrams,
                                                 query_first_code,
                                                 query_start_range);
          } else
          {
            const size_t ref_len
              = match_pos_pair_determine_length<RefHashedQgramsClass>
                                               (ref_hashed_qgrams,
                                                ref_first_code,
                                                ref_start_range);
            qgram_store_extract_coords<DecodedEnvCoords,
                                       QueryHashedQgramsClass>
                                      (&query_dec,
                                       query_hashed_qgrams,
                                       query_first_code,
                                       query_start_range);
#ifdef WITH_COMBINATION_HISTOGRAM
            combination_histogram[std::make_pair(ref_len,query_dec.size())]++;
#endif
            pos_pair_buf.clear();
            for (size_t r_idx = ref_start_range;
                 r_idx < ref_start_range + ref_len; r_idx++)
            {
              const size_t ref_seqnum
                = ref_hashed_qgrams.sequence_number_get(r_idx);
              const size_t ref_startpos = ref_hashed_qgrams.startpos_get(r_idx);
              for (DecodedEnvCoords &query_dv : query_dec)
              {
                if constexpr (self_match == Regular)
                {
                  if (ref_first_seqnum + ref_seqnum
                        < query_first_seqnum + query_dv.seqnum ||
                      (ref_first_seqnum + ref_seqnum ==
                       query_first_seqnum + query_dv.seqnum &&
                       ref_startpos < static_cast<size_t>(query_dv.startpos)))
                  {
                    pos_pair_buf.emplace_back(
                      PositionPair<uint64_t>
                        (static_cast<uint64_t>(ref_first_seqnum + ref_seqnum),
                         static_cast<uint64_t>(query_first_seqnum +
                                               query_dv.seqnum),
                         static_cast<uint64_t>(query_dv.startpos),
                         static_cast<uint64_t>(ref_startpos)));
                  }
                } else
                {
                  if constexpr (self_match == None)
                  {
                    pos_pair_buf.emplace_back(
                      PositionPair<uint64_t>
                        (static_cast<uint64_t>(ref_first_seqnum + ref_seqnum),
                         static_cast<uint64_t>(query_first_seqnum +
                                               query_dv.seqnum),
                         static_cast<uint64_t>(query_dv.startpos),
                         static_cast<uint64_t>(ref_startpos)));
                  } else
                  {
                    static_assert(self_match == WithReverseComplement);
                    bool accept_match;
#undef SKDEBUG
#ifdef SKDEBUG
                    std::cout << "# seed: s_seqnum = " << ref_seqnum
                              << ", q_seqnum = " << query_dv.seqnum
                              << ", s_start = " << ref_startpos
                              << ", q_start = " << query_dv.startpos
                              << ", diag = "
                              << static_cast<int>(query_dv.startpos -
                                                  ref_startpos)
                              << ", case ";
#endif
                    if ((ref_first_seqnum + ref_seqnum) % 2 == 0)
                    {
                      if ((query_first_seqnum + query_dv.seqnum) % 2 == 0)
                      {
                        /* both on forward strand */
                        accept_match
                          = ref_first_seqnum + ref_seqnum
                              < query_first_seqnum + query_dv.seqnum ||
                            (ref_first_seqnum + ref_seqnum ==
                             query_first_seqnum + query_dv.seqnum &&
                             ref_startpos <
                             static_cast<size_t>(query_dv.startpos));
#ifdef SKDEBUG
                        std::cout << "+ + ";
                        if (accept_match)
                        {
                          std::cout << "accept";
                        }
#endif
                      } else
                      {
                        /* ref on forward strand and query on reverse strand
                           transform both to opposite strands as
                           reverse complement is always on ref. */
                        accept_match
                          = ref_first_seqnum + ref_seqnum + 1
                              < query_first_seqnum + query_dv.seqnum ||
                            (ref_first_seqnum + ref_seqnum + 1 ==
                             query_first_seqnum + query_dv.seqnum) /* &&
                             ref_startpos >
                             static_cast<size_t>(query_dv.startpos)*/;
#ifdef SKDEBUG
                        std::cout << "+ - ";
                        if (accept_match)
                        {
                          std::cout << "accept";
                        }
#endif
                        if (accept_match)
                        {
                          const size_t ref_seqlen
                            = ref_hashed_qgrams.multiseq.sequence_length_get(
                                      ref_first_seqnum + ref_seqnum);
                          const size_t rev_ref_startpos
                            = ref_seqlen - (ref_startpos +
                                          ref_hashed_qgrams.qgram_length_get());
                          const size_t query_seqlen
                            = query_hashed_qgrams.multiseq.sequence_length_get(
                                      query_first_seqnum + query_dv.seqnum);
                          const size_t fwd_query_startpos
                            = query_seqlen - (query_dv.startpos +
                                        query_hashed_qgrams.qgram_length_get());
                          assert(query_first_seqnum + query_dv.seqnum > 0 &&
                                 ((query_first_seqnum + query_dv.seqnum - 1)
                                  % 2 == 0));
                          if (ref_first_seqnum + ref_seqnum + 1 ==
                              query_first_seqnum + query_dv.seqnum &&
                              fwd_query_startpos < rev_ref_startpos)
                          {
                            pos_pair_buf.emplace_back(
                              PositionPair<uint64_t>(
                                 static_cast<uint64_t>(ref_first_seqnum +
                                                       ref_seqnum + 1),
                                 static_cast<uint64_t>(query_first_seqnum +
                                                       query_dv.seqnum - 1),
                                 static_cast<uint64_t>(rev_ref_startpos),
                                 static_cast<uint64_t>(fwd_query_startpos)));
                          } else
                          {
                            pos_pair_buf.emplace_back(
                              PositionPair<uint64_t>(
                                 static_cast<uint64_t>(ref_first_seqnum +
                                                       ref_seqnum + 1),
                                 static_cast<uint64_t>(query_first_seqnum +
                                                       query_dv.seqnum - 1),
                                 static_cast<uint64_t>(fwd_query_startpos),
                                 static_cast<uint64_t>(rev_ref_startpos)));
                          }
                          accept_match = false;
                        }
                      }
                    } else
                    {
                      if ((query_first_seqnum + query_dv.seqnum) % 2 == 0)
                      {
                        /* ref on reverse strand and query on forward strand */
                        accept_match
                          = ref_first_seqnum + ref_seqnum
                              < query_first_seqnum + query_dv.seqnum ||
                            (ref_first_seqnum + ref_seqnum ==
                             query_first_seqnum + query_dv.seqnum + 1 /*&&
                             ref_startpos <=
                             static_cast<size_t>(query_dv.startpos)*/);
#ifdef SKDEBUG
                        std::cout << "- + ";
                        if (accept_match)
                        {
                          std::cout << "accept";
                        }
#endif
                      } else
                      {
                        accept_match
                          = ref_first_seqnum + ref_seqnum
                              < query_first_seqnum + query_dv.seqnum ||
                            (ref_first_seqnum + ref_seqnum ==
                             query_first_seqnum + query_dv.seqnum &&
                             ref_startpos >
                             static_cast<size_t>(query_dv.startpos));
#ifdef SKDEBUG
                        std::cout << "- - ";
                        if (accept_match)
                        {
                          std::cout << "accept";
                        }
#endif
                        if (accept_match)
                        {
                          const size_t ref_seqlen
                            = ref_hashed_qgrams.multiseq.sequence_length_get(
                                      ref_first_seqnum + ref_seqnum-1);
                          const size_t fwd_ref_startpos
                            = ref_seqlen - (ref_startpos +
                                          ref_hashed_qgrams.qgram_length_get());
                          const size_t query_seqlen
                            = query_hashed_qgrams.multiseq.sequence_length_get(
                                      query_first_seqnum + query_dv.seqnum - 1);
                          const size_t fwd_query_startpos
                            = query_seqlen - (query_dv.startpos +
                                        query_hashed_qgrams.qgram_length_get());
                          pos_pair_buf.emplace_back(
                            PositionPair<uint64_t>
                              (
                               static_cast<uint64_t>(ref_first_seqnum +
                                                     ref_seqnum - 1),
                               static_cast<uint64_t>(query_first_seqnum +
                                                     query_dv.seqnum - 1),
                               static_cast<uint64_t>(fwd_query_startpos),
                               static_cast<uint64_t>(fwd_ref_startpos)));
                          accept_match = false;
                        }
                      }
                    }
                    if (accept_match)
                    {
                      pos_pair_buf.emplace_back(
                        PositionPair<uint64_t>
                        (static_cast<uint64_t>(ref_first_seqnum + ref_seqnum),
                         static_cast<uint64_t>(query_first_seqnum +
                                               query_dv.seqnum),
                         static_cast<uint64_t>(query_dv.startpos),
                         static_cast<uint64_t>(ref_startpos)));
                    }
#ifdef SKDEBUG
                    std::cout << '\n';
#endif
                  }
                }
              }
            }
            buf_index = 0;
            ref_start_range += ref_len;
            query_start_range += query_dec.size();
          }
        }
      }
      pos_pair_buf_max_size = std::max(pos_pair_buf_max_size,
                                       pos_pair_buf.size());
    }

    IteratorStreamed(const RefHashedQgramsClass &_ref_hashed_qgrams,
                     const QueryHashedQgramsClass &_query_hashed_qgrams,
                     bool _exhausted,
                     FILE *_out_fp)
      : pos_pair_buf_max_size(0)
      , ref_hashed_qgrams(_ref_hashed_qgrams)
      , query_hashed_qgrams(_query_hashed_qgrams)
      , ref_start_range(0)
      , query_start_range(0)
      , buf_index(0)
      , exhausted(_exhausted)
      , out_fp(_out_fp)
    { }
#ifdef WITH_COMBINATION_HISTOGRAM
    ~IteratorStreamed(void)
    {
      for (auto &[p, c] : combination_histogram)
      {
        fprintf(out_fp,"# CC\t%zu\t%zu\t%zu\n",std::get<0>(p),std::get<1>(p),c);
      }
    }
#endif
    const PositionPair<uint64_t> &operator*() const
    {
      return position_pair;
    }
    auto& operator++() /* prefix increment*/
    {
      assign_next_match_to_pp();
      return *this;
    }
    bool operator != (const IteratorStreamed& other) const noexcept
    {
      return exhausted != other.exhausted;
    }

  };

  class IteratorStored
  {
    PositionPair<uint64_t> position_pair;
    const PositionPairStore &pos_pair_store;
    const PosPairPacker &pos_pair_packer;
    size_t current_idx;
    const DiagonalEncoder &diagonal_encoder;
    public:
    IteratorStored(const PositionPairStore &_pos_pair_store,
                   const PosPairPacker &_pos_pair_packer,
                   size_t _current_idx,
                   const DiagonalEncoder &_diagonal_encoder)
      : pos_pair_store(_pos_pair_store)
      , pos_pair_packer(_pos_pair_packer)
      , current_idx(_current_idx)
      , diagonal_encoder(_diagonal_encoder)
    { }
    const PositionPair<uint64_t> &operator*()
    {
      const PositionPairEncoding &bu = pos_pair_store[current_idx];
            /* XXX: switch 0/1 and 2/3 */
      position_pair.seqnum0 = bu.template decode_at<0>(pos_pair_packer);
      position_pair.seqnum1 = bu.template decode_at<1>(pos_pair_packer);
      position_pair.startpos1 = bu.template decode_at<3>(pos_pair_packer);
      if constexpr (use_diagonal_encoder)
      {
        position_pair.startpos0
          = diagonal_encoder.decode(bu.template decode_at<2>(pos_pair_packer),
                                    position_pair.startpos1);
      } else
      {
        position_pair.startpos0 = bu.template decode_at<2>(pos_pair_packer);
      }
      return position_pair;
    }
    IteratorStored& operator++() /* prefix increment*/
    {
      current_idx++;
      return *this;
    }
    bool operator != (const IteratorStored& other) const noexcept
    {
      return current_idx != other.current_idx;
    }
  };

  struct IteratorMEMs
  {
    private:
    const PositionPairStore &pos_pair_store;
    size_t pos_pair_store_idx;
    const PosPairPacker &pos_pair_packer;
    const DiagonalEncoder &diagonal_encoder;
    const size_t qgram_length;
    const std::vector<uint32_t> &consecutive_vec;
    size_t consecutive_vec_idx;
    PositionPair<uint64_t> position_pair;
    public:
    IteratorMEMs(const PositionPairStore &_pos_pair_store,
                 const PosPairPacker &_pos_pair_packer,
                 const DiagonalEncoder &_diagonal_encoder,
                 size_t _qgram_length,
                 const std::vector<uint32_t> &_consecutive_vec,
                 size_t _consecutive_vec_idx)
      : pos_pair_store(_pos_pair_store)
      , pos_pair_store_idx(0)
      , pos_pair_packer(_pos_pair_packer)
      , diagonal_encoder(_diagonal_encoder)
      , qgram_length(_qgram_length)
      , consecutive_vec(_consecutive_vec)
      , consecutive_vec_idx(_consecutive_vec_idx)
    {
      static_assert(use_diagonal_encoder);
    }
    const PositionPair<uint64_t> operator*()
    {
      const PositionPairEncoding &bu = pos_pair_store[pos_pair_store_idx];
      const uint64_t current_seqpos0
        = bu.template decode_at<0>(pos_pair_packer);
      const uint64_t current_seqpos1
        = bu.template decode_at<1>(pos_pair_packer);
      const uint64_t current_diag = bu.template decode_at<2>(pos_pair_packer);
      const uint64_t current_query_startpos
        = bu.template decode_at<3>(pos_pair_packer);
      return PositionPair<uint64_t>(
               current_seqpos0,
               current_seqpos1,
               current_query_startpos,
               diagonal_encoder.decode(current_diag,current_query_startpos),
               qgram_length + consecutive_vec[consecutive_vec_idx]);
    }
    auto& operator++() /* prefix increment*/
    {
      assert(consecutive_vec_idx < consecutive_vec.size());
      pos_pair_store_idx += 1 + consecutive_vec[consecutive_vec_idx];
      consecutive_vec_idx++;
      return *this;
    }
    bool operator != (const IteratorMEMs& other) const noexcept
    {
      return consecutive_vec_idx != other.consecutive_vec_idx;
    }
  };

  private:
    const RefHashedQgramsClass &ref_hashed_qgrams;
    const QueryHashedQgramsClass &query_hashed_qgrams;
    PosPairPacker pos_pair_packer;
    PositionPairStore pos_pair_store;
    const DiagonalEncoder diagonal_encoder;
    uint32_t max_consecutive;
    FILE *out_fp;
    std::vector<uint32_t> consecutive_vec;

    uint64_t refpos_encoding(size_t ref_startpos,
                             [[maybe_unused]] uint32_t query_startpos)
      const noexcept
    {
      if constexpr (use_diagonal_encoder)
      {
        return diagonal_encoder.encode(ref_startpos,query_startpos);
      } else
      {
        return ref_startpos;
      }
    }

    void store_seeds_in_vector(void) noexcept
    {
      size_t ref_start_range_local = 0, query_start_range_local = 0;
      std::vector<DecodedEnvCoords> query_dec_local;

      while (ref_start_range_local < ref_hashed_qgrams.size() &&
             query_start_range_local < query_hashed_qgrams.size())
      {
        const uint64_t ref_first_code
          = ref_hashed_qgrams.hash_value_get(ref_start_range_local);
        const uint64_t query_first_code
          = query_hashed_qgrams.hash_value_get(query_start_range_local);

        if (ref_first_code < query_first_code)
        {
          ref_start_range_local
            += match_pos_pair_determine_length<RefHashedQgramsClass>
                                              (ref_hashed_qgrams,
                                               ref_first_code,
                                               ref_start_range_local);
        } else
        {
          if (ref_first_code > query_first_code)
          {
            query_start_range_local
              += match_pos_pair_determine_length<QueryHashedQgramsClass>
                                                (query_hashed_qgrams,
                                                 query_first_code,
                                                 query_start_range_local);
          } else
          {
            const size_t ref_len
              = match_pos_pair_determine_length<RefHashedQgramsClass>
                                               (ref_hashed_qgrams,
                                                ref_first_code,
                                                ref_start_range_local);
            qgram_store_extract_coords<DecodedEnvCoords,
                                       QueryHashedQgramsClass>
                                      (&query_dec_local,
                                       query_hashed_qgrams,
                                       query_first_code,
                                       query_start_range_local);
            for (size_t r_idx = ref_start_range_local;
                 r_idx < ref_start_range_local + ref_len; r_idx++)
            {
              const size_t ref_seqnum
                = ref_hashed_qgrams.sequence_number_get(r_idx);
              const size_t ref_startpos
                = ref_hashed_qgrams.startpos_get(r_idx);
              for (DecodedEnvCoords &query_dv : query_dec_local)
              {
                bool accept_match;
                if constexpr (self_match == Regular)
                {
                  accept_match = ref_first_seqnum + ref_seqnum
                                   < query_first_seqnum + query_dv.seqnum ||
                                 (ref_first_seqnum + ref_seqnum ==
                                    query_first_seqnum + query_dv.seqnum &&
                                  ref_startpos <
                                  static_cast<size_t>(query_dv.startpos));
                } else
                {
                  if constexpr (self_match == WithReverseComplement)
                  {
                    accept_match = ref_first_seqnum + ref_seqnum
                                     < query_first_seqnum + query_dv.seqnum ||
                                   (ref_first_seqnum + ref_seqnum ==
                                      query_first_seqnum + query_dv.seqnum &&
                                    ref_startpos <
                                    static_cast<size_t>(query_dv.startpos));
                  } else
                  {
                    accept_match = true;
                  }
                }
                if (accept_match)
                {
                  pos_pair_store.emplace_back(
                           /* here a BytesUnit Object of 4 values
                              is created and appended to the
                              end of the pos_pair_store. */
                           PositionPairEncoding
                                    (pos_pair_packer,
                                     {static_cast<uint64_t>(ref_first_seqnum +
                                                            ref_seqnum),
                                      static_cast<uint64_t>(query_first_seqnum +
                                                            query_dv.seqnum),
                                      refpos_encoding(ref_startpos,
                                                      query_dv.startpos),
                                      static_cast<uint64_t>(query_dv.startpos)
                                      }));
                }
              }
            }
            ref_start_range_local += ref_len;
            query_start_range_local += query_dec_local.size();
          }
        }
      }
    }

    int sequences_bits_get(int num,int for_length) const noexcept
    {
      if (num == 0)
      {
        return ref_hashed_qgrams.packer_bit_group_size_get(1+for_length);
      }
      return query_hashed_qgrams.packer_bit_group_size_get(1+for_length);
    }

    int sequences_bits_ref_seqnum_get(void) const noexcept
    {
      return sequences_bits_get(0,0);
    }

    int sequences_bits_ref_length_get(void) const noexcept
    {
      return sequences_bits_get(0,1);
    }

    int sequences_bits_query_seqnum_get() const noexcept
    {
      return sequences_bits_get(1,0);
    }

    int sequences_bits_query_length_get() const noexcept
    {
      return sequences_bits_get(1,1);
    }

    int refpos_number_of_bits(void) const noexcept
    {
      if constexpr (use_diagonal_encoder)
      {
        const size_t ref_max_len
          = ref_hashed_qgrams.multiseq.sequences_maximum_length_get();
        const size_t query_max_len
          = query_hashed_qgrams.multiseq.sequences_maximum_length_get();
        return static_cast<int>(std::bit_width(ref_max_len +
                                               query_max_len + 1));
      } else
      {
        return ref_hashed_qgrams.packer_bit_group_size_get(2);
      }
    }

  public:
    static constexpr const bool possible_false_positive_matches
      = RefHashedQgramsClass::possible_false_positive_matches;
    static constexpr const bool delivers_length_value = use_diagonal_encoder;
    EnumPositionPairMatches(
      const RefHashedQgramsClass &_ref_hashed_qgrams,
      const QueryHashedQgramsClass &_query_hashed_qgrams,
      FILE *_out_fp)
      : ref_hashed_qgrams(_ref_hashed_qgrams)
      , query_hashed_qgrams(_query_hashed_qgrams)
      , pos_pair_packer(/* create pos_pair_packer, by providing the
                           number of bits required for each bitgroup. */
           /* XXX exchange order of 1/1 and 2/2 */
           {_ref_hashed_qgrams.packer_bit_group_size_get(1), /*rseq*/
            _query_hashed_qgrams.packer_bit_group_size_get(1),/*qseq*/
            refpos_number_of_bits(),
            _query_hashed_qgrams.packer_bit_group_size_get(2)})/*qpos*/
      , pos_pair_store(pos_pair_packer,
                       _ref_hashed_qgrams.multiseq.sequences_number_get() == 1
                       &&
                       _query_hashed_qgrams.multiseq.sequences_number_get()==1,
                       _ref_hashed_qgrams.multiseq
                         .has_reverse_complement_is_set(),
                       _query_hashed_qgrams.multiseq.has_read_pairs_is_set())
      , diagonal_encoder(DiagonalEncoder(_ref_hashed_qgrams.multiseq
                                           .sequences_maximum_length_get()))
      , max_consecutive(0)
      , out_fp(_out_fp)
    {
      assert(ref_hashed_qgrams.qgram_length_get() ==
             query_hashed_qgrams.qgram_length_get());
      if constexpr (store_seeds)
      {
        store_seeds_in_vector();
        if constexpr (use_diagonal_encoder)
        {
          if (pos_pair_store.size() >= size_t(1))
          {
            const bool reversed_byte_order = is_big_endian() ? false : true;
            ska_large_lsb_small_radix_sort(sizeof_unit_pospair,
                                           sizeof_unit_pospair * CHAR_BIT,
                                           pos_pair_store.byte_array(),
                                           pos_pair_store.size(),
                                           reversed_byte_order);
            Segmentation segmentation(pos_pair_store);
            for (auto &&[segment_start,segment_length] : segmentation)
            {
              const PositionPairEncoding &bu0 = pos_pair_store[segment_start];
              uint32_t consecutive = 0;
              uint64_t previous_diag
                = bu0.template decode_at<2>(pos_pair_packer);
              uint64_t previous_query_startpos
                = bu0.template decode_at<3>(pos_pair_packer);
              for (size_t idx = segment_start + 1;
                   idx < segment_start + segment_length; idx++)
              {
                const PositionPairEncoding &bu = pos_pair_store[idx];
                const uint64_t current_diag
                  = bu.template decode_at<2>(pos_pair_packer);
                const uint64_t current_query_startpos
                  = bu.template decode_at<3>(pos_pair_packer);
                if (previous_diag == current_diag &&
                    previous_query_startpos + 1 == current_query_startpos)
                {
                  consecutive++;
                } else
                {
                  max_consecutive = std::max(max_consecutive,consecutive);
                  consecutive_vec.push_back(consecutive);
                  consecutive = 0;
                }
                previous_query_startpos = current_query_startpos;
                previous_diag = current_diag;
              }
              consecutive_vec.push_back(consecutive);
              max_consecutive = std::max(max_consecutive,consecutive);
            }
          }
        }
      }
    }

    int sequences_bits_sum_get(void) const noexcept
    {
      int sum = 0;
      for (int num = 0; num <= 1; num++)
      {
        for (int for_length = 0; for_length <= 1; for_length++)
        {
          sum += sequences_bits_get(num,for_length);
        }
      }
      return sum;
    }

    auto begin(void) const
    {
      if constexpr (store_seeds)
      {
        if constexpr (use_diagonal_encoder)
        {
          return IteratorMEMs(pos_pair_store,
                              pos_pair_packer,
                              diagonal_encoder,
                              ref_hashed_qgrams.qgram_length_get(),
                              consecutive_vec,
                              0);
        } else
        {
          return IteratorStored(pos_pair_store,
                                pos_pair_packer,
                                0,
                                diagonal_encoder);
        }
      } else
      {
        IteratorStreamed iterator(ref_hashed_qgrams,query_hashed_qgrams,false,
                                  out_fp);
        iterator.assign_next_match_to_pp();
        return iterator;
      }
    }

    auto end(void) const
    {
      if constexpr (store_seeds)
      {
        if constexpr (use_diagonal_encoder)
        {
          return IteratorMEMs(pos_pair_store,
                              pos_pair_packer,
                              diagonal_encoder,
                              ref_hashed_qgrams.qgram_length_get(),
                              consecutive_vec,
                              consecutive_vec.size());
        } else
        {
          return IteratorStored(pos_pair_store,
                                pos_pair_packer,
                                pos_pair_store.size(),
                                diagonal_encoder);
        }
      } else
      {
        return IteratorStreamed(ref_hashed_qgrams, query_hashed_qgrams, true,
                                out_fp);
      }
    }

    size_t size(void) const noexcept
    {
      if constexpr (store_seeds)
      {
        return pos_pair_store.size();
      } else
      {
        return 0;
      }
    }

    size_t size_in_RAM(void) const noexcept
    {
      if constexpr (store_seeds)
      {
        return pos_pair_store.size_in_RAM();
      } else
      {
        return 0;
      }
    }

    template<int ref_idx>
    std::array<int,5> match_packer_order_units(int remaining_bits_for_length)
           const noexcept
    {
      if constexpr (ref_idx == 0)
      {
        return {
                sequences_bits_ref_seqnum_get(),
                sequences_bits_query_seqnum_get(),
                sequences_bits_ref_length_get(),
                sequences_bits_query_length_get(),
                remaining_bits_for_length
               };
      } else
      {
        return {
                sequences_bits_query_seqnum_get(),
                sequences_bits_ref_seqnum_get(),
                sequences_bits_query_length_get(),
                sequences_bits_ref_length_get(),
                remaining_bits_for_length
               };
      }
    }

    size_t maximum_match_length_get(void) const noexcept
    {
      if constexpr (use_diagonal_encoder)
      {
        return ref_hashed_qgrams.qgram_length_get() + max_consecutive;
      } else
      {
        return 0;
      }
    }

    size_t number_of_MEMs_get(void) const noexcept
    {
      return consecutive_vec.size();
    }
};
#endif
