#ifndef GUIDE_PARAMS_HPP
#define GUIDE_PARAMS_HPP
#include <cstddef>
#ifndef NDEBUG
#include <numeric>
#endif
#include <cmath>
#include "display_options.hpp"

class GuideParams
{
  const bool at_constant_distance,
             final_polishing,
             group_by_query,
             verbose_option;
  const size_t qgram_length;
  int hash_bits;
  const size_t minimum_mem_length;
  const size_t minimum_als_length;
  const size_t number_of_threads;
  const size_t thread_pool_max_size;
  const size_t ref_window_size;
  const size_t query_window_size;
  const size_t max_previous;
  const size_t query_split_size;
  const size_t max_replicate_ref;
  const size_t max_replicate_qry;
  const double maximum_error_percentage;
  const double polishing_error_percentage;
  const DisplayOptions &display_options;
  size_t cop_sample_distance_get(void) const noexcept
  {
    assert(qgram_length <= minimum_mem_length);
    if (minimum_mem_length == qgram_length)
    {
      return size_t(1);
    }
    size_t max_window = minimum_mem_length - qgram_length + 1,
           k1 = 1 + std::sqrt(max_window),
           k2 = k1 - 1;
    if (k1 * k2 > max_window)
    {
      assert(k2 > 0);
      k1--;
      k2--;
    }
    assert(k1 * k2 <= max_window);
    assert(std::gcd(k1,k2) == 1);
    return k1;
  }
  public:

  GuideParams(size_t _qgram_length,
              int _hash_bits,
              size_t _minimum_mem_length,
              size_t _minimum_als_length,
              size_t _number_of_threads,
              size_t _thread_pool_max_size,
              size_t _query_split_size,
              size_t _max_replicate_ref,
              size_t _max_replicate_qry,
              double _maximum_error_percentage,
              double _polishing_error_percentage,
              bool cop_option,
              bool _final_polishing,
              bool _group_by_query,
              bool _verbose_option,
              const DisplayOptions &_display_options)
    : at_constant_distance(cop_option)
    , final_polishing(_final_polishing or
                      _display_options.cigar_string_display())
    , group_by_query(_group_by_query)
    , verbose_option(_verbose_option)
    , qgram_length(_qgram_length)
    , hash_bits(_hash_bits)
    , minimum_mem_length(_minimum_mem_length)
    , minimum_als_length(_minimum_als_length)
    , number_of_threads(_number_of_threads)
    , thread_pool_max_size(_thread_pool_max_size)
    , ref_window_size(cop_option ? cop_sample_distance_get()
                                 : minimum_mem_length - qgram_length + 1)
    , query_window_size(cop_option
                         ? (minimum_mem_length == qgram_length
                              ? ref_window_size
                              : (ref_window_size - 1))
                         : ref_window_size)
    , max_previous(size_t(50))
    , query_split_size(_query_split_size)
    , max_replicate_ref(_max_replicate_ref)
    , max_replicate_qry(_max_replicate_qry)
    , maximum_error_percentage(_maximum_error_percentage)
    , polishing_error_percentage(_polishing_error_percentage)
    , display_options(_display_options)
  {
    assert(_minimum_mem_length >= _qgram_length);
  }
  void hash_bits_set_if_undefined(bool with_nt_hash,int sequences_bits) noexcept
  {
    if (hash_bits == -1)
    {
      if (with_nt_hash)
      {
        assert(sequences_bits <= 64);
        hash_bits = 64 - sequences_bits;
      } else
      {
        hash_bits = 2 * qgram_length;
      }
    }
  }
  bool perform_qgram_sampling(void) const noexcept
  {
    return qgram_length < minimum_mem_length;
  }
#define GUIDE_PARAMS_DECLARE_ACCESSOR(RETURN_TYPE,KEY) \
        RETURN_TYPE KEY##_get(void) const noexcept {return KEY;}
  GUIDE_PARAMS_DECLARE_ACCESSOR(bool,at_constant_distance)
  GUIDE_PARAMS_DECLARE_ACCESSOR(bool,final_polishing)
  GUIDE_PARAMS_DECLARE_ACCESSOR(bool,group_by_query);
  GUIDE_PARAMS_DECLARE_ACCESSOR(bool,verbose_option);
  GUIDE_PARAMS_DECLARE_ACCESSOR(int,hash_bits);
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,qgram_length)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,minimum_mem_length)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,minimum_als_length)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,number_of_threads)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,thread_pool_max_size)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,ref_window_size)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,query_window_size)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,max_previous)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,query_split_size)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,max_replicate_ref)
  GUIDE_PARAMS_DECLARE_ACCESSOR(size_t,max_replicate_qry)
  GUIDE_PARAMS_DECLARE_ACCESSOR(double,maximum_error_percentage)
  GUIDE_PARAMS_DECLARE_ACCESSOR(double,polishing_error_percentage)
  GUIDE_PARAMS_DECLARE_ACCESSOR(const DisplayOptions &,display_options)
};
#endif
