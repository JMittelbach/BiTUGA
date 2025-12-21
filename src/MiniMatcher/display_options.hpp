#ifndef DISPLAY_OPTIONS_HPP
#define DISPLAY_OPTIONS_HPP

#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include "utilities/string_tokenizer.hpp"
#include "utilities/find_lit_string.hpp"

class DisplayOptions
{
  static constexpr const GttlLitStringInitializerList display_keys
  {
    "swallow",
    "s.seq",
    "q.seq",
    "minimizers",
    "strand",
    "distance",
    "error_percentage",
    "q.read_pairs",
    "q.seqlen",
    "chain_elements",
    "chain_closed",
    "chain",
    "s.seqid",
    "q.seqid",
    "cigar_string",
    "evalue",
    "non_redundant"
  };
  uint64_t flags;
  public:
  DisplayOptions(void)
    : flags(0)
  { }
  std::string display_help_string(void) noexcept
  {
    std::string help_string{};
    for (auto &&arg : display_keys)
    {
      if (help_string.size() == 0)
      {
        help_string += std::string(arg);
      } else
      {
        assert(arg != nullptr);
        help_string += ", " + std::string(arg);
      }
    }
    return help_string;
  }

  void set_flags(const std::string &argstring)
  {
    StringTokenizer st(argstring);
    std::vector<std::string> display_keys_vector;
    for (auto &&dk : display_keys)
    {
      assert(dk != nullptr);
      display_keys_vector.push_back(std::string(dk));
    }
    for (auto &&arg : st.token_list_get())
    {
      auto found = std::find(display_keys_vector.begin(),
                             display_keys_vector.end(),
                             arg);
      if (found == display_keys_vector.end())
      {
        throw std::invalid_argument(std::string("illegal argument \"") + arg +
                                    std::string("\" to option -d/--display, "
                                                "possible values: ") +
                                    display_help_string());
      }
      this->flags |= (uint64_t(1) <<
                      static_cast<int>(found - display_keys_vector.begin()));
    }
  }
#define DISPLAY_FUNC_GENERATE2(KEY,KEY2DISPLAY)\
  bool KEY##_display(void) const noexcept\
  {\
    return flags & (uint64_t(1) << \
                    gttl_find_lit_string_at_compile_time(display_keys,\
                                                         #KEY2DISPLAY));\
  }
#define DISPLAY_FUNC_GENERATE1(KEY) DISPLAY_FUNC_GENERATE2(KEY,KEY)
  DISPLAY_FUNC_GENERATE1(swallow)
  DISPLAY_FUNC_GENERATE2(s_seq,s.seq)
  DISPLAY_FUNC_GENERATE2(q_seq,q.seq)
  DISPLAY_FUNC_GENERATE1(minimizers)
  DISPLAY_FUNC_GENERATE1(strand)
  DISPLAY_FUNC_GENERATE1(distance)
  DISPLAY_FUNC_GENERATE1(error_percentage)
  DISPLAY_FUNC_GENERATE2(q_read_pairs,q.read_pairs)
  DISPLAY_FUNC_GENERATE2(q_seqlen,q.seqlen)
  DISPLAY_FUNC_GENERATE1(chain_elements)
  DISPLAY_FUNC_GENERATE1(chain_closed)
  DISPLAY_FUNC_GENERATE2(s_seqid,s.seqid)
  DISPLAY_FUNC_GENERATE2(q_seqid,q.seqid)
  DISPLAY_FUNC_GENERATE1(cigar_string)
  DISPLAY_FUNC_GENERATE1(evalue)
  DISPLAY_FUNC_GENERATE1(non_redundant);

  bool chain(void) const noexcept
  {
    /* chain_elements implies chain */
    return (flags & (uint64_t(1)
                     << gttl_find_lit_string_at_compile_time(display_keys,
                                                             "chain"))) or
           chain_elements_display() or
           chain_closed_display();
  }

  std::string to_string(void) const noexcept
  {
    std::string fields("Fields: ");
    if (non_redundant_display())
    {
      fields += "non_red, ";
    }
    if (s_seqid_display())
    {
      fields += "s.seqid, ";
    } else
    {
      fields += "s.seqnum, ";
    }
    fields += "s.start, s.len, ";
    if (strand_display())
    {
      fields += "strand, ";
    }
    if (q_seqid_display())
    {
      fields += "q.seqid, ";
    } else
    {
      fields += "q.seqnum, ";
    }
    if (q_read_pairs_display())
    {
      fields += "q.readinst, ";
    }
#ifdef WITH_DIAGONAL
    fields += "q.start, diag, q.len";
#else
    fields += "q.start, q.len";
#endif
    if (distance_display())
    {
      fields += ", editdist";
    }
    if (error_percentage_display())
    {
      fields += ", err_perc";
    }
    if (evalue_display())
    {
      fields += ", evalue";
    }
    if (q_seqlen_display())
    {
      fields += ", q.seqlen";
    }
    if (s_seq_display())
    {
      fields += ", s.seq";
    }
    if (q_seq_display())
    {
      fields += ", q.seq";
    }
    if (cigar_string_display())
    {
      fields += ", cigar";
    }
    return fields;
  }
};
#endif
