#ifndef SEGMENTATION_HPP
#define SEGMENTATION_HPP
#include <cstddef>
#include <vector>
#include <cassert>
#include <tuple>

template<class Table>
class Segmentation
{
  struct Iterator
  {
    private:
    const Table &table;
    size_t segment_start,
           current_idx;
    bool exhausted;
    size_t find_end(size_t i,size_t j) const noexcept
    {
      assert(i < j);
      if (table.all_same_segment())
      {
        return table.size();
      }
      while (j < table.size() && table.same_segment(i,j))
      {
        j++;
      }
      return j;
    }
    public:
    Iterator(const Table &_table,
             bool _exhausted)
      : table(_table)
      , segment_start(0)
      , current_idx(find_end(0,1))
      , exhausted(_exhausted)
    {}
    std::pair<size_t,size_t> operator*(void) const noexcept
    {
      assert(segment_start < current_idx);
      return std::make_pair(segment_start,current_idx - segment_start);
    }
    Iterator& operator++(void) /* prefix increment*/
    {
      assert(!exhausted && segment_start < current_idx);
      if (current_idx < table.size())
      {
        segment_start = current_idx;
        current_idx = find_end(segment_start,current_idx+1);
      } else
      {
        exhausted = true;
      }
      return *this;
    }
    bool operator != (const Iterator& other) const noexcept
    {
      return exhausted != other.exhausted;
    }
  };
  const Table &table;
  public:
  Segmentation(const Table &_table)
    : table(_table)
  {}
  Iterator begin(void) const
  {
    return Iterator(table,false);
  }
  Iterator end(void) const
  {
    return Iterator(table,true);
  }
};
#endif
