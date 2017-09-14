#pragma once

#include <matrix.h>

#include <algorithm>
#include <initializer_list>

namespace mx
{
    class ArraySize 
    {
    public:
        ArraySize(mxArray const * m)
        :   dim_(mxGetNumberOfDimensions(m))
        ,   size_(mxGetDimensions(m))
        {
        }

        using const_iterator = mwSize const *;

        const_iterator begin() const
        {
            return size_;
        }

        const_iterator end() const
        {
            return size_ + dim_;
        }

        mwSize dim() const
        {
            return dim_;
        }

        mwSize operator[] (mwSize i) const
        {
            return i < dim_ ? size_[i] : 1;
        }

        template <typename T>
        bool isEqual(std::initializer_list<T> l) const
        {
            return dim() == l.size() && std::equal(begin(), end(), l.begin());
        }

    private:
        mwSize dim_;
        mwSize const * size_;
    };

    ArraySize GetSize(mxArray const * m)
    {
        return ArraySize(m);
    }

    inline bool operator==(ArraySize const& s1, ArraySize const& s2)
    {
        return s1.dim() == s2.dim() && std::equal(s1.begin(), s1.end(), s2.begin());
    }

    inline bool operator!=(ArraySize const& s1, ArraySize const& s2)
    {
        return !(s1 == s2);
    }
}