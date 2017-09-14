#pragma once

#include <matrix.h>

#include <cstdint>

namespace mx
{
	template<class T>
	struct TypeTraits;
	
	template<>
	struct TypeTraits<std::uint32_t>
	{
		static const auto mxClassID = mxUINT32_CLASS;
    };

    template<>
	struct TypeTraits<std::int32_t>
	{
		static const auto mxClassID = mxINT32_CLASS;
    };

    template<>
	struct TypeTraits<double>
	{
		static const auto mxClassID = mxDOUBLE_CLASS;
    };
}