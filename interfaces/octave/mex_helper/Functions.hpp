#pragma once

#include "TypeTraits.hpp"

#include <matrix.h>

#include <initializer_list>
#include <stdexcept>

namespace mx
{
	inline bool IsScalar(mxArray const * m)
	{
		return mxIsScalar(m);
    }
    
	// Check that a variable is 1xN or Nx1.
	inline bool IsVector(mxArray const * m)
	{
		return mxGetNumberOfDimensions(m) == 2 && (mxGetM(m) == 1 || mxGetN(m) == 1);
	}

	inline bool IsStruct(mxArray const * m)
	{
		return mxIsStruct(m);
    }

	template <typename T>
	inline bool IsOfType(mxArray const * m)
	{
		return mxGetClassID(m) == TypeTraits<T>::mxClassID;
	}

	inline size_t GetM(mxArray const * m)
	{
		return mxGetM(m);
	}

	inline void SetM(mxArray *pm, mwSize m)
	{
		mxSetM(pm, m);
	}

	inline size_t GetN(mxArray const * m)
	{
		return mxGetN(m);
	}

	inline void SetN(mxArray *pm, mwSize n)
	{
		mxSetN(pm, n);
	}

	inline mwSize GetNumberOfDimensions(const mxArray *pm)
	{
		return mxGetNumberOfDimensions(pm);
	}

	inline size_t GetNumberOfElements(const mxArray *pm)
	{
		return mxGetNumberOfElements(pm);
	}
	
	template<class T>
	inline T GetScalar(mxArray const * ptr_value)
	{
		if (!IsOfType<T>(ptr_value))
			throw std::invalid_argument("GetScalar(): MATLAB array class does not match the requested type");

		if (!IsScalar(ptr_value))
			throw std::invalid_argument("GetScalr(): MATLAB array is not scalar");
	
		return *reinterpret_cast<T const *>(mxGetData(ptr_value));
	};

	inline double const * GetPr(mxArray const * m)
	{
		double const * pr = mxGetPr(m);
		
		if (!pr && GetNumberOfElements(m) > 0)
			throw std::invalid_argument("GetPr(): Unable to get pointer to real data");
		
		return pr;
	}

	inline double * GetPr(mxArray * m)
	{
		double * pr = mxGetPr(m);
		
		if (!pr && GetNumberOfElements(m) > 0)
			throw std::invalid_argument("GetPr(): Unable to get pointer to real data");
		
		return pr;
	}

	template <typename T>
	inline T const * GetData(mxArray const * m)
	{
		if (!IsOfType<T>(m))
			throw std::invalid_argument("GetData(): MATLAB array class does not match the requested type");

		void const * data = mxGetData(m);
		
		if (!data && GetNumberOfElements(m) > 0)
			throw std::invalid_argument("GetData(): Unable to get pointer to data");
		
		return reinterpret_cast<T const *>(data);
	}

	template <typename T>
	inline T * GetData(mxArray * m)
	{
		if (!IsOfType<T>(m))
			throw std::invalid_argument("GetData(): MATLAB array class does not match the requested type");

		void * data = mxGetData(m);
		
		if (!data && GetNumberOfElements(m) > 0)
			throw std::invalid_argument("GetData(): Unable to get pointer to data");
		
		return reinterpret_cast<T *>(data);
	}

	inline mxArray const * GetField(mxArray const * m, mwIndex index, char const * fieldname)
	{
		mxArray const * field = mxGetField(m, index, fieldname);
		
		if (!field)
			throw std::invalid_argument("GetField(): Field " + std::string(fieldname) + " not found");
		
		return field;
	}

	void SetField(mxArray *pm, mwIndex index, const char *fieldname, mxArray *pvalue)
	{
		mxSetField(pm, index, fieldname, pvalue);
	}
	
	template<class T>
	inline mxArray * CreateNumericScalar(T val, mxComplexity complexity_flag = mxREAL)
	{
		mxArray * mat = mxCreateNumericMatrix(1, 1, TypeTraits<T>::mxClassID, complexity_flag);

		if (!mat)
			throw std::runtime_error("CreateNumericScalar(): mxCreateNumericMatrix() failed");

		*GetData<T>(mat) = val;
	
		return mat;
	};

	inline mxArray * CreateDoubleMatrix(mwSize m, mwSize n, mxComplexity ComplexFlag = mxREAL)
	{
		mxArray * mat = mxCreateDoubleMatrix(m, n, ComplexFlag);

		if (!mat)
			throw std::runtime_error("CreateDoubleMatrix(): mxCreateDoubleMatrix() failed");

		return mat;
	}

	template <typename T>
	inline mxArray * CreateNumericArray(std::initializer_list<mwSize> dims, mxComplexity ComplexFlag = mxREAL)
	{
		mxArray * mat = mxCreateNumericArray(dims.size(), dims.begin(), TypeTraits<T>::mxClassID, ComplexFlag);

		if (!mat)
			throw std::runtime_error("CreateNumericArray(): mxCreateNumericArray() failed");

		return mat;
	}

	inline mxArray * CreateStructMatrix(mwSize m, mwSize n, std::initializer_list<char const *> fieldnames)
	{
		mxArray * mat = mxCreateStructMatrix(m, n, fieldnames.size(), const_cast<char const **>(fieldnames.begin()));

		if (!mat)
			throw std::runtime_error("CreateStructMatrix(): mxCreateStructMatrix() failed");

		return mat;
	}
}