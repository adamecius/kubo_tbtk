/* Copyright 2019 Kristofer Björnson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/** @package TBTKcalc
 *  @file ComplexVector.h
 *  @brief N-dimensional vector.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_TBTK_COMPLEX_VECTOR
#define COM_DAFER45_TBTK_COMPLEX_VECTOR

#include "TBTK/TBTKMacros.h"

#include <cmath>
#include <complex>
#include <initializer_list>
#include <ostream>
#include <vector>

namespace TBTK{

class ComplexVector{
public:
	/** Constructor. */
	ComplexVector(unsigned int size);

	/** Copy constructor. */
	ComplexVector(const ComplexVector &complexVector);

	/** Move constructor. */
	ComplexVector(ComplexVector &&complexVector);

	/** Constructor. */
	ComplexVector(const std::vector<std::complex<double>> &components);

	/** Destructor. */
	~ComplexVector();

	/** Assignment operator. */
	ComplexVector& operator=(const ComplexVector &complexVector);

	/** Move assignment operator. */
	ComplexVector& operator=(ComplexVector &&complexVector);

	/** Accessor operator. */
	std::complex<double>& operator[](unsigned int n);

	/** Accessor operator. */
	const std::complex<double>& operator[](unsigned int n) const;

	/** Addition operator. */
	const ComplexVector operator+(const ComplexVector &rhs) const;

	/** Subtraction operator. */
	const ComplexVector operator-(const ComplexVector &rhs) const;

	/** Inversion operator. */
	const ComplexVector operator-() const;

	/** Multiplication operator (vector*scalar). */
	const ComplexVector operator*(std::complex<double> rhs) const;

	/** Multiplication operator (scalar*vector). */
	friend const ComplexVector operator*(
		std::complex<double> lhs,
		const ComplexVector &rhs
	);

	/** Division operator. */
	const ComplexVector operator/(std::complex<double> rhs) const;

	/** Returns a unit vector pointing in the same direction as the
	 *  original vector. */
	const ComplexVector unit() const;

	/** Returns a vector that is the component of the vector that is
	 *  parallel to the argument. */
	const ComplexVector parallel(const ComplexVector &v) const;

	/** Norm. */
	double norm() const;

	/** Dot product. */
	static std::complex<double> dotProduct(
		const ComplexVector &lhs,
		const ComplexVector &rhs
	);

	/** Get the size of the ComplexVector.
	 *
	 *  @return The size of the ComplexVector. */
	unsigned int getSize() const;

	/** Get a std::vector<double> representation of the vector. */
	const std::vector<std::complex<double>> getStdVector() const;

	/** operator<< for ostream. */
	friend std::ostream& operator<<(std::ostream &stream, const ComplexVector &v);
private:
	/** Data size. */
	unsigned int size;

	/** Data. */
	std::complex<double> *data;
};

inline std::complex<double>& ComplexVector::operator[](unsigned int n){
	return data[n];
}

inline const std::complex<double>& ComplexVector::operator[](unsigned int n) const{
	return data[n];
}

inline const ComplexVector ComplexVector::operator+(const ComplexVector &rhs) const{
	TBTKAssert(
		size == rhs.size,
		"ComplexVector::operator+()",
		"Incompatible dimensions. Left hand side has " << size
		<< " components, while the right hand side has " << rhs.size
		<< " components.",
		""
	);

	ComplexVector result(size);
	for(unsigned int n = 0; n < size; n++)
		result.data[n] = data[n] + rhs.data[n];

	return result;
}

inline const ComplexVector ComplexVector::operator-(const ComplexVector &rhs) const{
	TBTKAssert(
		size == rhs.size,
		"ComplexVector::operator-()",
		"Incompatible dimensions. Left hand side has " << size
		<< " components, while the right hand side has " << rhs.size
		<< " components.",
		""
	);

	ComplexVector result(size);
	for(unsigned int n = 0; n < size; n++)
		result.data[n] = data[n] - rhs.data[n];

	return result;
}

inline const ComplexVector ComplexVector::operator-() const{
	ComplexVector result(size);
	for(unsigned int n = 0; n < size; n++)
		result.data[n] = -data[n];

	return result;
}

inline const ComplexVector ComplexVector::operator*(std::complex<double> rhs) const{
	ComplexVector result(size);
	for(unsigned int n = 0; n < size; n++)
		result.data[n] = data[n]*rhs;

	return result;
}

inline const ComplexVector operator*(
	std::complex<double> lhs,
	const ComplexVector &rhs
){
	ComplexVector result(rhs.size);
	for(unsigned int n = 0; n < rhs.size; n++)
		result.data[n] = lhs*rhs.data[n];

	return result;
}

inline const ComplexVector ComplexVector::operator/(
	std::complex<double> rhs
) const{
	ComplexVector result(size);
	for(unsigned int n = 0; n < size; n++)
		result.data[n] = data[n]/rhs;

	return result;
}

inline const ComplexVector ComplexVector::unit() const{
	return (*this)/norm();
}

inline const ComplexVector ComplexVector::parallel(
	const ComplexVector &v
) const{
	TBTKAssert(
		size == v.size,
		"ComplexVector::parallel()",
		"Incompatible dimensions.",
		""
	);

	return dotProduct(*this, v.unit())*v.unit();
}

inline double ComplexVector::norm() const{
	return real(sqrt(dotProduct(*this, *this)));
}

inline std::complex<double> ComplexVector::dotProduct(
	const ComplexVector &lhs,
	const ComplexVector &rhs
){
	TBTKAssert(
		lhs.size == rhs.size,
		"ComplexVector::dotProduct()",
		"Incompatible dimensions. Left hand side has " << lhs.size
		<< " components, while the right hand side has " << rhs.size
		<< " components.",
		""
	);

	std::complex<double> dp = 0;
	for(unsigned int n = 0; n < lhs.size; n++)
		dp += conj(lhs.data[n])*rhs.data[n];

	return dp;
}

inline unsigned int ComplexVector::getSize() const{
	return size;
}

inline const std::vector<std::complex<double>> ComplexVector::getStdVector() const{
	std::vector<std::complex<double>> result;
	for(unsigned int n = 0; n < size; n++)
		result.push_back(data[n]);

	return result;
}

inline std::ostream& operator<<(std::ostream &stream, const ComplexVector &v){
	stream << "(";
	for(unsigned int n = 0; n < v.size; n++){
		if(n != 0)
			stream << ", ";
		stream << v.data[n];
	}

	return stream;
}

};	//End namespace TBTK

#endif
