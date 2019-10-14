/* Copyright 2019 Jose H. Garcia and Kristofer Björnson
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

/** @package Kubo
 *  @file main.cpp
 *  @brief Kubo
 *
 *  @author Jose H. Garcia
 *  @author Kristofer Björnson
 */

#include "TBTK/SparseMatrix.h"

#include "ComplexVector.h"

#include <complex>
#include <vector>

namespace TBTK{

//A sparse matrix class needed to perform matrix-vector multiplication
class KuboSparseMatrix{
	typedef std::complex<double> my_complex;
	typedef std::vector<std::complex<double> > my_complex_vec;
public:
	KuboSparseMatrix();

	KuboSparseMatrix(
		const SparseMatrix<std::complex<double>> &sparseMatrix
	);

	ComplexVector operator*(const ComplexVector &rhs);

	void self_rescaling(const my_complex a){
	}; // perform the transformation A ---> a A

	void self_rescaling(const double a){
	}; // perform the transformation A ---> a A
private:
	unsigned int numRows;
	unsigned int numColumns;
	std::vector<unsigned int> columnPointers;
	std::vector<unsigned int> rows;
	std::vector<std::complex<double>> values;
};

};	//End of namespace TBTK
