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

/** @file ComplexVector.cpp
 *
 *  @author Kristofer Björnson
 */

#include "KuboSparseMatrix.h"

#include "TBTK/TBTKMacros.h"

using namespace std;

namespace TBTK{

KuboSparseMatrix::KuboSparseMatrix(){
}

KuboSparseMatrix::KuboSparseMatrix(
	const SparseMatrix<complex<double>> &sparseMatrix
){
	numRows = sparseMatrix.getNumRows();
	numColumns = sparseMatrix.getNumColumns();

	const unsigned int *columnPointers
		= sparseMatrix.getCSCColumnPointers();
	for(unsigned int n = 0; n < numColumns+1; n++)
		this->columnPointers.push_back(columnPointers[n]);

	const unsigned int *rows = sparseMatrix.getCSCRows();
	for(unsigned int n = 0; n < sparseMatrix.getCSCNumMatrixElements(); n++)
		this->rows.push_back(rows[n]);

	const complex<double> *values = sparseMatrix.getCSCValues();
	for(unsigned int n = 0; n < sparseMatrix.getCSCNumMatrixElements(); n++)
		this->values.push_back(values[n]);
}

ComplexVector KuboSparseMatrix::operator*(const ComplexVector &rhs){
	TBTKAssert(
		numColumns == rhs.getSize(),
		"KuboSparseMatrix::operator*()",
		"Incompatible dimensions. The matrix has '" << numColumns
		<< "' columns, but the vector has '" << rhs.getSize()
		<< "' rows.",
		""
	);

	ComplexVector result(numRows);
	for(unsigned int n = 0; n < numRows; n++)
		result[n] = 0;

	for(unsigned int column = 0; column < columnPointers.size(); column++){
		unsigned int start = columnPointers[column];
		unsigned int end = columnPointers[column+1];
		for(unsigned int n = start; n < end; n++){
			unsigned int row = rows[n];
			result[row] += values[n]*rhs[column];
		}
	}

	return result;
}

};	//End of namespace TBTK
