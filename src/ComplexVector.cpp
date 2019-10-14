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

#include "ComplexVector.h"

#include "TBTK/TBTKMacros.h"

using namespace std;

namespace TBTK{

ComplexVector::ComplexVector(unsigned int size){
	this->size = size;
	data = new complex<double>[size];
}

ComplexVector::ComplexVector(const vector<complex<double>> &components){
	size = components.size();
	data = new complex<double>[size];
	for(unsigned int n = 0; n < size; n++)
		data[n] = components.at(n);
}

ComplexVector::ComplexVector(const ComplexVector &complexVector){
	size = complexVector.size;
	data = new complex<double>[size];
	for(unsigned int n = 0; n < size; n++)
		data[n] = complexVector.data[n];
}

ComplexVector::ComplexVector(ComplexVector &&complexVector){
	size = complexVector.size;
	data = complexVector.data;
	complexVector.data = nullptr;
}

ComplexVector::~ComplexVector(){
	if(data != nullptr)
		delete [] data;
}

ComplexVector& ComplexVector::operator=(const ComplexVector &rhs){
	if(this != &rhs){
		size = rhs.size;
		data = new complex<double>[size];
		for(unsigned int n = 0; n < size; n++)
			data[n] = rhs.data[n];
	}

	return *this;
}

ComplexVector& ComplexVector::operator=(ComplexVector &&rhs){
	if(this != &rhs){
		size = rhs.size;
		data = rhs.data;
		rhs.data = nullptr;
	}

	return *this;
}

};	//End of namespace TBTK
