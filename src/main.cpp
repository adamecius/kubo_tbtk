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

#include "TBTK/FileWriter.h"
#include "TBTK/Model.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
#include "TBTK/TBTK.h"
#include "TBTK/Timer.h"
#include "TBTK/Vector3d.h"

#include "ComplexVector.h"
#include "KuboSparseMatrix.h"

#include <complex>

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);

void swap(ComplexVector& x, ComplexVector& y){
	ComplexVector temp = std::move(x);
	x = std::move(y);
	y = std::move(temp);
};

//IndexedDataTree<complex<double>> generateVelocityOperator(
vector<tuple<complex<double>, Index>> generateVelocityOperator(
	const Model &model,
	const Vector3d &direction
){
	//Get the Geometry from the model.
	const Geometry &geometry = model.getGeometry();
	TBTKAssert(
		geometry.getDimensions() == 3,
		"generateVelocityOperator()",
		"The geometry dimension must have dimension '3' but has"
		<< " dimension '" << geometry.getDimensions() << "'.",
		""
	);

	//Create the velocity operator.
//	IndexedDataTree<complex<double>> velocityOperator;
	vector<tuple<complex<double>, Index>> velocityOperator;

	//Iterate over all HoppingAmplitudes.
	const HoppingAmplitudeSet &hoppingAmplitudeSet
		= model.getHoppingAmplitudeSet();
	for( HoppingAmplitudeSet::ConstIterator iterator = hoppingAmplitudeSet.cbegin();
		 iterator != hoppingAmplitudeSet.cend();
		 ++iterator)
	{
		//Get the amplitude, to-, and from-Indices.
		complex<double> amplitude = (*iterator).getAmplitude();
		const Index to = (*iterator).getToIndex();
		const Index from = (*iterator).getFromIndex();

		//Get the coordinates corresponding to the to- and from-Indices
		//and vonvert them to Vector3d.
		Vector3d toR = Vector3d(geometry.getCoordinate(to));
		Vector3d fromR = Vector3d(geometry.getCoordinate(from));

		//Calculate H_{ij}(direction.(R_i - R_j)).
		complex<double> i(0,1);
		complex<double> operatorValue =-i*amplitude*Vector3d::dotProduct( //FIXED lack of i. 
			direction,
			toR - fromR
		);

		//Avoid adding elements that are identically zero.
		if(fabs(operatorValue) == 0)
			continue;

		//Add the value to the velocity operator.
/*		try{
			//Try to add it to an already existing entry.
			velocityOperator.get({to, from}) += operatorValue;
		}
		catch(ElementNotFoundException e){
			//Create a new entry if it fails.
			velocityOperator.add(operatorValue, {to, from});
		}*/
		velocityOperator.push_back(make_tuple(operatorValue, Index({to, from})));
	}

	return velocityOperator;
}

KuboSparseMatrix convertModelToKuboSparseMatrix(const Model &model){
	
	SparseMatrix<complex<double>> sparseHamiltonian
		= model.getHoppingAmplitudeSet().getSparseMatrix();
	sparseHamiltonian.setStorageFormat(
		SparseMatrix<complex<double>>::StorageFormat::CSC
	);

	return KuboSparseMatrix(sparseHamiltonian);
}

SparseMatrix<complex<double>> convertIndexedDataTreeToSparseMatrix(
//	const IndexedDataTree<complex<double>> &indexedDataTree,
	const vector<tuple<complex<double>, Index>> &indexedDataTree,
	const Model &model
){
	//Get the HoppingAmplitudeSet.
	const HoppingAmplitudeSet &hoppingAmplitudeSet
		= model.getHoppingAmplitudeSet();

	//Create the resulting SparseMatrix.
	SparseMatrix<complex<double>> sparseMatrix(
		SparseMatrix<complex<double>>::StorageFormat::CSC,
		model.getBasisSize(),
		model.getBasisSize()
	);

	//Iterate over all elements in the IndexedDataTree.
	for(
/*		IndexedDataTree<complex<double>>::ConstIterator iterator
			= indexedDataTree.cbegin();
		iterator != indexedDataTree.cend();
		++iterator*/
		unsigned int n = 0; n < indexedDataTree.size(); n++
	){
		//Get the Index of the current element.
//		const Index &index = iterator.getCurrentIndex();
		const Index &index = get<1>(indexedDataTree[n]);

		//Split Indices such as {{x, y, s}, {x', y', s'}} into the two
		//components {x, y, s} and {x', y', s'}. Assert that two
		//components are obtained.
		vector<Index> components = index.split();
		TBTKAssert(
			components.size() == 2,
			"convertOperator()",
			"Invalid operator index '" << index.toString() << "'."
			<< " The Index must have two component indices, but"
			<< " have '" << components.size() << "'.",
			""
		);

		//Convert the two components (the to- and fromIndices) to their
		//corresponding linear Hilbert space indices.
		int row = hoppingAmplitudeSet.getBasisIndex(components[0]);
		int column = hoppingAmplitudeSet.getBasisIndex(components[1]);

		//Add te element to the SparseMatrix.
		sparseMatrix.add(row, column, get<0>(indexedDataTree[n]));
		std::cout<<row<<" "<<column<<" "<< get<0>(indexedDataTree[n])<<std::endl;
	}

	//Construct the SparseMatrix internal CSC representation.
	sparseMatrix.construct();

	return sparseMatrix;
}

KuboSparseMatrix convertIndexedDataTreeToKuboSparseMatrix(
//	const IndexedDataTree<complex<double>> &indexedDataTree,
	const vector<tuple<complex<double>, Index>> &indexedDataTree,
	const Model &model
){
	SparseMatrix<complex<double>> sparseMatrix
		= convertIndexedDataTreeToSparseMatrix(
			indexedDataTree,
			model
		);

	return KuboSparseMatrix(sparseMatrix);
}

ComplexVector generateRandomPhaseVector(unsigned int basisSize){
	ComplexVector randomPhaseVector(basisSize);
	for(unsigned int n = 0; n < basisSize; n++){
		randomPhaseVector[n] = exp(
			2.0*M_PI*i * (double)rand()/(double)RAND_MAX
		)/sqrt(basisSize);
	}

	return randomPhaseVector;
}

void kuboCalculation(
	const Model &model,
//	const IndexedDataTree<complex<double>> targetOperator,
	const vector<tuple<complex<double>, Index>> &targetOperator,
	double bandWidth,
	double bandCenter,
	int M0,
	int M1
){
	//Convert the Model into a KuboSparseMatrix.
	KuboSparseMatrix H = convertModelToKuboSparseMatrix(model); H.print(); //Check this, very good.	
	H.rescale(2.0/bandWidth); H.print(); //Check this, very good.	

	//Extract the basis size.
	const int basisSize = model.getBasisSize();

	//Convert targetOperator from IndexedDataTree to KuboSparseMatrix.
	KuboSparseMatrix Op = convertIndexedDataTreeToKuboSparseMatrix(
		targetOperator,
		model
	);

	//TODO: Vx should be calculated independently of Op.
	KuboSparseMatrix Vx = Op; Vx.print(); //Check this, very good.	

	//Generate the random phase vector.
	ComplexVector randomPhaseVector = generateRandomPhaseVector(basisSize);

	//Create the result.
	std::vector<complex<double>> mu2D(M0*M1);

	//Need five vectors for performing the iteration
	ComplexVector jLm0(basisSize);
	ComplexVector jLm1(basisSize);
	ComplexVector jRm0(basisSize);
	ComplexVector jRm1(basisSize);
	ComplexVector jV(basisSize);

	//Starts chebyshev iteration
	jLm0 = Vx*randomPhaseVector;
	jLm1 = H*jLm0;
	const double shift = -2*bandCenter/bandWidth;
	for(int m0 = 0; m0 < M0; m0++)
	{
		const complex<double > a=2.0, b=-1.0;
		jRm0 = randomPhaseVector;
		jRm1 = H*jRm0;
		for(int m1 = 0; m1 < M1; m1++){
			jV = Op*jLm0;
			mu2D[m0*M1 + m1] = ComplexVector::dotProduct(jV, jRm0);
			std::cout<<mu2D[m0*M1 + m1]<<std::endl;
			jLm0 = a*(H*jLm1) + b*jLm0;
			jLm0 = shift*jLm1 + jLm0;
			swap(jLm1, jLm0);
		}
		jLm0 = a*(H*jLm1) + b*jLm0;
		jLm0 = shift*jLm1 + jLm0;
		swap(jLm1, jLm0);
	}

	//postprocess moments
};

//As an example, let's compute the chebyshev moments of the one dimensional linear chain.

int main(int argc, char **argv)
{
	//Initialize TBTK.
	Initialize();

	Timer::tick("Full calculation");
	//Lattice size
	const int SIZE_X = 1000000;

	//Parameters
	double t = 1.0;

	Timer::tick("Create Model");
	//Create model and set up hopping parameters
	Model model;
	for(int x = 0; x < SIZE_X; x++) 		//Add hopping parameters corresponding to t
		if(x+1 < SIZE_X)
			model << HoppingAmplitude(-t,{x+1,0,0},{x,0,0}) + HC;

	//Construct model
	model.construct();
	Timer::tock();

	Timer::tick("Create Geometry");
	//Set up the geometry.
	Geometry &geometry = model.getGeometry();
	for(int x = 0; x < SIZE_X; x++)
		geometry.setCoordinate({x,0,0}, {(double)x,0,0} );
	Timer::tock();

	Timer::tick("Create velocity operator");
	//Generate the operator. The second argument is the direction of
	//interest.
	std::cout<<"Computing velocity operator"<<std::endl;
//	IndexedDataTree<complex<double>> customOperator = generateVelocityOperator(model, {1, 0, 0});//CHECK AND WORK
	vector<tuple<complex<double>, Index>> customOperator = generateVelocityOperator(model, {1, 0, 0});//CHECK AND WORK
	Timer::tock();

	Timer::tick("Calculate");
	//Run the calculation.
	const double BAND_WIDTH  = 8.0*t;
	const double BAND_CENTER = 0.0*t;
	const int M0=10;
	const int M1=10;
	std::cout<<"Computing Chebyshev moments"<<std::endl;
	kuboCalculation(model, customOperator, BAND_WIDTH, BAND_CENTER, M0, M1);
	Timer::tock();

	return 0;
}
