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

IndexedDataTree<complex<double>> generateVelocityOperator(
	const Model &model,
	const Vector3d &direction
){
	//Get the Geometry from the model.
	const Geometry *geometry = model.getGeometry();
	TBTKAssert(
		geometry != nullptr,
		"generateVelocityOperator()",
		"The model contains no Geometry.",
		"Use model.createGeometry() to create a geometry and"
		<< " geometry->setCoordinates() to set the coordinates."
	);
	TBTKAssert(
		geometry->getDimensions() == 3,
		"generateVelocityOperator()",
		"The geometry dimension must have dimension '3' but has"
		<< " dimension '" << geometry->getDimensions() << "'.",
		""
	);

	//Create the velocity operator.
	IndexedDataTree<complex<double>> velocityOperator;

	//Iterate over all HoppingAmplitudes.
	const HoppingAmplitudeSet &hoppingAmplitudeSet
		= model.getHoppingAmplitudeSet();
	for(
		HoppingAmplitudeSet::ConstIterator iterator
			= hoppingAmplitudeSet.cbegin();
		iterator != hoppingAmplitudeSet.cend();
		++iterator
	){
		//Get the amplitude, to-, and from-Indices.
		complex<double> amplitude = (*iterator).getAmplitude();
		const Index to = (*iterator).getToIndex();
		const Index from = (*iterator).getFromIndex();

		//Get the coordinates corresponding to the to- and from-Indices
		//and vonvert them to Vector3d.
		const double *toCoordinates = geometry->getCoordinates(to);
		Vector3d toR = Vector3d(
			{toCoordinates[0], toCoordinates[1], toCoordinates[2]}
		);
		const double *fromCoordinates = geometry->getCoordinates(from);
		Vector3d fromR = Vector3d(
			{fromCoordinates[0], fromCoordinates[1], fromCoordinates[2]}
		);

		//Calculate H_{ij}(direction.(R_i - R_j)).
		complex<double> operatorValue = amplitude*Vector3d::dotProduct(
			direction,
			toR - fromR
		);

		//Avoid adding elements that are identically zero.
		if(fabs(operatorValue) == 0)
			continue;

		//Add the value to the velocity operator.
		try{
			//Try to add it to an already existing entry.
			velocityOperator.get({to, from}) += operatorValue;
		}
		catch(ElementNotFoundException e){
			//Create a new entry if it fails.
			velocityOperator.add(operatorValue, {to, from});
		}
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
	const IndexedDataTree<complex<double>> &indexedDataTree,
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
		IndexedDataTree<complex<double>>::ConstIterator iterator
			= indexedDataTree.cbegin();
		iterator != indexedDataTree.cend();
		++iterator
	){
		//Get the Index of the current element.
		const Index &index = iterator.getCurrentIndex();

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
		sparseMatrix.add(row, column, *iterator);
	}

	//Construct the SparseMatrix internal CSC representation.
	sparseMatrix.construct();

	return sparseMatrix;
}

KuboSparseMatrix convertIndexedDataTreeToKuboSparseMatrix(
	const IndexedDataTree<complex<double>> &indexedDataTree,
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
	const IndexedDataTree<complex<double>> targetOperator,
	double bandWidth,
	double bandCenter,
	int M0,
	int M1
){
	//Convert the Model into a KuboSparseMatrix.
	KuboSparseMatrix H = convertModelToKuboSparseMatrix(model);
	H.rescale(2.0/bandWidth);

	//Extract the basis size.
	const int basisSize = model.getBasisSize();

	//Convert targetOperator from IndexedDataTree to KuboSparseMatrix.
	KuboSparseMatrix Op = convertIndexedDataTreeToKuboSparseMatrix(
		targetOperator,
		model
	);

	//TODO: Vx should be calculated independently of Op.
	KuboSparseMatrix Vx = Op;

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
	for(int m0 = 0; m0 < M0; m0++){
		const complex<double > a=2.0, b=1.0;
		jRm0 = a*(H*jRm1) + b*jRm0;
		jRm0 = shift*jRm1 + jRm0;	//Can be optimized by implementing operator+=().
		swap(jRm1, jRm0);
		jRm0 = randomPhaseVector;
		jRm1 = H*jRm0;
		for(int m1 = 0; m1 < M1; m1++){
			jLm1 = H*jLm0;
			jLm0 = a*(H*jLm1) + b*jLm0;
			jLm0 = shift*jLm1 + jLm0;	//Can be optimized by implementing operator+=().
			swap(jLm1, jLm0);
			jV = Op*jLm0;
			mu2D[m0*M1 + m1] = ComplexVector::dotProduct(jV, jRm0);
		}
	}

	//postprocess moments
};

int main(int argc, char **argv){
	//Lattice size
	const int SIZE_X = 20;
	const int SIZE_Y = 20;

	//Parameters
	complex<double> mu = 0.0;
	complex<double> t = 1.0;

	//Create model and set up hopping parameters
	Model model;
	for(int x = 0; x < SIZE_X; x++){
		for(int y = 0; y < SIZE_Y; y++){
			for(int s = 0; s < 2; s++){
				//Add hopping amplitudes corresponding to chemical potential
				model << HoppingAmplitude(
					 mu,
					{x, y, s},
					{x, y, s}
				);

				//Add hopping parameters corresponding to t
				if(x+1 < SIZE_X){
					model << HoppingAmplitude(
						-t,
						{(x+1)%SIZE_X,	y,	s}, // periodic boundary conditions
						{x,		y,	s}
					) + HC;
				}
				if(y+1 < SIZE_Y){
					model << HoppingAmplitude(
						-t,
						{x,	(y+1)%SIZE_Y,	s}, // periodic boundary conditions
						{x,	y,		s}
					) + HC;
				}
			}
		}
	}
	//Construct model
	model.construct();

	//Set up the geometry.
	model.createGeometry(3);
	Geometry *geometry = model.getGeometry();
	for(int x = 0; x < SIZE_X; x++){
		for(int y = 0; y < SIZE_Y; y++){
			for(int s = 0; s < 2; s++){
				geometry->setCoordinates(
					{x, y, s},
					{(double)x, (double)y, 0}
				);
			}
		}
	}

	//Generate the operator. The second argument is the direction of
	//interest.
	IndexedDataTree<complex<double>> customOperator
		= generateVelocityOperator(model, {1, 0, 0});

	//Run the calculation.
	const double BAND_WIDTH = 5;
	const double BAND_CENTER = 0;
	const int M0=1;
	const int M1=1;
	kuboCalculation(model, customOperator, BAND_WIDTH, BAND_CENTER, M0, M1);

	return 0;
}
