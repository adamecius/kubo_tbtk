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

#include "ComplexVector.h"
#include "KuboSparseMatrix.h"

#include <complex>

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);

void swap(
	ComplexVector& x,
	ComplexVector& y
){
	ComplexVector temp = std::move(x);
	x = std::move(y);
	y = std::move(temp);
}; //swap vectors by reference and not by copy

//The velocity operator is essential to compute any nonequilibrium property.
//In a real-space tight-binding formulation it can be easily computed as V = i [ H, X] = H_{ij} ( R_i - R_j );
class velocity_operator
{
	std::vector<double> get_from_to_differences(){
	}; //Get the difference from sites inital sities to final sites. Without periodic boundary conditions

	std::vector<complex<double> > get_from_to_amplitudes(){
	}; //Get the amplitudes from sites inital sities to final sites. Without periodic boundary conditions

	std::vector<int> get_from_indexes(){
	};

	std::vector<int> get_to_indexes(){
	};

	int getBasisSize(){
	};

//	sparse_matrix compute_velocity_operator(){
	KuboSparseMatrix compute_velocity_operator(){
		int    basisSize =  getBasisSize();
		std::vector<int>    i_idx =  get_from_indexes();
		std::vector<int>    j_idx =  get_to_indexes();
		std::vector<double> Rij = get_from_to_differences();
		std::vector<complex<double> >  Vij = get_from_to_amplitudes();

		for( int n = 0; n < Vij.size(); n++ )
			Vij[n] = i*Vij[n]*Rij[n];

		//Convert this data into a SparseMatrix on CSC format
		SparseMatrix<complex<double>> VijmatSparseMatrix(
			SparseMatrix<complex<double>>::StorageFormat::CSC,
			basisSize,
			basisSize
		);
		//...
		KuboSparseMatrix Vijmat(VijmatSparseMatrix);

		return Vijmat;
	};
};

SparseMatrix<complex<double>> convertIndexedDataTreeToSparseMatrix(
	const IndexedDataTree<complex<double>> &indexedDataTree,
	const Model &model
){
	const HoppingAmplitudeSet &hoppingAmplitudeSet
		= model.getHoppingAmplitudeSet();

	SparseMatrix<complex<double>> sparseMatrix(
		SparseMatrix<complex<double>>::StorageFormat::CSC,
		model.getBasisSize(),
		model.getBasisSize()
	);
	for(
		IndexedDataTree<complex<double>>::ConstIterator iterator
			= indexedDataTree.cbegin();
		iterator != indexedDataTree.cend();
		++iterator
	){
		const Index &index = iterator.getCurrentIndex();
		vector<Index> components = index.split();
		TBTKAssert(
			components.size() == 2,
			"convertOperator()",
			"Invalid operator index '" << index.toString() << "'."
			<< " The Index must have two component indices, but"
			<< " have '" << components.size() << "'.",
			""
		);
		int row = hoppingAmplitudeSet.getBasisIndex(components[0]);
		int column = hoppingAmplitudeSet.getBasisIndex(components[1]);

		sparseMatrix.add(row, column, *iterator);
	}

	sparseMatrix.construct();

	return sparseMatrix;
}

KuboSparseMatrix convertIndexedDataTreeToKuboSparseMatrix(
	const IndexedDataTree<complex<double>> &indexedDataTree,
	const Model &model
){
	SparseMatrix<complex<double>> sparseMatrix = convertIndexedDataTreeToSparseMatrix(
		indexedDataTree,
		model
	);

	return KuboSparseMatrix(sparseMatrix);
}

ComplexVector generateRandomPhaseVector(unsigned int basisSize){
	ComplexVector randomPhaseVector(basisSize);
	for(int n = 0; n < basisSize; n++){
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
	SparseMatrix<complex<double>> sparseHamiltonian
		= model.getHoppingAmplitudeSet().getSparseMatrix();
	sparseHamiltonian.setStorageFormat(
		SparseMatrix<complex<double>>::StorageFormat::CSC
	);
	KuboSparseMatrix H(sparseHamiltonian);

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

	//With this random phase approximation, the code goes as follow
	H.self_rescaling( 2.0/bandWidth);

	//Set somewhere the number of coefficients. For this there are two coefficients so one defines a matrix M0xM1
	std::vector < complex<double> > mu2D(M0*M1);

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
	for(int m0= 0; m0 < M0; m0++)
	{
		const complex<double > a=2.0, b=1.0;
		jRm0 = a*(H*jRm1) + b*jRm0;
		jRm0 = shift*jRm1 + jRm0;	//Can be optimized by implementing operator+=().
		swap(jRm1, jRm0);
		jRm0 = randomPhaseVector;
		jRm1 = H*jRm0;
		for(int m1 = 0; m1 < M1; m1++)
		{
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
	const int SIZE_X = 5;
	const int SIZE_Y = 5;

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

	//TODO: Generate a proper operator here.
	IndexedDataTree<complex<double>> customOperator;
	customOperator.add(1, {{1, 1, 0}, {2, 1, 0}});

	const double BAND_WIDTH = 5;
	const double BAND_CENTER = 0;
	const int M0=1;
	const int M1=1;
	kuboCalculation(model, customOperator, BAND_WIDTH, BAND_CENTER, M0, M1);

	return 0;
}
