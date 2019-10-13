/* Copyright 2016 Kristofer Björnson
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

/** @package TBTKtemp
 *  @file main.cpp
 *  @brief Basic diagonalization example
 *
 *  Basic example of diagonalization of a 2D tight-binding model with t = 1 and
 *  mu = -1. Lattice with edges and a size of 20x20 sites.
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/FileWriter.h"
#include "TBTK/Model.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
#include <complex>

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);

complex<double> dot(const std::vector<complex<double>>& x,const std::vector<complex<double>>& y){ return 0;}; //swap vectors by reference and not by copy
void swap( std::vector<complex<double>>& x,std::vector<complex<double>>& y){}; //swap vectors by reference and not by copy
void axpy( const complex<double>& a ,std::vector<complex<double>>& x,std::vector<complex<double>>& y){}; //PErform the operation y= ax+ y

//A sparse matrix class needed to perform matrix-vector multiplication
class sparse_matrix{
	typedef complex<double> my_complex;
	typedef std::vector<complex<double> > my_complex_vec;
	
	public: 
	sparse_matrix(){};
	
	sparse_matrix(const  std::vector<int> rows,const  std::vector<int> cols, const my_complex_vec vals, const int dim0, const int dim1 ){};
	
	void spMVmult( const my_complex_vec& x , const my_complex_vec& y){}; //return y = A x 

	void axpby(const my_complex &a,const my_complex_vec& x ,const my_complex& b, my_complex_vec& y){}; //return y = a 1A x + b y 

	void self_rescaling(const my_complex a){}; // perform the transformation A ---> a A

	void self_rescaling(const double a){}; // perform the transformation A ---> a A
};

//The velocity operator is essential to compute any nonequilibrium property.
//In a real-space tight-binding formulation it can be easily computed as V = i [ H, X] = H_{ij} ( R_i - R_j );
class velocity_operator
{
	std::vector<double> get_from_to_differences(){}; //Get the difference from sites inital sities to final sites. Without periodic boundary conditions
	
	std::vector<complex<double> > get_from_to_amplitudes(){}; //Get the amplitudes from sites inital sities to final sites. Without periodic boundary conditions

	std::vector<int> get_from_indexes(){}; 

	std::vector<int> get_to_indexes(){}; 

	int get_num_orbitals(){}; 

	sparse_matrix compute_velocity_operator(){
		int    num_orbs =  get_num_orbitals();
		std::vector<int>    i_idx =  get_from_indexes();
		std::vector<int>    j_idx =  get_to_indexes();
		std::vector<double> Rij = get_from_to_differences();
		std::vector<complex<double> >  Vij = get_from_to_amplitudes();

		complex<double> I(0,1);
		for( int i = 0; i < Vij.size(); i++ )
			Vij[i] = I*Vij[i]*Rij[i];
		
		//Convert this data into a COO matrix format
		sparse_matrix Vijmat( i_idx, j_idx, Vij, num_orbs ,num_orbs);
		
		return Vijmat;}; 

};


void kuboCalculation(const Model &model)
{
	
	//Convert model's hamiltoniano to sparse matrix;
	sparse_matrix H, Vx, Op; //The electric field is assume for simplicity in x, therefore Vx
	//For the conductivity Op = Vx, for other quantities Op is an arbitrary quantum mechanical operator

	//Define here trace solver. Meaning Tr[ A]. For large systems, stochastic trace approximation is needed. 
	const int num_orbs=1; //The number of sites in the models (including spin and internal degrees of freedom)
	std::vector<complex<double> >rphase_vec( num_orbs);
	complex<double> I(0,1);
	for( int i = 0; i <   num_orbs; i++ )
		rphase_vec[i] = exp( 2.0*M_PI*I * (double)rand()/(double)RAND_MAX )/sqrt( num_orbs );//With this definition of |phi>, is easy to show that <phi|A|phi> = Tr[A] + random_noise which dissapear for num_orbs--> infty

	//With this random phase approximation, the code goes as follow
	double bandWidth,bandCenter; //Extract somehow the BandWidth and BandCenter. Or set it
	H.self_rescaling( 2.0/bandWidth);
	
	//Set somewhere the number of coefficients. For this there are two coefficients so one defines a matrix M0xM1
	const int M0=1, M1=1;
	std::vector < complex<double> > mu2D( M0*M1);


	//Need five vectors for performing the iteration
	std::vector<complex<double> >  jLm0( num_orbs),jLm1( num_orbs),jRm0( num_orbs ),jRm1( num_orbs ),  jV( num_orbs );
	
	//Starts chebyshev iteration
	Vx.spMVmult(rphase_vec, jLm0);
	H.spMVmult(jLm0, jLm1);
	const double shift= -2*bandCenter/bandWidth;
	for(int  m0= 0; m0 <M0 ; m0++ )
	{
		const complex<double > a=2.0, b  = 1.0;
		H.axpby(a,jRm1,b,jRm0);
		axpy(shift, jRm1 , jRm0);
		swap ( jRm1, jRm0 );
		
		jRm0 = rphase_vec;
		H.spMVmult(jRm0, jRm1);
		for(int  m1= 0; m1 <M1 ; m1++ )
		{
			H.spMVmult(jLm0, jLm1);
			H.axpby(a,jLm1,b,jLm0);
			axpy(shift, jLm1 , jLm0);
			swap ( jLm1, jLm0 );
			Op.spMVmult(jLm0, jV);
			mu2D[ m0*M1 + m1 ] = dot( jV,jRm0 );
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
	
	kuboCalculation(model);

    return 0;
}

