// ================================================================================================
//
// This main file only provides a small example of how to use M2S2, a Header-only library to use
// Matrices for Mechanices of Solids and Structures together with MKL Pardiso to solve a system of
// equations.
// 
// ================================================================================================

// Include the M2S2 library. This is the only required file that must be included.
#include "M2S2.h"

// Include your solver, in this case, Pardiso (intel copyrights)
#include "mkl_pardiso.h"

// Std libraries, if any
#include <iostream>

// Unit tests
#include "UnitTest.h"

int main()
{
	std::cout << INFO("Copyright(C) 2024 Dorival Piedade Neto & Rodrigo Ribeiro Paccola & Rogerio Carrazedo") << std::endl
		<< INFO("All Rights Reserved.") << std::endl
		<< INFO("Structural Engineering Department.") << std::endl
		<< INFO("University of Sao Paulo at Sao Carlos School of Engineering.") << std::endl << std::endl
		<< WARN("This program is free: you can redistribute it under the terms of the License.")
		<< WARN("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;")
		<< WARN("Without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.")
		<< WARN("It is provided \"AS IS\".") << std::endl
		<< INFO("In no event shall the authors be liable for any claim, damages or other liability,") << std::endl
		<< INFO("whether in an action of contract, tort or otherwire, arising from, out of or in ") << std::endl
		<< INFO("connection with the software or the use of other dealing in the software.") << std::endl << std::endl
		<< INFO("Neither the name of the copyright holder nor the names of any other contributors") << std::endl
		<< INFO("may be used to endorse or promote products derived from this software without") << std::endl
		<< INFO("specific prior written permission.") << std::endl << std::endl
		<< INFO("If you are using, please cite it in your reseach. We have a DOI: ") << std::endl
		<< INFO("Well, not yet :)") << std::endl << std::endl;

	std::cout << INPUT("\nClick a button to continue!") << std::endl;
	//std::cin.get();

	M2S2::unitTest();

	/*
	Here is what we are trying to solve (notice that it is symmetric):
	| 10  1   2   0 | a   | 1 |
	| 1	  9   11  8 | b = | 1 |
	| 2	  11  14  0 | c   | 1 |
	| 0	  8   0   1 | d   | 1 |

	Solution is given by:
	a = 0.0953811; b = 0.121478; c = -0.0376443; d = 0.0281755;
	*/

	// Two dyadics, summed accordingly, results in the matrix we seek
	M2S2::Dyadic2S mi_t1({ 10, 1, 2, 3, 4, 5 });
	M2S2::Dyadic2S mi_t2({ 6, 7, 8, 9, 0, 1 });

	// The number of DOF is given by
	const int size = 4;

	// We must create a sparseMatrix and then send to a solver
	// Here, we initialize with size 4 and symmetric
	M2S2::sparseMatrix mi_spMatrix(size, true);

	// An initial bandwidth is required (just a hint, used to initial allocation)
	mi_spMatrix.setAllLineSize(size - 1);

	// We can push the dyadic to the sparse matrix as triplets, as MatrixS, or MatrixX
	// But since we are have a symmetric system of equations, MatrixX is not an alternative

	// ------------------------------------------------------------------------------------
	// Pushing the dyadics as triplets
	// M2S2 stores the matrices in Voigt notation using row major order
	// Notice that each triplet is composed by row index, col index, and value
	// Also notice that we pushed only the upper triangular part of the symmetric matrix
	std::vector<M2S2::triplet> mi_t3;

	// First dyadic -> vector of triplets
	mi_t3.emplace_back(0, 0, mi_t1.getVector().at(0));		// t1(0,0)
	mi_t3.emplace_back(0, 1, mi_t1.getVector().at(1));		// t1(0,1)
	mi_t3.emplace_back(0, 2, mi_t1.getVector().at(2));		// t1(0,2)
	mi_t3.emplace_back(1, 1, mi_t1.getVector().at(3));		// t1(1,1)
	mi_t3.emplace_back(1, 2, mi_t1.getVector().at(4));		// t1(1,2)
	mi_t3.emplace_back(2, 2, mi_t1.getVector().at(5));		// t1(2,2)

	// Push triplet to the sparse matrix
	mi_spMatrix.push(mi_t3);
	// ------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------
	// Second dyadic -> MatrixS and corresponding indexes
	M2S2::MatrixS mi_t4(mi_t2.rows(), mi_t2.getVector());
	// Dof of t4 in the sparse Matrix
	std::vector<int> mi_t4_index({ 1, 2, 3 });

	// Push MatrixS to the sparse matrix
	mi_spMatrix.push(mi_t4, mi_t4_index);
	// ------------------------------------------------------------------------------------

	// After getting every contribution, we add equal terms, reducing the size of the matrix
	mi_spMatrix.addEqualTerms();

	// Pardiso uses CSR matrices, thus needs conversion
	M2S2::CSR mi_csrMatrix;
	mi_spMatrix.saveAsCSR(mi_csrMatrix);

	// Creating the Right and Left hand side
	std::vector<double> mi_RHS({ 1, 1, 1, 1 });
	std::vector<double> mi_LHS(4);

	// ------------------------------------------------------------------------------------
	// Solving using pardiso - this was copied from Intel documentation

	int mtype = -2;       /* Real symmetric matrix */
	int nrhs = 1;         /* Number of right hand sides. */

	int idum;             /* Integer dummy. */
	double ddum;          /* Double dummy */

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void* pt[64];

	/* Pardiso control parameters. */
	int iparm[64];
	int maxfct, mnum, phase, error, msglvl;

	/* -------------------------------------*/
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------*/
	for (int i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[34] = 1;        /* PARDISO use C-style indexing for ia and ja arrays */

	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;             /* Which factorization to use. */
	msglvl = 0;           /* Do not print statistical information */
	error = 0;            /* Initialize error flag */

	/* ----------------------------------------------------------------*/
	/* .. Initialize the internal solver memory pointer. This is only  */
	/*   necessary for the FIRST call of the PARDISO solver.           */
	/* ----------------------------------------------------------------*/
	for (int i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}

	/* --------------------------------------------------------------------*/
	/* .. Reordering and Symbolic Factorization. This step also allocates  */
	/*    all memory that is necessary for the factorization.              */
	/* --------------------------------------------------------------------*/
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, mi_csrMatrix.mv_value.data(), mi_csrMatrix.mv_rowIndex.data(), mi_csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		std::cout << "\nERROR during symbolic factorization: " << error;
		std::cout << "\nClick a button to exit!";
		std::cin.get();
		exit(1);
	}
	std::cout << "\nReordering completed ... ";
	std::cout << "\nNumber of nonzeros in factors = " << iparm[17];
	std::cout << "\nNumber of factorization MFLOPS = " << iparm[18];

	/* ----------------------------*/
	/* .. Numerical factorization. */
	/* ----------------------------*/
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, mi_csrMatrix.mv_value.data(), mi_csrMatrix.mv_rowIndex.data(), mi_csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		std::cout << "\nERROR during numerical factorization: " << error;
		std::cout << "\n\nClick a button to exit!";
		std::cin.get();
		exit(2);
	}
	std::cout << "\nFactorization completed ... ";

	/* -----------------------------------------------*/
	/* .. Back substitution and iterative refinement. */
	/* -----------------------------------------------*/
	phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */

	/* Set right hand side to one. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, mi_csrMatrix.mv_value.data(), mi_csrMatrix.mv_rowIndex.data(), mi_csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, mi_RHS.data(), mi_LHS.data(), &error);
	if (error != 0)
	{
		std::cout << "\nERROR during solution: " << error;
		std::cout << "\n\nClick a button to exit!";
		std::cin.get();
		exit(3);
	}
	std::cout << "\nSolve completed ... ";
	std::cout << "\nThe solution of the system is: ";
	for (int i = 0; i < size; i++)
	{
		std::cout << "\n LHS [" << i << "] = " << mi_LHS.at(i);
	}
	std::cout << "\n";

	/* --------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------*/
	phase = -1;           /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, &ddum, mi_csrMatrix.mv_rowIndex.data(), mi_csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	// ------------------------------------------------------------------------------------

	std::cout << "\nExpected resuts are:";
	std::cout << "\n LHS [0] = 0.0953811";
	std::cout << "\n LHS [1] = 0.121478";
	std::cout << "\n LHS [2] = -0.0376443";
	std::cout << "\n LHS [3] = 0.0281755";

	std::cout << "\n\nClick a button to exit!";
	std::cin.get();

	return 0;
}


// Just don't forget:
// 
// std::invalid_argument 
// an argument value has not been accepted.
//
// std::domain_error
// the inputs are outside of the domain on which an operation is defined
//
// std::length_error
// attempts to exceed implementation defined length limits for some object.
//
// std::out_of_range
// attempt to access elements out of defined range.
//
// std::runtime_error
// due to events beyond the scope of the program and can not be easily predicted.
//
// std::range_error
// used to report range errors (that is, situations where a result of a computation cannot be represented by the destination type).
//
// std::overflow_error
// used to report arithmetic overflow errors (that is, situations where a result of a computation is too large for the destination type).
//
// std::underflow_error
// used to report arithmetic underflow errors (that is, situations where the result of a computation is a subnormal floating - point value).