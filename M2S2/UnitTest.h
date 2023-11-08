// ================================================================================================
// 
// This file is part of M2S2 - Matrices for Mechanices of Solids and Structures
//
// Copyright(C) 2023 
//		Dorival Piedade Neto &
//		Rodrigo Ribeiro Paccola &
//		Rogério Carrazedo
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// ================================================================================================
#pragma once

// Include the M2S2 library.
#include "M2S2.h"

namespace M2S2 {
	// ================================================================================================
	//
	// Unit testing - Checks every function from all classes
	// 
	// ================================================================================================
	//
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
	//
	// ================================================================================================
	void unitTest_Dyadis2S() {

		bool mi_check = true;

		try {
			std::cout << std::endl;
			std::cout << WARN("Tests every aspects of M2S2 library!");

			// ----------------------------------------------------------------------------------------
			// Dyadic2S tests
			std::cout << WARN("Beggining Dyadic2S Constructors evaluation!");

			// Constructors
			{
				Dyadic2S mi_d2S;		// Basic constructor
				if (mi_d2S.rows() != 0) mi_check = false;
				if (mi_d2S.size() != 0) mi_check = false;
				if (mi_d2S.getVector().size() != 0)  mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S Basic Constructor Failed!");
			}

			{
				Dyadic2S mi_d2S(2);	// Constructor of 2D space (2x2);
				if (mi_d2S.rows() != 2) mi_check = false;
				if (mi_d2S.size() != 3) mi_check = false;
				if (mi_d2S.getVector().size() != 3) mi_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != 0.) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S 2D Constructor Failed!");
			}

			{
				Dyadic2S mi_d2S(3);	// Constructor of 3D space (3x3);
				if (mi_d2S.rows() != 3) mi_check = false;
				if (mi_d2S.size() != 6) mi_check = false;
				if (mi_d2S.getVector().size() != 6) mi_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != 0.) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S 3D Constructor Failed!");
			}

			{
				std::vector<double> mi_Voigt2D({ 0,  1., 2. });
				Dyadic2S mi_d2S(mi_Voigt2D);	// Constructor from vector (2D space);
				if (mi_d2S.rows() != 2) mi_check = false;
				if (mi_d2S.size() != 3) mi_check = false;
				if (mi_d2S.getVector().size() != 3) mi_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != mi_Voigt2D.at(i)) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S 2D Constructor from Vector Input Failed!");
			}

			{
				std::vector<double> mi_Voigt3D({ 0., 1., 2., 3., 4., 5. });
				Dyadic2S mi_d2S(mi_Voigt3D);	// Constructor from vector (3D space);
				if (mi_d2S.rows() != 3) mi_check = false;
				if (mi_d2S.size() != 6) mi_check = false;
				if (mi_d2S.getVector().size() != 6) mi_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != mi_Voigt3D.at(i)) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S 3D Constructor from Vector Input Failed!");
			}

			{
				// Check if throw in constructor is working 
				// The vector has wrong dimensions
				try {
					mi_check = false;
					std::vector<double> mi_NonVoigt({ 0., 1., 2., 3., 4. });
					Dyadic2S mi_d2S(mi_NonVoigt);
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}
			}

			{
				// Copy constructor
				Dyadic2S mi_d2S_1({ 0,  1., 2. });
				Dyadic2S mi_d2S_2(mi_d2S_1);
				if (mi_d2S_2.rows() != 2) mi_check = false;
				if (mi_d2S_2.size() != 3) mi_check = false;
				if (mi_d2S_2.getVector().size() != 3) mi_check = false;
				for (int i = 0; i < mi_d2S_2.getVector().size(); ++i) {
					if (mi_d2S_2.getVector().at(i) != (double)i) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S Copy Constructor Failed!");
			}

			{
				// Explicit move constructor
				Dyadic2S mi_d2S_1({ 0., 1., 2., 3., 4., 5. });
				Dyadic2S mi_d2S_2 = std::move(mi_d2S_1);
				if (mi_d2S_2.rows() != 3) mi_check = false;
				if (mi_d2S_2.size() != 6) mi_check = false;
				if (mi_d2S_2.getVector().size() != 6) mi_check = false;
				for (int i = 0; i < mi_d2S_2.getVector().size(); ++i) {
					if (mi_d2S_2.getVector().at(i) != (double)i) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S 3D Explicit Move Constructor Failed!");
			}

			{
				// Implicit move constructor
				Dyadic2S mi_d2S(Dyadic2S::identity(3));
				if (mi_d2S.rows() != 3) mi_check = false;
				if (mi_d2S.size() != 6) mi_check = false;
				if (mi_d2S.getVector().size() != 6) mi_check = false;
				for (int i = 0; i < mi_d2S.rows(); ++i) {
					for (int j = i; j < mi_d2S.cols(); ++j) {
						if (i == j) {
							if (mi_d2S.at(i, j) != 1.) mi_check = false;
						}
						if (i != j) {
							if (mi_d2S.at(i, j) != 0.) mi_check = false;
						}
					}
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S 3D Implicit Move Constructor Failed!");
			}

			{
				// Input by stream
				Dyadic2S mi_d2S(2);
				std::stringstream str("0. 1. 2.");
				str >> mi_d2S;
				if (mi_d2S.rows() != 2) mi_check = false;
				if (mi_d2S.size() != 3) mi_check = false;
				if (mi_d2S.getVector().size() != 3) mi_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != (double)i) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - Value input by stream Failed!");
			}

			{
				// Input by stream
				Dyadic2S mi_d2S(3);
				std::stringstream str("0. 1. 2. 3. 4. 5.");
				str >> mi_d2S;
				if (mi_d2S.rows() != 3) mi_check = false;
				if (mi_d2S.size() != 6) mi_check = false;
				if (mi_d2S.getVector().size() != 6) mi_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != (double)i) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - Value input by stream Failed!");
			}

			{
				// Check if throw in constructor is working
				// The stream is unable to fill the dyadic
				try {
					mi_check = false;

					Dyadic2S mi_d2S(3);
					std::stringstream str("0. 1. 2.");
					str >> mi_d2S;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}
			}

			{
				// Check if throw in constructor is working
				// The dyadic size is smaller than stream
				try {
					mi_check = false;

					Dyadic2S mi_d2S;
					std::stringstream str("0. 1. 2.");
					str >> mi_d2S;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_check) {
				std::cout << OK("Dyadic2S constructors checked!\n\n");
				std::cout << WARN("Beggining Dyadic2S Methods evaluation!");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2S constructor error!"));
			}

			// ----------------------------------------------------------------------------------------
			{
				// Swap method
				Dyadic2S mi_d2S_1({ 0,  1., 2. });
				Dyadic2S mi_d2S_2({ 0., 1., 2., 3., 4., 5. });

				mi_d2S_2.swap(mi_d2S_1);
				if (mi_d2S_1.rows() != 3) mi_check = false;
				if (mi_d2S_1.size() != 6) mi_check = false;
				if (mi_d2S_1.getVector().size() != 6) mi_check = false;
				for (int i = 0; i < mi_d2S_1.getVector().size(); ++i) {
					if (mi_d2S_1.getVector().at(i) != (double)i) mi_check = false;
				}

				if (mi_d2S_2.rows() != 2) mi_check = false;
				if (mi_d2S_2.size() != 3) mi_check = false;
				if (mi_d2S_2.getVector().size() != 3) mi_check = false;
				for (int i = 0; i < mi_d2S_2.getVector().size(); ++i) {
					if (mi_d2S_2.getVector().at(i) != (double)i) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - Swap Method Failed!!");
			}

			{
				// Clear method
				Dyadic2S mi_d2S({ 0,  1., 2. });

				mi_d2S.clear();
				if (mi_d2S.rows() != 2) mi_check = false;
				if (mi_d2S.size() != 3) mi_check = false;
				if (mi_d2S.getVector().size() != 3) mi_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != 0.) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - Clear Method Failed!!");
			}

			{
				// getVoigtMnemonics method
				Dyadic2S mi_d2S_1({ 0,  1., 2. });
				Dyadic2S mi_d2S_2({ 0., 1., 2., 3., 4., 5. });

				auto mi_Voigt = mi_d2S_1.getVoigtMnemonics();
				if (mi_Voigt.size() != mi_d2S_1.size()) mi_check = false;
				for (int i = 0; i < mi_d2S_1.rows(); ++i) {
					if ((int)(mi_Voigt.at(i) - mi_d2S_1.at(i, i)) * 1000000 != 0) mi_check = false;

					for (unsigned int j = i + 1; j < mi_d2S_1.cols(); ++j)
						if ((int)(mi_Voigt.at(mi_Voigt.size() - i - j) - mi_d2S_1.at(i, j)) * 1000000 != 0) mi_check = false;
				}

				mi_Voigt = mi_d2S_2.getVoigtMnemonics();
				if (mi_Voigt.size() != mi_d2S_2.size()) mi_check = false;
				for (int i = 0; i < mi_d2S_2.rows(); ++i) {
					if ((int)(mi_Voigt.at(i) - mi_d2S_2.at(i, i)) * 1000000 != 0) mi_check = false;

					for (unsigned int j = i + 1; j < mi_d2S_2.cols(); ++j)
						if ((int)(mi_Voigt.at(mi_Voigt.size() - i - j) - mi_d2S_2.at(i, j)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - getVoigtMnemonics Method Failed!!");
			}

			{
				// Random symmetric matrices
				Dyadic2S mi_d2S_1 = Dyadic2S({ 2.5, 4.8, -3.6 });
				Dyadic2S mi_d2S_2 = Dyadic2S({ 2.5, 4.8, -3.6, 6.7, -1.4, 8.1 });

				// Determinant
				if ((int)(mi_d2S_1.determinant() - (-32.04)) * 1000000 != 0) mi_check = false;
				if ((int)(mi_d2S_2.determinant() - (-94.297)) * 1000000 != 0) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S - determinant Method Failed!!");

				// Norm
				if ((int)(mi_d2S_1.norm() - (7.430343195)) * 1000000 != 0) mi_check = false;
				if ((int)(mi_d2S_2.norm() - (16.34533572)) * 1000000 != 0) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S - norm Method Failed!!");

				// Trace
				if ((int)(mi_d2S_1.trace() - (-1.1)) * 1000000 != 0) mi_check = false;
				if ((int)(mi_d2S_2.trace() - (17.3)) * 1000000 != 0) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S - trace Method Failed!!");
			}

			{
				// Check if throws in previous methods are working
				// The dyadic is empty
				try {
					mi_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.determinant();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}

				try {
					mi_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.norm();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}

				try {
					mi_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.trace();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}
			}

			{
				// Random symmetric matrices
				Dyadic2S mi_d2S_1 = Dyadic2S({ 2.5, 4.8, -3.6 });
				Dyadic2S mi_d2S_2 = Dyadic2S({ 2.5, 4.8, -3.6, 6.7, -1.4, 8.1 });

				// Eigenvalues
				std::vector<double> eval2D = mi_d2S_1.eigenvalues();
				std::vector<double> eval3D = mi_d2S_2.eigenvalues();
				std::vector<double> expected2D({ -6.237046685231272,  5.1370466852312724 });
				std::vector<double> expected3D({ 12.377474154257952, -1.2368787192025401, 6.1594045649445857 });

				for (int i = 0; i < eval2D.size(); ++i) {
					if ((int)(eval2D.at(i) - expected2D.at(i)) * 1000000 != 0) mi_check = false;
				}
				for (int i = 0; i < eval3D.size(); ++i) {
					if ((int)(eval3D.at(i) - expected3D.at(i)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - eigenvalue Method Failed!!");

				// Invariants
				eval2D = mi_d2S_1.invariants();
				eval3D = mi_d2S_2.invariants();
				expected2D = std::vector<double>{ -1.1, -32.04 };
				expected3D = std::vector<double>{ 17.3, 53.31, -94.297 };

				for (int i = 0; i < eval2D.size(); ++i) {
					if ((int)(eval2D.at(i) - expected2D.at(i)) * 1000000 != 0) mi_check = false;
				}
				for (int i = 0; i < eval3D.size(); ++i) {
					if ((int)(eval3D.at(i) - expected3D.at(i)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - invariants Method Failed!!");

				// Inverse
				Dyadic2S mi_d2S_3 = mi_d2S_1.inverse();
				Dyadic2S mi_d2S_4 = mi_d2S_2.inverse();

				Dyadic2S mi_d2S_5 = Dyadic2S({ 0.1123595506, 0.1498127341, -0.0780274657 });
				Dyadic2S mi_d2S_6 = Dyadic2S({ -0.5547366300, 0.3588661357, -0.1845233677, -0.0773089282, 0.1461340234, 0.0667041369 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if ((int)(mi_d2S_3.getVector().at(i) - mi_d2S_5.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if ((int)(mi_d2S_4.getVector().at(i) - mi_d2S_6.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - inverse Method Failed!!");

				// A^transpose * A
				mi_d2S_3 = mi_d2S_1.getATA();
				mi_d2S_4 = mi_d2S_2.getATA();
				mi_d2S_5 = M2S2::Dyadic2S({ 29.29, -5.28,  36.00 });
				mi_d2S_6 = M2S2::Dyadic2S({ 42.25, 49.20, -44.88, 69.89, -38.0, 80.53 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if ((int)(mi_d2S_3.getVector().at(i) - mi_d2S_5.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if ((int)(mi_d2S_4.getVector().at(i) - mi_d2S_6.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - getATA Method Failed!!");

				// Spheric part
				mi_d2S_3 = mi_d2S_1.getSpheric();
				mi_d2S_4 = mi_d2S_2.getSpheric();
				mi_d2S_5 = M2S2::Dyadic2S({ -0.366666667, 0., -0.366666667 });
				mi_d2S_6 = M2S2::Dyadic2S({ 5.766666667, 0., 0., 5.766666667, 0., 5.766666667 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if ((int)(mi_d2S_3.getVector().at(i) - mi_d2S_5.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if ((int)(mi_d2S_4.getVector().at(i) - mi_d2S_6.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - getSpheric Method Failed!!");

				// Deviatoric part
				mi_d2S_3 = mi_d2S_1.getDeviator();
				mi_d2S_4 = mi_d2S_2.getDeviator();
				mi_d2S_5 = M2S2::Dyadic2S({ 2.866666667, 4.8, -3.233333333 });
				mi_d2S_6 = M2S2::Dyadic2S({ -3.266666667, 4.8, -3.6, 0.9333333333, -1.4, 2.3333333333 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if ((int)(mi_d2S_3.getVector().at(i) - mi_d2S_5.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if ((int)(mi_d2S_4.getVector().at(i) - mi_d2S_6.getVector().at(i)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S - getDeviator Method Failed!!");

				// Double Dot product
				if ((int)(mi_d2S_1.contraction(mi_d2S_5) - 53.4866666651) * 1000000 != 0) mi_check = false;
				if ((int)(mi_d2S_2.contraction(mi_d2S_6) - 82.0599999982) * 1000000 != 0) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S - contraction Method Failed!!");
			}

			{
				// Check if throws in previous methods are working
				// The dyadic is empty
				try {
					mi_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.eigenvalues();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}

				// The dyadic is empty
				try {
					mi_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.inverse();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}

				// Singular Dyadic (non-invertible)
				try {
					mi_check = false;
					Dyadic2S mi_d2S({ 1.0, -4.8, 23.04 });
					mi_d2S.inverse();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}

				// The dyadic is empty
				try {
					mi_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.invariants();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}

				// The dyadic is empty
				try {
					mi_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.getATA();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_check) {
				std::cout << OK("Dyadic2S Methods Checked!\n\n");
				std::cout << WARN("Beggining Dyadic2S Operators evaluation!");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2S Methods error!"));
			}

			// ----------------------------------------------------------------------------------------
			{
				// Comparison operators == and !=
				Dyadic2S mi_d2S_1({ 0,  1., 2. });
				Dyadic2S mi_d2S_2(mi_d2S_1);

				Dyadic2S mi_d2S_3({ 0., 1., 2., 3., 4., 5. });
				Dyadic2S mi_d2S_4(mi_d2S_3);

				if (!(mi_d2S_1 == mi_d2S_2)) mi_check = false;
				if (mi_d2S_3 != mi_d2S_4) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S comparison operators ( == and != ) Failed!!");
			}

			{
				// Operators + and +=
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_4;

				Dyadic2S mi_d2S_2({ 30., -5., 36. });
				Dyadic2S mi_d2S_3({ 12., -4., 14. });

				Dyadic2S mi_d2S_5({ 40., 49., -4., 12., -38.,  81. });
				Dyadic2S mi_d2S_6({ 23., -9., -3., 21.,  18., -60. });

				mi_d2S_1 = mi_d2S_2 + mi_d2S_3;
				mi_d2S_2 += mi_d2S_3;
				if (mi_d2S_1 != mi_d2S_2) mi_check = false;

				mi_d2S_4 = mi_d2S_5 + mi_d2S_6;
				mi_d2S_5 += mi_d2S_6;
				if (mi_d2S_4 != mi_d2S_5) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S addition operators ( + and += ) Failed!!");
			}

			{
				// Operators - and -=
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_4;

				Dyadic2S mi_d2S_2({ 30., -5., 36. });
				Dyadic2S mi_d2S_3({ 12., -4., 14. });

				Dyadic2S mi_d2S_5({ 40., 49., -4., 12., -38.,  81. });
				Dyadic2S mi_d2S_6({ 23., -9., -3., 21.,  18., -60. });

				mi_d2S_1 = mi_d2S_2 - mi_d2S_3;
				mi_d2S_2 -= mi_d2S_3;
				if (mi_d2S_1 != mi_d2S_2) mi_check = false;

				mi_d2S_4 = mi_d2S_5 - mi_d2S_6;
				mi_d2S_5 -= mi_d2S_6;
				if (mi_d2S_4 != mi_d2S_5) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S substraction operators ( - and -= ) Failed!!");
			}

			{
				// Operators * and *=
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_4;

				Dyadic2S mi_d2S_2({ 30., -5., 36. });
				Dyadic2S mi_d2S_3({ 12., -4., 14. });

				Dyadic2S mi_d2S_5({ 40., 49., -4., 12., -38.,  81. });
				Dyadic2S mi_d2S_6({ 23., -9., -3., 21.,  18., -60. });

				mi_d2S_1 = mi_d2S_2 * mi_d2S_3;
				mi_d2S_2 *= mi_d2S_3;
				if (mi_d2S_1 != mi_d2S_2) mi_check = false;

				mi_d2S_4 = mi_d2S_5 * mi_d2S_6;
				mi_d2S_5 *= mi_d2S_6;
				if (mi_d2S_4 != mi_d2S_5) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S multiplication operators ( * and *= ) Failed!!");
			}

			{
				// Operator * for dot product (with vector)
				std::vector<double> input1{ -1.1, -32.4 };
				std::vector<double> input2{ 17.3, 53.3, -94.3 };
				Dyadic2S mi_d2S_1({ 30., -5., 36. });
				Dyadic2S mi_d2S_2({ 40., 49., -4., 12., -38.,  81. });

				std::vector<double> expected1{ 129. , -1160.9 };
				std::vector<double> expected2{ 3680.9 , 5070.7, -9732.9 };

				auto output1 = mi_d2S_1 * input1;
				auto output2 = mi_d2S_2 * input2;

				for (int i = 0; i < expected1.size(); ++i) {
					if ((int)(output1.at(i) - expected1.at(i)) * 1000000 != 0) mi_check = false;
				}
				for (int i = 0; i < expected2.size(); ++i) {
					if ((int)(output2.at(i) - expected2.at(i)) * 1000000 != 0) mi_check = false;
				}
				if (!mi_check) std::cout << ERROR("Dyadic2S dot product operator ( Dyadis2S * vector ) Failed!!");
			}

			{
				// Operators + and += with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 + 2.5;
				mi_d2S_3 += 2.5;
				if (mi_d2S_1 != mi_d2S_3) mi_check = false;

				mi_d2S_2 = mi_d2S_4 + 4.8;
				mi_d2S_4 += 4.8;
				if (mi_d2S_2 != mi_d2S_4) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S addition operators with scalar ( + and += ) Failed!!");
			}

			{
				// Operators - and -= with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 - 2.5;
				mi_d2S_3 -= 2.5;
				if (mi_d2S_1 != mi_d2S_3) mi_check = false;

				mi_d2S_2 = mi_d2S_4 - 4.8;
				mi_d2S_4 -= 4.8;
				if (mi_d2S_2 != mi_d2S_4) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S substraction operators with scalar ( - and -= ) Failed!!");
			}

			{
				// Operators * and *= with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 * 0.25;
				mi_d2S_3 *= 0.25;
				if (mi_d2S_1 != mi_d2S_3) mi_check = false;

				mi_d2S_2 = mi_d2S_4 * 0.48;
				mi_d2S_4 *= 0.48;
				if (mi_d2S_2 != mi_d2S_4) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S multiplicative operators with scalar ( * and *= ) Failed!!");
			}

			{
				// Operators / and /= with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 / 0.25;
				mi_d2S_3 /= 0.25;
				if (mi_d2S_1 != mi_d2S_3) mi_check = false;

				mi_d2S_2 = mi_d2S_4 / 0.48;
				mi_d2S_4 /= 0.48;
				if (mi_d2S_2 != mi_d2S_4) mi_check = false;
				if (!mi_check) std::cout << ERROR("Dyadic2S division operators with scalar ( / and /= ) Failed!!");
			}

			{
				// Check if throws in previous methods are working
				// Division by zero
				try {
					mi_check = false;
					Dyadic2S mi_d2S({ 5., -5., 5. });
					mi_d2S /= 0.;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}

				// Check if throws in previous methods are working
				// Division by zero
				try {
					mi_check = false;
					Dyadic2S mi_d2S({ 5., -5., 5. });
					mi_d2S = mi_d2S / 0.;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_check = true;
				}
			}
			std::cout << WARN("Got here\n");

		}
		catch (std::runtime_error &error) {
			std::cerr << error.what() << std::endl;
		}
	}
}


