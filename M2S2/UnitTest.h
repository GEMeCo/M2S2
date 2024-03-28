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
	// Testing all features of Dyadic2S
	// 
	// ================================================================================================
	void unitTest_Dyadic2S() {

		bool mi_loc_check = true;
		bool mi_glob_check = true;

		try {
			// ----------------------------------------------------------------------------------------
			// Dyadic2S tests
			std::cout << std::endl << WARN("Beggining Dyadic2S Constructors evaluation!");

			// Constructors
			{
				Dyadic2S mi_d2S;		// Basic constructor
				if (mi_d2S.rows() != 0) mi_loc_check = false;
				if (mi_d2S.size() != 0) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 0)  mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S Basic Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				Dyadic2S mi_d2S(2);	// Constructor of 2D space (2x2);
				if (mi_d2S.rows() != 2) mi_loc_check = false;
				if (mi_d2S.size() != 3) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 3) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != 0.) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S 2D Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				Dyadic2S mi_d2S(3);	// Constructor of 3D space (3x3);
				if (mi_d2S.rows() != 3) mi_loc_check = false;
				if (mi_d2S.size() != 6) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 6) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != 0.) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S 3D Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				std::vector<double> mi_Voigt2D({ 0., 1., 2. });
				Dyadic2S mi_d2S(mi_Voigt2D);	// Constructor from vector (2D space);
				if (mi_d2S.rows() != 2) mi_loc_check = false;
				if (mi_d2S.size() != 3) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 3) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != mi_Voigt2D.at(i)) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S 2D Constructor from Vector Input Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				std::vector<double> mi_Voigt3D({ 0., 1., 2., 3., 4., 5. });
				Dyadic2S mi_d2S(mi_Voigt3D);	// Constructor from vector (3D space);
				if (mi_d2S.rows() != 3) mi_loc_check = false;
				if (mi_d2S.size() != 6) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 6) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != mi_Voigt3D.at(i)) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S 3D Constructor from Vector Input Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throw in constructor is working 
				// The vector has wrong dimensions
				try {
					mi_loc_check = false;
					std::vector<double> mi_NonVoigt({ 0., 1., 2., 3., 4. });
					Dyadic2S mi_d2S(mi_NonVoigt);
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			{
				// Copy constructor
				Dyadic2S mi_d2S_1({ 0., 1., 2. });
				Dyadic2S mi_d2S_2(mi_d2S_1);
				if (mi_d2S_2.rows() != 2) mi_loc_check = false;
				if (mi_d2S_2.size() != 3) mi_loc_check = false;
				if (mi_d2S_2.getVector().size() != 3) mi_loc_check = false;
				for (int i = 0; i < mi_d2S_2.getVector().size(); ++i) {
					if (mi_d2S_2.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S Copy Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Explicit move constructor
				Dyadic2S mi_d2S_1({ 0., 1., 2., 3., 4., 5. });
				Dyadic2S mi_d2S_2 = std::move(mi_d2S_1);
				if (mi_d2S_2.rows() != 3) mi_loc_check = false;
				if (mi_d2S_2.size() != 6) mi_loc_check = false;
				if (mi_d2S_2.getVector().size() != 6) mi_loc_check = false;
				for (int i = 0; i < mi_d2S_2.getVector().size(); ++i) {
					if (mi_d2S_2.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S Explicit Move Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Implicit move constructor
				Dyadic2S mi_d2S(Dyadic2S::identity(3));
				if (mi_d2S.rows() != 3) mi_loc_check = false;
				if (mi_d2S.size() != 6) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 6) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.rows(); ++i) {
					for (int j = i; j < mi_d2S.cols(); ++j) {
						if (i == j) {
							if (mi_d2S.at(i, j) != 1.) mi_loc_check = false;
						}
						if (i != j) {
							if (mi_d2S.at(i, j) != 0.) mi_loc_check = false;
						}
					}
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S Implicit Move Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Input by stream
				Dyadic2S mi_d2S(2);
				std::stringstream str("0. 1. 2.");
				str >> mi_d2S;
				if (mi_d2S.rows() != 2) mi_loc_check = false;
				if (mi_d2S.size() != 3) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 3) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - Value input by stream Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Input by stream
				Dyadic2S mi_d2S(3);
				std::stringstream str("0. 1. 2. 3. 4. 5.");
				str >> mi_d2S;
				if (mi_d2S.rows() != 3) mi_loc_check = false;
				if (mi_d2S.size() != 6) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 6) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - Value input by stream Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throw in constructor is working
				// The stream is unable to fill the dyadic
				try {
					mi_loc_check = false;

					Dyadic2S mi_d2S(3);
					std::stringstream str("0. 1. 2.");
					str >> mi_d2S;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			{
				// Check if throw in constructor is working
				// The dyadic size is smaller than stream
				try {
					mi_loc_check = false;

					Dyadic2S mi_d2S;
					std::stringstream str("0. 1. 2.");
					str >> mi_d2S;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_glob_check) {
				std::cout << OK("Dyadic2S constructors checked!\n\n");
				std::cout << WARN("Beggining Dyadic2S Methods evaluation!");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2S constructor error!"));
			}

			// ----------------------------------------------------------------------------------------
			{
				// Swap method
				Dyadic2S mi_d2S_1({ 0., 1., 2. });
				Dyadic2S mi_d2S_2({ 0., 1., 2., 3., 4., 5. });

				mi_d2S_2.swap(mi_d2S_1);
				if (mi_d2S_1.rows() != 3) mi_loc_check = false;
				if (mi_d2S_1.size() != 6) mi_loc_check = false;
				if (mi_d2S_1.getVector().size() != 6) mi_loc_check = false;
				for (int i = 0; i < mi_d2S_1.getVector().size(); ++i) {
					if (mi_d2S_1.getVector().at(i) != (double)i) mi_loc_check = false;
				}

				if (mi_d2S_2.rows() != 2) mi_loc_check = false;
				if (mi_d2S_2.size() != 3) mi_loc_check = false;
				if (mi_d2S_2.getVector().size() != 3) mi_loc_check = false;
				for (int i = 0; i < mi_d2S_2.getVector().size(); ++i) {
					if (mi_d2S_2.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - Swap Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Clear method
				Dyadic2S mi_d2S({ 0., 1., 2. });

				mi_d2S.clear();
				if (mi_d2S.rows() != 2) mi_loc_check = false;
				if (mi_d2S.size() != 3) mi_loc_check = false;
				if (mi_d2S.getVector().size() != 3) mi_loc_check = false;
				for (int i = 0; i < mi_d2S.getVector().size(); ++i) {
					if (mi_d2S.getVector().at(i) != 0.) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - Clear Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// getVoigtMnemonics method
				Dyadic2S mi_d2S_1({ 0., 1., 2. });
				Dyadic2S mi_d2S_2({ 0., 1., 2., 3., 4., 5. });

				auto mi_Voigt = mi_d2S_1.getVoigtMnemonics();
				if (mi_Voigt.size() != mi_d2S_1.size()) mi_loc_check = false;
				for (int i = 0; i < mi_d2S_1.rows(); ++i) {
					if (!almost_equal(mi_Voigt.at(i), mi_d2S_1.at(i, i))) mi_loc_check = false;

					for (unsigned int j = i + 1; j < mi_d2S_1.cols(); ++j)
						if (!almost_equal(mi_Voigt.at(mi_Voigt.size() - i - j), mi_d2S_1.at(i, j))) mi_loc_check = false;
				}

				mi_Voigt = mi_d2S_2.getVoigtMnemonics();
				if (mi_Voigt.size() != mi_d2S_2.size()) mi_loc_check = false;
				for (int i = 0; i < mi_d2S_2.rows(); ++i) {
					if (!almost_equal(mi_Voigt.at(i), mi_d2S_2.at(i, i))) mi_loc_check = false;

					for (unsigned int j = i + 1; j < mi_d2S_2.cols(); ++j)
						if (!almost_equal(mi_Voigt.at(mi_Voigt.size() - i - j), mi_d2S_2.at(i, j))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - getVoigtMnemonics Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Random symmetric matrices
				Dyadic2S mi_d2S_1 = Dyadic2S({ 2.5, 4.8, -3.6 });
				Dyadic2S mi_d2S_2 = Dyadic2S({ 2.5, 4.8, -3.6, 6.7, -1.4, 8.1 });

				// Determinant
				if (!almost_equal(mi_d2S_1.determinant(), -32.04)) mi_loc_check = false;
				if (!almost_equal(mi_d2S_2.determinant(), -94.297)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - determinant Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Norm
				if (!almost_equal(mi_d2S_1.norm(), 8.0802227692063045)) mi_loc_check = false;
				if (!almost_equal(mi_d2S_2.norm(), 13.880561948278607)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - norm Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Trace
				if (!almost_equal(mi_d2S_1.trace(), -1.1)) mi_loc_check = false;
				if (!almost_equal(mi_d2S_2.trace(), 17.3)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - trace Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throws in previous methods are working
				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.determinant();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.norm();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.trace();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
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
					if (!almost_equal(eval2D.at(i), expected2D.at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < eval3D.size(); ++i) {
					if (!almost_equal(eval3D.at(i), expected3D.at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - eigenvalue Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Invariants
				eval2D = mi_d2S_1.invariants();
				eval3D = mi_d2S_2.invariants();
				expected2D = std::vector<double>{ -1.1, -32.04 };
				expected3D = std::vector<double>{ 17.3, 53.31, -94.297 };

				for (int i = 0; i < eval2D.size(); ++i) {
					if (!almost_equal(eval2D.at(i), expected2D.at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < eval3D.size(); ++i) {
					if (!almost_equal(eval3D.at(i), expected3D.at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - invariants Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Inverse
				Dyadic2S mi_d2S_3 = mi_d2S_1.inverse();
				Dyadic2S mi_d2S_4 = mi_d2S_2.inverse();

				Dyadic2S mi_d2S_5 = Dyadic2S({ 0.112359550561798, 0.149812734082397, -0.0780274656679151 });
				Dyadic2S mi_d2S_6 = Dyadic2S({ -0.554736630009438, 0.358866135720118, -0.184523367657508, -0.0773089281737489, 0.14613402335175, 0.0667041369290645 });
				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_3.getVector().at(i), mi_d2S_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_4.getVector().at(i), mi_d2S_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - inverse Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// A^transpose * A
				mi_d2S_3 = mi_d2S_1.getATA();
				mi_d2S_4 = mi_d2S_2.getATA();
				mi_d2S_5 = Dyadic2S({ 29.29, -5.28,  36.00 });
				mi_d2S_6 = Dyadic2S({ 42.25, 49.20, -44.88, 69.89, -38.0, 80.53 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_3.getVector().at(i), mi_d2S_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_4.getVector().at(i), mi_d2S_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - getATA Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Spheric part
				mi_d2S_3 = mi_d2S_1.getSpheric();
				mi_d2S_4 = mi_d2S_2.getSpheric();
				mi_d2S_5 = Dyadic2S({ -0.366666666666667, 0., -0.366666666666667 });
				mi_d2S_6 = Dyadic2S({ 5.766666666666667, 0., 0., 5.766666666666667, 0., 5.766666666666667 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_3.getVector().at(i), mi_d2S_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_4.getVector().at(i), mi_d2S_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - getSpheric Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Deviatoric part
				mi_d2S_3 = mi_d2S_1.getDeviator();
				mi_d2S_4 = mi_d2S_2.getDeviator();
				mi_d2S_5 = Dyadic2S({ 2.8666666666667, 4.8, -3.2333333333333 });
				mi_d2S_6 = Dyadic2S({ -3.2666666666667, 4.8, -3.6, 0.93333333333333, -1.4, 2.33333333333333 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_3.getVector().at(i), mi_d2S_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_4.getVector().at(i), mi_d2S_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - getDeviator Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Double Dot product
				mi_d2S_5 = Dyadic2S({ 2.9, 5.8, -3.2 });
				mi_d2S_6 = Dyadic2S({ -3.3, 8.8, -2.6, 0.9, -1.8, 2.3 });

				if (!almost_equal(mi_d2S_1.contraction(mi_d2S_5), 74.45)) mi_loc_check = false;
				if (!almost_equal(mi_d2S_2.contraction(mi_d2S_6), 124.65)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S - contraction Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throws in previous methods are working
				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.eigenvalues();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.inverse();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// Singular Dyadic (non-invertible)
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S({ 1.0, -4.8, 23.04 });
					mi_d2S.inverse();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.invariants();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S;
					mi_d2S.getATA();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_glob_check) {
				std::cout << OK("Dyadic2S Methods checked!\n\n");
				std::cout << WARN("Beggining Dyadic2S Operators evaluation!");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2S Methods error!"));
			}

			// ----------------------------------------------------------------------------------------
			{
				// Comparison operators == and !=
				Dyadic2S mi_d2S_1({ 0., 1., 2. });
				Dyadic2S mi_d2S_2({ 0., 1., 2. });
					//mi_d2S_1);

				Dyadic2S mi_d2S_3({ 0., 1., 2., 3., 4., 5. });
				Dyadic2S mi_d2S_4(mi_d2S_3);

				if (!(mi_d2S_1 == mi_d2S_2)) mi_loc_check = false;
				if (mi_d2S_3 != mi_d2S_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S comparison operators ( == and != ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
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
				if (mi_d2S_1 != mi_d2S_2) mi_loc_check = false;

				mi_d2S_4 = mi_d2S_5 + mi_d2S_6;
				mi_d2S_5 += mi_d2S_6;
				if (mi_d2S_4 != mi_d2S_5) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S addition operators ( + and += ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
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
				if (mi_d2S_1 != mi_d2S_2) mi_loc_check = false;

				mi_d2S_4 = mi_d2S_5 - mi_d2S_6;
				mi_d2S_5 -= mi_d2S_6;
				if (mi_d2S_4 != mi_d2S_5) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S substraction operators ( - and -= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
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
				if (mi_d2S_1 != mi_d2S_2) mi_loc_check = false;

				mi_d2S_4 = mi_d2S_5 * mi_d2S_6;
				mi_d2S_5 *= mi_d2S_6;
				if (mi_d2S_4 != mi_d2S_5) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S multiplication operators ( * and *= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
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
					if (!almost_equal(output1.at(i), expected1.at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < expected2.size(); ++i) {
					if (!almost_equal(output2.at(i), expected2.at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S dot product operator ( Dyadis2S * vector ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators + and += with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 + 2.5;
				mi_d2S_3 += 2.5;
				if (mi_d2S_1 != mi_d2S_3) mi_loc_check = false;

				mi_d2S_2 = mi_d2S_4 + 4.8;
				mi_d2S_4 += 4.8;
				if (mi_d2S_2 != mi_d2S_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S addition operators with scalar ( + and += ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators - and -= with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 - 2.5;
				mi_d2S_3 -= 2.5;
				if (mi_d2S_1 != mi_d2S_3) mi_loc_check = false;

				mi_d2S_2 = mi_d2S_4 - 4.8;
				mi_d2S_4 -= 4.8;
				if (mi_d2S_2 != mi_d2S_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S substraction operators with scalar ( - and -= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators * and *= with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 * 0.25;
				mi_d2S_3 *= 0.25;
				if (mi_d2S_1 != mi_d2S_3) mi_loc_check = false;

				mi_d2S_2 = mi_d2S_4 * 0.48;
				mi_d2S_4 *= 0.48;
				if (mi_d2S_2 != mi_d2S_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S multiplicative operators with scalar ( * and *= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators / and /= with scalars
				Dyadic2S mi_d2S_1;
				Dyadic2S mi_d2S_2;

				Dyadic2S mi_d2S_3({ 30., -5., 36. });
				Dyadic2S mi_d2S_4({ 40., 49., -4., 12., -38.,  81. });

				mi_d2S_1 = mi_d2S_3 / 0.25;
				mi_d2S_3 /= 0.25;
				if (mi_d2S_1 != mi_d2S_3) mi_loc_check = false;

				mi_d2S_2 = mi_d2S_4 / 0.48;
				mi_d2S_4 /= 0.48;
				if (mi_d2S_2 != mi_d2S_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2S division operators with scalar ( / and /= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throws in previous methods are working
				// Division by zero
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S({ 5., -5., 5. });
					mi_d2S /= 0.;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// Check if throws in previous methods are working
				// Division by zero
				try {
					mi_loc_check = false;
					Dyadic2S mi_d2S({ 5., -5., 5. });
					mi_d2S = mi_d2S / 0.;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_glob_check) {
				std::cout << OK("Dyadic2S Operators checked!\n\n");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2S Operators error!"));
			}
		}
		catch (std::runtime_error &error) {
			std::cerr << error.what() << std::endl;
		}
	}

	// ================================================================================================
	//
	// Testing all features of Dyadic2N
	//
	// ================================================================================================
	void unitTest_Dyadic2N() {

		bool mi_loc_check = true;
		bool mi_glob_check = true;

		try {
			// ----------------------------------------------------------------------------------------
			// Dyadic2N tests
			std::cout << std::endl << WARN("Beggining Dyadic2N Constructors evaluation!");

			// Constructors
			{
				Dyadic2N mi_d2N;		// Basic constructor
				if (mi_d2N.rows() != 0) mi_loc_check = false;
				if (mi_d2N.size() != 0) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 0)  mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N Basic Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				Dyadic2N mi_d2N(2);	// Constructor of 2D space (2x2);
				if (mi_d2N.rows() != 2) mi_loc_check = false;
				if (mi_d2N.size() != 4) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 4) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.getVector().size(); ++i) {
					if (mi_d2N.getVector().at(i) != 0.) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N 2D Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				Dyadic2N mi_d2N(3);	// Constructor of 3D space (3x3);
				if (mi_d2N.rows() != 3) mi_loc_check = false;
				if (mi_d2N.size() != 9) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 9) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.getVector().size(); ++i) {
					if (mi_d2N.getVector().at(i) != 0.) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N 3D Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				std::vector<double> mi_Voigt2D({ 0., 1., 2. , 3.});
				Dyadic2N mi_d2N(mi_Voigt2D);	// Constructor from vector (2D space);
				if (mi_d2N.rows() != 2) mi_loc_check = false;
				if (mi_d2N.size() != 4) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 4) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.getVector().size(); ++i) {
					if (mi_d2N.getVector().at(i) != mi_Voigt2D.at(i)) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N 2D Constructor from Vector Input Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				std::vector<double> mi_Voigt3D({ 0., 1., 2., 3., 4., 5., 6., 7., 8. });
				Dyadic2N mi_d2N(mi_Voigt3D);	// Constructor from vector (3D space);
				if (mi_d2N.rows() != 3) mi_loc_check = false;
				if (mi_d2N.size() != 9) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 9) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.getVector().size(); ++i) {
					if (mi_d2N.getVector().at(i) != mi_Voigt3D.at(i)) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N 3D Constructor from Vector Input Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throw in constructor is working 
				// The vector has wrong dimensions
				try {
					mi_loc_check = false;
					std::vector<double> mi_NonVoigt({ 0., 1., 2., 3., 4. });
					Dyadic2N mi_d2N(mi_NonVoigt);
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			{
				// Copy constructor
				Dyadic2N mi_d2N_1({ 0., 1., 2., 3. });
				Dyadic2N mi_d2N_2(mi_d2N_1);
				if (mi_d2N_2.rows() != 2) mi_loc_check = false;
				if (mi_d2N_2.size() != 4) mi_loc_check = false;
				if (mi_d2N_2.getVector().size() != 4) mi_loc_check = false;
				for (int i = 0; i < mi_d2N_2.getVector().size(); ++i) {
					if (mi_d2N_2.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N Copy Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Explicit move constructor
				Dyadic2N mi_d2N_1({ 0., 1., 2., 3., 4., 5., 6., 7., 8. });
				Dyadic2N mi_d2N_2 = std::move(mi_d2N_1);
				if (mi_d2N_2.rows() != 3) mi_loc_check = false;
				if (mi_d2N_2.size() != 9) mi_loc_check = false;
				if (mi_d2N_2.getVector().size() != 9) mi_loc_check = false;
				for (int i = 0; i < mi_d2N_2.getVector().size(); ++i) {
					if (mi_d2N_2.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N Explicit Move Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Implicit move constructor
				Dyadic2N mi_d2N(Dyadic2N::identity(3));
				if (mi_d2N.rows() != 3) mi_loc_check = false;
				if (mi_d2N.size() != 9) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 9) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.rows(); ++i) {
					for (int j = i; j < mi_d2N.cols(); ++j) {
						if (i == j) {
							if (mi_d2N.at(i, j) != 1.) mi_loc_check = false;
						}
						if (i != j) {
							if (mi_d2N.at(i, j) != 0.) mi_loc_check = false;
						}
					}
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N Implicit Move Constructor Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Input by stream
				Dyadic2N mi_d2N(2);
				std::stringstream str("0. 1. 2. 3.");
				str >> mi_d2N;
				if (mi_d2N.rows() != 2) mi_loc_check = false;
				if (mi_d2N.size() != 4) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 4) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.getVector().size(); ++i) {
					if (mi_d2N.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - Value input by stream Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Input by stream
				Dyadic2N mi_d2N(3);
				std::stringstream str("0. 1. 2. 3. 4. 5. 6. 7. 8.");
				str >> mi_d2N;
				if (mi_d2N.rows() != 3) mi_loc_check = false;
				if (mi_d2N.size() != 9) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 9) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.getVector().size(); ++i) {
					if (mi_d2N.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - Value input by stream Failed!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throw in constructor is working
				// The stream is unable to fill the dyadic
				try {
					mi_loc_check = false;

					Dyadic2N mi_d2N(3);
					std::stringstream str("0. 1. 2.");
					str >> mi_d2N;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			{
				// Check if throw in constructor is working
				// The dyadic size is smaller than stream
				try {
					mi_loc_check = false;

					Dyadic2N mi_d2N;
					std::stringstream str("0. 1. 2. 3.");
					str >> mi_d2N;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_glob_check) {
				std::cout << OK("Dyadic2N constructors checked!\n\n");
				std::cout << WARN("Beggining Dyadic2N Methods evaluation!");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2N constructor error!"));
			}

			// ----------------------------------------------------------------------------------------
			{
				// Swap method
				Dyadic2N mi_d2N_1({ 0., 1., 2., 3. });
				Dyadic2N mi_d2N_2({ 0., 1., 2., 3., 4., 5., 6., 7., 8. });

				mi_d2N_2.swap(mi_d2N_1);

				if (mi_d2N_1.rows() != 3) mi_loc_check = false;
				if (mi_d2N_1.size() != 9) mi_loc_check = false;
				if (mi_d2N_1.getVector().size() != 9) mi_loc_check = false;
				for (int i = 0; i < mi_d2N_1.getVector().size(); ++i) {
					if (mi_d2N_1.getVector().at(i) != (double)i) mi_loc_check = false;
				}

				if (mi_d2N_2.rows() != 2) mi_loc_check = false;
				if (mi_d2N_2.size() != 4) mi_loc_check = false;
				if (mi_d2N_2.getVector().size() != 4) mi_loc_check = false;
				for (int i = 0; i < mi_d2N_2.getVector().size(); ++i) {
					if (mi_d2N_2.getVector().at(i) != (double)i) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - Swap Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Clear method
				Dyadic2N mi_d2N({ 0., 1., 2., 3. });

				mi_d2N.clear();
				if (mi_d2N.rows() != 2) mi_loc_check = false;
				if (mi_d2N.size() != 4) mi_loc_check = false;
				if (mi_d2N.getVector().size() != 4) mi_loc_check = false;
				for (int i = 0; i < mi_d2N.getVector().size(); ++i) {
					if (mi_d2N.getVector().at(i) != 0.) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - Clear Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// getVoigtMnemonics method
				Dyadic2N mi_d2N_1({ 0., 1., 2., 3. });
				Dyadic2N mi_d2N_2({ 0., 1., 2., 3., 4., 5., 6., 7., 8. });

				auto mi_Voigt = mi_d2N_1.getVoigtMnemonics();
				unsigned int mi_nVoigt = 3 * mi_d2N_1.rows() - 3;
				if (mi_Voigt.size() != mi_d2N_1.size()) mi_loc_check = false;
				for (int i = 0; i < mi_d2N_1.rows(); ++i) {
					if (!almost_equal(mi_Voigt.at(i), mi_d2N_1.at(i, i))) mi_loc_check = false;

					for (unsigned int j = i + 1; j < mi_d2N_1.cols(); ++j) {
						if (!almost_equal(mi_Voigt.at(mi_nVoigt - i - j), mi_d2N_1.at(i, j))) mi_loc_check = false;
						if (!almost_equal(mi_Voigt.at(mi_d2N_1.size() - i - j), mi_d2N_1.at(j, i))) mi_loc_check = false;
					}
				}

				mi_Voigt = mi_d2N_2.getVoigtMnemonics();
				mi_nVoigt = 3 * mi_d2N_2.rows() - 3;
				if (mi_Voigt.size() != mi_d2N_2.size()) mi_loc_check = false;
				for (int i = 0; i < mi_d2N_2.rows(); ++i) {
					if (!almost_equal(mi_Voigt.at(i), mi_d2N_2.at(i, i))) mi_loc_check = false;

					for (unsigned int j = i + 1; j < mi_d2N_2.cols(); ++j) {
						if (!almost_equal(mi_Voigt.at(mi_nVoigt - i - j), mi_d2N_2.at(i, j))) mi_loc_check = false;
						if (!almost_equal(mi_Voigt.at(mi_d2N_2.size() - i - j), mi_d2N_2.at(j, i))) mi_loc_check = false;
					}
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - getVoigtMnemonics Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Random symmetric matrices
				Dyadic2N mi_d2N_1 = Dyadic2N({ 2.5, 4.8, 4.8, -3.6 });
				Dyadic2N mi_d2N_2 = Dyadic2N({ 2.5, 4.8, -3.6, 4.8, 6.7, -1.4, -3.6, -1.4, 8.1 });

				// Determinant
				if (!almost_equal(mi_d2N_1.determinant(), -32.04)) mi_loc_check = false;
				if (!almost_equal(mi_d2N_2.determinant(), -94.297)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - determinant Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Norm
				if (!almost_equal(mi_d2N_1.norm(), 8.0802227692063045)) mi_loc_check = false;
				if (!almost_equal(mi_d2N_2.norm(), 13.880561948278607)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - norm Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Trace
				if (!almost_equal(mi_d2N_1.trace(), -1.1)) mi_loc_check = false;
				if (!almost_equal(mi_d2N_2.trace(), 17.3)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - trace Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throws in previous methods are working
				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.determinant();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.norm();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.trace();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			{
				// Random matrices
				Dyadic2N mi_d2N_1 = Dyadic2N({ 2.5, 4.8, 4.8, -3.6 });
				Dyadic2N mi_d2N_2 = Dyadic2N({ 2.5, 4.8, -3.6, 4.8, 6.7, -1.4, -3.6, -1.4, 8.1 });

				// Eigenvalues
				std::vector<double> eval2D = mi_d2N_1.eigenvalues();
				std::vector<double> eval3D = mi_d2N_2.eigenvalues();
				std::vector<double> expected2D({ -6.237046685231272,  5.1370466852312724 });
				std::vector<double> expected3D({ 12.377474154257952, -1.2368787192025401, 6.1594045649445857 });

				for (int i = 0; i < eval2D.size(); ++i) {
					if (!almost_equal(eval2D.at(i), expected2D.at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < eval3D.size(); ++i) {
					if (!almost_equal(eval3D.at(i), expected3D.at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - eigenvalue Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Invariants
				eval2D = mi_d2N_1.invariants();
				eval3D = mi_d2N_2.invariants();
				expected2D = std::vector<double>{ -1.1, -32.04 };
				expected3D = std::vector<double>{ 17.3, 53.31, -94.297 };

				for (int i = 0; i < eval2D.size(); ++i) {
					if (!almost_equal(eval2D.at(i), expected2D.at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < eval3D.size(); ++i) {
					if (!almost_equal(eval3D.at(i), expected3D.at(i))) mi_loc_check = false;

				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - invariants Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Inverse
				Dyadic2N mi_d2N_3 = mi_d2N_1.inverse();
				Dyadic2N mi_d2N_4 = mi_d2N_2.inverse();

				Dyadic2N mi_d2N_5 = Dyadic2N({ 0.112359550561798, 0.149812734082397, 0.149812734082397, -0.0780274656679151 });
				Dyadic2N mi_d2N_6 = Dyadic2N({ -0.554736630009438, 0.358866135720118, -0.184523367657508, 0.358866135720118, -0.0773089281737489, 0.14613402335175, -0.184523367657508, 0.14613402335175, 0.0667041369290645 });

				for (int i = 0; i < mi_d2N_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_3.getVector().at(i), mi_d2N_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2N_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_4.getVector().at(i), mi_d2N_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - inverse Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// A^transpose * A
				Dyadic2S mi_d2S_3 = mi_d2N_1.getATA();
				Dyadic2S mi_d2S_4 = mi_d2N_2.getATA();
				Dyadic2S mi_d2S_5 = Dyadic2S({ 29.29, -5.28,  36.00 });
				Dyadic2S mi_d2S_6 = Dyadic2S({ 42.25, 49.20, -44.88, 69.89, -38.0, 80.53 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_3.getVector().at(i), mi_d2S_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_4.getVector().at(i), mi_d2S_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - getATA Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Spheric part
				mi_d2S_3 = mi_d2N_1.getSpheric();
				mi_d2S_4 = mi_d2N_2.getSpheric();
				mi_d2S_5 = Dyadic2S({ -0.366666666666667, 0., -0.366666666666667 });
				mi_d2S_6 = Dyadic2S({ 5.766666666666667, 0., 0., 5.766666666666667, 0., 5.766666666666667 });

				for (int i = 0; i < mi_d2S_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_3.getVector().at(i), mi_d2S_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2S_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_4.getVector().at(i), mi_d2S_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - getSpheric Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Deviatoric part
				mi_d2N_3 = mi_d2N_1.getDeviator();
				mi_d2N_4 = mi_d2N_2.getDeviator();
				mi_d2N_5 = Dyadic2N({ 2.8666666666667, 4.8, 4.8, -3.2333333333333 });
				mi_d2N_6 = Dyadic2N({ -3.2666666666667, 4.8, -3.6, 4.8, 0.93333333333333, -1.4, -3.6, -1.4, 2.33333333333333 });

				for (int i = 0; i < mi_d2N_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_3.getVector().at(i), mi_d2N_5.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2N_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_4.getVector().at(i), mi_d2N_6.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - getDeviator Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// Double Dot product
				mi_d2N_5 = Dyadic2N({ 2.9, 5.8, 5.8, -3.2 });
				mi_d2N_6 = Dyadic2N({ -3.3, 8.8, -2.6, 8.8, 0.9, -1.8, -2.6, -1.8, 2.3 });

				if (!almost_equal(mi_d2N_1.contraction(mi_d2N_5), 74.45)) mi_loc_check = false;
				if (!almost_equal(mi_d2N_2.contraction(mi_d2N_6), 124.65)) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - contraction Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Methods only for nonsymmetric tensors
				// Transpose
				Dyadic2N mi_d2N_1 = Dyadic2N({ 2.5, 4.8, -3.6, 4.8 });
				Dyadic2N mi_d2N_2 = Dyadic2N({ 2.5, 4.8, -3.6, 6.7, -1.4, 4.8, 8.1, -3.6, -1.4 });

				Dyadic2N expected2D = Dyadic2N({ 2.5, -3.6, 4.8, 4.8 });
				Dyadic2N expected3D = Dyadic2N({ 2.5, 6.7, 8.1, 4.8, -1.4, -3.6, -3.6, 4.8, -1.4 });

				Dyadic2N mi_d2N_3 = mi_d2N_1.transpose();
				Dyadic2N mi_d2N_4 = mi_d2N_2.transpose();

				for (int i = 0; i < mi_d2N_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_3.getVector().at(i), expected2D.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2N_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_4.getVector().at(i), expected3D.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - transpose Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// getSymmetric
				Dyadic2S mi_d2S_1 = mi_d2N_1.getSymmetric();
				Dyadic2S mi_d2S_2 = mi_d2N_2.getSymmetric();

				Dyadic2S expectedS2D = Dyadic2S({ 2.5, 0.6, 4.8 });
				Dyadic2S expectedS3D = Dyadic2S({ 2.5, 5.75, 2.25, -1.4, 0.6, -1.4 });

				for (int i = 0; i < mi_d2S_1.getVector().size(); ++i) {
					if (!almost_equal(mi_d2S_1.getVector().at(i), expectedS2D.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2N_2.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_2.getVector().at(i), expectedS3D.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - getSymmetric Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}

				// getAsymmetric
				mi_d2N_3 = mi_d2N_1.getAsymmetric();
				mi_d2N_4 = mi_d2N_2.getAsymmetric();

				expected2D = mi_d2N_1 - expectedS2D;
				expected3D = mi_d2N_2 - expectedS3D;

				for (int i = 0; i < mi_d2N_3.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_3.getVector().at(i), expected2D.getVector().at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < mi_d2N_4.getVector().size(); ++i) {
					if (!almost_equal(mi_d2N_4.getVector().at(i), expected3D.getVector().at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N - getAsymmetric Method Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throws in previous methods are working
				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.eigenvalues();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.inverse();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// Singular Dyadic (non-invertible)
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N({ 1.0, -4.8, 23.04 });
					mi_d2N.inverse();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.invariants();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.getATA();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.transpose();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.getSymmetric();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// The dyadic is empty
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N;
					mi_d2N.getAsymmetric();
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_glob_check) {
				std::cout << OK("Dyadic2N Methods checked!\n\n");
				std::cout << WARN("Beggining Dyadic2N Operators evaluation!");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2N Methods error!"));
			}

			// ----------------------------------------------------------------------------------------
			{
				// Comparison operators == and !=
				Dyadic2N mi_d2N_1({ 0., 1., 2., 3. });
				Dyadic2N mi_d2N_2(mi_d2N_1);

				Dyadic2N mi_d2N_3({ 0., 1., 2., 3., 4., 5., 6., 7., 8. });
				Dyadic2N mi_d2N_4(mi_d2N_3);

				if (!(mi_d2N_1 == mi_d2N_2)) mi_loc_check = false;
				if (mi_d2N_3 != mi_d2N_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N comparison operators ( == and != ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators + and +=
				Dyadic2N mi_d2N_1;
				Dyadic2N mi_d2N_4;

				Dyadic2N mi_d2N_2({ 30., -5., 36., -22. });
				Dyadic2N mi_d2N_3({ 12., -4., 14., -17. });

				Dyadic2N mi_d2N_5({ 40., 49., -4., 12., -38.,  81., 18., -72., 32.8 });
				Dyadic2N mi_d2N_6({ 23., -9., -3., 21.,  18., -60., 51., -24., 26.5 });

				mi_d2N_1 = mi_d2N_2 + mi_d2N_3;
				mi_d2N_2 += mi_d2N_3;
				if (mi_d2N_1 != mi_d2N_2) mi_loc_check = false;

				mi_d2N_4 = mi_d2N_5 + mi_d2N_6;
				mi_d2N_5 += mi_d2N_6;
				if (mi_d2N_4 != mi_d2N_5) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N addition operators ( + and += ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators - and -=
				Dyadic2N mi_d2N_1;
				Dyadic2N mi_d2N_4;

				Dyadic2N mi_d2N_2({ 30., -5., 36., -22. });
				Dyadic2N mi_d2N_3({ 12., -4., 14., -17. });

				Dyadic2N mi_d2N_5({ 40., 49., -4., 12., -38.,  81., 18., -72., 32.8 });
				Dyadic2N mi_d2N_6({ 23., -9., -3., 21.,  18., -60., 51., -24., 26.5 });

				mi_d2N_1 = mi_d2N_2 - mi_d2N_3;
				mi_d2N_2 -= mi_d2N_3;
				if (mi_d2N_1 != mi_d2N_2) mi_loc_check = false;

				mi_d2N_4 = mi_d2N_5 - mi_d2N_6;
				mi_d2N_5 -= mi_d2N_6;
				if (mi_d2N_4 != mi_d2N_5) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N substraction operators ( - and -= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators * and *=
				Dyadic2N mi_d2N_1;
				Dyadic2N mi_d2N_4;

				Dyadic2N mi_d2N_2({ 30., -5., 36., -22. });
				Dyadic2N mi_d2N_3({ 12., -4., 14., -17. });

				Dyadic2N mi_d2N_5({ 40., 49., -4., 12., -38.,  81., 18., -72., 32.8 });
				Dyadic2N mi_d2N_6({ 23., -9., -3., 21.,  18., -60., 51., -24., 26.5 });

				mi_d2N_1 = mi_d2N_2 * mi_d2N_3;
				mi_d2N_2 *= mi_d2N_3;
				if (mi_d2N_1 != mi_d2N_2) mi_loc_check = false;

				mi_d2N_4 = mi_d2N_5 * mi_d2N_6;
				mi_d2N_5 *= mi_d2N_6;
				if (mi_d2N_4 != mi_d2N_5) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N multiplication operators ( * and *= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operator * for dot product (with vector)
				std::vector<double> input1{ -1.1, -32.4 };
				std::vector<double> input2{ 17.3, 53.3, -94.3 };
				Dyadic2N mi_d2N_1({ 30., -5., -5, 36. });
				Dyadic2N mi_d2N_2({ 40., 49., -4., 49., 12., -38., -4., -38., 81. });

				std::vector<double> expected1{ 129. , -1160.9 };
				std::vector<double> expected2{ 3680.9 , 5070.7, -9732.9 };

				auto output1 = mi_d2N_1 * input1;
				auto output2 = mi_d2N_2 * input2;

				for (int i = 0; i < expected1.size(); ++i) {
					if (!almost_equal(output1.at(i), expected1.at(i))) mi_loc_check = false;
				}
				for (int i = 0; i < expected2.size(); ++i) {
					if (!almost_equal(output2.at(i), expected2.at(i))) mi_loc_check = false;
				}
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N dot product operator ( Dyadis2N * vector ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators + and += with scalars
				Dyadic2N mi_d2N_1;
				Dyadic2N mi_d2N_2;

				Dyadic2N mi_d2N_3({ 30., -5., 36., -22. });
				Dyadic2N mi_d2N_4({ 40., 49., -4., 12., -38.,  81., 18., -72., 32.8 });

				mi_d2N_1 = mi_d2N_3 + 2.5;
				mi_d2N_3 += 2.5;
				if (mi_d2N_1 != mi_d2N_3) mi_loc_check = false;

				mi_d2N_2 = mi_d2N_4 + 4.8;
				mi_d2N_4 += 4.8;
				if (mi_d2N_2 != mi_d2N_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N addition operators with scalar ( + and += ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators - and -= with scalars
				Dyadic2N mi_d2N_1;
				Dyadic2N mi_d2N_2;

				Dyadic2N mi_d2N_3({ 30., -5., 36., -22. });
				Dyadic2N mi_d2N_4({ 40., 49., -4., 12., -38.,  81., 18., -72., 32.8 });

				mi_d2N_1 = mi_d2N_3 - 2.5;
				mi_d2N_3 -= 2.5;
				if (mi_d2N_1 != mi_d2N_3) mi_loc_check = false;

				mi_d2N_2 = mi_d2N_4 - 4.8;
				mi_d2N_4 -= 4.8;
				if (mi_d2N_2 != mi_d2N_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N substraction operators with scalar ( - and -= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators * and *= with scalars
				Dyadic2N mi_d2N_1;
				Dyadic2N mi_d2N_2;

				Dyadic2N mi_d2N_3({ 30., -5., 36., -22. });
				Dyadic2N mi_d2N_4({ 40., 49., -4., 12., -38.,  81., 18., -72., 32.8 });

				mi_d2N_1 = mi_d2N_3 * 0.25;
				mi_d2N_3 *= 0.25;
				if (mi_d2N_1 != mi_d2N_3) mi_loc_check = false;

				mi_d2N_2 = mi_d2N_4 * 0.48;
				mi_d2N_4 *= 0.48;
				if (mi_d2N_2 != mi_d2N_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N multiplicative operators with scalar ( * and *= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Operators / and /= with scalars
				Dyadic2N mi_d2N_1;
				Dyadic2N mi_d2N_2;

				Dyadic2N mi_d2N_3({ 30., -5., 36., -22. });
				Dyadic2N mi_d2N_4({ 40., 49., -4., 12., -38.,  81., 18., -72., 32.8 });

				mi_d2N_1 = mi_d2N_3 / 0.25;
				mi_d2N_3 /= 0.25;
				if (mi_d2N_1 != mi_d2N_3) mi_loc_check = false;

				mi_d2N_2 = mi_d2N_4 / 0.48;
				mi_d2N_4 /= 0.48;
				if (mi_d2N_2 != mi_d2N_4) mi_loc_check = false;
				if (!mi_loc_check) {
					std::cout << ERROR("Dyadic2N division operators with scalar ( / and /= ) Failed!!");
					mi_glob_check = false;
					mi_loc_check = true;
				}
			}

			{
				// Check if throws in previous methods are working
				// Division by zero
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N({ 5., -5., -5., 5. });
					mi_d2N /= 0.;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}

				// Check if throws in previous methods are working
				// Division by zero
				try {
					mi_loc_check = false;
					Dyadic2N mi_d2N({ 5., -5., -5., 5. });
					mi_d2N = mi_d2N / 0.;
				}
				catch (std::exception& ierr) {
					//std::cout << ierr.what() << std::endl;
					mi_loc_check = true;
				}
			}

			// ----------------------------------------------------------------------------------------
			if (mi_glob_check) {
				std::cout << OK("Dyadic2N Operators checked!\n\n");
			}
			else {
				throw std::runtime_error(ERROR("Check Dyadic2N Operators error!"));
			}
		}
		catch (std::runtime_error& error) {
			std::cerr << error.what() << std::endl;
		}
	}

	// ================================================================================================
	//
	// Unit testing - Checks every function from all classes
	// 
	// ================================================================================================
	void unitTest_Global() {
		std::cout << std::endl << WARN("Testing every aspects of M2S2 library!");

		unitTest_Dyadic2S();
		unitTest_Dyadic2N();
	}
}
