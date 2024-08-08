// ================================================================================================
// 
// This file is part of M2S2 - Matrices for Mechanices of Solids and Structures
//
// Copyright(C) 2024 
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

// M2S2 libraries
#include "Common.h"
#include "Dyadic2S.h"
#include "Dyadic2N.h"
#include "Dyadic4S.h"
#include "Dyadic4C.h"
#include "MatrixS.h"
#include "MatrixX.h"
#include "MatrixSparse.h"

namespace M2S2 {
	// ================================================================================================
	//
	// Helper functions
	//
	// ================================================================================================
	//
	// Operators for Symmetric and Asymmetric dyadics
	//
	// ================================================================================================

	/** Overloads operator + for addition -> R_ij = F_ij + S_ij
	  * @param first First dyadic to be added.
	  * @param second Second dyadic to be added.
	  * @return an asymmetric dyadic with the addition of components.
	  */
	inline Dyadic2N operator + (const Dyadic2S& first, const Dyadic2N& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on add operator: Empty tensor!"));
		assert(first.rows() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(second.rows());
		for (unsigned int i = 0; i < result.rows(); ++i) {
			for (unsigned int j = 0; j < result.cols(); ++j) {
				result.at(i, j) = first.at(i, j) + second.at(i, j);
			}
		}
		return result;
	}

	/** Overloads operator - for substraction -> R_ij = F_ij - S_ij
	  * @param first First dyadic to be substracted.
	  * @param second Second dyadic to be substracted.
	  * @return an asymmetric dyadic with the substraction of components.
	  */
	inline Dyadic2N operator - (const Dyadic2S& first, const Dyadic2N& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on sub operator: Empty tensor!"));
		assert(first.rows() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(second.rows());
		for (unsigned int i = 0; i < result.rows(); ++i) {
			for (unsigned int j = 0; j < result.cols(); ++j) {
				result.at(i, j) = first.at(i, j) - second.at(i, j);
			}
		}
		return result;
	}

	/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
	  * @param first Dyadic to multiply.
	  * @param second Dyadic to be multiplied with.
	  * @return an asymmetric dyadic with the multiplication of components.
	  */
	inline Dyadic2N operator * (const Dyadic2S& first, const Dyadic2N& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on mult operator: Empty tensor!"));
		assert(first.cols() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(second.rows());
		for (unsigned int i = 0; i < result.rows(); ++i) {
			for (unsigned int j = 0; j < first.cols(); ++j) {
				for (unsigned int k = 0; k < second.rows(); ++k) {
					result.at(i, j) += first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	/** Overloads operator + for addition -> R_ij = F_ij + S_ij
	  * @param first First dyadic to be added.
	  * @param second Second dyadic to be added.
	  * @return an asymmetric dyadic with the addition of components.
	  */
	inline Dyadic2N operator + (const Dyadic2N& first, const Dyadic2S& second)
	{
		return (second + first);
	}

	/** Overloads operator - for substraction -> R_ij = F_ij - S_ij
	  * @param first First dyadic to be substracted.
	  * @param second Second dyadic to be substracted.
	  * @return an asymmetric dyadic with the substraction of components.
	  */
	inline Dyadic2N operator - (const Dyadic2N& first, const Dyadic2S& second)
	{
		return (second - first) * (-1.);
	}

	/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
	  * @param first Dyadic to multiply.
	  * @param second Dyadic to be multiplied with.
	  * @return an asymmetric dyadic with the multiplication of components.
	  */
	inline Dyadic2N operator * (const Dyadic2N& first, const Dyadic2S& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on mult operator: Empty tensor!"));
		assert(first.cols() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(first.rows());
		for (unsigned int i = 0; i < first.rows(); ++i) {
			for (unsigned int j = 0; j < second.cols(); ++j) {
				for (unsigned int k = 0; k < first.cols(); ++k) {
					result.at(i, j) += first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
	  * @param first Dyadic to multiply.
	  * @param second Dyadic to be multiplied with.
	  * @return an asymmetric dyadic with the multiplication of components.
	  */
	inline Dyadic2N operator * (const Dyadic2S& first, const Dyadic2S& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on mult operator: Empty tensor!"));
		assert(first.cols() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(first.rows());
		for (unsigned int i = 0; i < first.rows(); i++) {
			for (unsigned int j = 0; j < second.cols(); j++) {
				for (unsigned int k = 0; k < first.cols(); k++) {
					result.at(i, j) += first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Dyadics and Scalars
	//
	// ================================================================================================

	/** Overloads operator + to sum a scalar with a dyadic -> T_ij = T_ij + alfa
	  * @param alfa Scalar to be added.
	  * @param dyadic
	  * @return a dyadic with the result.
	  */
	inline Dyadic2S operator + (const double alfa, const Dyadic2S& dyadic)
	{
		return dyadic + alfa;
	}

	/** Overloads operator - to substract a scalar with a dyadic -> T_ij = alfa - T_ij
	  * @param alfa Scalar to be substracted.
	  * @param dyadic Dyadic to be substracted with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic2S operator - (const double alfa, const Dyadic2S& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
	  * @param alfa Scalar to be multiplied with.
	  * @param dyadic Dyadic to be multiplied with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic2S operator * (const double alfa, const Dyadic2S& dyadic)
	{
		return dyadic * alfa;
	}

	/** Overloads operator + to sum a scalar with a dyadic -> T_ij = T_ij + alfa
	  * @param alfa Scalar to be added.
	  * @param dyadic Dyadic to be added with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic2N operator + (const double alfa, const Dyadic2N& dyadic)
	{
		return dyadic + alfa;
	}

	/** Overloads operator - to substract a scalar with a dyadic -> T_ij = alfa - T_ij
	  * @param alfa Scalar to be substracted.
	  * @param dyadic Dyadic to be substracted with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic2N operator - (const double alfa, const Dyadic2N& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
	  * @param alfa Scalar to be multiplied with.
	  * @param dyadic Dyadic to be multiplied with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic2N operator * (const double alfa, const Dyadic2N& dyadic)
	{
		return dyadic * alfa;
	}

	/** Overloads operator + to sum a scalar with a dyadic -> T_ij = T_ij + alfa
	  * @param alfa Scalar to be added.
	  * @param dyadic Dyadic to be added with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic4C operator + (const double alfa, const Dyadic4C& dyadic)
	{
		return dyadic + alfa;
	}

	/** Overloads operator - to substract a scalar with a dyadic -> T_ij = alfa - T_ij
	  * @param alfa Scalar to be substracted.
	  * @param dyadic Dyadic to be substracted with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic4C operator - (const double alfa, const Dyadic4C& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
	  * @param alfa Scalar to be multiplied with.
	  * @param dyadic Dyadic to be multiplied with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic4C operator * (const double alfa, const Dyadic4C& dyadic)
	{
		return dyadic * alfa;
	}

	/** Overloads operator + to sum a scalar with a dyadic -> T_ij = T_ij + alfa
	  * @param alfa Scalar to be added.
	  * @param dyadic Dyadic to be added with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic4S operator + (const double alfa, const Dyadic4S& dyadic)
	{
		return dyadic + alfa;
	}

	/** Overloads operator - to substract a scalar with a dyadic -> T_ij = alfa - T_ij
	  * @param alfa Scalar to be substracted.
	  * @param dyadic Dyadic to be substracted with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic4S operator - (const double alfa, const Dyadic4S& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
	  * @param alfa Scalar to be multiplied with.
	  * @param dyadic Dyadic to be multiplied with.
	  * @return a dyadic with the result.
	  */
	inline Dyadic4S operator * (const double alfa, const Dyadic4S& dyadic)
	{
		return dyadic * alfa;
	}

	// ================================================================================================
	//
	// Operators for Symmetric and Asymmetric matrices
	//
	// ================================================================================================
	//
	// Symmetric + Asymmetric matrices
	//
	// ================================================================================================

	/** Overloads operator + for addition -> R_ij = F_ij + S_ij
	  * @param first First matrix to be added.
	  * @param second Second matrix to be added.
	  * @return an asymmetric matrix with the addition of components.
	  */
	inline MatrixX operator + (const MatrixS& first, const MatrixX& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on add operator: Empty tensor!"));
		assert(first.rows() == second.rows());		// Size of dyadics does not correspond!
		assert(first.cols() == second.cols());		// Size of dyadics does not correspond!

		MatrixX result(first.rows(), first.cols());
		for (unsigned int i = 0; i < first.rows(); ++i) {
			for (unsigned int j = 0; j < first.cols(); ++j) {
				result.at(i, j) = first.at(i, j) + second.at(i, j);
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Symmetric - Asymmetric matrices
	//
	// ================================================================================================

	/** Overloads operator - for substraction -> R_ij = F_ij - S_ij
	  * @param first First matrix to be substracted.
	  * @param second Second matrix to be substracted.
	  * @return an asymmetric matrix with the substraction of components.
	  */
	inline MatrixX operator - (const MatrixS& first, const MatrixX& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on substraction operator: Empty tensor!"));
		assert(first.rows() == second.rows());		// Size of matrices does not correspond!
		assert(first.cols() == second.cols());		// Size of matrices does not correspond!

		MatrixX result(first.rows(), first.cols());
		for (unsigned int i = 0; i < first.rows(); ++i) {
			for (unsigned int j = 0; j < first.cols(); ++j) {
				result.at(i, j) = first.at(i, j) - second.at(i, j);
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Symmetric * Asymmetric matrices
	//
	// ================================================================================================

	/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
	  * @param first Matrix to multiply.
	  * @param second Matrix to be multiplied with.
	  * @return an asymmetric matrix with the multiplication of components.
	  */
	inline MatrixX operator * (const MatrixS& first, const MatrixX& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on multiplication operator: Empty tensor!"));
		assert(first.cols() == second.rows());		// Size of matrices does not correspond!

		MatrixX result(first.rows(), second.cols());
		for (unsigned int i = 0; i < first.rows(); ++i) {
			for (unsigned int j = 0; j < second.cols(); ++j) {
				for (unsigned int k = 0; k < first.cols(); ++k) {
					result.at(i, j) += first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Asymmetric + Symmetric matrices
	//
	// ================================================================================================

	/** Overloads operator + for addition -> R_ij = F_ij + S_ij
	  * @param first First matrix to be added.
	  * @param second Second matrix to be added.
	  * @return an asymmetric matrix with the addition of components.
	  */
	inline MatrixX operator + (const MatrixX& first, const MatrixS& second)
	{
		return (second + first);
	}

	// ================================================================================================
	//
	// Asymmetric - Symmetric matrices
	//
	// ================================================================================================

	/** Overloads operator - for substraction -> R_ij = F_ij - S_ij
	  * @param first First matrix to be substracted.
	  * @param second Second matrix to be substracted.
	  * @return an asymmetric matrix with the substraction of components.
	  */
	inline MatrixX operator - (const MatrixX& first, const MatrixS& second)
	{
		return (second - first) * (-1.);
	}

	// ================================================================================================
	//
	// Asymmetric * Symmetric matrices
	//
	// ================================================================================================

	/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
	  * @param first Matrix to multiply.
	  * @param second Matrix to be multiplied with.
	  * @return an asymmetric matrix with the multiplication of components.
	  */
	inline MatrixX operator * (const MatrixX& first, const MatrixS& second)
	{
		assert(first.cols() == second.rows());		// Size of dyadics does not correspond!

		MatrixX result(first.rows(), second.cols());
		for (unsigned int i = 0; i < first.rows(); ++i) {
			for (unsigned int j = 0; j < second.cols(); ++j) {
				for (unsigned int k = 0; k < first.cols(); ++k) {
					result.at(i, j) += first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
	  * @param first Matrix to multiply.
	  * @param second Matrix to be multiplied with.
	  * @return an asymmetric matrix with the multiplication of components.
	  */
	inline MatrixX operator * (const MatrixS& first, const MatrixS& second)
	{
		if (first.getVector().empty() || second.getVector().empty()) throw std::runtime_error(ERROR("Invalid request on mult operator: Empty tensor!"));
		assert(first.cols() == second.rows());		// Size of matrices does not correspond!

		MatrixX result(first.rows(), second.cols());
		for (unsigned int i = 0; i < first.rows(); i++) {
			for (unsigned int j = 0; j < second.cols(); j++) {
				for (unsigned int k = 0; k < first.cols(); k++) {
					result.at(i, j) += first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Matrix and Scalars
	//
	// ================================================================================================
	inline MatrixX operator + (const double alfa, const MatrixX& dyadic)
	{
		return dyadic + alfa;
	}

	inline MatrixX operator - (const double alfa, const MatrixX& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	inline MatrixX operator * (const double alfa, const MatrixX& dyadic)
	{
		return dyadic * alfa;
	}

	inline MatrixS operator + (const double alfa, const MatrixS& dyadic)
	{
		return dyadic + alfa;
	}

	inline MatrixS operator - (const double alfa, const MatrixS& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	inline MatrixS operator * (const double alfa, const MatrixS& dyadic)
	{
		return dyadic * alfa;
	}

	// ================================================================================================
	//
	// Functions
	//
	// ================================================================================================

	/** @return the determinant of the dyadic.
	  */
	static inline double determinant(const Dyadic2N& T) {
		return T.determinant();
	}

	/** @return the determinant of the dyadic.
	  */
	static inline double determinant(const Dyadic2S& T) {
		return T.determinant();
	}

	/** @return the norm of the dyadic.
	  */
	static inline double norm(const Dyadic2N& T) {
		return T.norm();
	}

	/** @return the norm of the dyadic.
	  */
	static inline double norm(const Dyadic2S& T) {
		return T.norm();
	}

	/** @return the trace of the dyadic.
	  */
	static inline double trace(const Dyadic2N& T) {
		return T.trace();
	}

	/** @return the trace of the dyadic.
	  */
	static inline double trace(const Dyadic2S& T) {
		return T.trace();
	}

	/** @return a vector with the dyadic' eigenvalues.
	  */
	static inline std::vector<double> eigenvalues(const Dyadic2N& T) {
		return T.eigenvalues();
	}

	/** @return a vector with the dyadic' eigenvalues.
	  */
	static inline std::vector<double> eigenvalues(const Dyadic2S& T) {
		return T.eigenvalues();
	}

	/** @return a vector with the dyadic' inverse.
	  */
	static inline std::vector<double> invariants(const Dyadic2N& T) {
		return T.invariants();
	}

	/** @return a vector with the dyadic' inverse.
	  */
	static inline std::vector<double> invariants(const Dyadic2S& T) {
		return T.invariants();
	}

	/** @return the transpose of the dyadic.
	  */
	static inline Dyadic2N transpose(const Dyadic2N& T) {
		return T.transpose();
	}

	/** @return the inverse of the dyadic.
	  */
	static inline Dyadic2N inverse(const Dyadic2N& T) {
		return T.inverse();
	}

	/** @return the inverse of the dyadic.
	  */
	static inline Dyadic2S inverse(const Dyadic2S& T) {
		return T.inverse();
	}

	/** @return the symmetric part of Tensor.
	  */
	static inline Dyadic2S getSymmetric(const Dyadic2N& T) {
		return T.getSymmetric();
	}

	/** @return the antisymmetric part of Tensor.
	  */
	static inline Dyadic2N getAsymmetric(const Dyadic2N& T) {
		return T.getAsymmetric();
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2N& T1, const Dyadic2N& T2) {
		return T1.contraction(T2);
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2S& T1, const Dyadic2S& T2) {
		return T1.contraction(T2);
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2S& T1, const Dyadic2N& T2) {
		assert(T1.rows() == T2.rows());		// Size of dyadics does not correspond!

		double result = 0.;

		for (unsigned int i = 0; i < T1.rows(); ++i) {
			for (unsigned int j = 0; j < T1.cols(); ++j) {
				result += T1.at(i, j) * T2.at(i, j);
			}
		}
		return result;
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2N& T1, const Dyadic2S& T2) {
		return contraction(T2, T1);
	}

	/** @return the transpose of a CSR sparse matrix. Notice that, for symmetric matrices, a row major matrix will become column major matrix.
	  */
	static inline CSR transpose(const CSR& input) {
		CSR output;

		output.mv_rows = input.mv_cols;
		output.mv_cols = input.mv_rows;
		output.mv_sym = input.mv_sym;
		output.mv_value.resize(input.mv_value.size(), 0.0);
		output.mv_colIndex.resize(input.mv_value.size(), 0);
		output.mv_rowIndex.resize(input.mv_cols + 2, 0);

		// Number of itens per column
		for (unsigned int i = 0; i < input.mv_value.size(); ++i) {
			++output.mv_rowIndex.at(input.mv_colIndex.at(i) + 2);
		}

		// Generate new rowIndex
		for (unsigned int i = 2; i < output.mv_rowIndex.size(); ++i) {
			output.mv_rowIndex.at(i) += output.mv_rowIndex.at(i - 1);
		}

		// Index of transposed matrix
		for (unsigned int i = 0; i < input.mv_rows; ++i) {
			for (unsigned int j = input.mv_rowIndex.at(i); j < input.mv_rowIndex.at(i + 1); ++j) {
				const unsigned int newIndex = output.mv_rowIndex.at(input.mv_colIndex.at(j) + 1)++;
				output.mv_value.at(newIndex) = input.mv_value.at(j);
				output.mv_colIndex.at(newIndex) = i;
			}
		}
		output.mv_rowIndex.pop_back(); // Exclude the extra index
		return output;
	}
}  // End of M2S2 namespace
