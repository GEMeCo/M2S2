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

// Libraries
#include "Common.h"
#include "Dyadic2S.h"

// ================================================================================================
//
// Dyadic2N class
//
// ================================================================================================
namespace M2S2 {
	/** @class Dyadic2N
	 * @brief NON-Symmetric 2nd order tensor.
	 * @details Asymmetric 2nd order (rank) tensors of 2 or 3 dimensional vector space (saved in Voigt notation using row major).
	 */
	class Dyadic2N {

	// Member Variables
	private:
		unsigned int mv_nDim;		// Dimensionality of space (2D or 3D)
		unsigned int mv_nSize;		// Number of components in Voigt notation
		std::vector<double> mv_Values;

	public:
		/** Asymmetric 2nd order tensor.
		  */
		Dyadic2N()
		{
			mv_nDim = 0;
			mv_nSize = 0;
			mv_Values.clear();
		}

		/** Asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space
		  */
		Dyadic2N(const unsigned int nDim) : mv_nDim(nDim)
		{
			if (mv_nDim != 2 && mv_nDim != 3) throw std::invalid_argument(ERROR("Invalid argument on Dyadic2N constructor: Wrong input in the size of vector space!"));

			mv_nSize = nDim * nDim;
			mv_Values.resize(mv_nSize);
		}

		/** Asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space.
		  * @param value Value to initiate the entire dyadic.
		  */
		Dyadic2N(const unsigned int nDim, const double value) : mv_nDim(nDim)
		{
			if (mv_nDim != 2 && mv_nDim != 3) throw std::invalid_argument(ERROR("Invalid argument on Dyadic2N constructor: Wrong input in the size of vector space!"));

			mv_nSize = nDim * nDim;
			mv_Values.resize(mv_nSize, value);
		}

		/** Asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param value Vector with either 4 or 9 values.
		  */
		Dyadic2N(const std::vector<double>& value)
		{
			if (value.size() != 4 && value.size() != 9) throw std::invalid_argument(ERROR("Invalid argument on Dyadic2N constructor: Size of input vector is not from 2D or 3D spaces!"));

			mv_nDim = (int)(std::sqrt(value.size()));
			mv_nSize = mv_nDim * mv_nDim;
			mv_Values = { value.begin(), value.begin() + mv_nSize };
		}

		/** Copy constructor for asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param other Asymmetric dyadic to be copied.
		  */
		Dyadic2N(const Dyadic2N& other) {
			mv_nDim = other.mv_nDim;
			mv_nSize = other.mv_nSize;
			mv_Values = other.mv_Values;
		}

		/** Move constructor for asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param other Dyadic to be moved.
		  */
		Dyadic2N(Dyadic2N&& other) noexcept : mv_nDim(other.mv_nDim), mv_nSize(other.mv_nSize), mv_Values(std::move(other.mv_Values)) {
			other.mv_nDim = 0;
			other.mv_nSize = 0;
		}

		/** Destructor.
		  */
		~Dyadic2N() = default;

		/** Generate an identity with the required dimensionality.
		  * @param nDim Dyadic dimensionality.
		  * @param value Diagonal value. Default is 1.
		  */
		static Dyadic2N identity(const unsigned int nDim, const double value = double(1))
		{
			Dyadic2N result(nDim);
			for (unsigned int i = 0; i < result.mv_nDim; ++i) {
				result.at(i, i) = value;
			}
			return result;
		}

		/** Overloads operator << to stream the dyadic. */
		friend std::ostream& operator<<(std::ostream& output, const Dyadic2N& tensor)
		{
			output << tensor.print();
			return output;
		}

		/** Overloads operator >> to stream the dyadic. */
		friend std::istream& operator>>(std::istream& input, Dyadic2N& tensor)
		{
			if (tensor.size() == 0) throw std::out_of_range(ERROR("Out of range in stream operator >>: Unable to fill the Dyadic!"));
			for (unsigned int i = 0; i < tensor.rows(); i++) {
				for (unsigned int j = 0; j < tensor.cols(); j++) {
					if (input.eof()) throw std::length_error(ERROR("Length error in stream operator >>: Unable to fill the Dyadic!"));
					input >> tensor.at(i, j);
				}
			}
			return input;
		}

		/** Prepare a string to print (to file or screen)
		  * @param precision Number of decimal digits after the decimal point (default is 4)
		  * @param width Minimum number of characters to be written (default is 8)
		  */
		const std::string print(const unsigned int precision = 4, const unsigned int width = 8) const
		{
			std::ostringstream output;
			output << std::endl;
			for (unsigned int i = 0; i < mv_nDim; i++) {
				output << "\t" << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, 0);

				for (unsigned int j = 1; j < mv_nDim; j++) {
					output << " " << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, j);
				}
				output << "\n";
			}
			output << std::endl;
			return output.str();
		}

		/** Swap dyadics
		  * @param other Dyadic to be swaped.
		  */
		void swap(Dyadic2N& other) noexcept
		{
			std::swap(mv_nDim, other.mv_nDim);
			std::swap(mv_nSize, other.mv_nSize);
			mv_Values.swap(other.mv_Values);
		}

		/** Set all values to zero. Size remains unchanged.
		  */
		void clear()
		{
			memset(&mv_Values[0], 0., mv_Values.size() * sizeof(double));
		}

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& at(const unsigned int i, const unsigned int j) {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2N method .at"));

			unsigned int pos = i * mv_nDim + j;
			if (pos > mv_nSize) throw std::out_of_range(ERROR("Out of range in Dyadic2N method .at: Unable to access required element!"));
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& at(const unsigned int i, const unsigned int j) const {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2N method .at"));

			unsigned int pos = i * mv_nDim + j;
			if (pos > mv_nSize) throw std::out_of_range(ERROR("Out of range in Dyadic2N method .at: Unable to access required element!"));
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& operator () (const unsigned int i, const unsigned int j) {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2N method .at"));

			unsigned int pos = i * mv_nDim + j;
			if (pos > mv_nSize) throw std::out_of_range(ERROR("Out of range in Dyadic2N method .at: Unable to access required element!"));
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& operator () (const unsigned int i, const unsigned int j) const {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2N method .at"));

			unsigned int pos = i * mv_nDim + j;
			if (pos > mv_nSize) throw std::out_of_range(ERROR("Out of range in Dyadic2N method .at: Unable to access required element!"));
			return mv_Values.at(pos);
		}

		/** @return the row size.
		  */
		unsigned int rows() const noexcept { return mv_nDim; }

		/** @return the column size.
		  */
		unsigned int cols() const noexcept { return mv_nDim; }

		/** @return the number of components in Voigt notation.
		  */
		unsigned int size() const noexcept { return mv_nSize; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		std::vector<double>& getVector() { return mv_Values; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		const std::vector<double>& getVector() const { return mv_Values; }

		/** @return current tensor written as vector in Voigt notation with mnemonics rule.
		  * V(1) = T(1,1); V(2) = T(2,2); V(3) = T(3,3); V(4) = T(2,3); V(5) = T(1,3); V(6) = T(1,2); V(7) = T(3,2); V(8) = T(3,1); V(9) = T(2,1)
		  */
		std::vector<double> getVoigtMnemonics() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .getVoigtMnemonics: Empty tensor!"));

			unsigned int mi_nVoigt = 3 * mv_nDim - 3;
			std::vector<double> result(mv_nSize);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.at(i) = at(i, i);
				for (unsigned int j = i + 1; j < mv_nDim; ++j) {
					result.at(mi_nVoigt - i - j) = at(i, j);
					result.at(mv_nSize - i - j) = at(j, i);
				}
			}
			return result;
		}

		/** @return a iterator to the beginning of the storing vector, providing direct access to the values.
		  */
		std::vector<double>::const_iterator begin() const { return mv_Values.begin(); }

		/** @return a iterator to the end of the storing vector, providing direct access to the values.
		  */
		std::vector<double>::const_iterator end() const { return mv_Values.end(); }

		/** Overloads operator = for copy assignment operations -> T = O
		  * @param other Dyadic to be copied.
		  */
		Dyadic2N& operator = (const Dyadic2N& other)
		{
			if (this == &other) {
				return *this;
			}

			if (mv_nDim == other.mv_nDim)
			{
				mv_Values = other.mv_Values;
			}
			else
			{
				mv_nDim = other.mv_nDim;
				mv_nSize = other.mv_nSize;
				mv_Values.resize(mv_nSize);
				mv_Values = other.mv_Values;
			}
			return *this;
		}

		/** Overloads operator = for move assignment operations -> T = &O
		  * @param other Dyadic to be moved.
		  */
		Dyadic2N& operator = (Dyadic2N&& other) noexcept
		{
			if (this != &other) {
				mv_nDim = other.mv_nDim;
				mv_nSize = other.mv_nSize;
				mv_Values = std::move(other.mv_Values);

				other.mv_nDim = 0;
				other.mv_nSize = 0;
			}
			return *this;
		}

		/** Overloads operator == for comparison operations
		  * @param other Dyadic to be compared.
		  */
		bool operator == (const Dyadic2N& other) const
		{
			bool mi_check = true;
			if (mv_nDim == other.mv_nDim && mv_nSize == other.mv_nSize) {
				for (int i = 0; i < mv_nSize; ++i) {
					if (!almost_equal(mv_Values.at(i), other.mv_Values.at(i))) {
						mi_check = false;
						break;
					}
				}
				return mi_check;
			}
			return false;
		}

		/** Overloads operator == for comparison operations
		  * @param other Dyadic to be compared.
		  */
		bool operator == (const Dyadic2S& other) const
		{
			bool mi_check = true;
			if (rows() == other.rows() && cols() == other.cols()) {
				for (int i = 0; i < rows(); ++i) {
					for (int j = 0; j < cols(); ++j) {
						if (!almost_equal(at(i, j), other.at(i, j))) {
							mi_check = false;
							break;
						}
					}
				}
				return mi_check;
			}
			return false;
		}

		/** Overloads operator != for comparison operations
		  * @param other Dyadic to be compared.
		  */
		bool operator != (const Dyadic2N& other) const
		{
			if (*this == other) return false;
			return true;
		}

		/** Overloads operator != for comparison operations
		  * @param other Dyadic to be compared.
		  */
		bool operator != (const Dyadic2S& other) const
		{
			if (*this == other) return false;
			return true;
		}

		/** Overloads operator += for cumulative addition -> T += O -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic2N& operator += (const Dyadic2N& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) += other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator -= for cumulative substraction -> T -= O -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2N& operator -= (const Dyadic2N& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator *= for cumulative multiplication -> T *= O -> T_ij = T_ik * O_kj
		  * @param other Dyadic to be multiplied with.
		  */
		Dyadic2N& operator *= (const Dyadic2N& other)
		{
			check_order(other);
			Dyadic2N result(mv_nDim);
			auto& values = result.getVector();
			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					for (unsigned int k = 0; k < mv_nDim; k++) {
						values.at(i * mv_nDim + j) += mv_Values.at(i * mv_nDim + k) * other.mv_Values.at(k * mv_nDim + j);
					}
				}
			}
			swap(result);
			return *this;
		}

		/** Overloads operator + for addition -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic2N operator + (const Dyadic2N& other) const
		{
			check_order(other);
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) += other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for substraction -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2N operator - (const Dyadic2N& other) const
		{
			check_order(other);
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for unary negation -> T = -T
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2N operator - () const
		{
			Dyadic2N result(this->mv_nDim);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= this->mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
		  * @param other Dyadic to be multiplied with.
		  */
		Dyadic2N operator * (const Dyadic2N& other) const
		{
			check_order(other);
			Dyadic2N result(mv_nDim);
			auto& values = result.getVector();
			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					for (unsigned int k = 0; k < mv_nDim; k++) {
						values.at(i * mv_nDim + j) += mv_Values.at(i * mv_nDim + k) * other.mv_Values.at(k * mv_nDim + j);
					}
				}
			}
			return result;
		}

		/** Overloads operator * for Dot product -> V_i = T_ij * Vo_j
		  * @param vec Vector to be multiplied with.
		  * @return a vector with the inner product.
		  */
		std::vector<double> operator * (const std::vector<double>& vec)
		{
			if (vec.size() != mv_nDim) throw std::range_error(ERROR("Invalid request in Dyadic2N dot product: Order of tensor and vector differ!"));

			std::vector<double> result(mv_nDim);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					result.at(i) += at(i, j) * vec.at(j);
				}
			}
			return result;
		}

		/** Overloads operator + to sum a scalar -> T_ij = T_ij + alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator + (const double alfa) const
		{
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.mv_Values.at(i) += alfa;
			}
			return result;
		}

		/** Overloads operator - to substract a scalar -> T_ij = T_ij - alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator - (const double alfa) const
		{
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.mv_Values.at(i) -= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator * (const double alfa) const
		{
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator / to divide with a non-zero scalar -> T_ij = T_ij / alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator / (const double alfa) const
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic2N division operator: Division by zero!"));
			}

			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				result.mv_Values.at(i) /= alfa;
			}
			return result;
		}

		/** Overloads operator += to sum a scalar -> T_ij += alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator += (const double alfa)
		{
			for (unsigned int i = 0; i < mv_nDim; i++) {
				mv_Values.at(i) += alfa;
			}
			return *this;
		}

		/** Overloads operator -= to substract a scalar -> T_ij -= alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator -= (const double alfa)
		{
			for (unsigned int i = 0; i < mv_nDim; i++) {
				mv_Values.at(i) -= alfa;
			}
			return *this;
		}

		/** Overloads operator *= to multiply with a scalar -> T_ij *= alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator *= (const double alfa)
		{
			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) *= alfa;
			}
			return *this;
		}

		/** Overloads operator /= to divide with a non-zero scalar -> T_ij /= alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator /= (const double alfa)
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic2N division operator: Division by zero!"));
			}

			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) /= alfa;
			}
			return *this;
		}

		/** @return the determinant of the dyadic.
		  */
		double determinant() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .determinant: Empty tensor!"));

			if (mv_nDim == 2) {
				return mv_Values.at(0) * mv_Values.at(3) - mv_Values.at(1) * mv_Values.at(2);
			}
			else {
				return mv_Values.at(0) * mv_Values.at(4) * mv_Values.at(8) +
					mv_Values.at(1) * mv_Values.at(5) * mv_Values.at(6) +
					mv_Values.at(2) * mv_Values.at(3) * mv_Values.at(7) -
					mv_Values.at(0) * mv_Values.at(5) * mv_Values.at(7) -
					mv_Values.at(1) * mv_Values.at(3) * mv_Values.at(8) -
					mv_Values.at(2) * mv_Values.at(4) * mv_Values.at(6);
			}
		}

		/** @return the Frobenius norm of the dyadic, the square root of the sum of absolute square of all components -> square root of T:T.
		  */
		double norm() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .norm: Empty tensor!"));

			return std::sqrt(contraction(*this));
		}

		/** @return the trace of the dyadic -> T_ii
		  */
		double trace() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .trace: Empty tensor!"));

			double result = 0;
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				result += at(i, i);
			}
			return result;
		}

		/** @return a vector with the dyadic' eigenvalues.
		  */
		std::vector<double> eigenvalues() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .eigenvalues: Empty tensor!"));

			std::cout << WARN("Eigenvalues of symmetric part only.\tOnly those will be always real.");
			auto TS = getSymmetric();
			std::vector<double> result = TS.eigenvalues();
			return result;
		}

		/** @return a vector with the dyadic' invariants.
		  */
		std::vector<double> invariants() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .invariants: Empty tensor!"));

			std::vector<double> result(mv_nDim);
			if (mv_nDim == 2) {
				result.at(0) = trace();
				result.at(1) = determinant();
			}
			else {
				result.at(0) = trace();
				result.at(1) = mv_Values.at(0) * mv_Values.at(4) + mv_Values.at(4) * mv_Values.at(8) + mv_Values.at(8) * mv_Values.at(0)
					- mv_Values.at(1) * mv_Values.at(3) - mv_Values.at(5) * mv_Values.at(7) - mv_Values.at(2) * mv_Values.at(6);
				result.at(2) = determinant();
			}
			return result;
		}

		/** @return the transpose of the dyadic -> T^t
		  */
		Dyadic2N transpose() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .transpose: Empty tensor!"));

			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				for (unsigned int j = i + 1; j < mv_nDim; ++j) {
					result.at(i, j) = at(j, i);
					result.at(j, i) = at(i, j);
				}
			}
			return result;
		}

		/** @return the inverse of the dyadic -> T^-1
		  */
		Dyadic2N inverse() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .invariants: Empty tensor!"));

			Dyadic2N result(mv_nDim);
			double det = determinant();
			if ((int)(det * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic2N method .inverse: Tensor is singular!"));
			}
			det = 1. / det;

			if (mv_nDim == 2) {
				result.at(0, 0) = det * at(1, 1);
				result.at(0, 1) = -1. * det * at(0, 1);
				result.at(1, 0) = -1. * det * at(1, 0);
				result.at(1, 1) = det * at(0, 0);
			}
			else
			{
				// For non-symmetric tensors
				result.at(0, 0) = at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1);
				result.at(0, 1) = at(0, 2) * at(2, 1) - at(0, 1) * at(2, 2);
				result.at(0, 2) = at(0, 1) * at(1, 2) - at(0, 2) * at(1, 1);
				result.at(1, 0) = at(1, 2) * at(2, 0) - at(1, 0) * at(2, 2);
				result.at(1, 1) = at(0, 0) * at(2, 2) - at(0, 2) * at(2, 0);
				result.at(1, 2) = at(0, 2) * at(1, 0) - at(0, 0) * at(1, 2);
				result.at(2, 0) = at(1, 0) * at(2, 1) - at(1, 1) * at(2, 0);
				result.at(2, 1) = at(0, 1) * at(2, 0) - at(0, 0) * at(2, 1);
				result.at(2, 2) = at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);

				result *= det;
			}
			return result;
		}

		/** @return the symmetric part of Tensor -> Sym(T)
		  */
		Dyadic2S getSymmetric() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .getSymmetric: Empty tensor!"));

			Dyadic2S result(mv_nDim);
			auto& values = result.getVector();

			if (mv_nDim == 2) {
				values.at(0) = mv_Values.at(0);
				values.at(1) = 0.5 * (mv_Values.at(1) + mv_Values.at(2));
				values.at(2) = mv_Values.at(3);
			}
			else
			{
				values.at(0) = mv_Values.at(0);
				values.at(1) = 0.5 * (mv_Values.at(1) + mv_Values.at(3));
				values.at(2) = 0.5 * (mv_Values.at(2) + mv_Values.at(6));
				values.at(3) = mv_Values.at(4);
				values.at(4) = 0.5 * (mv_Values.at(5) + mv_Values.at(7));
				values.at(5) = mv_Values.at(8);
			}
			return result;
		}

		/** @return the antisymmetric (skew) part of Tensor -> Asym(T)
		  */
		Dyadic2N getAsymmetric() const
		{
			Dyadic2S Sym = this->getSymmetric();
			Dyadic2N result(this->mv_nDim);
			for (unsigned int i = 0; i < result.rows(); ++i) {
				for (unsigned int j = 0; j < result.cols(); ++j) {
					result.at(i, j) = at(i, j) - Sym.at(i, j);
				}
			}
			return result;
		}

		/** @return the product of the transposed tensor and itself (always symmetric) -> T1^T * T1
		  */
		Dyadic2S getATA() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .getATA: Empty tensor!"));

			Dyadic2S result(mv_nDim);
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				for (unsigned int j = i; j < mv_nDim; ++j) {
					for (unsigned int k = 0; k < mv_nDim; ++k) {
						result.at(i, j) += mv_Values.at(k * mv_nDim + i) * mv_Values.at(k * mv_nDim + j);
					}
				}
			}
			return result;
		}

		/** @return the spherical part of Tensor
		  */
		Dyadic2S getSpheric() const
		{
			double mi_diag = trace() / 3.;

			return Dyadic2S::identity(mv_nDim, mi_diag);
		}

		/** @return the deviator part of Tensor
		  */
		Dyadic2N getDeviator() const
		{
			double mi_diag = trace() / 3.;
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < result.rows(); ++i) {
				result.at(i, i) -= mi_diag;
			}
			return result;
		}

		/** @return the double-dot product of the dyadics (contraction) -> T1:T2.
		  * @param other Dyadic to be used on the double-dot product.
		  */
		double contraction(const Dyadic2N& other) const
		{
			check_order(other);
			double result = 0.;
			for (unsigned int i = 0; i < mv_nSize; ++i) {
				result += mv_Values.at(i) * other.mv_Values.at(i);
			}
			return result;
		}

	private:

		inline void check_order(const Dyadic2N& other) const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2N method .check_order: Empty tensor!"));
			assert(mv_nDim == other.mv_nDim);	// Order of tensors differ
		}
	};
}  // End of M2S2 namespace
