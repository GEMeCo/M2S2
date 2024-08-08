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

// ================================================================================================
//
// Dyadic2S class
//
// ================================================================================================
namespace M2S2 {
	/** @class Dyadic2S
	 * @brief Symmetric 2nd order tensor.
	 * @details Symmetric 2nd order (rank) tensors of 2 or 3 dimensional vector space (saved in Voigt notation using row major).
	 */
	class Dyadic2S {

	// Member Variables
	private:
		unsigned int mv_nDim;		// Dimensionality of space (2D or 3D)
		unsigned int mv_nVoigt;		// Number of components in Voigt notation
		std::vector<double> mv_Values;

	public:
		/** Symmetric 2nd order tensor.
		  */
		Dyadic2S()
		{
			mv_nDim = 0;
			mv_nVoigt = 0;
			mv_Values.clear();
		}

		/** Symmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space
		  */
		Dyadic2S(const unsigned int nDim) : mv_nDim(nDim)
		{
			if (mv_nDim != 2 && mv_nDim != 3) throw std::invalid_argument(ERROR("Invalid argument on Dyadic2S constructor: Wrong input in the size of vector space!"));

			mv_nVoigt = 3 * nDim - 3;
			mv_Values.resize(mv_nVoigt);
		}

		/** Symmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space.
		  * @param value Value to initiate the entire dyadic.
		  */
		Dyadic2S(const unsigned int nDim, const double value) : mv_nDim(nDim)
		{
			if (mv_nDim != 2 && mv_nDim != 3) throw std::invalid_argument(ERROR("Invalid argument on Dyadic2S constructor: Wrong input in the size of vector space!"));

			mv_nVoigt = 3 * nDim - 3;
			mv_Values.resize(mv_nVoigt, value);
		}

		/** Symmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param value Vector with either 3 or 6 values.
		  */
		Dyadic2S(const std::vector<double>& value)
		{
			if(value.size() != 3 && value.size() != 6) throw std::invalid_argument(ERROR("Invalid argument on Dyadic2S constructor: Size of input vector is not from 2D or 3D spaces!"));

			mv_nDim = (int)((value.size() + 3) / 3);
			mv_nVoigt = 3 * mv_nDim - 3;
			mv_Values = { value.begin(), value.begin() + mv_nVoigt };
		}

		/** Copy constructor for symmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param other Dyadic to be copied.
		  */
		Dyadic2S(const Dyadic2S& other) {
			mv_nDim = other.mv_nDim;
			mv_nVoigt = other.mv_nVoigt;
			mv_Values = other.mv_Values;
		}

		/** Move constructor for symmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param other Dyadic to be moved.
		  */
		Dyadic2S(Dyadic2S&& other) noexcept : mv_nDim(other.mv_nDim), mv_nVoigt(other.mv_nVoigt), mv_Values(std::move(other.mv_Values)) {
			other.mv_nDim = 0;
			other.mv_nVoigt = 0;
		}

		/** Destructor.
		  */
		~Dyadic2S() = default;

		/** Generate an identity with the required dimensionality.
		  * @param nDim Dyadic dimensionality.
		  * @param value Diagonal value. Default is 1.
		  */
		static Dyadic2S identity(const unsigned int nDim, const double value = double(1))
		{
			Dyadic2S result(nDim);
			for (unsigned int i = 0; i < result.mv_nDim; ++i) {
				result.at(i, i) = value;
			}
			return result;
		}

		/** Overloads operator << to stream the dyadic. */
		friend std::ostream& operator<<(std::ostream& output, const Dyadic2S& tensor)
		{
			output << tensor.print();
			return output;
		}
		
		/** Overloads operator >> to stream the dyadic. */
		friend std::istream& operator>>(std::istream& input, Dyadic2S& tensor)
		{
			if(tensor.size() == 0) throw std::out_of_range(ERROR("Out of range in stream operator >>: Unable to fill the Dyadic!"));
			for (unsigned int i = 0; i < tensor.rows(); i++) {
				for (unsigned int j = i; j < tensor.cols(); j++) {
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
		void swap(Dyadic2S& other) noexcept
		{
			std::swap(mv_nDim, other.mv_nDim);
			std::swap(mv_nVoigt, other.mv_nVoigt);
			mv_Values.swap(other.mv_Values);
		}

		/** Set all values to zero. Size remains unchanged.
		  */
		void clear()
		{
			memset(&mv_Values[0], 0., mv_Values.size() * sizeof(double));
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& at(const unsigned int i, const unsigned int j) {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2S method .at"));

			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
			if(pos > mv_nVoigt) throw std::out_of_range(ERROR("Out of range in Dyadic2S method .at: Unable to access required element!"));
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& at(const unsigned int i, const unsigned int j) const {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2S method .at"));

			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
			if (pos > mv_nVoigt) throw std::out_of_range(ERROR("Out of range in Dyadic2S method .at: Unable to access required element!"));
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& operator () (const unsigned int i, const unsigned int j) {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2S method .at"));

			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
			if (pos > mv_nVoigt) throw std::out_of_range(ERROR("Out of range in Dyadic2S method .at: Unable to access required element!"));
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& operator () (const unsigned int i, const unsigned int j) const {
			if (i > mv_nDim || j > mv_nDim) throw std::out_of_range(ERROR("One of components is out of range in Dyadic2S method .at"));

			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
			if (pos > mv_nVoigt) throw std::out_of_range(ERROR("Out of range in Dyadic2S method .at: Unable to access required element!"));
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
		unsigned int size() const noexcept { return mv_nVoigt; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		std::vector<double>& getVector() { return mv_Values; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		const std::vector<double>& getVector() const { return mv_Values; }

		/** @return current tensor written as vector in Voigt notation with mnemonics rule.
		  * V(1) = T(1,1); V(2) = T(2,2); V(3) = T(3,3); V(4) = T(2,3); V(5) = T(1,3); V(6) = T(1,2)
		  */
		std::vector<double> getVoigtMnemonics() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .getVoigtMnemonics: Empty tensor!"));

			std::vector<double> result(mv_nVoigt);
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				result.at(i) = at(i, i);
				for (unsigned int j = i + 1; j < mv_nDim; ++j) {
					result.at(mv_nVoigt - i - j) = at(i, j);
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

		/** Overloads operator = for copy assignment operations
		  * @param other Dyadic to be copied.
		  */
		Dyadic2S& operator = (const Dyadic2S& other)
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
				mv_nVoigt = other.mv_nVoigt;
				mv_Values.resize(mv_nVoigt);
				mv_Values = other.mv_Values;
			}
			return *this;
		}

		/** Overloads operator = for move assignment operations
		  * @param other Dyadic to be moved.
		  */
		Dyadic2S& operator = (Dyadic2S&& other) noexcept
		{
			if (this != &other) {
				mv_nDim = other.mv_nDim;
				mv_nVoigt = other.mv_nVoigt;
				mv_Values = std::move(other.mv_Values);

				other.mv_nDim = 1;
				other.mv_nVoigt = 1;
			}
			return *this;
		}

		/** Overloads operator == for comparison operations
		  * @param other Dyadic to be compared.
		  */
		bool operator == (const Dyadic2S& other) const
		{
			bool mi_check = true;
			if(mv_nDim == other.mv_nDim && mv_nVoigt == other.mv_nVoigt) {
				for (int i = 0; i < mv_nVoigt; ++i) {

					if (!almost_equal(mv_Values.at(i), other.mv_Values.at(i))) {
						mi_check = false;
						break;
					}
				}
				return mi_check;
			}
			return false;
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
		Dyadic2S& operator += (const Dyadic2S& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				mv_Values.at(i) += other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator -= for cumulative substraction -> T -= O -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2S& operator -= (const Dyadic2S& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return *this;
		}

		// Overload operator *= for cumulative multiplication is defined elsewhere
		// The return value may be asymetric
		// Dyadic2N& operator *= (const Dyadic2S& other);

		/** Overloads operator + for addition -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic2S operator + (const Dyadic2S& other) const
		{
			check_order(other);
			Dyadic2S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) += other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for substraction -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2S operator - (const Dyadic2S& other) const
		{
			check_order(other);
			Dyadic2S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for unary negation -> T = -T
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2S operator - () const
		{
			Dyadic2S result(this->mv_nDim);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= this->mv_Values.at(i);
			}
			return result;
		}

		// Overload operator * for multiplication is defined elsewhere
		// The return value may be asymetric
		// Dyadic2N operator * (const Dyadic2S& other) const;

		/** Overloads operator * for Dot product -> V_i = T_ij * Vo_j
		  * @param vec Vector to be multiplied with.
		  * @return a vector with the inner product.
		  */
		std::vector<double> operator * (const std::vector<double>& vec)
		{
			if (vec.size() != mv_nDim) throw std::range_error(ERROR("Invalid request in Dyadic2S dot product: Order of tensor and vector differ!"));

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
		Dyadic2S operator + (const double alfa) const
		{
			Dyadic2S result(*this);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.mv_Values.at(i) += alfa;
			}
			return result;
		}

		/** Overloads operator - to substract a scalar -> T_ij = T_ij - alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic2S operator - (const double alfa) const
		{
			Dyadic2S result(*this);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.mv_Values.at(i) -= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic2S operator * (const double alfa) const
		{
			Dyadic2S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator / to divide with a non-zero scalar -> T_ij = T_ij / alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic2S operator / (const double alfa) const
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic2S division operator: Division by zero!"));
			}

			Dyadic2S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) /= alfa;
			}
			return result;
		}

		/** Overloads operator += to sum a scalar -> T_ij += alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic2S& operator += (const double alfa)
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
		Dyadic2S& operator -= (const double alfa)
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
		Dyadic2S& operator *= (const double alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) *= alfa;
			}
			return *this;
		}

		/** Overloads operator /= to divide with a non-zero scalar -> T_ij /= alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic2S& operator /= (const double alfa)
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic2S division operator: Division by zero!"));
			}

			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) /= alfa;
			}
			return *this;
		}

		/** @return the determinant of the dyadic.
		  */
		double determinant() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .determinant: Empty tensor!"));

			if (mv_nDim == 2) {
				return mv_Values.at(0) * mv_Values.at(2) - mv_Values.at(1) * mv_Values.at(1);
			}
			else {
				return mv_Values.at(0) * mv_Values.at(3) * mv_Values.at(5) +
					2. * mv_Values.at(1) * mv_Values.at(2) * mv_Values.at(4) -
					mv_Values.at(0) * mv_Values.at(4) * mv_Values.at(4) -
					mv_Values.at(1) * mv_Values.at(1) * mv_Values.at(5) -
					mv_Values.at(2) * mv_Values.at(2) * mv_Values.at(3);
			}
		}

		/** @return the Frobenius norm of the dyadic, the square root of the sum of absolute square of all components -> square root of T:T.
		  */
		double norm() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .norm: Empty tensor!"));

			return std::sqrt(contraction(*this));
		}

		/** @return the trace of the dyadic -> T_ii.
		  */
		double trace() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .trace: Empty tensor!"));

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
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .eigenvalues: Empty tensor!"));

			std::vector<double> result(mv_nDim);
			if (mv_nDim == 2) {
				double delta = std::sqrt(
					(mv_Values.at(0) + mv_Values.at(2)) * (mv_Values.at(0) + mv_Values.at(2))
					- 4 * (mv_Values.at(0) * mv_Values.at(2) - mv_Values.at(1) * mv_Values.at(1)));
				result.at(0) = 0.5 * (mv_Values.at(0) + mv_Values.at(2) - delta);
				result.at(1) = 0.5 * (mv_Values.at(0) + mv_Values.at(2) + delta);
			}
			else {
				result = invariants();
				double P, Q, R, S, T, A;

				P = result.at(0) / 3.;
				Q = P * result.at(1) - result.at(2) - 2. * P * P * P;
				R = 3. * P * P - result.at(1);
				S = std::sqrt(R / 3.);
				T = R * S / 3.;
				A = std::acos(-0.5 * Q / T);

				result.at(0) = P + 2. * S * std::cos(A / 3.);
				result.at(1) = P + 2. * S * std::cos(A / 3. + 2. * MV_PI / 3.);
				result.at(2) = P + 2. * S * std::cos(A / 3. + 4. * MV_PI / 3.);
			}
			return result;
		}

		/** @return a vector with the dyadic' invariants.
		  */
		std::vector<double> invariants() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .invariants: Empty tensor!"));

			std::vector<double> result(mv_nDim);
			if (mv_nDim == 2) {
				result.at(0) = trace();
				result.at(1) = determinant();
			}
			else {
				result.at(0) = trace();
				result.at(1) = mv_Values.at(0) * mv_Values.at(3) + mv_Values.at(3) * mv_Values.at(5)
					+ mv_Values.at(5) * mv_Values.at(0) - mv_Values.at(1) * mv_Values.at(1)
					- mv_Values.at(2) * mv_Values.at(2) - mv_Values.at(4) * mv_Values.at(4);
				result.at(2) = determinant();
			}
			return result;
		}

		/** @return the inverse of the dyadic -> T^-1.
		  */
		Dyadic2S inverse() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .inverse: Empty tensor!"));

			Dyadic2S result(mv_nDim);
			double det = determinant();
			if ((int)(det * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic2S method .inverse: Tensor is singular!"));
			}
			det = 1. / det;

			if (mv_nDim == 2) {
				result.at(0, 0) = det * at(1, 1);
				result.at(0, 1) = -1. * det * at(0, 1);
				result.at(1, 1) = det * at(0, 0);
			}
			else
			{
				// For symmetric tensors
				result.at(0, 0) = at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1);
				result.at(0, 1) = at(0, 2) * at(2, 1) - at(0, 1) * at(2, 2);
				result.at(0, 2) = at(0, 1) * at(1, 2) - at(0, 2) * at(1, 1);
				//result.at(1, 0) = at(1, 2) * at(2, 0) - at(1, 0) * at(2, 2);
				result.at(1, 1) = at(0, 0) * at(2, 2) - at(0, 2) * at(2, 0);
				result.at(1, 2) = at(0, 2) * at(1, 0) - at(0, 0) * at(1, 2);
				//result.at(2, 0) = at(1, 0) * at(2, 1) - at(1, 1) * at(2, 0);
				//result.at(2, 1) = at(0, 1) * at(2, 0) - at(0, 0) * at(2, 1);
				result.at(2, 2) = at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);

				result *= det;
			}
			return result;
		}

		/** @return the product of the transposed tensor and itself -> T1^T * T1.
		  */
		Dyadic2S getATA() const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .getATA: Empty tensor!"));

			Dyadic2S result(mv_nDim);
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				for (unsigned int j = i; j < mv_nDim; ++j) {
					for (unsigned int k = 0; k < mv_nDim; ++k) {
						result.at(i, j) += at(i, k) * at(k, j);
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
		Dyadic2S getDeviator() const
		{
			double mi_diag = trace() / 3.;
			Dyadic2S result(*this);
			for (unsigned int i = 0; i < result.rows(); ++i) {
				result.at(i, i) -= mi_diag;
			}
			return result;
		}

		/** @return the double-dot product of the dyadics (contraction - sum of product of all entries) -> T1:T2.
		  * @param other Dyadic to be used on the double-dot product.
		  */
		double contraction(const Dyadic2S& other) const
		{
			check_order(other);

			double result = 0.;
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				for (unsigned int j = 0; j < mv_nDim; ++j) {
					result += at(i, j) * other.at(i, j);
				}
			}
			return result;
		}

	private:

		inline void check_order(const Dyadic2S& other) const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic2S method .check_order: Empty tensor!"));
			assert(mv_nDim == other.mv_nDim);	// Order of tensors differ
		}
	};
}  // End of M2S2 namespace
