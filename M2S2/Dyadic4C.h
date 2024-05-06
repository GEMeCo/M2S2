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
#include "Dyadic2N.h"

// ================================================================================================
//
// Dyadic4C class
//
// ================================================================================================
namespace M2S2 {
	/** @class Dyadic4C
	 * @brief Symmetric 4th order tensor for orthotropic constitutive matrices.
	 * @details Symmetric 4th order (rank) tensors of 2 or 3 dimensional vector space (saved in Voigt notation using row major), used for orthotropic constitutive matrices.
	 */
	class Dyadic4C {
	public:
		/** Symmetric 4th order tensor for orthotropic constitutive matrices.
		  */
		Dyadic4C()
		{
			mv_nDim = 0;
			mv_nSize = 0;
			mv_nVoigt = 0;
			mv_Values.clear();
		}

		/** Symmetric 4th order tensor for orthotropic constitutive matrices, in 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space
		  */
		Dyadic4C(unsigned int nDim) : mv_nDim(nDim)
		{
			if (mv_nDim != 2 && mv_nDim != 3) throw std::invalid_argument(ERROR("Invalid argument on Dyadic4C constructor: Wrong input in the size of vector space!"));

			mv_nVoigt = 3 * nDim - 3;
			mv_nSize = mv_nDim * mv_nDim;
			mv_Values.resize(mv_nSize);
		}

		/** Symmetric 4th order tensor for orthotropic constitutive matrices, in 2 or 3 dimensional vector space.
		  * @param value Vector with either 4 or 9 values.
		  */
		Dyadic4C(const std::vector<double>& value)
		{
			mv_nSize = (unsigned int)(value.size());
			if (mv_nSize != 4 && mv_nSize != 9) throw std::invalid_argument(ERROR("Invalid argument on Dyadic4C constructor: Size of input vector is not from 2D or 3D spaces!"));

			mv_nDim = (mv_nSize == 4) ? 2 : 3;
			mv_nVoigt = 3 * mv_nDim - 3;
			mv_Values = { value.begin(), value.begin() + mv_nSize };
		}

		/** Copy constructor for symmetric 4th order tensor for orthotropic constitutive matrices, in 2 or 3 dimensional vector space.
		  * @param other Dyadic to be copied.
		  */
		Dyadic4C(const Dyadic4C& other) {
			mv_nDim = other.mv_nDim;
			mv_nSize = other.mv_nSize;
			mv_nVoigt = other.mv_nVoigt;
			mv_Values = other.mv_Values;
		}

		/** Move constructor for symmetric 4th order tensor for orthotropic constitutive matrices, in 2 or 3 dimensional vector space.
		  * @param other Dyadic to be moved.
		  */
		Dyadic4C(Dyadic4C&& other) noexcept : mv_nDim(other.mv_nDim), mv_nSize(other.mv_nSize), mv_nVoigt(other.mv_nVoigt), mv_Values(std::move(other.mv_Values)) { }

		/** Destructor.
		  */
		~Dyadic4C() { }

		/** Generate an identity with the required dimensionality.
		  * @param nDim Dyadic dimensionality.
		  * @param value Diagonal value. Default is 1.
		  */
		static Dyadic4C identity(unsigned int nDim, const double& value = double(1))
		{
			Dyadic4C result(nDim);
			for (unsigned int i = 0; i < result.mv_nVoigt; ++i) {
				result.at(i, i) = value;
			}
			return result;
		}

		/** Overloads operator << to stream the dyadic. */
		friend std::ostream& operator<<(std::ostream& output, const Dyadic4C& tensor)
		{
			output << tensor.print();
			return output;
		}

		/** Overloads operator >> to stream the dyadic. */
		friend std::istream& operator>>(std::istream& input, Dyadic4C& tensor)
		{
			if (tensor.size() == 0) throw std::out_of_range(ERROR("Out of range in stream operator >>: Unable to fill the Dyadic!"));
			for (unsigned int i = 0; i < tensor.mv_nSize; i++) {
				if (input.eof()) throw std::length_error(ERROR("Length error in stream operator >>: Unable to fill the Dyadic!"));
				input >> tensor.getVector().at(i);
			}
			return input;
		}

		/** Prepare a string to print (to file or screen)
		  * @param precision Number of decimal digits after the decimal point (default is 4)
		  * @param width Minimum number of characters to be written (default is 8)
		  */
		const std::string print(const int precision = 4, const int width = 8) const
		{
			std::ostringstream output;
			output << std::endl;
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				output << "\t" << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, 0);

				for (unsigned int j = 1; j < mv_nVoigt; j++) {
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
		void swap(Dyadic4C& other) noexcept
		{
			std::swap(mv_nDim, other.mv_nDim);
			std::swap(mv_nSize, other.mv_nSize);
			std::swap(mv_nVoigt, other.mv_nVoigt);
			mv_Values.swap(other.mv_Values);
		}

		/** Set all values to zero. Size remains unchanged.
		  */
		void clear()
		{
			memset(&mv_Values[0], 0., mv_Values.size() * sizeof(double));
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& at(unsigned int i, unsigned int j) {
			if (i > mv_nVoigt || j > mv_nVoigt) throw std::out_of_range(ERROR("One of components is out of range in Dyadic4C method .at"));

			if (i < mv_nDim && j < mv_nDim) {
				unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
				return mv_Values.at(pos);
			}
			else if (i < mv_nVoigt && j < mv_nVoigt) {
				if (i == j) { return mv_Values.at(mv_nVoigt - mv_nDim + i); }
				else { return mv_trash; }
			}
			else throw std::out_of_range(ERROR("Out of range in Dyadic4C method .at: Unable to access required element!"));
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& at(unsigned int i, unsigned int j) const {
			if (i > mv_nVoigt || j > mv_nVoigt) throw std::out_of_range(ERROR("One of components is out of range in Dyadic4C method .at"));

			if (i < mv_nDim && j < mv_nDim) {
				unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
				return mv_Values.at(pos);
			}
			else if (i < mv_nVoigt && j < mv_nVoigt) {
				if (i == j) { return mv_Values.at(mv_nVoigt - mv_nDim + i); }
				else { return mg_zero; }
			}
			else throw std::out_of_range(ERROR("Out of range in Dyadic4C method .at: Unable to access required element!"));
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& operator()(unsigned int i, unsigned int j) {
			if (i > mv_nVoigt || j > mv_nVoigt) throw std::out_of_range(ERROR("One of components is out of range in Dyadic4C method .at"));

			if (i < mv_nDim && j < mv_nDim) {
				unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
				return mv_Values.at(pos);
			}
			else if (i < mv_nVoigt && j < mv_nVoigt) {
				if (i == j) { return mv_Values.at(mv_nVoigt - mv_nDim + i); }
				else { return mv_trash; }
			}
			else throw std::out_of_range(ERROR("Out of range in Dyadic4C method .at: Unable to access required element!"));
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& operator()(unsigned int i, unsigned int j) const {
			if (i > mv_nVoigt || j > mv_nVoigt) throw std::out_of_range(ERROR("One of components is out of range in Dyadic4C method .at"));

			if (i < mv_nDim && j < mv_nDim) {
				unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nDim - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nDim - i * 0.5 - 0.5) + j);
				return mv_Values.at(pos);
			}
			else if (i < mv_nVoigt && j < mv_nVoigt) {
				if (i == j) { return mv_Values.at(mv_nVoigt - mv_nDim + i); }
				else { return mg_zero; }
			}
			else throw std::out_of_range(ERROR("Out of range in Dyadic4C method .at: Unable to access required element!"));
		}

		/** @return the row size.
		  */
		unsigned int rows() const noexcept { return mv_nVoigt; }

		/** @return the column size.
		  */
		unsigned int cols() const noexcept { return mv_nVoigt; }

		/** @return the number of itens.
		  */
		unsigned int size() const noexcept { return mv_nSize; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		std::vector<double>& getVector() { return mv_Values; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		const std::vector<double>& getVector() const { return mv_Values; }

		/** @return a iterator to the beginning of the storing vector, providing direct access to the values.
		  */
		std::vector<double>::const_iterator begin() const { return mv_Values.begin(); }

		/** Overloads operator = for copy assignment operations -> T = O
		  * @param other Dyadic to be copied.
		  */
		Dyadic4C& operator=(const Dyadic4C& other)
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
				mv_nVoigt = other.mv_nVoigt;
				mv_Values.resize(mv_nSize);
				mv_Values = other.mv_Values;
			}
			return *this;
		}

		/** Overloads operator = for move assignment operations -> T = &O
		  * @param other Dyadic to be moved.
		  */
		Dyadic4C& operator=(Dyadic4C&& other) noexcept
		{
			if (this != &other) {
				mv_nDim = other.mv_nDim;
				mv_nSize = other.mv_nSize;
				mv_nVoigt = other.mv_nVoigt;
				mv_Values = std::move(other.mv_Values);

				other.mv_nDim = 0;
				other.mv_nSize = 0;
				other.mv_nVoigt = 0;
			}
			return *this;
		}

		/** Overloads operator == for comparison operations
		  * @param other Dyadic to be compared.
		  */
		bool operator==(const Dyadic4C& other)
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

		/** Overloads operator != for comparison operations
		  * @param other Dyadic to be compared.
		  */
		bool operator!=(const Dyadic4C& other)
		{
			if (*this == other) return false;
			return true;
		}


		/** Overloads operator += for cumulative addition -> T += O -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic4C& operator+=(const Dyadic4C& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) += other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator -= for cumulative substraction -> T -= O -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic4C& operator-=(const Dyadic4C& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return *this;
		}

		// Overloading operator *= for cumulative multiplication is not defined
		// Dyadic4C& operator*=(const Dyadic4C& other);

		/** Overloads operator + for addition -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic4C operator+(const Dyadic4C& other) const
		{
			check_order(other);
			Dyadic4C result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) += other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for substraction -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic4C operator-(const Dyadic4C& other) const
		{
			check_order(other);
			Dyadic4C result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return result;
		}

		// Overloading operator * for multiplication is not defined
		// Dyadic4C operator*(const Dyadic4C& other) const;

		/** Overloads operator * for Dot product -> L_ij = T_ijkl * T_kl
		  * @param T Symmetric dyadic to be multiplied with.
		  * @return a symmetric dyadic with the inner product.
		  */
		Dyadic2S operator*(const Dyadic2S& T)
		{
			if (T.size() != mv_nVoigt) throw std::range_error(ERROR("Invalid request in Dyadic4C dot product: Order of tensor are not suitable for this operation!"));

			Dyadic2S result(mv_nDim);
			std::vector<double> R = T.getVoigtMnemonics();
			std::vector<double> L(R.size());

			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					L.at(i) += at(i, j) * R.at(j);
				}
			}
			for (unsigned int i = mv_nDim; i < mv_nVoigt; i++) {
				L.at(i) = at(i, i) * R.at(i);
			}

			// L is written in Voigt Mnemonics. Have to save the other way around to result
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.at(i, i) = L.at(i);
				
				for (unsigned int j = i + 1; j < mv_nDim; j++) {
					result.at(i, j) = L.at(mv_nVoigt - i - j);
				}
			}

			return result;
		}

		/** Overloads operator * for Dot product -> L_ij = T_ijkl * T_kl
		  * @param T Asymmetric dyadic to be multiplied with.
		  * @return a asymmetric dyadic with the inner product.
		  */
		Dyadic2N operator*(const Dyadic2N& T)
		{
			if (T.rows() != mv_nDim) throw std::range_error(ERROR("Invalid request in Dyadic4C dot product: Order of tensor are not suitable for this operation!"));

			Dyadic2N result(mv_nDim);
			std::vector<double> R = T.getVoigtMnemonics();
			std::vector<double> L(R.size());

			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					L.at(i) += at(i, j) * R.at(j);
				}
			}
			for (unsigned int i = mv_nDim; i < mv_nVoigt; i++) {
				L.at(i) = at(i, i) * R.at(i);
				L.at(i + mv_nVoigt - mv_nDim) = at(i, i) * R.at(i + mv_nVoigt - mv_nDim);
			}

			// L is written in Voigt Mnemonics. Have to save the other way around to result
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.at(i, i) = L.at(i);

				for (unsigned int j = i + 1; j < mv_nDim; j++) {
					result.at(i, j) = L.at(mv_nVoigt - i - j);
					result.at(j, i) = L.at(L.size() - i - j);
				}
			}

			return result;
		}

		/** Overloads operator + to sum a scalar -> T_ij = T_ij + alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic4C operator+(const double& alfa) const
		{
			Dyadic4C result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) += alfa;
			}
			return result;
		}

		/** Overloads operator - to substract a scalar -> T_ij = T_ij - alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic4C operator-(const double& alfa) const
		{
			Dyadic4C result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) -= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic4C operator*(const double& alfa)
		{
			Dyadic4C result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic4C operator*(const double& alfa) const
		{
			Dyadic4C result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator / to divide with a non-zero scalar -> T_ij = T_ij / alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic4C operator/(const double& alfa) const
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic4C division operator: Division by zero!"));
			}

			Dyadic4C result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) /= alfa;
			}
			return result;
		}

		/** Overloads operator += to sum a scalar -> T_ij += alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic4C& operator+=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) += alfa;
			}
			return *this;
		}

		/** Overloads operator -= to substract a scalar -> T_ij -= alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic4C& operator-=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) -= alfa;
			}
			return *this;
		}

		/** Overloads operator *= to multiply with a scalar -> T_ij *= alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic4C& operator*=(const double& alfa)
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
		Dyadic4C& operator/=(const double& alfa)
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in Dyadic4C division operator: Division by zero!"));
			}

			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) /= alfa;
			}
			return *this;
		}

	private:
		unsigned int mv_nDim;		// Dimensionality of space (2D or 3D)
		unsigned int mv_nSize;		// Total number of itens
		unsigned int mv_nVoigt;		// Number of rows and columns (always symmetric)
		inline static double mv_trash = 0.;
		std::vector<double> mv_Values;

		inline void check_order(const Dyadic4C& other) const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in Dyadic4C method .check_order: Empty tensor!"));
			assert(mv_nVoigt == other.mv_nVoigt);	// Order of tensors differ
		}
	};
}  // End of M2S2 namespace
