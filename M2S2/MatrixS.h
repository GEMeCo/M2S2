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
// MatrixS class
//
// ================================================================================================
namespace M2S2 {
	/** @class MatrixS
	 * @brief Symmetric square matrix of any order.
	 * @details Symmetric square matrix of any order (saved as row major).
	 */
	class MatrixS {
	public:
		/** Symmetric square matrix of any order.
		  */
		MatrixS()
		{
			mv_nSize = 0;
			mv_numItens = 0;
			mv_Values.clear();
		}

		/** Symmetric square matrix of any order.
		  * @param nSize Number of columns / rows of the square matrix
		  */
		MatrixS(unsigned int nSize) : mv_nSize(nSize)
		{
			mv_numItens = (unsigned int)(0.5 * mv_nSize * mv_nSize + 0.5 * mv_nSize);
			mv_Values.resize(mv_numItens);
		}

		/** Symmetric square matrix of any order.
		  * @param nSize Number of columns / rows of the square matrix
		  * @param value Value to initiate the entire matrix.
		  */
		MatrixS(unsigned int nSize, const double& value) : mv_nSize(nSize)
		{
			mv_numItens = (unsigned int)(0.5 * mv_nSize * mv_nSize + 0.5 * mv_nSize);
			mv_Values.resize(mv_numItens, value);
		}

		/** Symmetric square matrix of any order.
		  * @param nSize Number of columns / rows of the square matrix
		  * @param value Vector with enough values to initiate the matrix.
		  */
		MatrixS(unsigned int nSize, const std::vector<double>& value) : mv_nSize(nSize)
		{
			mv_numItens = (unsigned int)(0.5 * mv_nSize * mv_nSize + 0.5 * mv_nSize);
			if (value.size() != mv_numItens) throw std::invalid_argument(ERROR("Invalid argument on MatrixS constructor: Size of input vector does not correspond to required matrix!"));
			mv_Values = { value.begin(), value.begin() + mv_numItens };
		}

		/** Copy constructor for symmetric square matrix of any order.
		  * @param other Matrix to be copied.
		  */
		MatrixS(const MatrixS& other) {
			mv_nSize = other.mv_nSize;
			mv_numItens = other.mv_numItens;
			mv_Values = other.mv_Values;
		}

		/** Move constructor for symmetric square matrix of any order.
		  * @param other Matrix to be moved.
		  */
		MatrixS(MatrixS&& other) noexcept : mv_nSize(other.mv_nSize), mv_numItens(other.mv_numItens), mv_Values(std::move(other.mv_Values)) { }

		/** Destructor.
		  */
		~MatrixS() { }

		/** Generate an identity with the required dimensionality.
		  * @param nOrder Matrix dimensionality.
		  * @param value Diagonal value. Default is 1.
		  */
		static MatrixS identity(unsigned int nOrder, const double& value = double(1))
		{
			MatrixS result(nOrder);
			auto& R = result.getVector();

			for (unsigned int i = 0; i < result.mv_nSize; ++i) {
				R.at((unsigned int)(i * (result.mv_nSize - i * 0.5 + 0.5))) = value;
			}
			return result;
		}

		/** Overloads operator << to stream the matrix. */
		friend std::ostream& operator<<(std::ostream& output, const MatrixS& matrix)
		{
			output << matrix.print();
			return output;
		}

		/** Overloads operator >> to stream the matrix. */
		friend std::istream& operator>>(std::istream& input, MatrixS& matrix)
		{
			if (matrix.size() == 0) throw std::out_of_range(ERROR("Out of range in stream operator >>: Unable to fill the Dyadic!"));
			for (unsigned int i = 0; i < matrix.rows(); i++) {
				for (unsigned int j = i; j < matrix.cols(); j++) {
					if (input.eof()) throw std::length_error(ERROR("Length error in stream operator >>: Unable to fill the Dyadic!"));
					input >> matrix.at(i, j);
				}
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
			for (unsigned int i = 0; i < mv_nSize; i++) {
				output << "\t" << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, 0);

				for (unsigned int j = 1; j < mv_nSize; j++) {
					output << " " << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, j);
				}
				output << "\n";
			}
			output << std::endl;
			return output.str();
		}

		/** Swap matrices
		  * @param other Matrix to be swaped.
		  */
		void swap(MatrixS& other) noexcept
		{
			std::swap(mv_nSize, other.mv_nSize);
			std::swap(mv_numItens, other.mv_numItens);
			mv_Values.swap(other.mv_Values);
		}

		/** Set all values to zero. Size remains unchanged.
		  */
		void clear()
		{
			memset(&mv_Values[0], 0., mv_Values.size() * sizeof(double));
		}

		/** Resize the matrix.
		  * @param nSize Number of columns / rows of the square matrix
		  */
		void resize(const unsigned int& nSize)
		{
			mv_nSize = nSize;
			mv_numItens = (unsigned int)(0.5 * mv_nSize * mv_nSize + 0.5 * mv_nSize);
			mv_Values.resize(mv_numItens);
			clear();
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& at(unsigned int i, unsigned int j) {
			if (i > mv_nSize || j > mv_nSize) throw std::out_of_range(ERROR("One of components is out of range in MatrixS method .at"));
			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nSize - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nSize - i * 0.5 - 0.5) + j);
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& at(unsigned int i, unsigned int j) const {
			if (i > mv_nSize || j > mv_nSize) throw std::out_of_range(ERROR("One of components is out of range in MatrixS method .at"));
			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nSize - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nSize - i * 0.5 - 0.5) + j);
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& operator()(unsigned int i, unsigned int j) {
			if (i > mv_nSize || j > mv_nSize) throw std::out_of_range(ERROR("One of components is out of range in MatrixS method .at"));
			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nSize - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nSize - i * 0.5 - 0.5) + j);
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& operator()(unsigned int i, unsigned int j) const {
			if (i > mv_nSize || j > mv_nSize) throw std::out_of_range(ERROR("One of components is out of range in MatrixS method .at"));
			unsigned int pos = (i > j) ? (unsigned int)(j * (mv_nSize - j * 0.5 - 0.5) + i) : (unsigned int)(i * (mv_nSize - i * 0.5 - 0.5) + j);
			return mv_Values.at(pos);
		}

		/** @return the row size.
		  */
		unsigned int rows() const noexcept { return mv_nSize; }

		/** @return the column size.
		  */
		unsigned int cols() const noexcept { return mv_nSize; }

		/** @return the number of itens.
		  */
		unsigned int size() const noexcept { return mv_numItens; }

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
		  * @param other Matrix to be copied.
		  */
		MatrixS& operator=(const MatrixS& other)
		{
			if (this == &other) {
				return *this;
			}

			if (mv_nSize == other.mv_nSize)
			{
				mv_Values = other.mv_Values;
			}
			else
			{
				mv_nSize = other.mv_nSize;
				mv_numItens = other.mv_numItens;
				mv_Values.resize(mv_numItens);
				mv_Values = other.mv_Values;
			}
			return *this;
		}

		/** Overloads operator = for move assignment operations -> T = &O
		  * @param other Matrix to be moved.
		  */
		MatrixS& operator=(MatrixS&& other) noexcept
		{
			if (this != &other) {
				mv_nSize = other.mv_nSize;
				mv_numItens = other.mv_numItens;
				mv_Values = std::move(other.mv_Values);

				other.mv_nSize = 0;
				other.mv_numItens = 0;
			}
			return *this;
		}

		/** Overloads operator == for comparison operations
		  * @param other Matrix to be compared.
		  */
		bool operator==(const MatrixS& other)
		{
			bool mi_check = true;
			if (mv_nSize == other.mv_nSize && mv_numItens == other.mv_numItens) {
				for (int i = 0; i < mv_numItens; ++i) {
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
		  * @param other Matrix to be compared.
		  */
		bool operator!=(const MatrixS& other)
		{
			if (*this == other) return false;
			return true;
		}

		/** Overloads operator += for cumulative addition -> T += O -> T = T + O
		  * @param other Matrix to be added.
		  */
		MatrixS& operator+=(const MatrixS& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) += other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator -= for cumulative substraction -> T -= O -> T = T - O
		  * @param other Matrix to be substracted.
		  */
		MatrixS& operator-=(const MatrixS& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator *= for cumulative multiplication -> T *= O -> T_ij = T_ik * O_kj
		  * @param other Matrix to be multiplied with.
		  */
		MatrixS& operator*=(const MatrixS& other)
		{
			check_order(other);
			MatrixS result(mv_nSize);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				for (unsigned int j = i; j < mv_nSize; j++) {
					for (unsigned int k = 0; k < mv_nSize; k++) {
						result.at(i, j) += at(i, k) * other.at(k, j);
					}
				}
			}
			swap(result);
			return *this;
		}

		/** Overloads operator + for addition -> T = T + O
		  * @param other Matrix to be added.
		  */
		MatrixS operator+(const MatrixS& other) const
		{
			check_order(other);
			MatrixS result(mv_nSize);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) = mv_Values.at(i) + other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for substraction -> T = T - O
		  * @param other Matrix to be substracted.
		  */
		MatrixS operator-(const MatrixS& other) const
		{
			check_order(other);
			MatrixS result(mv_nSize);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) = mv_Values.at(i) - other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
		  * @param other Matrix to be multiplied with.
		  */
		MatrixS operator*(const MatrixS& other) const
		{
			check_order(other);
			MatrixS result(mv_nSize);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				for (unsigned int j = i; j < mv_nSize; j++) {
					for (unsigned int k = 0; k < mv_nSize; k++) {
						result.at(i, j) += at(i, k) * other.at(k, j);
					}
				}
			}
			return result;
		}

		/** Overloads operator * for Dot product -> V_i = T_ij * vec_j
		  * @param vec Vector to be multiplied with.
		  * @return a vector with the inner product.
		  */
		std::vector<double> operator*(const std::vector<double>& vec)
		{
			if (vec.size() != mv_nSize) throw std::range_error(ERROR("Invalid request in MatrixS dot product: Order of tensor and vector differ!"));

			std::vector<double> L(mv_nSize);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				for (unsigned int j = 0; j < mv_nSize; j++) {
					L.at(i) += at(i, j) * vec.at(j);
				}
			}
			return L;
		}

		/** Overloads operator + to sum a scalar -> T_ij = T_ij + alfa
		  * @param alfa Scalar to be added.
		  * @return a matrix with the result.
		  */
		MatrixS operator+(const double& alfa) const
		{
			MatrixS result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) += alfa;
			}
			return result;
		}

		/** Overloads operator - to substract a scalar -> T_ij = T_ij - alfa
		  * @param alfa Scalar to be substracted.
		  * @return a matrix with the result.
		  */
		MatrixS operator-(const double& alfa) const
		{
			MatrixS result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) -= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a matrix with the result.
		  */
		MatrixS operator*(const double& alfa) const
		{
			MatrixS result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator / to divide with a non-zero scalar -> T_ij = T_ij / alfa
		  * @param alfa Scalar to be divided with.
		  * @return a matrix with the result.
		  */
		MatrixS operator/(const double& alfa) const
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in MatrixS division operator: Division by zero!"));
			}

			MatrixS result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) /= alfa;
			}
			return result;
		}

		/** Overloads operator += to sum a scalar -> T_ij += alfa
		  * @param alfa Scalar to be added.
		  * @return a matrix with the result.
		  */
		MatrixS& operator+=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) += alfa;
			}
			return *this;
		}

		/** Overloads operator -= to substract a scalar -> T_ij -= alfa
		  * @param alfa Scalar to be substracted.
		  * @return a matrix with the result.
		  */
		MatrixS& operator-=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) -= alfa;
			}
			return *this;
		}

		/** Overloads operator *= to multiply with a scalar -> T_ij *= alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a matrix with the result.
		  */
		MatrixS& operator*=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) *= alfa;
			}
			return *this;
		}

		/** Overloads operator /= to divide with a non-zero scalar -> T_ij /= alfa
		  * @param alfa Scalar to be divided with.
		  * @return a matrix with the result.
		  */
		MatrixS& operator/=(const double& alfa)
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in MatrixS division operator: Division by zero!"));
			}

			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) /= alfa;
			}
			return *this;
		}

	private:
		unsigned int mv_nSize;		// Number of columns / rows
		unsigned int mv_numItens;	// Number of stored values
		std::vector<double> mv_Values;

		inline void check_order(const MatrixS& other) const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in MatrixS method .check_order: Empty matrix!"));
			assert(mv_nSize == other.mv_nSize);	// Order of matrices differ
		}
	};
}  // End of M2S2 namespace
