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
// MatrixX class
//
// ================================================================================================
namespace M2S2 {
	/** @class MatrixX
	 * @brief 2nd rank matrices of any size, even not square.
	 * @details 2nd rank matrices of any size, even not square (saved as row major).
	 */
	class MatrixX {

	public:
		/** 2nd rank matrices of any size. Values are saved in a single vector, as row major.
		  */
		MatrixX()
		{
			mv_nCol = 0;
			mv_nRow = 0;
			mv_nSize = 0;
			mv_Values.clear();
		}

		/** 2nd rank matrices of any size. Values are saved in a single vector, as row major.
		  * @param nRow Number of rows
		  * @param nCol Number of columns
		  */
		MatrixX(unsigned int nRow, unsigned int nCol) : mv_nRow(nRow), mv_nCol(nCol)
		{
			mv_nSize = mv_nRow * mv_nCol;
			mv_Values.resize(mv_nSize);
		}

		/** 2nd rank matrices of any size. Values are saved in a single vector, as row major.
		  * @param nRow Number of rows
		  * @param nCol Number of columns
		  * @param value Value to initiate the entire matrix.
		  */
		MatrixX(unsigned int nRow, unsigned int nCol, const double& value) : mv_nRow(nRow), mv_nCol(nCol)
		{
			mv_nSize = mv_nRow * mv_nCol;
			mv_Values.resize(mv_nSize, value);
		}

		/** 2nd rank matrices of any size. Values are saved in a single vector, as row major.
		  * @param nRow Number of rows
		  * @param nCol Number of columns
		  * @param value Vector with enough values to initiate the matrix.
		  */
		MatrixX(unsigned int nRow, unsigned int nCol, const std::vector<double>& value) : mv_nRow(nRow), mv_nCol(nCol)
		{
			mv_nSize = mv_nRow * mv_nCol;

			if (value.size() != mv_nSize) throw std::invalid_argument(ERROR("Invalid argument on MatrixX constructor: Size of input vector does not correspond to required matrix!"));
			mv_Values = { value.begin(), value.begin() + mv_nSize };
		}

		/** Copy constructor for 2nd rank matrices of any size.
		  * @param other Matrix to be copied.
		  */
		MatrixX(const MatrixX& other) {
			mv_nRow = other.mv_nRow;
			mv_nCol = other.mv_nCol;
			mv_nSize = other.mv_nSize;
			mv_Values = other.mv_Values;
		}

		/** Move constructor for symmetric square matrix of any order.
		  * @param other Matrix to be moved.
		  */
		MatrixX(MatrixX&& other) noexcept
			: mv_nRow(other.mv_nRow), mv_nCol(other.mv_nCol), mv_nSize(other.mv_nSize), mv_Values(std::move(other.mv_Values)) { }

		/** Destructor.
		  */
		~MatrixX() { }

		/** Overloads operator << to stream the matrix. */
		friend std::ostream& operator<<(std::ostream& output, const MatrixX& matrix)
		{
			output << matrix.print();
			return output;
		}

		/** Overloads operator >> to stream the matrix. */
		friend std::istream& operator>>(std::istream& input, MatrixX& matrix)
		{
			if (matrix.size() == 0) throw std::out_of_range(ERROR("Out of range in stream operator >>: Unable to fill the Dyadic!"));
			for (unsigned int i = 0; i < matrix.rows(); i++) {
				for (unsigned int j = 0; j < matrix.cols(); j++) {
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
			for (unsigned int i = 0; i < mv_nRow; i++) {
				output << "\t" << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, 0);

				for (unsigned int j = 1; j < mv_nCol; j++) {
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
		void swap(MatrixX& other) noexcept
		{
			std::swap(mv_nRow, other.mv_nRow);
			std::swap(mv_nCol, other.mv_nCol);
			std::swap(mv_nSize, other.mv_nSize);
			mv_Values.swap(other.mv_Values);
		}

		/** Set all values to zero. Size remains unchanged.
		  */
		void clear()
		{
			memset(&mv_Values[0], 0., mv_Values.size() * sizeof(double));
		}

		/** Resize the matrix.
		  * @param nRow Number of rows
		  * @param nCol Number of columns
		  */
		void resize(const unsigned int& nRow, const unsigned int& nCol)
		{
			mv_nRow = nRow;
			mv_nCol = nCol;
			mv_nSize = mv_nRow * mv_nCol;
			mv_Values.resize(mv_nSize);
			clear();
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& at(unsigned int i, unsigned int j) { 
			if (i > mv_nRow || j > mv_nCol) throw std::out_of_range(ERROR("One of components is out of range in MatrixX method .at"));
			return mv_Values.at(i * mv_nCol + j);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& at(unsigned int i, unsigned int j) const {
			if (i > mv_nRow || j > mv_nCol) throw std::out_of_range(ERROR("One of components is out of range in MatrixX method .at"));
			return mv_Values.at(i * mv_nCol + j);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& operator()(unsigned int i, unsigned int j) {
			if (i > mv_nRow || j > mv_nCol) throw std::out_of_range(ERROR("One of components is out of range in MatrixX method .at"));
			return mv_Values.at(i * mv_nCol + j);
		}

		/** Access specified element - returns Matrix_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& operator()(unsigned int i, unsigned int j) const {
			if (i > mv_nRow || j > mv_nCol) throw std::out_of_range(ERROR("One of components is out of range in MatrixX method .at"));
			return mv_Values.at(i * mv_nCol + j);
		}

		/** @return the row size.
		  */
		unsigned int rows() const noexcept { return mv_nRow; }

		/** @return the column size.
		  */
		unsigned int cols() const noexcept { return mv_nCol; }

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
		  * @param other Matrix to be copied.
		  */
		MatrixX& operator=(const MatrixX& other)
		{
			if (this == &other) {
				return *this;
			}

			if (mv_nCol == other.mv_nCol && mv_nRow == other.mv_nRow)
			{
				mv_Values = other.mv_Values;
			}
			else
			{
				mv_nRow = other.mv_nRow;
				mv_nCol = other.mv_nCol;
				mv_nSize = other.mv_nSize;
				mv_Values.resize(mv_nSize);
				mv_Values = other.mv_Values;
			}
			return *this;
		}

		/** Overloads operator = for move assignment operations -> T = &O
		  * @param other Matrix to be moved.
		  */
		MatrixX& operator=(MatrixX&& other) noexcept
		{
			if (this != &other) {
				mv_nRow = other.mv_nRow;
				mv_nCol = other.mv_nCol;
				mv_nSize = other.mv_nSize;
				mv_Values = std::move(other.mv_Values);

				other.mv_nRow = 0;
				other.mv_nCol = 0;
				other.mv_nSize = 0;
			}
			return *this;
		}

		/** Overloads operator == for comparison operations
		  * @param other Matrix to be compared.
		  */
		bool operator==(const MatrixX& other)
		{
			bool mi_check = true;
			if (mv_nSize == other.mv_nSize && mv_nRow == other.mv_nRow && mv_nCol == other.mv_nCol) {
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
		  * @param other Matrix to be compared.
		  */
		bool operator!=(const MatrixX& other)
		{
			if (*this == other) return false;
			return true;
		}

		/** Overloads operator += for cumulative addition -> T += O -> T = T + O
		  * @param other Matrix to be added.
		  */
		MatrixX& operator+=(const MatrixX& other)
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
		MatrixX& operator-=(const MatrixX& other)
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
		MatrixX& operator*=(const MatrixX& other)
		{
			if (mv_nCol != other.mv_nRow) throw std::range_error(ERROR("Invalid request in MatrixX dot product: Internal order of matrices differ!"));

			MatrixX result(mv_nRow, other.mv_nCol);
			for (unsigned int i = 0; i < mv_nRow; i++) {
				for (unsigned int j = 0; j < other.mv_nCol; j++) {
					for (unsigned int k = 0; k < mv_nCol; k++) {
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
		MatrixX operator+(const MatrixX& other) const
		{
			check_order(other);
			MatrixX result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) += other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for substraction -> T = T - O
		  * @param other Matrix to be substracted.
		  */
		MatrixX operator-(const MatrixX& other) const
		{
			check_order(other);
			MatrixX result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
		  * @param other Matrix to be multiplied with.
		  */
		MatrixX operator*(const MatrixX& other) const
		{
			if (mv_nCol != other.mv_nRow) throw std::range_error(ERROR("Invalid request in MatrixX dot product: Internal order of matrices differ!"));

			MatrixX result(mv_nRow, other.mv_nCol);
			for (unsigned int i = 0; i < mv_nRow; i++) {
				for (unsigned int j = 0; j < other.mv_nCol; j++) {
					for (unsigned int k = 0; k < mv_nCol; k++) {
						result.at(i, j) += at(i, k) * other.at(k, j);
					}
				}
			}
			return result;
		}

		/** Overloads operator * for Dot product -> V_i = T_ij * R_j
		  * @param vec Vector to be multiplied with.
		  * @return a vector with the inner product.
		  */
		std::vector<double> operator*(const std::vector<double>& vec)
		{
			if (vec.size() != mv_nCol) throw std::range_error(ERROR("Invalid request in MatrixX dot product: Order of tensor and vector differ!"));

			std::vector<double> L(mv_nRow);
			for (unsigned int i = 0; i < mv_nRow; i++) {
				for (unsigned int j = 0; j < mv_nCol; j++) {
					L.at(i) += at(i, j) * vec.at(j);
				}
			}
			return L;
		}

		/** Overloads operator + to sum a scalar -> T_ij = T_ij + alfa
		  * @param alfa Scalar to be added.
		  * @return a matrix with the result.
		  */
		MatrixX operator+(const double& alfa) const
		{
			MatrixX result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) += alfa;
			}
			return result;
		}

		/** Overloads operator - to substract a scalar -> T_ij = T_ij - alfa
		  * @param alfa Scalar to be substracted.
		  * @return a matrix with the result.
		  */
		MatrixX operator-(const double& alfa) const
		{
			MatrixX result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) -= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a matrix with the result.
		  */
		MatrixX operator*(const double& alfa) const
		{
			MatrixX result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator / to divide with a non-zero scalar -> T_ij = T_ij / alfa
		  * @param alfa Scalar to be divided with.
		  * @return a matrix with the result.
		  */
		MatrixX operator/(const double& alfa) const
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in MatrixX division operator: Division by zero!"));
			}

			MatrixX result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) /= alfa;
			}
			return result;
		}

		/** Overloads operator += to sum a scalar -> T_ij += alfa
		  * @param alfa Scalar to be added.
		  * @return a matrix with the result.
		  */
		MatrixX& operator+=(const double& alfa)
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
		MatrixX& operator-=(const double& alfa)
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
		MatrixX& operator*=(const double& alfa)
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
		MatrixX& operator/=(const double& alfa)
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument(ERROR("Invalid argument in MatrixX division operator: Division by zero!"));
			}

			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) /= alfa;
			}
			return *this;
		}

		/** @return the transpose of the matrix.
		  */
		MatrixX transpose() const
		{
			MatrixX result(mv_nCol, mv_nRow);
			for (unsigned int i = 0; i < mv_nCol; ++i) {
				for (unsigned int j = 0; j < mv_nRow; ++j) {
					result.at(i, j) = at(j, i);
				}
			}
			return result;
		}

		/** @return the matrix as symmetric (even if it is not). Must at least be square.
		  */
		MatrixS SaveAsMatrixS() const
		{
			assert(mv_nRow == mv_nCol);		// Number of rows and columns must at least be the same

			MatrixS result(mv_nCol);
			for (unsigned int i = 0; i < mv_nRow; ++i) {
				for (unsigned int j = i; j < mv_nCol; ++j) {
					result.at(i, j) = at(i, j);
				}
			}
			return result;
		}

	private:
		unsigned int mv_nRow;		// Number of rows
		unsigned int mv_nCol;		// Number of columns
		unsigned int mv_nSize;		// Number of stored values
		std::vector<double> mv_Values;

		inline void check_order(const MatrixX& other) const
		{
			if (mv_Values.empty()) throw std::runtime_error(ERROR("Invalid request in MatrixX method .check_order: Empty matrix!"));
			assert(mv_nCol == other.mv_nCol && mv_nRow == other.mv_nRow);	// Order of matrices differ
		}
	};
}  // End of M2S2 namespace
