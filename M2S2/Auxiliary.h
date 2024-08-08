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

namespace M2S2 {
    /** @class triplet
      * @brief A triplet is a record for each non-zero entry in the matrix (row number, column number, and value)
      * @details These triplet format is auxiliary only for sparseMatrix. Use with caution.
      */
    class triplet {
    public:
        /** M2S2 triplet implementation. A triplet corresponds to the non-zero values of a sparse matrix along with their row and column index values
          */
        triplet() {
            mv_row = 0;
            mv_col = 0;
            mv_value = 0.;
        }

        /** M2S2 triplet implementation. A triplet corresponds to the non-zero values of a sparse matrix along with their row and column index values
          * @param row Row number
          * @param col Column number
          * @param value Value of the entry
          */
        triplet(const unsigned int row, const unsigned int col, const double value) {
            mv_row = row;
            mv_col = col;
            mv_value = value;
        }

        /** Copy constructor for M2S2 triplet.
          * @param other Triplet to be copied.
          */
        triplet(const triplet& other) {
            mv_row = other.mv_row;
            mv_col = other.mv_col;
            mv_value = other.mv_value;
        }

        /** Move constructor for M2S2 triplet.
          * @param other Triplet to be moved.
          */
        triplet(triplet&& other) noexcept : mv_row(other.mv_row), mv_col(other.mv_col), mv_value(other.mv_value) { }

        /** Destructor.
          */
        ~triplet() = default;

    // Member Variables
    public:
        /** @brief row index */
        int mv_row;

        /** @brief column index */
        int mv_col;

        /** @brief Value of the matrix element. */
        double mv_value;
    };

    /** @class CSR
      * @brief CSR stands for Compressed Sparse Row, a Sparse BLAS format for matrices.
      * @details This CSR class was created only to use in Sparse packages. All operations must be performed in sparseMatrix class.
      */
    class CSR {

    public:
        /** M2S2 Compressed Sparse Row matrices.
          */
        CSR() {
            mv_sym = false;
            mv_rows = 0;
            mv_cols = 0;
        }

        /** Copy constructor for M2S2 CSR.
          * @param other CSR to be copied.
          */
        CSR(const CSR& other) {
            mv_sym = other.mv_sym;
            mv_rows = other.mv_rows;
            mv_cols = other.mv_cols;

            mv_rowIndex.assign(other.mv_rowIndex.begin(), other.mv_rowIndex.end());
            mv_colIndex.assign(other.mv_colIndex.begin(), other.mv_colIndex.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());
        }

        /** Move constructor for M2S2 CSR.
          * @param other CSR to be moved.
          */
        CSR(CSR&& other) noexcept 
            : mv_sym(other.mv_sym), mv_rows(other.mv_rows), mv_cols(other.mv_cols), mv_rowIndex(other.mv_rowIndex), mv_colIndex(other.mv_colIndex), mv_value(other.mv_value) { }

        /** Destructor.
          */
        ~CSR() = default;

        /** Removes all elements from the matrix, leaving the container empty
          */
        void destroy()
        {
            mv_sym = false;
            mv_rows = 0;
            mv_cols = 0;

            std::vector<int>().swap(mv_rowIndex);
            std::vector<int>().swap(mv_colIndex);
            std::vector<double>().swap(mv_value);
        }

        /** Overloads operator << to stream the matrix. */
        friend std::ostream& operator<<(std::ostream& output, const CSR& matrix)
        {
            output << matrix.print();
            return output;
        }

        /** Prepare a string to print (to file or screen)
          * @param precision Number of decimal digits after the decimal point (default is 3)
          * @param width Minimum number of characters to be written (default is 7)
          */
        const std::string print(const int precision = 3, const int width = 7) const
        {
            std::ostringstream output;
            unsigned int vecIndex;
            unsigned int rowStart;
            unsigned int rowEnd;

            output << std::endl;
            if (mv_sym) output << WARN("Printing only the symmetric part - the rest is filled by zeros!");

            for (unsigned int i = 0; i < mv_rows; ++i) {
                rowStart = mv_rowIndex.at(i);
                rowEnd = mv_rowIndex.at(i + 1) - 1;
                vecIndex = rowStart;

                for (unsigned int j = 0; j < mv_cols; ++j) {
                    if (j != mv_colIndex.at(vecIndex)) {
                        output << std::setw(width) << std::setprecision(precision) << mg_zero << " ";
                    }
                    else {
                        output << std::setw(width) << std::setprecision(precision) << mv_value.at(vecIndex) << " ";
                        if (vecIndex != rowEnd) ++vecIndex;
                    }
                }
                output << std::endl;
            }
            output << std::endl;
            return output.str();
        }

        /** Save current CSR sparse matrix as MatrixS.
          * @param other Dense symmetric matrix.
          */
        void saveAsMatrixS(M2S2::MatrixS& other) {
            assert(mv_sym); // Cannot save an asymmetric sparse matrix as a symmetric matrix!

            other.resize(mv_rows);
            other.clear();

            unsigned int vecIndex;
            unsigned int rowStart;
            unsigned int rowEnd;

            for (unsigned int i = 0; i < mv_rows; ++i) {
                rowStart = mv_rowIndex.at(i);
                rowEnd = mv_rowIndex.at(i + 1) - 1;
                vecIndex = rowStart;

                for (unsigned int j = 0; j < mv_cols; ++j) {
                    if (j == mv_colIndex.at(vecIndex)) {
                        other.at(i, j) = mv_value.at(vecIndex);

                        if (vecIndex != rowEnd) {
                            ++vecIndex;
                        }
                    }
                }
            }
        }

        /** Save current CSR sparse matrix as MatrixX.
          * @param other Dense asymmetric matrix.
          */
        void saveAsMatrixX(M2S2::MatrixX& other) {
            other.resize(mv_rows, mv_cols);
            other.clear();

            unsigned int vecIndex;
            unsigned int rowStart;
            unsigned int rowEnd;

            if(mv_sym) {
                for (unsigned int i = 0; i < mv_rows; ++i) {
                    rowStart = mv_rowIndex.at(i);
                    rowEnd = mv_rowIndex.at(i + 1) - 1;
                    vecIndex = rowStart;

                    for (unsigned int j = 0; j < mv_cols; ++j) {
                        if (j == mv_colIndex.at(vecIndex)) {
                            other.at(i, j) = mv_value.at(vecIndex);
                            other.at(j, i) = mv_value.at(vecIndex);

                            if (vecIndex != rowEnd) {
                                ++vecIndex;
                            }
                        }
                    }
                }
            }
            else {
                for (unsigned int i = 0; i < mv_rows; ++i) {
                    rowStart = mv_rowIndex.at(i);
                    rowEnd = mv_rowIndex.at(i + 1) - 1;
                    vecIndex = rowStart;

                    for (unsigned int j = 0; j < mv_cols; ++j) {
                        if (j == mv_colIndex.at(vecIndex)) {
                            other.at(i, j) = mv_value.at(vecIndex);

                            if (vecIndex != rowEnd) {
                                ++vecIndex;
                            }
                        }
                    }
                }
            }
        }

    public:
        /** @brief If not symmetric, sym = false */
        bool mv_sym = false;

        /** @brief Number of rows */
        int mv_rows;

        /** @brief Number of columns */
        int mv_cols;

        /** @brief rowIndex[i] - start of line i; rowIndex[i+1] - end of line i. */
        std::vector<int> mv_rowIndex;

        /** @brief Column index. */
        std::vector<int> mv_colIndex;

        /** @brief Value of the matrix elements. */
        std::vector<double> mv_value;
    };


    /** @class CSC
      * @brief CSC stands for Compressed Sparse Column, a Sparse BLAS format for matrices.
      * @details These CSC class was created only to output data. All operations must be performed in sparseMatrix class.
      */
    class CSC {

    public:
        /** M2S2 Compressed Sparse Row matrices.
          */
        CSC() {
            mv_sym = false;
            mv_rows = 0;
            mv_cols = 0;
        }

        /** Copy constructor for M2S2 CSC.
          * @param other CSC to be copied.
          */
        CSC(const CSC& other) {
            mv_sym = other.mv_sym;
            mv_rows = other.mv_rows;
            mv_cols = other.mv_cols;

            mv_colIndex.assign(other.mv_colIndex.begin(), other.mv_colIndex.end());
            mv_rowIndex.assign(other.mv_rowIndex.begin(), other.mv_rowIndex.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());
        }

        /** Move constructor for M2S2 CSC.
          * @param other CSC to be moved.
          */
        CSC(CSC&& other) noexcept
            : mv_sym(other.mv_sym), mv_rows(other.mv_rows), mv_cols(other.mv_cols), mv_colIndex(other.mv_colIndex), mv_rowIndex(other.mv_rowIndex), mv_value(other.mv_value) { }

        /** Destructor.
          */
        ~CSC() = default;

        /** Removes all elements from the matrix, leaving the container empty
          */
        void destroy()
        {
            mv_sym = false;
            mv_rows = 0;
            mv_cols = 0;

            std::vector<int>().swap(mv_colIndex);
            std::vector<int>().swap(mv_rowIndex);
            std::vector<double>().swap(mv_value);

            mv_colIndex.clear();
            mv_rowIndex.clear();
            mv_value.clear();
        }

        /** Overloads operator << to stream the matrix. */
        friend std::ostream& operator<<(std::ostream& output, const CSC& matrix)
        {
            output << matrix.print();
            return output;
        }

        /** Prepare a string to print (to file or screen)
          * @param precision Number of decimal digits after the decimal point (default is 3)
          * @param width Minimum number of characters to be written (default is 7)
          */
        const std::string print(const int precision = 3, const int width = 7) const
        {
            std::ostringstream output;
            int vecIndex;
            unsigned int colStart;
            unsigned int colEnd;

            output << std::endl;
            if (mv_sym) output << WARN("Printing only the symmetric part - the rest is filled by zeros!");

            for (unsigned int i = 0; i < mv_rows; ++i) {
                for (unsigned int j = 0; j < mv_cols; ++j) {
                    vecIndex = -1;
                    colStart = mv_colIndex.at(j);
                    colEnd = mv_colIndex.at(j + 1);

                    for (unsigned int k = colStart; k < colEnd; ++k) {
                        if (mv_rowIndex.at(k) == i) {
                            vecIndex = k;
                            break;
                        }
                        if (mv_rowIndex.at(k) > i) break;
                    }
                    if (vecIndex != -1) {
                        output << std::setw(width) << std::setprecision(precision) << mv_value.at(vecIndex) << " ";
                    }
                    else {
                        output << std::setw(width) << std::setprecision(precision) << mg_zero << " ";
                    }
                }
                output << std::endl;
            }
            output << std::endl;
            return output.str();
        }

        /** Save current CSC sparse matrix as MatrixS.
          * @param other Dense symmetric matrix.
          */
        void saveAsMatrixS(M2S2::MatrixS& other) {
            assert(mv_sym); // Cannot save an asymmetric sparse matrix as a symmetric matrix!

            other.resize(mv_cols);
            other.clear();

            unsigned int vecIndex;
            unsigned int colStart;
            unsigned int colEnd;

            for (unsigned int i = 0; i < mv_cols; ++i) {
                colStart = mv_colIndex.at(i);
                colEnd   = mv_colIndex.at(i + 1) - 1;
                vecIndex = colStart;

                for (unsigned int j = 0; j < mv_rows; ++j) {
                    if (j == mv_rowIndex.at(vecIndex)) {
                        other.at(j, i) = mv_value.at(vecIndex);

                        if (vecIndex != colEnd) {
                            ++vecIndex;
                        }
                    }
                }
            }
        }

        /** Save current CSC sparse matrix as MatrixX.
          * @param other Dense asymmetric matrix.
          */
        void saveAsMatrixX(M2S2::MatrixX& other) {
            other.resize(mv_rows, mv_cols);
            other.clear();

            unsigned int vecIndex;
            unsigned int colStart;
            unsigned int colEnd;

            if (mv_sym) {
                for (unsigned int i = 0; i < mv_cols; ++i) {
                    colStart = mv_colIndex.at(i);
                    colEnd = mv_colIndex.at(i + 1) - 1;
                    vecIndex = colStart;

                    for (unsigned int j = 0; j < mv_rows; ++j) {
                        if (j == mv_rowIndex.at(vecIndex)) {
                            other.at(j, i) = mv_value.at(vecIndex);
                            other.at(i, j) = mv_value.at(vecIndex);

                            if (vecIndex != colEnd) {
                                ++vecIndex;
                            }
                        }
                    }
                }
            }
            else {
                for (unsigned int i = 0; i < mv_cols; ++i) {
                    colStart = mv_colIndex.at(i);
                    colEnd = mv_colIndex.at(i + 1) - 1;
                    vecIndex = colStart;

                    for (unsigned int j = 0; j < mv_rows; ++j) {
                        if (j == mv_rowIndex.at(vecIndex)) {
                            other.at(j, i) = mv_value.at(vecIndex);
                            other.at(i, j) = mv_value.at(vecIndex);

                            if (vecIndex != colEnd) {
                                ++vecIndex;
                            }
                        }
                    }
                }
            }
        }

    public:
        /** @brief If not symmetric, sym = false */
        bool mv_sym = false;

        /** @brief Number of rows */
        int mv_rows;

        /** @brief Number of columns */
        int mv_cols;

        /** @brief colIndex[i] - start of line i; colIndex[i+1] - end of line i. */
        std::vector<int> mv_colIndex;

        /** @brief Row index. */
        std::vector<int> mv_rowIndex;

        /** @brief Value of the matrix elements. */
        std::vector<double> mv_value;
    };


    /** @class line
      * @brief Auxiliary class for sparseMatrix.
      * @details Holds information about each line of the sparseMatrix.
      *
      * Comment (1)
      * resize variable is used to control the resize, according to the
      * adopted policy
      * First attempt to create a policy:
      * 1st resize -> new size is 2x initial size;
      * 2nd resize -> new size is 3x initial size;
      * 3rd resize -> new size is 4x initial size;
      * 4th resize -> new size is 2x last size;
      * 5th resize -> new size is 2x last size;
      * 6th resize -> stop program and show message!
      *
      * Comment (2)
      * Each time new values are pushed in the line, assembled is set to 0, i.e. false.
      * Each time the line is sorted and equal terms are summed, assembled is set to 1.
      * In other words, this variables is used to verify if it is necessary to assemble
      * the line. It is not necessary to assemble it to push new terms, but to use the
      * sparse matrix, all lines must be assembled.
      */
    class line {
    public:
        /** M2S2 row / column constructor.
          */
        line() = default;

        /** M2S2 row / column constructor.
          * @param index vector of indexes of the values in the line.
          * @param value vector of values in the line.
          */
        line(const std::vector<int>& index, const std::vector<double>& value) {
            assert(index.size() == value.size()); // Size of input vectors are not of the same size

            mv_index.assign(index.begin(), index.end());
            mv_value.assign(value.begin(), value.end());
        }

        /** M2S2 row / column constructor.
          * @param size Size to be allocated.
          */
        line(const unsigned int size) {
            reserve(size);
        }

        /** Copy constructor for M2S2 line.
          * @param other Line to be copied.
          */
        line(const line& other) {
            mv_index.reserve(other.mv_index.capacity());
            mv_value.reserve(other.mv_value.capacity());

            mv_index.assign(other.mv_index.begin(), other.mv_index.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());

            mv_assembled = other.mv_assembled;
        }

        /** Move constructor for M2S2 line.
          * @param other Line to be moved.
          */
        line(line&& other) noexcept
        {
            mv_index = std::move(other.mv_index);
            mv_value = std::move(other.mv_value);

            mv_assembled = other.mv_assembled;
        }

        /** Destructor.
          */
        ~line() = default;

        /** Move assignment operator for M2S2 line.
          * @param other Line to be moved.
          */
        line& operator = (line&& other) noexcept
        {
            mv_index = std::move(other.mv_index);
            mv_value = std::move(other.mv_value);

            mv_assembled = other.mv_assembled;
            return *this;
        }

        /** Copy assignment operator for M2S2 line.
          * @param other Line to be copied.
          */
        line& operator = (const line& other)
        {
            mv_index.reserve(other.mv_index.capacity());
            mv_value.reserve(other.mv_value.capacity());

            mv_index.assign(other.mv_index.begin(), other.mv_index.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());

            mv_assembled = other.mv_assembled;

            return *this;
        }

    public:
        /** @brief index of the values in the line */
        std::vector<int> mv_index;

        /** @brief value of the line item */
        std::vector<double> mv_value;

        /** @brief number of times line was resized */
        int mv_resizeCount = 0;

        /** @brief false -> not assembled; true -> assembled */
        bool mv_assembled = false;

    public:
        /** Reserve capacity to the line. Notice that size is not checked (if it is smaller than current).
          * @param size New size to be allocated.
          */
        void reserve(const unsigned int& size)
        {
            mv_index.reserve(size);
            mv_value.reserve(size);

            mv_resizeCount = 0;
            mv_assembled = false;
        }

        /** Removes all elements from the line, leaving the line empty
          * Notice that std::vector::clear does not change capacity, but here it does.
          */
        void destroy()
        {
            if (mv_index.size()) {
                std::vector<int>().swap(mv_index);
                std::vector<double>().swap(mv_value);

                mv_resizeCount = 0;
                mv_assembled = false;
            }
        }

        /** Deletes all elements in the line, but does not change its capacity. Like std::vector::clear()
          */
        void clear()
        {
            mv_index.clear();
            mv_value.clear();
            mv_resizeCount = 0;
            mv_assembled = false;
        }

        /** Overloads operator << to stream the line. */
        friend std::ostream& operator<<(std::ostream& output, const line& line)
        {
            output << line.print();
            return output;
        }

        /** Overloads operator >> to stream the line. */
        friend std::istream& operator>>(std::istream& input, line& line)
        {
            for (int i = 0; i < line.mv_index.size(); ++i) {
                input >> line.mv_index.at(i);
            }
            for (int i = 0; i < line.mv_index.size(); ++i) {
                input >> line.mv_value.at(i);
            }
            return input;
        }

        /** Prepare a string to print (to file or screen) line
          * @param precision Number of decimal digits after the decimal point (default is 4)
          * @param width Minimum number of characters to be written (default is 8)
          */
        const std::string print(const int precision = 4, const int width = 8) const
        {
            std::ostringstream output;
            output << "Line size: " << mv_index.size() << "\t" << std::fixed << std::setprecision(precision) << std::setw(width);
            for (int i = 0; i < mv_index.size(); ++i) {
                output << mv_index.at(i) << ":" << mv_value.at(i) << " ";
            }
            std::cout << std::endl;
            return output.str();
        }

        /** If there are terms with repeated index, these are summed
          */
        bool addEqualTerms() {
            if (mv_index.size() && !mv_assembled) {
                sort();
                /* Adding equal terms will only work if line is alredy sorted */
                int pos = 0;
                for (int i = 1; i < mv_index.size(); i++) {
                    if (mv_index.at(i) == mv_index.at(pos))
                        mv_value.at(pos) += mv_value.at(i);
                    else {
                        pos++;
                        mv_index.at(pos) = mv_index.at(i);
                        mv_value.at(pos) = mv_value.at(i);
                    }
                }
                mv_index.resize(++pos); /* new size, but keeping memory (reserve) */
                mv_value.resize(pos); /* new size, but keeping memory (reserve) */
                mv_assembled = true;
            }
            return mv_assembled;
        }

        /** @return the position in line of index. If index is not found return -1.
          * @param index Index to be found.
          */
        unsigned int search(const unsigned int index)
        {
            // Trying to enhance performance by using a bissection search algorithm
            if (mv_index.size()) {
                if (!mv_assembled) addEqualTerms();
                int ib, im, ie;
                ib = 0;
                ie = mv_index.size() - 1;
                if ((mv_index.at(ib) > index) || (mv_index.at(ie) < index))
                    return -1;
                if (mv_index.at(ib) == index)
                    return ib;
                if (mv_index.at(ie) == index)
                    return ie;

                /* int iter = 0;*/
                while ((ie - ib) > 1) {
                    if (mv_index.at(ib) == index)
                        return ib;
                    if (mv_index.at(ie) == index)
                        return ie;
                    im = ib + (ie - ib) / 2;
                    if (mv_index.at(im) == index)
                        return im;
                    if (index <= mv_index.at(im)) {
                        ie = im;
                    }
                    else {
                        ib = im;
                    }
                }
                return -1;
            }
            else {
                return -1;
            }
        }

        /** @return the position in line of index. If index is not found return -1.
          * @param index Index to be found.
          */
        unsigned int search(const unsigned int index) const
        {
            assert(mv_assembled); // Const version can only be used with already sorted line.

            // Trying to enhance performance by using a bissection search algorithm
            if (mv_index.size()) {
                int ib, im, ie;
                ib = 0;
                ie = mv_index.size() - 1;
                if ((mv_index.at(ib) > index) || (mv_index.at(ie) < index))
                    return -1;
                if (mv_index.at(ib) == index)
                    return ib;
                if (mv_index.at(ie) == index)
                    return ie;

                /* int iter = 0;*/
                while ((ie - ib) > 1) {
                    if (mv_index.at(ib) == index)
                        return ib;
                    if (mv_index.at(ie) == index)
                        return ie;
                    im = ib + (ie - ib) / 2;
                    if (mv_index.at(im) == index)
                        return im;
                    if (index <= mv_index.at(im)) {
                        ie = im;
                    }
                    else {
                        ib = im;
                    }
                }
                return -1;
            }
            else {
                return -1;
            }
        }

        /* Increase initial capacity up to five times */
        void resize() {
            /* To be used for automatic resize */
            int initSize, lastSize, newSize;
            switch (mv_resizeCount) {
            case 0:
                /* First resize */
                initSize = mv_index.capacity();
                newSize = 2 * initSize;
                break;
            case 1:
                /* Second resize */
                initSize = mv_index.capacity() / 2;
                newSize = 3 * initSize;
                break;
            case 2:
                /* Third resize */
                initSize = mv_index.capacity() / 3;
                newSize = 4 * initSize;
                break;
            case 3:
                /* Fourth resize */
                lastSize = mv_index.capacity();
                newSize = 2 * lastSize;
                break;
            case 4:
                /* Fifth resize */
                lastSize = mv_index.capacity();
                newSize = 2 * lastSize;
                break;
            default:
                throw std::invalid_argument("\n\nLine resize failed!\nVerify initial size!");
            }
            mv_resizeCount++;
            reserve(newSize);
        }

    private:

        /* Intended just to test search */
        const unsigned int naiveSearch(const unsigned int index)
        {
            if (mv_index.size()) {
                if (!mv_assembled) addEqualTerms();
                if ((mv_index.at(0) > index) || (mv_index.at(mv_index.size() - 1) < index))
                    return -1;
                for (int i = 0; i < mv_index.size(); ++i) {
                    if (index == mv_index.at(i))
                        return i;
                    else if (index < mv_index.at(i))
                        return -1;
                }
            }
            else {
                return -1;
            }
            return -1;
        }

        /* Sort line by indexes */
        void quicksort(int* index, double* value, const int& begin, const int& end)
        {
            int i, j, pivo, aux;
            double daux;
            i = begin;
            j = end - 1;
            pivo = index[(begin + end) / 2];
            while (i <= j)
            {
                while (index[i] < pivo && i < end)
                {
                    i++;
                }
                while (index[j] > pivo && j > begin)
                {
                    j--;
                }
                if (i <= j)
                {
                    aux = index[i];
                    index[i] = index[j];
                    index[j] = aux;

                    daux = value[i];
                    value[i] = value[j];
                    value[j] = daux;

                    i++;
                    j--;
                }
            }
            if (j > begin)
                quicksort(index, value, begin, j + 1);
            if (i < end)
                quicksort(index, value, i, end);
        }

        /* Sort line by indexes */
        void quicksort(std::vector<int>& index, std::vector<double>& value, const int& begin, const int& end)
        {
            int i, j, pivo, aux;
            double daux;
            i = begin;
            j = end - 1;
            pivo = index[(begin + end) / 2];
            while (i <= j)
            {
                while (index[i] < pivo && i < end)
                {
                    i++;
                }
                while (index[j] > pivo && j > begin)
                {
                    j--;
                }
                if (i <= j)
                {
                    aux = index[i];
                    index[i] = index[j];
                    index[j] = aux;

                    daux = value[i];
                    value[i] = value[j];
                    value[j] = daux;

                    i++;
                    j--;
                }
            }
            if (j > begin)
                quicksort(index, value, begin, j + 1);
            if (i < end)
                quicksort(index, value, i, end);
        }

        /* Sort line by indexes */
        void sort() {
            quicksort(mv_index, mv_value, 0, mv_index.size());
        }
    };
}  // End of M2S2 namespace
