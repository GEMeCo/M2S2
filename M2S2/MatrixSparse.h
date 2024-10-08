// ================================================================================================
// 
// This file is part of M2S2 - Matrices for Mechanices of Solids and Structures
//
// Copyright(C) 2024 
//		Dorival Piedade Neto &
//		Rodrigo Ribeiro Paccola &
//		Rog�rio Carrazedo
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// ================================================================================================
#ifndef _SPARSESET_
#define _SPARSESET_

// M2S2 libraries
#include "Common.h"
#include "Auxiliary.h"
#include "MatrixS.h"
#include "MatrixX.h"

// ================================================================================================
//
// sparseMatrix class
//
// ================================================================================================
namespace M2S2 {
    /** @class sparseMatrix
      * @brief Sparse matrix is a matrix in which most of the elements are zero. Thus, only non-zero values are stored.
      */
    class sparseMatrix {

    // Member Variables
    public:
        /** @brief If not symmetric, sym = false */
        bool mv_sym = false;

        /** @brief May be 'R' or 'C' (Row major or Column major) */
        char mv_type = 'R';

        /** @brief Lines of the sparse matrix. */
        std::vector<M2S2::line> mv_line;

    public:
        /** M2S2 sparse matrix implementation.
          */
        sparseMatrix() = default;

        /** M2S2 sparse matrix implementation.
          * @param nLines Number of lines of the sparse matrix.
          * @param sym If not symmetric, sym = false.
          * @param type Row major or Column major ('R' or 'C'). Default is row major.
          */
        sparseMatrix(const unsigned int nLines, const bool sym = false, const char type = 'R') {
            mv_type = toupper(type);
            assert(mv_type == 'R' || mv_type == 'C'); // Wrong type

            mv_line.resize(nLines);
            mv_sym = sym;
        }

        /** Move constructor for M2S2 sparse matrix.
          * @param other sparseMatrix to be moved.
          */
        sparseMatrix(sparseMatrix&& other) noexcept :
            mv_sym(other.mv_sym), mv_type(other.mv_type), mv_line(std::move(other.mv_line)) {
            other.mv_sym = false;
            other.mv_type = 'R';
        }

        /** Destructor.
          */
        ~sparseMatrix() = default;

    private:
        // Copy constructor is deleted for obvious reason (2 sparse matrices in memory?).
        // Use copy function instead
        sparseMatrix(const sparseMatrix& other) = delete;

    public:
        /** Overloads operator << to stream the matrix. */
        friend std::ostream& operator<<(std::ostream& output, const sparseMatrix& matrix)
        {
            output << matrix.print();
            return output;
        }

        /** Print the non-zero values of the sparse matrix
          * @param precision Number of decimal digits after the decimal point (default is 4)
          * @param width Minimum number of characters to be written (default is 8)
          */
        const std::string print(const unsigned int precision = 4, const unsigned int width = 8) const
        {
            std::ostringstream output;
            output << std::endl << std::fixed;

            // Print non-zero values
            for (unsigned int i = 0; i < mv_line.size(); i++) {
                if (mv_line.at(i).mv_assembled) output << "Assembled Line: " << i << " ";
                else output << "Non assembled Line: " << i << " ";

                for (unsigned int j = 0; j < mv_line.at(i).mv_index.size(); j++) {
                    output << mv_line.at(i).mv_index.at(j) << " : "
                        << std::setw(width) << std::setprecision(precision) << mv_line.at(i).mv_value.at(j) << " ";
                }
                output << std::endl;
            }
            output << std::endl;
            return output.str();
        }

        /** Make a copy of other sparse matrix
          * @param other sparse matrix to be copied.
          * @param swap Swap between row major and column major.
          */
        void copy(M2S2::sparseMatrix& other, const bool swap)
        {
            destroy();
            resize(other.mv_line.size());
            mv_sym = other.mv_sym;

            if (swap) {
                if (other.mv_type == 'C') mv_type = 'R';
                else mv_type = 'C';

                if (mv_sym) {
                    for (int i = 0; i < other.mv_line.size(); i++) {
                        mv_line.assign(other.mv_line.begin(), other.mv_line.end());
                    }
                }
                else {
                    // meanSize: used to set the initial line sizes
                    int meanSize = 0;
                    for (int i = 0; i < other.mv_line.size(); i++) {
                        meanSize += other.mv_line.at(i).mv_index.size();
                    }
                    meanSize /= other.mv_line.size();
                    setAllLineSize(meanSize);

                    // Now, let's copy all information, one by one
                    int currentIndex;
                    double currentValue;
                    for (int i = 0; i < other.mv_line.size(); i++) {
                        for (int j = 0; j < other.mv_line.at(i).mv_index.size(); j++) {
                            currentIndex = other.mv_line.at(i).mv_index.at(j);
                            currentValue = other.mv_line.at(i).mv_value.at(j);

                            mv_line.at(currentIndex).mv_index.push_back(i);
                            mv_line.at(currentIndex).mv_value.push_back(currentValue);
                        }
                    }
                    addEqualTerms();
                }
            }
            else {
                mv_type = other.mv_type;

                for (int i = 0; i < other.mv_line.size(); i++) {
                    mv_line.assign(other.mv_line.begin(), other.mv_line.end());
                }
            }
        }

        /** Save current sparse matrix in CSR format. Notice that this will destroy current CSR matrix.
          * @param other sparse matrix in CSR format.
          */
        void saveAsCSR(M2S2::CSR& other)
        {
            // we must first add equal terms, reducing the size of the matrix
            addEqualTerms();

            // Clean up current CSR matrix
            other.destroy();

            if (mv_type == 'C' && !mv_sym) { // Not symmetric and column major 
                // Swap and then save in CSR format
                M2S2::sparseMatrix swaped;
                swaped.copy(*this, true);
                swaped.saveAsCSR(other);
            }
            else {
                // Compute number of itens in CSR
                int nnz = 0;
                int col = 0;
                other.mv_sym = mv_sym;
                other.mv_rows = mv_line.size();
                other.mv_rowIndex.reserve(mv_line.size() + 1);
                // nnz -> number of the first component in line
                for (int i = 0; i < mv_line.size(); ++i) {
                    other.mv_rowIndex.push_back(nnz);
                    nnz += mv_line.at(i).mv_index.size();
                }
                other.mv_rowIndex.push_back(nnz);
                other.mv_colIndex.reserve(nnz);
                other.mv_value.reserve(nnz);
                // copy terms
                for (int i = 0; i < mv_line.size(); ++i) {
                    for (int j = 0; j < mv_line.at(i).mv_index.size(); ++j) {
                        other.mv_colIndex.push_back(mv_line.at(i).mv_index.at(j));
                        other.mv_value.push_back(mv_line.at(i).mv_value.at(j));

                        col = (col < other.mv_colIndex.back()) ? col : other.mv_colIndex.back();
                    }
                }
                // If not sym, the number of columns is the greatest value from colIndex
                other.mv_cols = mv_sym ? mv_line.size() : col + 1;
            }
        }

        /** Save current sparse matrix in CSC format. Notice that this will destroy current CSC matrix.
          * @param other sparse matrix in CSC format.
          */
        void saveAsCSC(M2S2::CSC& other)
        {
            // we must first add equal terms, reducing the size of the matrix
            addEqualTerms();

            // Clean up current CSC matrix
            other.destroy();

            if (mv_type == 'R' && !mv_sym) { // Not symmetric and column major 
                // Swap and then save in CSC format
                M2S2::sparseMatrix swaped;
                swaped.copy(*this, true);
                swaped.saveAsCSC(other);
            }
            else {
                // Compute number of itens in CSR
                int nnz = 0;
                other.mv_sym = mv_sym;
                other.mv_rows = mv_line.size();
                other.mv_cols = mv_line.size();
                other.mv_colIndex.reserve(mv_line.size() + 1);
                // nnz -> number of the first component in line
                for (int i = 0; i < mv_line.size(); ++i) {
                    other.mv_colIndex.push_back(nnz);
                    nnz += mv_line.at(i).mv_index.size();
                }
                other.mv_colIndex.push_back(nnz);
                other.mv_rowIndex.reserve(nnz);
                other.mv_value.reserve(nnz);
                // copy terms
                for (int i = 0; i < mv_line.size(); ++i) {
                    for (int j = 0; j < mv_line.at(i).mv_index.size(); ++j) {
                        other.mv_rowIndex.push_back(mv_line.at(i).mv_index.at(j));
                        other.mv_value.push_back(mv_line.at(i).mv_value.at(j));
                    }
                }
            }
        }

        /** Resize the matrix. Notice that size is not checked (if it is smaller than current).
          * @param size Size to be allocated and initiated (usually, the number of DOF).
          */
        void resize(const unsigned int size)
        {
            if (mv_line.size()) destroy();
            mv_line.resize(size);
        }

        /** Removes all elements from the matrix, leaving the container empty
          */
        void destroy()
        {
            std::vector<M2S2::line>().swap(mv_line);
            mv_sym = false;
            mv_type = 'R';
        }

        /** Deletes all elements in every lines, but does not change the sparse matrix size, line capacity, type or symmetry status.
          * It is just like setting lines to zero.
          */
        void clear()
        {
            for (int i = 0; i < mv_line.size(); ++i) {
                mv_line.at(i).clear();
            }
        }

        /** Reserve capacity to a specific line. Notice that size is not checked (if it is smaller than current).
          * @param lineIndex Index of the line to be modified.
          * @param size New size to be allocated.
          */
        void setLineSize(const unsigned int lineIndex, const unsigned int size) {
            mv_line.at(lineIndex).reserve(size);
        }

        /** Reserve capacity to all lines of the sparse matrix. Notice that size is not checked (if it is smaller than current).
          * @param size New size to be allocated.
          */
        void setAllLineSize(const unsigned int size) {
            for (int i = 0; i < mv_line.size(); i++)
                setLineSize(i, size);
        }

        /** If there are terms with repeated index, these are summed
          */
        void addEqualTerms() {
            bool mi_assembled;
            for (int i = 0; i < mv_line.size(); i++) {
                /* NOW sortLine IS CALLED INSIDE addEqualTerms*/
                /* sort(); */
                mi_assembled = mv_line.at(i).addEqualTerms();
                assert(mi_assembled); // ERROR: Probably, a line is empty!
            }
        }

        /** Push a M2S2::MatrixX to the asymmetric sparse matrix
          * @param matrix M2S2::MatrixX to be included (must be Square).
          * @param indexes Vector with global indexes (dof).
          */
        void push(M2S2::MatrixX& matrix, const std::vector<int>& indexes)
        {
            assert(!mv_sym); // Cannot push an asymmetric matrix to a symmetric sparse matrix!
            assert(matrix.cols() == indexes.size() && matrix.rows() == indexes.size()); // Matrix and indexing are not of the same size

            // Reallocate sparse matrix lines when necessary
            reallocateLines(indexes);

            // Auxiliary
            int row, col;

            if (mv_type == 'R') {
                for (int i = 0; i < matrix.rows(); i++) {
                    row = indexes.at(i);
                    for (int j = 0; j < matrix.cols(); j++) {
                        col = indexes.at(j);

                        mv_line.at(row).mv_index.push_back(col);
                        mv_line.at(row).mv_value.push_back(matrix.at(i, j));
                    }
                    mv_line.at(row).mv_assembled = false;
                }
            }
            else {
                for (int i = 0; i < matrix.rows(); i++) {
                    row = indexes.at(i);
                    for (int j = 0; j < matrix.cols(); j++) {
                        col = indexes.at(j);

                        mv_line.at(col).mv_index.push_back(row);
                        mv_line.at(col).mv_value.push_back(matrix.at(i, j));
                    }
                    mv_line.at(row).mv_assembled = false;
                }
            }
        }

        /** Push a M2S2::MatrixS to the sparse matrix.
          * @param matrix M2S2::MatrixS to be included.
          * @param indexes Vector with global indexes (dof).
          */
        void push(M2S2::MatrixS& matrix, const std::vector<int>& indexes)
        {
            assert(matrix.rows() == indexes.size()); // Matrix and indexing are not of the same size

            // Reallocate sparse matrix lines when necessary
            reallocateLines(indexes);

            // Auxiliary
            int row, col;

            if (mv_sym) {
                if (mv_type == 'R') {   // col must be >= row
                    for (int i = 0; i < matrix.rows(); i++) {
                        row = indexes.at(i);
                        for (int j = i; j < matrix.cols(); j++) {
                            col = indexes.at(j);

                            if (row < col) {
                                mv_line.at(row).mv_index.push_back(col);
                                mv_line.at(row).mv_value.push_back(matrix.at(i, j));
                            }
                            else {
                                mv_line.at(col).mv_index.push_back(row);
                                mv_line.at(col).mv_value.push_back(matrix.at(i, j));
                            }

                            // FIX: If there are repetead indexes, must also push row > col
                            // Only happens for symmetric matrices
                            if (row == col && i != j) {
                                mv_line.at(row).mv_index.push_back(col);
                                mv_line.at(row).mv_value.push_back(matrix.at(i, j));
                            }
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
                else { // row must be >= col
                    for (int i = 0; i < matrix.rows(); i++) {
                        row = indexes.at(i);
                        for (int j = i; j < matrix.cols(); j++) {
                            col = indexes.at(j);

                            if (col < row) {
                                mv_line.at(col).mv_index.push_back(row);
                                mv_line.at(col).mv_value.push_back(matrix.at(i, j));
                            }
                            else {
                                mv_line.at(row).mv_index.push_back(col);
                                mv_line.at(row).mv_value.push_back(matrix.at(i, j));
                            }

                            // FIX: If there are repetead indexes, must also push row > col
                            // Only happens for symmetric matrices
                            if (row == col && i != j) {
                                mv_line.at(col).mv_index.push_back(row);
                                mv_line.at(col).mv_value.push_back(matrix.at(i, j));
                            }
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
            }
            else {
                if (mv_type == 'R') {
                    for (int i = 0; i < matrix.rows(); i++) {
                        row = indexes.at(i);
                        for (int j = 0; j < matrix.cols(); j++) {
                            col = indexes.at(j);

                            mv_line.at(row).mv_index.push_back(col);
                            mv_line.at(row).mv_value.push_back(matrix.at(i, j));
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
                else {
                    for (int i = 0; i < matrix.rows(); i++) {
                        row = indexes.at(i);
                        for (int j = 0; j < matrix.cols(); j++) {
                            col = indexes.at(j);

                            mv_line.at(col).mv_index.push_back(row);
                            mv_line.at(col).mv_value.push_back(matrix.at(i, j));
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
            }
        }

        /** Push a double** to the sparse matrix (C-like matrices)
          * @param matrix double** to be included (must be Square).
          * @param indexes Vector with global indexes (dof).
          * @param nItens Size of matrix (number of row/columns)
          */
        void push(const double** matrix, const int* indexes, const unsigned int nItens)
        {
            // Reallocate sparse matrix lines when necessary
            reallocateLines(indexes, nItens);

            // Auxiliary
            int row, col;

            if (mv_sym) {
                if (mv_type == 'R') {
                    for (int i = 0; i < nItens; i++) {
                        row = indexes[i];
                        for (int j = i; j < nItens; j++) {
                            col = indexes[j];

                            // FIX: Repeated indexes are now pushed correctly
                            if (row < col) {
                                mv_line.at(row).mv_index.push_back(col);
                                mv_line.at(row).mv_value.push_back(matrix[i][j]);
                            }
                            else {
                                mv_line.at(col).mv_index.push_back(row);
                                mv_line.at(col).mv_value.push_back(matrix[i][j]);
                            }

                            // FIX: If there are repetead indexes, must also push row > col
                            // Only happens for symmetric matrices
                            if (row == col && i != j) {
                                mv_line.at(row).mv_index.push_back(col);
                                mv_line.at(row).mv_value.push_back(matrix[i][j]);
                            }
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
                else {
                    for (int i = 0; i < nItens; i++) {
                        row = indexes[i];
                        for (int j = i; j < nItens; j++) {
                            col = indexes[j];

                            // FIX: Repeated indexes are now pushed correctly
                            if (col < row) {
                                mv_line.at(col).mv_index.push_back(row);
                                mv_line.at(col).mv_value.push_back(matrix[i][j]);
                            }
                            else {
                                mv_line.at(row).mv_index.push_back(col);
                                mv_line.at(row).mv_value.push_back(matrix[i][j]);
                            }

                            // FIX: If there are repetead indexes, must also push row > col
                            // Only happens for symmetric matrices
                            if (row == col && i != j) {
                                mv_line.at(col).mv_index.push_back(row);
                                mv_line.at(col).mv_value.push_back(matrix[i][j]);
                            }
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
            }
            else {
                if (mv_type == 'R') {
                    for (int i = 0; i < nItens; i++) {
                        row = indexes[i];
                        for (int j = 0; j < nItens; j++) {
                            col = indexes[j];

                            mv_line.at(row).mv_index.push_back(col);
                            mv_line.at(row).mv_value.push_back(matrix[i][j]);
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
                else {
                    for (int i = 0; i < nItens; i++) {
                        row = indexes[i];
                        for (int j = 0; j < nItens; j++) {
                            col = indexes[j];

                            mv_line.at(col).mv_index.push_back(row);
                            mv_line.at(col).mv_value.push_back(matrix[i][j]);
                        }
                        mv_line.at(row).mv_assembled = false;
                    }
                }
            }
        }

        /** Push a vector of triplets to the sparse matrix
          * @param triplets values to be pushed. 
          */
        void push(const std::vector<M2S2::triplet>& triplets)
        {
            int row_ind, col_ind;
            double val;
            for (int i = 0; i < triplets.size(); i++) {
                row_ind = triplets.at(i).mv_row;
                col_ind = triplets.at(i).mv_col;
                val = triplets.at(i).mv_value;

                if (mv_type == 'R') {
                    if (mv_line.at(row_ind).mv_index.capacity() <= mv_line.at(row_ind).mv_index.size()) {
                        mv_line.at(row_ind).resize();
                    }
                    mv_line.at(row_ind).mv_index.push_back(col_ind);
                    mv_line.at(row_ind).mv_value.push_back(val);
                    mv_line.at(row_ind).mv_assembled = 0;
                }
                else {
                    if (mv_line.at(col_ind).mv_index.capacity() <= mv_line.at(col_ind).mv_index.size()) {
                        mv_line.at(col_ind).resize();
                    }
                    mv_line.at(col_ind).mv_index.push_back(row_ind);
                    mv_line.at(col_ind).mv_value.push_back(val);
                    mv_line.at(col_ind).mv_assembled = 0;
                }
            }
        }

        /** Push triplets to the sparse matrix (C-like triplets)
          * @param row C-like vector with nItens containing row indexes.
          * @param col C-like vector with nItens containing column indexes.
          * @param value C-like vector with nItens containing values.
          * @param nItens Number of itens to be pushed.
          */
        void push(const int* row, const int* col, const double* value, const unsigned int nItens)
        {
            int row_ind, col_ind;
            double val;
            for (int i = 0; i < nItens; i++) {
                row_ind = row[i];
                col_ind = col[i];
                val = value[i];

                if (mv_type == 'R') {
                    if (mv_line.at(row_ind).mv_index.capacity() <= mv_line.at(row_ind).mv_index.size()) {
                        mv_line.at(row_ind).resize();
                    }
                    mv_line.at(row_ind).mv_index.push_back(col_ind);
                    mv_line.at(row_ind).mv_value.push_back(val);
                    mv_line.at(row_ind).mv_assembled = 0;
                }
                else {
                    if (mv_line.at(col_ind).mv_index.capacity() <= mv_line.at(col_ind).mv_index.size()) {
                        mv_line.at(col_ind).resize();
                    }
                    mv_line.at(col_ind).mv_index.push_back(row_ind);
                    mv_line.at(col_ind).mv_value.push_back(val);
                    mv_line.at(col_ind).mv_assembled = 0;
                }
            }
        }

        /** Retrieve a copy of a row from the sparse matrix
          * @param index Container index of the sparse matrix (begins with 0)
          * @return M2S2::line with row information
          */
        M2S2::line getRow(const unsigned int index)
        {
            // if the sparse matrix was saved as column major, and it is not symmetric, this search will take a while...
            if (mv_type == 'C' && !mv_sym) {
                int size = mv_line.size();
                M2S2::line mi_line(size);

                int ind;
                for (int i = 0; i < mv_line.size(); i++) {
                    ind = mv_line.at(i).search(index);

                    if (ind >= 0) {
                        mi_line.mv_index.push_back(i);
                        mi_line.mv_value.push_back(mv_line.at(i).mv_value.at(ind));
                    }
                }
                return mi_line;
            }
            else {
                return mv_line.at(index);
            }
        }

        /** Retrieve a copy of a column from the sparse matrix
          * @param index Container index of the sparse matrix (begins with 0)
          * @return M2S2::line with column information
          */
        M2S2::line getColumn(const unsigned int index)
        {
            // if the sparse matrix was saved as row major, and it is not symmetric, this search will take a while...
            if (mv_type == 'R' && !mv_sym) {
                int size = mv_line.size();
                M2S2::line mi_line(size);

                int ind;
                for (int i = 0; i < mv_line.size(); i++) {
                    ind = mv_line.at(i).search(index);

                    if (ind >= 0) {
                        mi_line.mv_index.push_back(i);
                        mi_line.mv_value.push_back(mv_line.at(i).mv_value.at(ind));
                    }
                }
                return mi_line;
            }
            else {
                return mv_line.at(index);
            }
        }

        /** Access a term from the sparse matrix
          * @param row Row index component
          * @param col Column index component
          * @return Required term, if any
          */
        inline double& at(const unsigned int row, const unsigned int col)
        {
            int ind;

            if (mv_type == 'C') {
                ind = mv_line.at(col).search(row);

                if (ind >= 0)
                    return mv_line.at(col).mv_value.at(ind);
                else {
                    throw(std::out_of_range("Value not found!"));
                    //return *mv_line.at(col).mv_value.end();
                }
            }
            else { // type is 'R'
                ind = mv_line.at(row).search(col);

                if (ind >= 0)
                    return mv_line.at(row).mv_value.at(ind);
                else {
                    throw(std::out_of_range("Value not found!"));
                    //return *mv_line.at(col).mv_value.end();
                }
            }
        }

        /** Access a term from the sparse matrix
          * @param row Row index component
          * @param col Column index component
          * @return Required term, if any
          */
        inline double& operator () (const unsigned int row, const unsigned int col) {
            return at(row, col);
        }

        /** Set a single value to all row components
          * @param index Row index component
          * @param value valeu to be imposed
          */
        void setValueToRow(const unsigned int index, const double value)
        {
            if (mv_type == 'C') {
                int ind;

                for (int i = 0; i < mv_line.size(); i++) {
                    ind = mv_line.at(i).search(index);
                    if (ind >= 0)
                        mv_line.at(i).mv_value.at(ind) = value;
                }
            }
            else {
                for (int i = 0; i < mv_line.at(index).mv_index.size(); i++) {
                    mv_line.at(index).mv_value.at(i) = value;
                }
            }
        }

        /** Set a single value to all column components
          * @param index Row index component
          * @param value valeu to be imposed
          */
        void setValueToColumn(const unsigned int index, const double value)
        {
            if (mv_type == 'R') {
                int ind;
                for (int i = 0; i < mv_line.size(); i++) {
                    ind = mv_line.at(i).search(index);
                    if (ind >= 0)
                        mv_line.at(i).mv_value.at(ind) = value;
                }
            }
            else {
                for (int i = 0; i < mv_line.at(index).mv_index.size(); i++) {
                    mv_line.at(index).mv_value.at(i) = value;
                }
            }
        }

        /** Set row and column to zero, main diagonal to one
          * @param index Row / Column index component
          */
        void zeroRowAndColumn(const unsigned int index) {
            setValueToRow(index, 0.0);
            setValueToColumn(index, 0.0);
            at(index, index) = 1.0;
        }

    private:
        /** Reallocate lines of the sparse matrix following indexes.
          * @param indexes Vector with global indexes (dof).
          */
        void reallocateLines(const std::vector<int>& indexes) {
            int mi_curIndex;        // Current index

            for (int i = 0; i < indexes.size(); ++i) {
                mi_curIndex = indexes.at(i);

                if (mv_line.at(mi_curIndex).mv_index.capacity() - mv_line.at(mi_curIndex).mv_index.size() < indexes.size()) {
                    mv_line.at(mi_curIndex).resize();
                }
            }
        }

        /** Reallocate lines of the sparse matrix following indexes.
          * @param indexes Vector with global indexes (dof).
          * @param nItens Number of itens to be pushed
          */
        void reallocateLines(const int* indexes, const unsigned int nItens) {
            int mi_curIndex;        // Current index

            for (int i = 0; i < nItens; ++i) {
                mi_curIndex = indexes[i];

                if (mv_line.at(mi_curIndex).mv_index.capacity() - mv_line.at(mi_curIndex).mv_index.size() < nItens) {
                    mv_line.at(mi_curIndex).resize();
                }
            }
        }

    };
}  // End of M2S2 namespace

#endif
