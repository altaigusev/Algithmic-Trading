#pragma once

#include <cmath>
#include <cstring>
#include <exception>
#include <iostream>

//___________________EXCEPTION CLASS FOR MATRIX__________________
class MatrixException : public std::exception {
 public:
  explicit MatrixException(const std::string& strExcep) : str(strExcep) {}
  ~MatrixException() throw() {}
  const char* what() const throw() { return str.c_str(); }

 private:
  std::string str;
};

template <class T>
class S21Matrix {
 public:
  //___________________CONSTRUCTORS__________________
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix<T>& other);
  S21Matrix(S21Matrix<T>&& other);

  //___________________DESTRUCTOR__________________
  ~S21Matrix();

  //___________________FUNCTIONS___________________
  bool eq_matrix(const S21Matrix<T>& other);
  void sum_matrix(const S21Matrix<T>& other);
  void sub_matrix(const S21Matrix<T>& other);
  void mul_number(const T& val);
  void mul_matrix(const S21Matrix<T>& other);
  S21Matrix<T> transpose();
  void fillMatrix(const T& val);
  void printMatrix();

  //___________________OPERATOR___________________
  S21Matrix& operator=(const S21Matrix<T>& other);
  S21Matrix& operator=(S21Matrix<T>&& other);
  S21Matrix<T>& operator+=(const S21Matrix& other);
  S21Matrix<T>& operator-=(const S21Matrix<T>& other);
  S21Matrix<T>& operator*=(const S21Matrix<T>& other);
  S21Matrix<T>& operator*=(const T& val);
  S21Matrix<T> operator+(const S21Matrix<T>& other);
  S21Matrix<T> operator-(const S21Matrix<T>& other);
  S21Matrix<T> operator*(const T& val);
  S21Matrix<T> operator*(const S21Matrix<T>& other);
  bool operator==(const S21Matrix<T>& other);
  T& operator()(int row, int col);
  T operator()(int row, int col) const;

  //___________________GET SET___________________
  int get_rows();
  int get_rows() const;
  void set_rows(int rows);
  int get_cols();
  int get_cols() const;
  void set_cols(int cols);

  //___________________TEMPLATE SPECIALIZATION___________________
  S21Matrix<T> getReducedMatrix(int rowsReduced, int colsReduced);
  double determinant();
  S21Matrix<T> calc_complements();
  S21Matrix<T> inverse_matrix();

 private:
  //___________________HELPER FOO___________________
  void memoryAllocation();
  void memoryClean();
  void sumOrSub(const S21Matrix<T>& other, int sign);

  //_________________________________________________________
  int _rows = 0, _cols = 0;
  T** _matrix = nullptr;
};

#include "s21_matrix_oop.inl"
