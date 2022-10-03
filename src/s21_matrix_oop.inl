
//______________________Constructor______________________
template <typename T>
S21Matrix<T>::S21Matrix() {
  _rows = 2;
  _cols = 2;
  _matrix = nullptr;
  memoryAllocation();
}

template <typename T>
S21Matrix<T>::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    MatrixException nonExist(
        "The number of rows/columns of the matrix is less than 1");
    throw nonExist;
  }
  _rows = rows;
  _cols = cols;
  _matrix = nullptr;
  memoryAllocation();
}

template <typename T>
S21Matrix<T>::S21Matrix(const S21Matrix<T>& other) {
  _matrix = nullptr;
  _rows = 0;
  _cols = 0;
  *this = other;
}

template <typename T>
S21Matrix<T>::S21Matrix(S21Matrix<T>&& other) {
  _matrix = nullptr;
  _rows = 0;
  _cols = 0;
  *this = std::move(other);
}

//______________________Destructor______________________
template <typename T>
S21Matrix<T>::~S21Matrix() {
  memoryClean();
}

//______________________Function______________________
template <typename T>
void S21Matrix<T>::memoryAllocation() {
  _matrix = new T*[_rows];
  for (int i = 0; i < _rows; i++) {
    _matrix[i] = new T[_cols]{};
  }
}

template <typename T>
void S21Matrix<T>::memoryClean() {
  if (_matrix != nullptr) {
    for (int i = 0; i < _rows; i++) {
      delete[] _matrix[i];
    }
    delete[] _matrix;
    _matrix = nullptr;
    _rows = 0;
    _cols = 0;
  }
}

template <typename T>
void S21Matrix<T>::sumOrSub(const S21Matrix<T>& other, int sign) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = _matrix[i][j] + (sign * other._matrix[i][j]);
    }
  }
}

template <typename T>
bool S21Matrix<T>::eq_matrix(const S21Matrix<T>& other) {
  bool rez = false;
  if (_rows == other._rows && _cols == other._cols) {
    rez = true;
    for (int i = 0; i < _rows && rez == true; i++) {
      for (int j = 0; j < _cols && rez == true; j++) {
        if (fabs(_matrix[i][j] - other._matrix[i][j]) >= 1e-7) {
          rez = false;
        }
      }
    }
  }
  return rez;
}

template <typename T>
void S21Matrix<T>::sum_matrix(const S21Matrix<T>& other) {
  if (_rows != other._rows || _cols != other._cols) {
    MatrixException inconsistency("Rows/columns of matrices are not equal");
    throw inconsistency;
  }
  sumOrSub(other, 1);
}

template <typename T>
void S21Matrix<T>::sub_matrix(const S21Matrix<T>& other) {
  if (_rows != other._rows || _cols != other._cols) {
    MatrixException inconsistency("Rows/columns of matrices are not equal");
    throw inconsistency;
  }
  sumOrSub(other, -1);
}

template <typename T>
void S21Matrix<T>::mul_number(const T& val) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = _matrix[i][j] * val;
    }
  }
}

template <typename T>
void S21Matrix<T>::mul_matrix(const S21Matrix<T>& other) {
  if (_cols != other._rows) {
    MatrixException nonExist("Columns are not equal to rows");
    throw nonExist;
  }
  S21Matrix<T> buf(_rows, other._cols);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      for (int n = 0; n < other._rows; n++) {
        buf._matrix[i][j] += _matrix[i][n] * other._matrix[n][j];
      }
    }
  }
  *this = buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::transpose() {
  S21Matrix<T> buf(_cols, _rows);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      buf._matrix[j][i] = _matrix[i][j];
    }
  }
  return buf;
}

template <typename T>
void S21Matrix<T>::fillMatrix(const T& val) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = val;
    }
  }
}

template <typename T>
void S21Matrix<T>::printMatrix() {
  std::cout << "\u001b[1;32;1;1m";
  for (int i = 0; i < get_rows(); ++i) {
    for (int j = 0; j < get_cols(); ++j) {
      std::cout << this->operator()(i, j);
      if (this->operator()(i, j) < 10)
        std::cout << "   ";
      else if (this->operator()(i, j) >= 10 && this->operator()(i, j) < 100)
        std::cout << "  ";
      else
        std::cout << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "\n"
            << "\u001b[0m";
}

//________________________Operator________________________
template <typename T>
S21Matrix<T>& S21Matrix<T>::operator=(const S21Matrix<T>& other) {
  if (this == &other) {
    MatrixException nonExist("Cannot assign an object's value to itself");
    throw nonExist;
  }
  memoryClean();
  _rows = other._rows;
  _cols = other._cols;
  memoryAllocation();

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = other._matrix[i][j];
    }
  }
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator=(S21Matrix<T>&& other) {
  if (this == &other) {
    MatrixException nonExist("Cannot assign an object's value to itself");
    throw nonExist;
  }
  memoryClean();
  std::swap(_rows, other._rows);
  std::swap(_cols, other._cols);
  std::swap(_matrix, other._matrix);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator+=(const S21Matrix<T>& other) {
  sum_matrix(other);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator-=(const S21Matrix<T>& other) {
  sub_matrix(other);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator*=(const T& val) {
  mul_number(val);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator*=(const S21Matrix<T>& other) {
  mul_matrix(other);
  return *this;
}

template <typename T>
bool S21Matrix<T>::operator==(const S21Matrix<T>& other) {
  return eq_matrix(other);
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator+(const S21Matrix<T>& other) {
  S21Matrix<T> buf(*this);
  buf.sum_matrix(other);
  return buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator-(const S21Matrix<T>& other) {
  S21Matrix<T> buf(*this);
  buf.sub_matrix(other);
  return buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator*(const T& val) {
  S21Matrix<T> buf(*this);
  buf.mul_number(val);
  return buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator*(const S21Matrix<T>& other) {
  S21Matrix<T> buf(*this);
  buf.mul_matrix(other);
  return buf;
}

template <typename T>
T& S21Matrix<T>::operator()(int row, int col) {
  if (row < 0 || col < 0 || row >= _rows || col >= _cols) {
    MatrixException non("Incorrect row/col");
    throw non;
  }
  return _matrix[row][col];
}

template <typename T>
T S21Matrix<T>::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row >= _rows || col >= _cols) {
    MatrixException non("Incorrect row/col");
    throw non;
  }
  return _matrix[row][col];
}

//________________________GET SET__________________________
template <typename T>
int S21Matrix<T>::get_rows() {
  return _rows;
}

template <typename T>
int S21Matrix<T>::get_rows() const {
  return _rows;
}

template <typename T>
int S21Matrix<T>::get_cols() {
  return _cols;
}

template <typename T>
int S21Matrix<T>::get_cols() const {
  return _cols;
}

template <typename T>
void S21Matrix<T>::set_rows(int rows) {
  S21Matrix<T> buf(rows, _cols);
  int minRow = _rows > rows ? rows : _rows;
  for (int i = 0; i < minRow; i++) {
    for (int j = 0; j < _cols; j++) {
      buf._matrix[i][j] = _matrix[i][j];
    }
  }
  *this = buf;
}

template <typename T>
void S21Matrix<T>::set_cols(int cols) {
  S21Matrix<T> buf(_rows, cols);
  int minCol = _cols > cols ? cols : _cols;
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < minCol; j++) {
      buf._matrix[i][j] = _matrix[i][j];
    }
  }
  *this = buf;
}

//___________________Template specialization______________________
template <typename T>
S21Matrix<T> S21Matrix<T>::getReducedMatrix(int rowsReduced, int colsReduced) {
  S21Matrix<T> buf(_rows - 1, _cols - 1);
  int newRows = 0, newCols = 0;
  for (int i = 0; i < _rows; i++) {
    if (i != rowsReduced) {
      for (int j = 0; j < _rows; j++) {
        if (j != colsReduced) {
          buf._matrix[newRows][newCols] = _matrix[i][j];
          newCols++;
        }
      }
      newRows++;
      newCols = 0;
    }
  }
  return buf;
}

template <typename T>
double S21Matrix<T>::determinant() {
  double det = 0.0;
  if (_rows != _cols) {
    MatrixException nonExist("Error, matrix is not square");
    throw nonExist;
  }
  if (_rows == 1) {
    det = _matrix[0][0];
  } else if (_rows == 2) {
    det = _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
  } else {
    for (int i = 0; i < _rows; i++) {
      S21Matrix<T> buf = getReducedMatrix(i, 0);
      det += _matrix[i][0] * pow(-1, i) * buf.determinant();
    }
  }
  return det;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::calc_complements() {
  if (_rows != _cols || _rows == 1) {
    MatrixException nonExist("Error, matrix is not square");
    throw nonExist;
  }
  S21Matrix<T> rez(_rows, _cols);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      S21Matrix<T> buf = getReducedMatrix(i, j);
      rez._matrix[i][j] = pow(-1, i + j) * buf.determinant();
    }
  }
  return rez;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::inverse_matrix() {
  double det = determinant();
  if (_rows != _cols || det == 0.0 || _rows == 1) {
    MatrixException nonExist("Determinant is zero, or matrix is not square");
    throw nonExist;
  }
  S21Matrix<T> rez(_rows, _cols);
  rez = calc_complements();
  rez = rez.transpose();
  rez.mul_number(1.0 / det);
  return rez;
}
