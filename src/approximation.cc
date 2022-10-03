#include "approximation.h"

#include <cmath>
#include <stdexcept>

#include "s21_matrix_oop.h"

namespace s21 {

std::vector<double> Approximation::getX(
    const std::vector<std::pair<double, double>>& data, const int& pow) {
  if (pow < 1 || pow > 6 || data.size() < 1) {
    throw std::invalid_argument("Error incorrect pow or size datasets");
  }
  std::vector<double> rezX;
  for (int i = pow * 2; i >= 1; --i) {
    double tmp = 0.0;
    for (auto j : data) {
      tmp += std::pow(j.second, i);
    }
    rezX.push_back(tmp);
  }
  return rezX;
}

std::vector<double> Approximation::getY(
    const std::vector<std::pair<double, double>>& data, const int& pow) {
  std::vector<double> rezY;
  for (int i = pow; i >= 1; --i) {
    double tmp = 0.0;
    for (auto j : data) {
      tmp += j.first * std::pow(j.second, i);
    }
    rezY.push_back(tmp);
  }
  double tmpSumY = 0.0;
  for (auto i : data) {
    tmpSumY += i.first;
  }
  rezY.push_back(tmpSumY);
  return rezY;
}

double Approximation::getResult(
    const int& pow, const std::vector<std::pair<double, double>>& data,
    const double& x) {
  std::vector<double> rezX = getX(data, pow);
  std::vector<double> rezY = getY(data, pow);
  S21Matrix<double> firstMatrix =
      calculateFirstMatrix(pow + 1, data.size(), rezX);
  std::vector<double> coefficients = calculateCoefficients(firstMatrix, rezY);
  double newY = calcNewY(x, coefficients);
  return newY;
}

S21Matrix<double> Approximation::calculateFirstMatrix(
    const int& sizeMatrix, const int& dataSize,
    const std::vector<double>& powXSum) {
  S21Matrix<double> result(sizeMatrix, sizeMatrix);
  result.fillMatrix(0.0);
  for (int row = 0; row < sizeMatrix; ++row) {
    for (int col = 0; col < sizeMatrix; ++col) {
      if (row == sizeMatrix - 1 && col == sizeMatrix - 1) {
        result.operator()(row, col) = dataSize;
      } else {
        result.operator()(row, col) = powXSum[col + row];
      }
    }
  }
  return result;
}

std::vector<double> Approximation::calculateCoefficients(
    const S21Matrix<double>& first, const std::vector<double>& tmpSumY) {
  std::vector<double> result;
  std::vector<S21Matrix<double>> tmp;
  for (int countK = 0; countK < (int)tmpSumY.size(); countK++) {
    S21Matrix<double> buf(first);
    for (int row = 0; row < first.get_rows(); ++row) {
      buf.operator()(row, countK) = tmpSumY[row];
    }
    tmp.push_back(buf);
  }
  for (auto x : tmp) {
    S21Matrix<double> buf(first);
    result.push_back((x.determinant() / buf.determinant()));
  }
  return result;
}

// double Approximation::calcNewY(const unsigned long& x,
//                                const std::vector<double>& coefficients) {
//   double y = 0.0;
//   for (int i = 0; i < (int)coefficients.size(); ++i) {
//     y += coefficients[i] * (std::pow(x, (coefficients.size() - 1) - i));
//   }
//   return y;
// }

double Approximation::calcNewY(const double& x,
                               const std::vector<double>& coefficients) {
  double y = 0.0;
  for (int i = 0; i < (int)coefficients.size(); ++i) {
    y += coefficients[i] * (std::pow(x, (coefficients.size() - 1) - i));
  }
  return y;
}

std::vector<double> Approximation::calcALLNewY(
    const int& pow, const std::vector<std::pair<double, double>>& data,
    int countPoints) {
  std::vector<double> rezX = getX(data, pow);
  std::vector<double> rezY = getY(data, pow);
  S21Matrix<double> firstMatrix =
      calculateFirstMatrix(pow + 1, data.size(), rezX);
  std::vector<double> coefficients = calculateCoefficients(firstMatrix, rezY);
  std::vector<double> allNewY;
  std::vector<double> newX;
  double step = data[data.size() - 1].second / countPoints;
  for (double i = 1.0; i < data[data.size() - 1].second; i += step) {
    newX.push_back(i);
  }
  for (int i = 0; i < (int)newX.size(); ++i) {
    allNewY.push_back(calcNewY(newX[i], coefficients));
  }
  return allNewY;
}

double Approximation::calcNewYWithW(const double& newY, const double& k) {
  if (k == 0 || k > 1 || k < -1) {
    throw std::invalid_argument("Incorrect weight");
  }
  return newY * k;
}

std::vector<double> Approximation::calcAllNewYWithW(
    const std::vector<double>& k, const int& pow,
    const std::vector<std::pair<double, double>>& data, int countPoints) {
  if (k.size() != data.size()) {
    throw std::invalid_argument("Count weight != count data");
  }
  std::vector<std::pair<double, double>> data2 = getNewData(k, pow, data);
  std::vector<double> rez = calcALLNewY(pow, data2, countPoints);
  return rez;
}

std::vector<std::pair<double, double>> Approximation::getNewData(
    const std::vector<double>& k, const int& pow,
    const std::vector<std::pair<double, double>>& data) {
  if (k.size() != data.size()) {
    throw std::invalid_argument("Count weight != count data");
  }
  std::vector<double> allNewY = calcALLNewY(pow, data, data.size());
  std::vector<double> allNewYWithK;
  for (int i = 0; i < (int)allNewY.size(); ++i) {
    allNewYWithK.push_back(calcNewYWithW(allNewY[i], k[i]));
  }
  std::vector<std::pair<double, double>> data2;
  for (int i = 0; i < (int)allNewYWithK.size(); ++i) {
    data2.push_back({allNewYWithK[i], data[i].second});
  }
  return data2;
}

}  // namespace s21
