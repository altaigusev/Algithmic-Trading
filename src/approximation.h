#pragma once

#include <string>
#include <vector>

template <class T>
class S21Matrix;

namespace s21 {
class Approximation {
 public:
  static std::vector<double> getX(
      const std::vector<std::pair<double, double>>& data, const int& pow);

  static std::vector<double> getY(
      const std::vector<std::pair<double, double>>& data, const int& pow);

  static double getResult(const int& pow,
                          const std::vector<std::pair<double, double>>& data,
                          const double& x);

  //  double calcNewY(const unsigned long& x,
  //                  const std::vector<double>& coefficients);

  static double calcNewY(const double& x,
                         const std::vector<double>& coefficients);

  static double calcNewYWithW(const double& newY, const double& k);

  static std::vector<double> calcALLNewY(
      const int& pow, const std::vector<std::pair<double, double>>& data,
      int countPoints);

  static std::vector<double> calcAllNewYWithW(
      const std::vector<double>& k, const int& pow,
      const std::vector<std::pair<double, double>>& data, int countPoints);

  static S21Matrix<double> calculateFirstMatrix(
      const int& sizeMatrix, const int& dataSize,
      const std::vector<double>& powXSum);

  static std::vector<double> calculateCoefficients(
      const S21Matrix<double>& first, const std::vector<double>& tmpSumY);

  static std::vector<std::pair<double, double>> getNewData(
      const std::vector<double>& k, const int& pow,
      const std::vector<std::pair<double, double>>& data);
};
//
}  // namespace s21
