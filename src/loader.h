#pragma once

#include <string>
#include <vector>

namespace s21 {

class Loader {
 public:
  static std::string getFirstDate(const std::string& filePath);
  static std::vector<std::pair<double, double>> getDataNorm(
      const std::string& filePath);
  static std::vector<std::pair<std::string, std::string>> loadFile(
      const std::string& filePath);
  static void validation(const std::string& date);
  static double convertDateTime(std::string buff);
  static std::vector<double> convertAllDAte(
      const std::vector<std::string>& date);
  static std::vector<double> loadWeights(std::string filePath);
  static std::vector<std::string> splitString(const std::string& line,
                                              const std::string& splitters);
};
}  // namespace s21
