#include "loader.h"

#include <fstream>
#include <set>

namespace s21 {

std::string Loader::getFirstDate(const std::string& filePath) {
  std::vector<std::pair<std::string, std::string>> result = loadFile(filePath);
  std::set<std::string> tmp;
  for (auto x : result) {
    tmp.insert(x.second);
  }
  return *(tmp.begin());
}

std::vector<std::pair<double, double>> Loader::getDataNorm(
    const std::string& filePath) {
  //  первоя строка - У(стоимость), вторая - Х(дата)
  std::vector<std::pair<std::string, std::string>> result = loadFile(filePath);
  std::vector<std::pair<double, double>> data;
  std::vector<std::string> tmpDate;
  for (auto x : result) {
    tmpDate.push_back(x.second);
  }
  std::vector<double> normalDate = convertAllDAte(tmpDate);
  for (int i = 0; i < (int)result.size(); ++i) {
    validation(result[i].second);
    double cost = std::stod(result[i].first);
    double date = normalDate[i];
    std::pair<double, double> tmp = {cost, date};
    data.push_back(tmp);
  }
  return data;
}

std::vector<double> Loader::loadWeights(std::string filePath) {
  std::ifstream ifs(filePath);
  if (!ifs.is_open()) {
    throw std::runtime_error("Can`t open file");
  }
  if (ifs.eof() == true) {
    throw std::runtime_error("Empty file");
  }
  std::string line;
  getline(ifs, line);
  if (line.compare("Date,Close") != 0) {
    throw std::runtime_error("Error data on file");
  }
  std::vector<double> rez;

  while (getline(ifs, line)) {
    std::vector<std::string> result = splitString(line, ",");
    if (result.size() < 2 || result.size() > 3) {
      throw std::runtime_error("Error data on file");
    }
    if (result.size() == 3) {
      rez.push_back(std::stod(result[2]));
    } else {
      rez.push_back(1.0);
    }
  }
  ifs.close();
  return rez;
}

std::vector<std::pair<std::string, std::string>> Loader::loadFile(
    const std::string& filePath) {
  std::ifstream ifs(filePath);
  if (!ifs.is_open()) {
    throw std::runtime_error("Can`t open file");
  }
  if (ifs.eof() == true) {
    throw std::runtime_error("Empty file");
  }
  std::string line;
  getline(ifs, line);
  if (line.compare("Date,Close") != 0) {
    throw std::runtime_error("Error data on file");
  }
  std::vector<std::pair<std::string, std::string>> data;
  while (getline(ifs, line)) {
    std::vector<std::string> result = splitString(line, ",");
    if (result.size() < 2 || result.size() > 3) {
      throw std::runtime_error("Error data on file");
    }
    data.push_back({result[1], result[0]});
  }
  ifs.close();
  return data;
}

void Loader::validation(const std::string& date) {
  std::vector<std::string> tmp = splitString(date, "-");
  if (tmp.size() != 3 || std::stod(tmp[0]) < 1970 || std::stod(tmp[0]) > 2099 ||
      std::stoi(tmp[1]) < 1 || std::stod(tmp[1]) > 12 ||
      std::stoi(tmp[2]) < 1 || std::stod(tmp[2]) > 31) {
    throw std::invalid_argument("Error invalid date");
  }
}

std::vector<std::string> Loader::splitString(const std::string& line,
                                             const std::string& splitters) {
  std::vector<std::string> result;
  std::string::size_type start(0), finish(0);
  while (finish != std::string::npos) {
    while (splitters.find_first_of(line[start]) != std::string::npos) ++start;
    if (start == std::string::npos) break;
    finish = line.find_first_of(splitters, start);
    result.emplace_back(line.substr(start, finish - start));
    if (finish == std::string::npos) break;
    start = finish;
  }
  return result;
}

double Loader::convertDateTime(std::string buff) {
  int yy, mm, dd;
  struct tm when;
  long tme;
  std::vector<std::string> tmp = splitString(buff, "-");
  yy = std::stoi(tmp[0]);
  mm = std::stoi(tmp[1]);
  dd = std::stoi(tmp[2]);
  time(&tme);
  when = *localtime(&tme);
  when.tm_year = yy;
  when.tm_mon = mm - 1;
  when.tm_mday = dd;
  double rez = (mktime(&when));
  return rez;
}

std::vector<double> Loader::convertAllDAte(
    const std::vector<std::string>& date) {
  std::vector<double> rez;
  for (int i = 0; i < (int)date.size(); ++i) {
    double tmp = convertDateTime(date[i]);
    rez.push_back(tmp);
  }
  double tmp = rez[0];
  for (int i = 0; i < (int)rez.size(); ++i) {
    rez[i] = 1 + (rez[i] - tmp) / 86400;
  }
  return rez;
}

}  // namespace s21
