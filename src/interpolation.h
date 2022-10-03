#pragma once

#include <QColor>
#include <QString>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

class Interpolation {
 public:
  double GetValueNewtonPolynome(const std::vector<double> &x,
                                const std::vector<double> &data,
                                const double &current_x);
  double GetValueCubicSpline(const std::vector<double> &x,
                             const std::vector<double> &data,
                             const double &current_x);

 private:
  double InterpolateValueNewtonPolynome(const std::vector<double> &x,
                                        const std::vector<double> &b,
                                        const double &current_x);
  std::vector<double> GetBVector(const std::vector<double> &x_vec,
                                 const std::vector<double> &data);
  double GetInternal(const int p, const int q,
                     const std::vector<double> &x_vec);
  double GetInternal2(const int p, const double inner_x,
                      const std::vector<double> &x_vec);
  std::vector<double> GetDividers(const std::vector<double> &x_vec);

  // Cubic Interpolation
  // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
  class Spline {
   public:
    Spline(const std::vector<double> &X, const std::vector<double> &Y) {
      this->set_points(X, Y);
    }
    Spline() = delete;
    Spline(const Spline &other) = delete;
    Spline &operator=(const Spline &other) = delete;
    Spline(Spline &&other) = delete;
    Spline &operator=(Spline &&other) = delete;
    ~Spline(){};
    void set_points(const std::vector<double> &x, const std::vector<double> &y);
    double operator()(double x) const;
    std::vector<double> get_x() const { return m_x; }
    std::vector<double> get_y() const { return m_y; }
    double get_x_min() const;
    double get_x_max() const;

   private:
    std::vector<double> m_x, m_y;
    std::vector<double> m_b, m_c, m_d;
    size_t find_closest(double x) const;
    double deriv(double x) const;
    std::vector<double> solve(double y) const;
  };

  class BandMatrix {
   public:
    BandMatrix() = delete;
    BandMatrix(int dim, int n_u, int n_l);
    BandMatrix(const BandMatrix &other) = delete;
    BandMatrix &operator=(const BandMatrix &other) = delete;
    BandMatrix(BandMatrix &&other) = delete;
    BandMatrix &operator=(BandMatrix &&other) = delete;
    ~BandMatrix(){};
    void resize(int dim, int n_u, int n_l);
    int dim() const;
    int num_upper() const { return (int)m_upper.size() - 1; }
    int num_lower() const { return (int)m_lower.size() - 1; }
    double &operator()(int i, int j);
    double operator()(int i, int j) const;
    double &saved_diag(int i);
    double saved_diag(int i) const;
    void lu_decompose();
    std::vector<double> r_solve(const std::vector<double> &b) const;
    std::vector<double> l_solve(const std::vector<double> &b) const;
    std::vector<double> lu_solve(const std::vector<double> &b,
                                 bool is_lu_decomposed = false);

   private:
    std::vector<std::vector<double> > m_upper;
    std::vector<std::vector<double> > m_lower;
  };
};
// Common methods
double get_eps();
std::vector<double> solve_cubic(double a, double b, double c, double d,
                                int newton_iter = 0);
