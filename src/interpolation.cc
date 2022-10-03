#include "interpolation.h"

// Newton polynome
std::vector<double> Interpolation::GetDividers(
    const std::vector<double> &x_vec) {
  std::vector<double> divider(x_vec.size(), 1);
  for (size_t i = 0; i < x_vec.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      divider[i] *= x_vec[i] - x_vec[j];
    }
  }
  return divider;
}

double Interpolation::GetInternal(const int p, const int q,
                                  const std::vector<double> &x_vec) {
  double result = 1;
  for (int i = 0; i < p; ++i) {
    result *= x_vec[q] - x_vec[i];
  }
  return result;
}

double Interpolation::GetInternal2(const int p, const double inner_x,
                                   const std::vector<double> &x_vec) {
  double result = 1;
  for (int i = 0; i < p; ++i) {
    result *= inner_x - x_vec[i];
  }
  return result;
}

std::vector<double> Interpolation::GetBVector(const std::vector<double> &x_vec,
                                              const std::vector<double> &data) {
  std::vector<double> divider = GetDividers(x_vec);
  std::vector<double> b(x_vec.size());
  b[0] = data[0];
  for (size_t i = 0; i < x_vec.size(); ++i) {
    std::vector<double> tmp(x_vec.size());
    tmp[0] = 1;
    double ttt = 0;
    for (size_t j = 0; j < i; ++j) {
      tmp[j] = GetInternal(j, i, x_vec);
    }
    for (size_t j = 0; j < i; ++j) {
      ttt += b[j] * tmp[j];
    }
    b[i] = (data[i] - ttt) / divider[i];
  }
  return b;
}

double Interpolation::InterpolateValueNewtonPolynome(
    const std::vector<double> &x, const std::vector<double> &b,
    const double &current_x) {
  double value = 0;
  for (size_t j = 0; j < x.size(); j++) {
    value += b[j] * GetInternal2(j, current_x, x);
  }
  return value;
}

double Interpolation::GetValueNewtonPolynome(const std::vector<double> &x,
                                             const std::vector<double> &data,
                                             const double &current_x) {
  auto max = std::max_element(x.begin(), x.end());
  auto min = std::min_element(x.begin(), x.end());
  if (current_x > *max || current_x < *min) {
    throw std::invalid_argument(
        "ERROR: in GetValueNewtonPolynome x not inside the data");
  }
  std::vector<double> b = GetBVector(x, data);
  double value = InterpolateValueNewtonPolynome(x, b, current_x);
  return value;
}

// Cubic Polynome
double Interpolation::GetValueCubicSpline(const std::vector<double> &x,
                                          const std::vector<double> &data,
                                          const double &current_x) {
  Spline s(x, data);
  double result = s(current_x);
  return result;
}

void Interpolation::Spline::set_points(const std::vector<double> &x,
                                       const std::vector<double> &y) {
  if (x.size() != y.size()) {
    throw std::out_of_range(
        "ERROR: number of x values must be equal number of y values");
  }
  if (x.size() < 3) {
    throw std::invalid_argument(
        "At least 3 points needed to calculate cubic spline");
  }
  m_x = x;
  m_y = y;
  int n = (int)x.size();
  for (int i = 0; i < n - 1; i++) {
    if (m_x[i] >= m_x[i + 1]) {
      throw std::invalid_argument("Dataset not sorted, x values not in order");
    }
  }

  BandMatrix A(n, 1, 1);
  std::vector<double> rhs(n);
  for (int i = 1; i < n - 1; i++) {
    A(i, i - 1) = 1.0 / 3.0 * (x[i] - x[i - 1]);
    A(i, i) = 2.0 / 3.0 * (x[i + 1] - x[i - 1]);
    A(i, i + 1) = 1.0 / 3.0 * (x[i + 1] - x[i]);
    rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
             (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
  }
  // boundary conditions
  A(0, 0) = 2.0;
  A(0, 1) = 0.0;
  rhs[0] = 0;

  A(n - 1, n - 1) = 2.0;
  A(n - 1, n - 2) = 0.0;
  rhs[n - 1] = 0;

  // solve the equation system to obtain the parameters c[]
  m_c = A.lu_solve(rhs);

  m_d.resize(n);
  m_b.resize(n);
  for (int i = 0; i < n - 1; i++) {
    m_d[i] = 1.0 / 3.0 * (m_c[i + 1] - m_c[i]) / (x[i + 1] - x[i]);
    m_b[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
             1.0 / 3.0 * (2.0 * m_c[i] + m_c[i + 1]) * (x[i + 1] - x[i]);
  }
  double h = x[n - 1] - x[n - 2];
  m_d[n - 1] = 0.0;
  m_b[n - 1] = 3.0 * m_d[n - 2] * h * h + 2.0 * m_c[n - 2] * h + m_b[n - 2];
}

size_t Interpolation::Spline::find_closest(double x) const {
  std::vector<double>::const_iterator it;
  it = std::upper_bound(m_x.begin(), m_x.end(), x);     // *it > x
  size_t idx = std::max(int(it - m_x.begin()) - 1, 0);  // m_x[idx] <= x
  return idx;
}

double Interpolation::Spline::operator()(double x) const {
  size_t idx = find_closest(x);
  double h = x - m_x[idx];
  double interpol = ((m_d[idx] * h + m_c[idx]) * h + m_b[idx]) * h + m_y[idx];
  return interpol;
}

double Interpolation::Spline::deriv(double x) const {
  size_t idx = find_closest(x);
  double h = x - m_x[idx];
  double interpol;
  interpol = (3.0 * m_d[idx] * h + 2.0 * m_c[idx]) * h + m_b[idx];
  return interpol;
}

std::vector<double> Interpolation::Spline::solve(double y) const {
  std::vector<double> x;
  std::vector<double> root;
  const size_t n = m_x.size();
  for (size_t i = 0; i < n - 1; i++) {
    root = solve_cubic(m_y[i] - y, m_b[i], m_c[i], m_d[i], 1);
    for (size_t j = 0; j < root.size(); j++) {
      double h = (i > 0) ? (m_x[i] - m_x[i - 1]) : 0.0;
      double eps = get_eps() * 512.0 * std::min(h, 1.0);
      if ((-eps <= root[j]) && (root[j] < m_x[i + 1] - m_x[i])) {
        double new_root = m_x[i] + root[j];
        if (x.size() > 0 && x.back() + eps > new_root) {
          x.back() = new_root;
        } else {
          x.push_back(new_root);
        }
      }
    }
  }
  return x;
}

double Interpolation::Spline::get_x_min() const {
  if (m_x.empty()) {
    throw std::invalid_argument("ERROR: no x-values");
  }
  return m_x.front();
}
double Interpolation::Spline::get_x_max() const {
  if (m_x.empty()) {
    throw std::invalid_argument("ERROR: no x-values");
  }
  return m_x.back();
}

Interpolation::BandMatrix::BandMatrix(int dim, int n_u, int n_l) {
  resize(dim, n_u, n_l);
}
void Interpolation::BandMatrix::resize(int dim, int n_u, int n_l) {
  if (dim <= 0 && n_u < 0 && n_l < 0) {
    throw std::out_of_range("Matrices cannot be less than 0");
  }
  m_upper.resize(n_u + 1);
  m_lower.resize(n_l + 1);
  for (size_t i = 0; i < m_upper.size(); i++) {
    m_upper[i].resize(dim);
  }
  for (size_t i = 0; i < m_lower.size(); i++) {
    m_lower[i].resize(dim);
  }
}
int Interpolation::BandMatrix::dim() const {
  if (m_upper.size() > 0) {
    return m_upper[0].size();
  } else {
    return 0;
  }
}

double &Interpolation::BandMatrix::operator()(int i, int j) {
  int k = j - i;
  if ((i < 0) && (i >= dim()) && (j < 0) && (j >= dim())) {
    throw std::out_of_range("Index out of range");
  }
  if (k >= 0)
    return m_upper[k][i];
  else
    return m_lower[-k][i];
}
double Interpolation::BandMatrix::operator()(int i, int j) const {
  int k = j - i;
  if ((i < 0) && (i >= dim()) && (j < 0) && (j >= dim())) {
    throw std::out_of_range("Index out of range");
  }
  if (k >= 0)
    return m_upper[k][i];
  else
    return m_lower[-k][i];
}

double Interpolation::BandMatrix::saved_diag(int i) const {
  if ((i < 0) && (i >= dim())) {
    throw std::out_of_range("Index out of range");
  }
  return m_lower[0][i];
}
double &Interpolation::BandMatrix::saved_diag(int i) {
  if ((i < 0) && (i >= dim())) {
    throw std::out_of_range("Index out of range");
  }
  return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void Interpolation::BandMatrix::lu_decompose() {
  int i_max, j_max;
  int j_min;
  double x;

  for (int i = 0; i < this->dim(); i++) {
    if (this->operator()(i, i) == 0.0) {
      std::invalid_argument("Division by zero");
    }
    this->saved_diag(i) = 1.0 / this->operator()(i, i);
    j_min = std::max(0, i - this->num_lower());
    j_max = std::min(this->dim() - 1, i + this->num_upper());
    for (int j = j_min; j <= j_max; j++) {
      this->operator()(i, j) *= this->saved_diag(i);
    }
    this->operator()(i, i) = 1.0;  // prevents rounding errors
  }

  // Gauss LR-Decomposition
  for (int k = 0; k < this->dim(); k++) {
    i_max = std::min(this->dim() - 1,
                     k + this->num_lower());  // num_lower not a mistake!
    for (int i = k + 1; i <= i_max; i++) {
      if (this->operator()(i, i) == 0.0) {
        std::invalid_argument("Dividion by zero");
      }
      x = -this->operator()(i, k) / this->operator()(k, k);
      this->operator()(i, k) = -x;  // assembly part of L
      j_max = std::min(this->dim() - 1, k + this->num_upper());
      for (int j = k + 1; j <= j_max; j++) {
        // assembly part of R
        this->operator()(i, j) =
            this->operator()(i, j) + x * this->operator()(k, j);
      }
    }
  }
}
// solves Ly=b
std::vector<double> Interpolation::BandMatrix::l_solve(
    const std::vector<double> &b) const {
  if (this->dim() != (int)b.size()) {
    throw std::invalid_argument("Invalid vector b");
  }
  std::vector<double> x(this->dim());
  int j_start;
  double sum;
  for (int i = 0; i < this->dim(); i++) {
    sum = 0;
    j_start = std::max(0, i - this->num_lower());
    for (int j = j_start; j < i; j++) sum += this->operator()(i, j) * x[j];
    x[i] = (b[i] * this->saved_diag(i)) - sum;
  }
  return x;
}
// solves Rx=y
std::vector<double> Interpolation::BandMatrix::r_solve(
    const std::vector<double> &b) const {
  if (this->dim() != (int)b.size()) {
    throw std::invalid_argument("Invalid vector b");
  }
  std::vector<double> x(this->dim());
  int j_stop;
  double sum;
  for (int i = this->dim() - 1; i >= 0; i--) {
    sum = 0;
    j_stop = std::min(this->dim() - 1, i + this->num_upper());
    for (int j = i + 1; j <= j_stop; j++) sum += this->operator()(i, j) * x[j];
    x[i] = (b[i] - sum) / this->operator()(i, i);
  }
  return x;
}

std::vector<double> Interpolation::BandMatrix::lu_solve(
    const std::vector<double> &b, bool is_lu_decomposed) {
  if (this->dim() != (int)b.size()) {
    throw std::invalid_argument("Invalid vector b");
  }
  std::vector<double> x, y;
  if (is_lu_decomposed == false) {
    this->lu_decompose();
  }
  y = this->l_solve(b);
  x = this->r_solve(y);
  return x;
}

// machine precision of a double, i.e. the successor of 1 is 1+eps
double get_eps() { return std::numeric_limits<double>::epsilon(); }

std::vector<double> solve_cubic(double a, double b, double c, double d,
                                int newton_iter) {
  if (d != 1.0) {
    a /= d;
    b /= d;
    c /= d;
  }

  std::vector<double> z;  // roots of the depressed cubic
  double p = -(1.0 / 3.0) * b + (1.0 / 9.0) * (c * c);
  double r = 2.0 * (c * c) - 9.0 * b;
  double q = -0.5 * a - (1.0 / 54.0) * (c * r);
  double discr = p * p * p - q * q;  // discriminant
  if (discr > 0) {
    z.resize(3);
    double ac = (1.0 / 3.0) * acos(q / (p * sqrt(p)));
    double sq = 2.0 * sqrt(p);
    z[0] = sq * cos(ac);
    z[1] = sq * cos(ac - 2.0 * M_PI / 3.0);
    z[2] = sq * cos(ac - 4.0 * M_PI / 3.0);
  } else if (discr < 0.0) {
    z.resize(1);
    double sgnq = (q >= 0 ? 1 : -1);
    double basis = fabs(q) + sqrt(-discr);
    double C = sgnq * pow(basis, 1.0 / 3.0);  // c++11 has std::cbrt()
    z[0] = C + p / C;
  }
  for (size_t i = 0; i < z.size(); i++) {
    z[i] -= (1.0 / 3.0) * c;
    for (int k = 0; k < newton_iter; k++) {
      double f = ((z[i] + c) * z[i] + b) * z[i] + a;
      double f1 = (3.0 * z[i] + 2.0 * c) * z[i] + b;
      if (fabs(f1) > 1e-8) {
        z[i] -= f / f1;
      }
    }
  }
  if (a == 0.0) {
    if (z.size() == 0) {
      throw std::invalid_argument("Roots not found");
    }
    double xmin = fabs(z[0]);
    size_t imin = 0;
    for (size_t i = 1; i < z.size(); i++) {
      if (xmin > fabs(z[i])) {
        xmin = fabs(z[i]);
        imin = i;
      }
    }
    z[imin] = 0.0;  // replace the smallest absolute value with 0
  }
  std::sort(z.begin(), z.end());
  return z;
}
