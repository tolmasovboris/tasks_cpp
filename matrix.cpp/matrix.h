#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

class BigInteger {
public:
  static const long long chunk_len = 9;
  static const long long chunk_count = 1e9;
  size_t number_length(std::string str) {
    size_t sz = str.size();
    if (str[0] == '-')
      --sz;
    return (sz + chunk_len - 1) / chunk_len;
  }
  size_t number_length(long long x) {
    if (x == 0) {
      return 1;
    }
    return std::ceil(std::log(std::abs(x) + 1) / std::log(chunk_count));
  }
  static std::pair<BigInteger, BigInteger> Division(const BigInteger &n,
                                                    const BigInteger &d) {
    if (d.is_negative) {
      std::pair<BigInteger, BigInteger> ans = Division(n, -d);
      return {-ans.first, ans.second};
    }
    if (n.is_negative) {
      std::pair<BigInteger, BigInteger> ans = Division(-n, d);
      return {-ans.first, -ans.second};
    }
    if (n < d) {
      return {0, n};
    }
    BigInteger ceill = 0;
    BigInteger ost = 0;
    for (size_t i = n.numbers.size(); i > 0; i--) {
      ost.push_top();
      ost += n.numbers[i - 1];
      if (ost < d) {
        ceill.push_top();
        continue;
      }
      long long j = BaseDivider(ost, d);
      ost -= d * j;
      ceill.push_top();
      ceill += j;
    }
    return {ceill, ost};
  }
  static BigInteger gcd(const BigInteger &p, const BigInteger &q) {
    if (q == 0) {
      return p;
    }
    return gcd(q, p % q);
  }

  BigInteger() : is_negative(false) {}
  BigInteger(int x) : is_negative(x < 0), numbers(number_length(x), 0) {
    x = std::abs(x);
    for (size_t i = 0; i < numbers.size(); ++i) {
      numbers[i] = x % chunk_count;
      x /= chunk_count;
    }
  }
  BigInteger(const BigInteger &other)
      : is_negative(other.is_negative), numbers(other.numbers) {}
  BigInteger(std::string str)
      : is_negative(str[0] == '-'), numbers(number_length(str), 0) {
    size_t k = 0;
    if (str == "-0") {
      numbers[0] = 0;
      return;
    }
    if (is_negative) {
      str.erase(0, 1);
    }
    int start = str.size();
    int end = 0;
    for (int i = start; i > end; i -= chunk_len) {
      if (i >= chunk_len) {

        numbers[k++] = stoi(str.substr(i - chunk_len, chunk_len));

      } else {

        numbers[k++] = stoi(str.substr(end, i - end));
      }
    }
  }
  std::string toString(bool zeros = true) const {
    std::stringstream ans;
    if (is_negative) {
      ans << '-';
    }
    auto iter = numbers.rbegin();
    ans << std::setfill('0');
    if (!zeros) {
      ans << std::setw(chunk_len);
    }
    ans << *iter++;

    for (; iter != numbers.rend(); ++iter) {
      ans << std::setw(chunk_len) << *iter;
    }
    return ans.str();
  }
  int compare(const BigInteger &other) const {
    if (is_negative != other.is_negative) {
      return is_negative ? -1 : 1;
    }
    size_t n = numbers.size();
    size_t m = other.numbers.size();
    if (n != m) {
      return (is_negative != (n < m)) != 0 ? -1 : 1;
    }
    for (size_t i = n; i > 0; i--) {
      if (numbers[i - 1] != other.numbers[i - 1]) {
        return (is_negative != (numbers[i - 1] < other.numbers[i - 1])) != 0
                   ? -1
                   : 1;
      }
    }
    return 0;
  }
  void push_top() {
    if (*this != 0) {
      numbers.insert(numbers.begin(), 0);
    }
  }

  BigInteger &operator=(BigInteger other) {
    swap(other);
    return *this;
  }
  BigInteger &operator+=(const BigInteger &other) {
    if (is_negative != other.is_negative) {
      if (is_negative) {
        return *this = other - (-*this);
      }
      return *this -= -other;
    }

    long long carry = 0;
    size_t n = numbers.size();
    size_t m = other.numbers.size();
    size_t max_size = std::max(n, m);
    for (size_t i = 0; i < max_size || carry != 0; ++i) {
      long long x = carry;
      if (i < n) {
        x += numbers[i];
      }
      if (i < m) {
        x += other.numbers[i];
      }
      carry = x >= chunk_count ? 1 : 0;
      if (carry != 0) {
        x -= chunk_count;
      }

      if (i >= n) {
        numbers.push_back(0);
      }
      numbers[i] = x;
    }

    return *this;
  }
  BigInteger &operator-=(const BigInteger &other) {
    if (other.is_negative) {
      return *this += -other;
    }
    if (*this == other) {
      return *this = 0;
    }
    if (*this < other) {
      *this = other - *this;
      is_negative = true;
      return *this;
    }

    size_t n = numbers.size();
    size_t m = other.numbers.size();
    long long carry = 0;
    for (size_t i = 0; i < n; ++i) {
      long long x = -carry;
      if (i < n) {
        x += numbers[i];
      }
      if (i < m) {
        x -= other.numbers[i];
      }
      carry = x < 0 ? 1 : 0;
      if (carry != 0) {
        x += chunk_count;
      }
      numbers[i] = x;
    }

    delete_zeros();
    return *this;
  }
  BigInteger &operator*=(const BigInteger &other) {
    size_t n = numbers.size();
    size_t m = other.numbers.size();
    if (n < m) {
      return *this = other * (*this);
    }

    std::vector<int> ans(n + m);
    for (size_t i = 0; i < n; ++i) {
      long long carry = 0;
      for (size_t j = 0; j < m || carry != 0; ++j) {
        long long x = ans[i + j] + carry;
        if (j < m) {
          x += numbers[i] * other.numbers[j];
        }
        ans[i + j] = x % chunk_count;
        carry = x / chunk_count;
      }
    }

    numbers = ans;
    delete_zeros();
    is_negative = numbers.back() != 0 && (is_negative != other.is_negative);
    return *this;
  }

  BigInteger &operator/=(const BigInteger &other) {
    return *this = Division(*this, other).first;
  }

  BigInteger &operator%=(const BigInteger &other) {
    return *this = Division(*this, other).second;
  }
  BigInteger operator-() const {
    BigInteger ans(*this);
    if (!ans.is_negative && ans.toString() != "0") {
      ans.is_negative = true;
    } else if (ans.is_negative && ans.toString() != "0") {
      ans.is_negative = false;
    }
    return ans;
  }
  BigInteger &operator++() { return *this += 1; }

  BigInteger operator++(int) {
    BigInteger copy = *this;
    ++*this;
    return copy;
  }

  BigInteger &operator--() { return *this -= 1; }

  BigInteger operator--(int) {
    BigInteger copy = *this;
    --*this;
    return copy;
  }

  explicit operator bool() const { return numbers.back() != 0; }

  friend BigInteger operator+(const BigInteger &lhs, const BigInteger &rhs) {
    BigInteger ans(lhs);
    ans += rhs;
    return ans;
  }
  friend BigInteger operator-(const BigInteger &lhs, const BigInteger &rhs) {
    BigInteger ans(lhs);
    ans -= rhs;
    return ans;
  }

  friend BigInteger operator*(const BigInteger &lhs, const BigInteger &rhs) {
    BigInteger ans(lhs);
    ans *= rhs;
    return ans;
  }

  friend BigInteger operator/(const BigInteger &lhs, const BigInteger &rhs) {
    BigInteger ans(lhs);
    ans /= rhs;
    return ans;
  }
  friend BigInteger operator%(const BigInteger &lhs, const BigInteger &rhs) {
    BigInteger ans(lhs);
    ans %= rhs;
    return ans;
  }
  bool operator==(const BigInteger &other) const {

    return this->compare(other) == 0;
  }

  bool operator!=(const BigInteger &other) const { return !(*this == other); }
  bool operator<(const BigInteger &other) const {
    return this->compare(other) == -1;
  }
  bool operator>(const BigInteger &other) const { return other < *this; }

  bool operator<=(const BigInteger &other) const { return !(other < *this); }

  bool operator>=(const BigInteger &other) const { return !(*this < other); }

private:
  bool is_negative = false;
  std::vector<int> numbers;
  static long long BaseDivider(const BigInteger &n, const BigInteger &d) {
    long long left = 1;
    long long right = chunk_count;
    while (right - left > 1) {
      long long mid = left + (right - left) / 2;
      (mid * d <= n ? left : right) = mid;
    }
    return left;
  }

  void swap(BigInteger &other) {
    std::swap(is_negative, other.is_negative);
    std::swap(numbers, other.numbers);
  }
  void delete_zeros() {
    while (numbers.size() > 1 && numbers.back() == 0) {
      numbers.pop_back();
    }
  }
};

BigInteger operator""_bi(unsigned long long x) {

  return BigInteger(std::to_string(x));
}
BigInteger operator""_bi(const char *str, size_t) {

  std::string str_ = std::string(str);
  if (str_ == "-0")
    str_ = "0";
  return BigInteger(str_);
}
std::ostream &operator<<(std::ostream &os, const BigInteger &other) {

  return os << other.toString();
}
std::istream &operator>>(std::istream &in, BigInteger &other) {
  std::string str;
  in >> str;
  if (str == "-0")
    str = "0";

  other = str;

  return in;
}
class Rational {
public:
  Rational() : numerator(0), denominator(1) {}
  Rational(const Rational &other)
      : numerator(other.numerator), denominator(other.denominator) {
    reduce();
  }
  Rational(const BigInteger &other) : numerator(other), denominator(1) {}
  Rational(int other) : numerator(other), denominator(1) {}
  Rational(const BigInteger &numerator, const BigInteger &denominator)
      : numerator(numerator), denominator(denominator) {
    reduce();
  }
  Rational &operator=(Rational other) {
    swap(other);
    return *this;
  }

  std::string toString() const {

    if (denominator == 1 || numerator == 0)
      return numerator.toString();
    return numerator.toString() + "/" + denominator.toString();
  }
  std::string asDecimal(size_t precision) const {
    std::pair<BigInteger, BigInteger> div =
        BigInteger::Division(numerator, denominator);
    std::string ceil;
    if (numerator < 0) {
      div.second = -div.second;
      ceil += "-";
    }
    ceil += div.first.toString();
    if (precision == 0) {
      return ceil;
    }

    std::string floor;
    std::size_t end =
        (precision + BigInteger::chunk_len - 1) / BigInteger::chunk_len;
    for (size_t i = 0; i < end; ++i) {
      div.second *= BigInteger::chunk_count;
      div = BigInteger::Division(div.second, denominator);
      floor += div.first.toString(false);
    }
    floor = floor.substr(0, precision);
    return ceil + "." + floor;
  }
  Rational &operator+=(const Rational &other) {
    numerator = numerator * other.denominator + other.numerator * denominator;
    denominator *= other.denominator;
    //reduce();
    return *this;
  }
  Rational &operator-=(const Rational &other) {
    numerator = numerator * other.denominator - other.numerator * denominator;
    denominator *= other.denominator;
    reduce();
    return *this;
  }
  Rational &operator*=(const Rational &other) {
    numerator *= other.numerator;
    denominator *= other.denominator;
    reduce();
    return *this;
  }

  Rational &operator/=(const Rational &other) {
    numerator *= other.denominator;
    denominator *= other.numerator;
    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }
    reduce();
    return *this;
  }
  Rational operator-() { return Rational(-numerator, denominator); }
  friend Rational operator+(const Rational &lhs, const Rational &rhs) {
    Rational ans(lhs);
    ans += rhs;
    return ans;
  }
  friend Rational operator-(const Rational &lhs, const Rational &rhs) {
    Rational ans(lhs);
    ans -= rhs;
    return ans;
  }
  friend Rational operator*(const Rational &lhs, const Rational &rhs) {
    Rational ans(lhs);
    ans *= rhs;
    return ans;
  }
  friend Rational operator/(const Rational &lhs, const Rational &rhs) {
    Rational ans(lhs);
    ans /= rhs;
    return ans;
  }
  bool operator==(const Rational &other) const {
    return (numerator == other.numerator && denominator == other.denominator);
  }
  bool operator!=(const Rational &other) const { return !(*this == other); }
  bool operator<(const Rational &other) const {
    return numerator * other.denominator < denominator * other.numerator;
  }
  bool operator>(const Rational &other) const { return (other < *this); }
  bool operator>=(const Rational &other) const { return !(*this < other); }
  bool operator<=(const Rational &other) const { return !(other < *this); }
  explicit operator double() const { return std::stod(asDecimal(16)); }

private:
  BigInteger numerator;
  BigInteger denominator;
  void swap(Rational &other) {
    std::swap(numerator, other.numerator);
    std::swap(denominator, other.denominator);
  }
  void reduce() {
    if (numerator == 0) {
      denominator = 1;
      return;
    }
    BigInteger gcd =
        BigInteger::gcd(numerator < 0 ? -numerator : numerator, denominator);
    numerator /= gcd;
    denominator /= gcd;
  }
};
std::istream &operator>>(std::istream &in, Rational &other) {
  BigInteger bigint;
    in >> bigint;
    other = Rational(bigint);
    return in;
}
//Проверка на простоту
template <size_t N, size_t M> struct IsDivisable {
  static const constexpr bool value = N % M == 0;
};

template <size_t N, size_t M>
const constexpr bool IsDivisableV = IsDivisable<N, M>::value;

template <size_t N, size_t M> struct LessSqrt {
  static const constexpr bool value = M * M < N;
};

template <size_t N, size_t M> const bool LessSqrtV = LessSqrt<N, M>::value;

template <size_t N, size_t M, bool Condition> struct is_prime_v {
  static const constexpr bool value =
      (IsDivisableV<N, M> || is_prime_v<N, M + 1, LessSqrtV<N, M + 1>>::value);
};

template <size_t N, size_t M> struct is_prime_v<N, M, false> {
  static const constexpr bool value = false;
};
template <size_t N, size_t M>
const bool is_prime_q = is_prime_v<N, M, LessSqrtV<N, M>>::value;
template <size_t N> struct is_prime {
  static const constexpr bool value = !is_prime_q<N, 2>;
};

template <size_t N> class Residue {
public:
  explicit Residue(int value)
      : value(((value % static_cast<int>(N)) + N) % N) {}
  Residue() : value(0) {}
  explicit operator int() { return static_cast<int>(value); }
  Residue<N> &operator+=(const Residue<N> &other) {
    value += other.value;
    value %= N;
    return *this;
  }
  Residue<N> &operator*=(const Residue<N> &other) {
    value *= other.value;
    value %= N;
    return *this;
  }
  Residue<N> &operator-=(const Residue<N> &other) {
    value -= other.value - N;
    value %= N;
    return *this;
  }
  Residue<N> &operator/=(const Residue<N> &other) {
    static_assert(is_prime<N>::value, "IS NOT PRIME!!!");
    Residue<N> inverse(1);
    for (size_t i = 0; i < N - 2; ++i) {
        inverse *= other;
    }
    value *= inverse.value;
    value %= N;
    return *this;
  }
  size_t GetValue() { return value; }
  size_t GetValue() const { return value; }
  bool operator==(const Residue<N> &other) const{
    return value == other.value;
  }
  bool operator!=(const Residue<N> &other) const{
    return value != other.value;
  }
  bool operator==(int other) const {
        return value == ((other % N + N) % N);
    }

    bool operator!=(int other) const {
        return value != ((other % N + N) % N);
    }
private:
  int value;
};
template <size_t N>
Residue<N> operator+(const Residue<N> &lhs, const Residue<N> &rhs) {
  Residue ans(lhs);
  ans += rhs;
  return ans;
}
template <size_t N>
Residue<N> operator-(const Residue<N> &lhs, const Residue<N> &rhs) {
  Residue ans(lhs);
  ans -= rhs;
  return ans;
}
template <size_t N>
Residue<N> operator*(const Residue<N> &lhs, const Residue<N> &rhs) {
  Residue ans(lhs);
  ans *= rhs;
  return ans;
}
template <size_t N>
Residue<N> operator/(const Residue<N> &lhs, const Residue<N> &rhs) {
  Residue ans(lhs);
  ans /= rhs;
  return ans;
}

template <size_t N>
std::ostream &operator<<(std::ostream &os, const Residue<N> &elem) {
  return os << elem.GetValue();
}
template <size_t M, size_t N, typename Field = Rational> class Matrix {
public:
  std::array<Field, N> getRow(size_t row) const { return table[row]; }
  std::array<Field, M> getColumn(size_t col) const {
    std::array<Field, M> ans;
    for (size_t i = 0; i < col; ++i) {
      ans[i] = table[col][i];
    }
    return ans;
  }
  Matrix() : table(std::array<std::array<Field, N>, M>()) {}
  Matrix(std::initializer_list<std::initializer_list<int>> matrix) {
    size_t i = 0;
    for (const auto &row : matrix) {
      size_t j = 0;
      for (const int &elem : row) {
        table[i][j++] = Field(elem);
      }
      ++i;
    }
  }
  const Matrix<N, M, Field> transposed() const {
    Matrix<N, M, Field> ans;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        ans[j][i] = table[i][j];
      }
    }
    return ans;
  }

  Field trace() {
    static_assert(M == N, "Differance in size");
    Field ans = Field(0);
    for (size_t i = 0; i < N; ++i) {
      ans += table[i][i];
    }
    return ans;
  }
  Matrix<M, N, Field> Gauss() const {
    Matrix<M, N, Field> gauss_matrix = *this;
    size_t row_now = 0;
    for (size_t col = 0; col < N; ++col) {
      bool permut = false;
      for (size_t row = row_now; row < M; ++row) {
        if (gauss_matrix[row][col] != Field(0)) {
          std::swap(gauss_matrix[row_now], gauss_matrix[row]);
          permut = true;
          if (row_now != row) {
            gauss_matrix.sign_permutation *= -1;
          }
          break;
        }
      }
      if (!permut) {
        continue;
      }
      Field main_elem = gauss_matrix[row_now][col];
      for (size_t row = row_now + 1; row < M; ++row) {
        Field cf = gauss_matrix[row][col] / main_elem;
        if (cf == Field(0)) {
          continue;
        }
        for (size_t i = 0; i < N; ++i) {
          gauss_matrix[row][i] -= cf * gauss_matrix[row_now][i];
        }
      }
      ++row_now;
    }
    return gauss_matrix;
  }
  Field det() const {
    static_assert(N == M, "not sqrted");
    Matrix<M, N, Field> gauss = Gauss();
    
    Field ans = Field(gauss.sign_permutation);
    for (size_t i = 0; i < M; ++i) {
      ans *= gauss[i][i];
    }
    return ans;
  }
  size_t rank() const {
    Matrix<M, N, Field> gauss = Gauss();
    size_t rg = M;
    for (size_t i = 0; i < M; ++i) {
      bool zeros = true;
      for (size_t j = 0; j < N; ++j) {
        if (gauss[i][j] != 0) {
          zeros = false;
          break;
        }
      }
      if (zeros) {
        rg -= 1;
      }
    }
    return rg;
  }
  Matrix<M, N, Field> inverted() {
    static_assert(N == M, "Differance in size");
    Matrix<N, N, Field> invertion;
    for (size_t i = 0; i < M; ++i) {
      invertion[i][i] = Field(1);
    }
     

    Matrix<N, N, Field> gauss = *this;
    size_t row_now = 0;
    for (size_t col = 0; col < M; ++col) {
      bool permut = false;
      for (size_t row = row_now; row < M; ++row) {
        if (gauss[row][col] != Field(0)) {
          std::swap(gauss[row_now], gauss[row]);
          std::swap(invertion[row_now], invertion[row]);
          permut = true;
          break;
        }
      }
      if (!permut) {
        continue;
      }
      Field main_elem = gauss[row_now][col];
      for (size_t row = row_now + 1; row < M; ++row) {
        Field cf = gauss[row][col] / main_elem;
        if (cf == Field(0)) {
          continue;
        }
        for (size_t i = 0; i < M; ++i) {
          gauss[row][i] -= cf * gauss[row_now][i];
          invertion[row][i] -= cf * invertion[row_now][i];
        }
      }
      ++row_now;
    }

    for (size_t i = M; i-- > 0; ) {
      Field main_elem = gauss[i][i];
      for (size_t row = 0; row < i; ++row) {
        Field cf = gauss[row][i] / main_elem;
        if (cf == Field(0)) {
          continue;
        }
        for (size_t j = 0; j < M; ++j) {
          gauss[row][j] -= cf * gauss[i][j];
          invertion[row][j] -= cf * invertion[i][j];
        }
      }
      for (size_t j = 0; j < M; ++j) {
        invertion[i][j] /= gauss[i][i];
      }
    }

    return invertion;
  }
  void invert() { *this = this->inverted(); }
  Matrix<M, N, Field> &operator+=(const Matrix<M, N, Field> &other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        table[i][j] += other[i][j];
      }
    }
    return *this;
  }
  Matrix<M, N, Field> &operator-=(const Matrix<M, N, Field> &other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        table[i][j] -= other[i][j];
      }
    }
    return *this;
  }
  Matrix<M, N, Field> &operator*=(const Field &elem) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        table[i][j] *= elem;
      }
    }
    return *this;
  }
  Matrix<M, N, Field> &operator*=(const Matrix<M, N, Field> &other) {
    static_assert(N == M, "Differance in size");
    *this = *this * other;
    return *this;
  }
  std::array<Field, N> &operator[](size_t row) { return table[row]; }
  const std::array<Field, N> &operator[](size_t row) const {
    return table[row];
  }

private:
  std::array<std::array<Field, N>, M> table;
  int sign_permutation = 1;
};
template <size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;
template <size_t M, size_t N, typename Field>
bool operator==(const Matrix<M, N, Field> &lhs,
                const Matrix<M, N, Field> &rhs) {
  bool equal = true;
  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < N; ++j) {
      if (lhs[i][j] != rhs[i][j]) {
        equal = false;
      }
    }
  }
  return equal;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> &lhs,
                              const Matrix<M, N, Field> &rhs) {
  Matrix<M, N, Field> ans(lhs);
  ans += rhs;
  return ans;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> &lhs,
                              const Matrix<M, N, Field> &rhs) {
  Matrix<M, N, Field> ans(lhs);
  ans -= rhs;
  return ans;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(Matrix<M, N, Field> lhs, const Field &rhs) {
  return lhs *= rhs;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field &lhs, Matrix<M, N, Field> rhs) {
  return rhs *= lhs;
}
template <size_t M, size_t N, size_t K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field> &lhs,
                              const Matrix<N, K, Field> &rhs) {
  Matrix<M, K, Field> res;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t k = 0; k < N; ++k) {
        res[i][j] += lhs[i][k] * rhs[k][j];
      }
    }
  }
  return res;
}
