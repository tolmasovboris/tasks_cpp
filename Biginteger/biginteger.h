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
    if (str[0] == '-') --sz;
    return (sz + chunk_len - 1) / chunk_len;
  }
  size_t number_length(long long x) {
    if (x == 0) {
      return 1;
    }
    return std::ceil(std::log(std::abs(x) + 1) / std::log(chunk_count));
  }
  static std::pair<BigInteger, BigInteger> Division(const BigInteger& n,
                                                    const BigInteger& d) {
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
    BigInteger ceil = 0;
    BigInteger remainder = 0;
    for (size_t i = n.numbers.size(); i > 0; i--) {
      remainder.push_top();
      remainder += n.numbers[i - 1];
      if (remainder < d) {
        ceil.push_top();
        continue;
      }
      long long j = BaseDivider(remainder, d);
      BigInteger red = d;
      red *= j;
      remainder -= red;
      ceil.push_top();
      ceil += j;
    }
    return {ceil, remainder};
  }
  static BigInteger gcd(const BigInteger& p, const BigInteger& q) {
    if (q == 0) {
      return p;
    }
    BigInteger red = p;
    red %= q;
    return gcd(q, red);
  }

  BigInteger() : is_negative(false) {}
  BigInteger(int x) : is_negative(x < 0), numbers(number_length(x), 0) {
    x = std::abs(x);
    for (size_t i = 0; i < numbers.size(); ++i) {
      numbers[i] = x % chunk_count;
      x /= chunk_count;
    }
  }
  BigInteger(const BigInteger& other)
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
  int compare(const BigInteger& other) const {
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

  BigInteger& operator=(const BigInteger& other) {
    if (this == &other) {
       return *this;
    }
    BigInteger other_copy(other);
    swap(other_copy);
    return *this;
  }
  void add(std::vector<int>& lhs, const std::vector<int>& rhs) {
    long long carry = 0;
    size_t n = lhs.size();
    size_t m = rhs.size();
    size_t max_size = std::max(n, m);
    for (size_t i = 0; i < max_size || carry != 0; ++i) {
      long long x = carry;
      if (i < n) {
        x += lhs[i];
      }
      if (i < m) {
        x += rhs[i];
      }
      carry = x >= chunk_count ? 1 : 0;
      if (carry != 0) {
        x -= chunk_count;
      }

      if (i >= n) {
        lhs.push_back(0);
      }
      lhs[i] = x;
    }
  }
  BigInteger& operator+=(const BigInteger& other) {
    if (is_negative != other.is_negative) {
      if (is_negative) {
        BigInteger red = other;
        red -= -*this;
        return *this = red;
      }
      return *this -= -other;
    }
    add(numbers, other.numbers);
    return *this;
  }
  void substract(std::vector<int>& lhs, const std::vector<int>& rhs) {
    size_t n = lhs.size();
    size_t m = rhs.size();
    long long carry = 0;
    for (size_t i = 0; i < n; ++i) {
      long long x = -carry;
      if (i < n) {
        x += lhs[i];
      }
      if (i < m) {
        x -= rhs[i];
      }
      carry = x < 0 ? 1 : 0;
      if (carry != 0) {
        x += chunk_count;
      }
      lhs[i] = x;
    }
  }
  BigInteger& operator-=(const BigInteger& other) {
    if (other.is_negative) {
      return *this += -other;
    }
    if (*this < other) {
      BigInteger red = other;
      red -= *this;
      *this = red;
      is_negative = true;
      return *this;
    }
    substract(numbers, other.numbers);
    delete_zeros();
    return *this;
  }
  BigInteger& operator*=(const BigInteger& other) {
    size_t n = numbers.size();
    size_t m = other.numbers.size();
    if (n < m) {
      BigInteger red = other;
      red *= *this;
      return *this = red;
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

  BigInteger& operator/=(const BigInteger& other) {
    return *this = Division(*this, other).first;
  }

  BigInteger& operator%=(const BigInteger& other) {
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
  BigInteger& operator++() { return *this += 1; }

  BigInteger operator++(int) {
    BigInteger copy = *this;
    ++*this;
    return copy;
  }

  BigInteger& operator--() { return *this -= 1; }

  BigInteger operator--(int) {
    BigInteger copy = *this;
    --*this;
    return copy;
  }

  explicit operator bool() const { return numbers.back() != 0; }

  bool operator==(const BigInteger& other) const {
    return this->compare(other) == 0;
  }

  bool operator!=(const BigInteger& other) const { return !(*this == other); }
  bool operator<(const BigInteger& other) const {
    return this->compare(other) == -1;
  }
  bool operator>(const BigInteger& other) const { return other < *this; }

  bool operator<=(const BigInteger& other) const { return !(other < *this); }

  bool operator>=(const BigInteger& other) const { return !(*this < other); }

 private:
  bool is_negative = false;
  std::vector<int> numbers;
  static long long BaseDivider(const BigInteger& n, const BigInteger& d) {
    long long left = 1;
    long long right = chunk_count;
    while (right - left > 1) {
      long long mid = left + (right - left) / 2;
      BigInteger red = mid;
      red *= d;
      (red <= n ? left : right) = mid;
    }
    return left;
  }

  void swap(BigInteger& other) {
    std::swap(is_negative, other.is_negative);
    std::swap(numbers, other.numbers);
  }
  void delete_zeros() {
    while (numbers.size() > 1 && numbers.back() == 0) {
      numbers.pop_back();
    }
  }
};
BigInteger operator+(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger ans(lhs);
  ans += rhs;
  return ans;
}
BigInteger operator-(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger ans(lhs);
  ans -= rhs;
  return ans;
}
BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger ans(lhs);
  ans *= rhs;
  return ans;
}
BigInteger operator%(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger ans(lhs);
  ans %= rhs;
  return ans;
}
BigInteger operator/(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger ans(lhs);
  ans /= rhs;
  return ans;
}
BigInteger operator""_bi(unsigned long long x) {
  return BigInteger(std::to_string(x));
}
BigInteger operator""_bi(const char* str, size_t) {
  std::string str_ = std::string(str);
  if (str_ == "-0") str_ = "0";
  return BigInteger(str_);
}
std::ostream& operator<<(std::ostream& os, const BigInteger& other) {
  return os << other.toString();
}
std::istream& operator>>(std::istream& in, BigInteger& other) {
  std::string str;
  in >> str;
  if (str == "-0") str = "0";

  other = str;

  return in;
}
class Rational {
 public:
  Rational() : numerator(0), denominator(1) {}
  Rational(const Rational& other)
      : numerator(other.numerator), denominator(other.denominator) {
    reduce();
  }
  Rational(const BigInteger& other) : numerator(other), denominator(1) {}
  Rational(int other) : numerator(other), denominator(1) {}
  Rational(const BigInteger& numerator, const BigInteger& denominator)
      : numerator(numerator), denominator(denominator) {
    reduce();
  }
  Rational& operator=(const Rational& other) {
    if (other == *this) {
       return *this;
    }
    Rational other_copy(other);
    swap(other_copy);
    return *this;
  }

  std::string toString() const {
    if (denominator == 1 || numerator == 0) return numerator.toString();
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
  Rational& operator+=(const Rational& other) {
    numerator = numerator * other.denominator + other.numerator * denominator;
    denominator *= other.denominator;
    reduce();
    return *this;
  }
  Rational& operator-=(const Rational& other) {
    numerator = numerator * other.denominator - other.numerator * denominator;
    denominator *= other.denominator;
    reduce();
    return *this;
  }
  Rational& operator*=(const Rational& other) {
    numerator *= other.numerator;
    denominator *= other.denominator;
    reduce();
    return *this;
  }

  Rational& operator/=(const Rational& other) {
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
  bool operator==(const Rational& other) const {
    return (numerator == other.numerator && denominator == other.denominator);
  }
  bool operator!=(const Rational& other) const { return !(*this == other); }
  bool operator<(const Rational& other) const {
    return numerator * other.denominator < denominator * other.numerator;
  }
  bool operator>(const Rational& other) const { return (other < *this); }
  bool operator>=(const Rational& other) const { return !(*this < other); }
  bool operator<=(const Rational& other) const { return !(other < *this); }
  explicit operator double() const { return std::stod(asDecimal(16)); }

 private:
  BigInteger numerator;
  BigInteger denominator;
  void swap(Rational& other) {
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
Rational operator+(const Rational& lhs, const Rational& rhs) {
  Rational ans(lhs);
  ans += rhs;
  return ans;
}
Rational operator-(const Rational& lhs, const Rational& rhs) {
  Rational ans(lhs);
  ans -= rhs;
  return ans;
}
Rational operator*(const Rational& lhs, const Rational& rhs) {
  Rational ans(lhs);
  ans *= rhs;
  return ans;
}
Rational operator/(const Rational& lhs, const Rational& rhs) {
  Rational ans(lhs);
  ans /= rhs;
  return ans;
} 
