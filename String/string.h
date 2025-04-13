#include <cassert>
#include <cstring>
#include <iostream>
class String {
private:
  char *str = nullptr;
  size_t sz;
  size_t cap;

public:
  String(char x) : str(new char[2]), sz(1), cap(2) {
    str[0] = x;
    str[1] = '\0';
  }
  String(const char *other = "")
      : str(new char[strlen(other) + 1]), sz(strlen(other)), cap(sz + 1) {
    memcpy(str, other, sz);
    str[sz] = '\0';
  }
  String(size_t count, char elem)
      : str(new char[count + 1]), sz(count), cap(count + 1) {
    memset(str, elem, count);
    str[count] = '\0';
  }

  String(const String &other)
      : str(new char[other.cap]), sz(other.sz), cap(other.cap) {
    memcpy(str, other.str, sz + 1);
  }
  String &operator=(String other) {
    swap(other);
    return *this;
  }
  String &operator+=(char x) {
    push_back(x);
    return *this;
  }

  String &operator+=(const String &other) {
    if (cap < other.sz + sz + 1) {
      cap = std::max(cap * 2, sz + other.sz + 1);
      char *new_data = new char[cap];
      memcpy(new_data, str, sz);
      delete[] str;
      str = new_data;
    }
    memcpy(str + sz, other.str, other.sz);
    sz += other.sz;
    str[sz] = '\0';
    return *this;
  }
  char &operator[](size_t index) { return str[index]; }
  const char &operator[](size_t index) const { return str[index]; }
  const char &back() const { return str[sz - 1]; }
  const char &front() const { return str[0]; }
  size_t length() const { return sz; }
  size_t size() { return sz; }
  size_t capacity() { return cap - 1; }
  char &back() { return str[sz - 1]; }
  char &front() { return str[0]; }
  char pop_back() {
    char elem = str[sz - 1];
    str[sz - 1] = '\0';
    sz--;
    return elem;
  }
  ~String() { delete[] str; }
  void swap(String &other) {
    std::swap(sz, other.sz);
    std::swap(cap, other.cap);
    std::swap(str, other.str);
  }
  void push_back(char x) {
    if (sz + 1 >= cap) {
      cap *= 2;
      char *new_data = new char[cap];
      memcpy(new_data, str, sz + 1);
      delete[] str;
      str = new_data;
    }
    str[sz++] = x;
    str[sz] = '\0';
  }
  size_t find(const String &substring) const {
    for (size_t i = 0; i <= length() - substring.length(); ++i) {
      if (strncmp(data() + i, substring.data(), substring.length()) == 0) {
        return i;
      }
    }
    return length();
  }
  size_t rfind(const String &substring) const {
    for (size_t i = length() - substring.length(); i > 0; --i) {
      if (strncmp(data() + i, substring.data(), substring.length()) == 0) {
        return i;
      }
    }
    return length();
  }
  bool empty() { return sz == 0; }
  bool empty() const { return sz == 0; }
  void clear() {
    str[0] = '\0';
    sz = 0;
  }
  void shrink_to_fit() {
    char *copy = str;
    str = new char[sz + 1];
    memcpy(str, copy, sz + 1);
    cap = sz + 1;
    delete[] copy;
  }
  String substr(size_t start, size_t count) const {
    String result(count, ' ');
    for (size_t i = start; i < start + count; ++i) {
      result[i - start] = str[i];
    }
    return result;
  }
  const char *data() const { return str; }
  char *data() { return str; }
};
bool operator==(const String &first, const String &second) {
  return first.length() == second.length() &&
         std::strcmp(first.data(), second.data()) == 0;
}
bool operator!=(const String &first, const String &second) {
  return !(first == second);
}

bool operator<(const String &first, const String &second) {
  return strcmp(first.data(), second.data()) < 0;
}
bool operator>=(const String &first, const String &second) {
  return !(first < second);
}
bool operator>(const String &first, const String &second) {
  return !(first < second) && (first != second);
}
bool operator<=(const String &first, const String &second) {
  return !(first > second);
}
String operator+(const String &lhs, const String &rhs) {
  String left(lhs.data());
  String right(rhs.data());
  left += right;
  return left;
}
std::ostream &operator<<(std::ostream &out, const String &c) {
  out << c.data();
  return out;
}
std::istream &operator>>(std::istream &in, String &c) {
  c.clear();
  char symbol = in.get();
  while (symbol != EOF && !std::isspace(symbol)) {
    c.push_back(symbol);
    symbol = in.get();
  }
  return in;
}
