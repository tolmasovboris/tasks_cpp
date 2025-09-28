#pragma once
#include <compare>
#include <limits>
#include <span>

namespace circular_buffer_details {

const std::size_t DYNAMIC_CAPACITY = std::numeric_limits<std::size_t>::max();

template <std::size_t Capacity>
class CapacityHelper {
 public:
  static constexpr std::size_t cp = Capacity;
};

template <>
class CapacityHelper<DYNAMIC_CAPACITY> {
 public:
  std::size_t cp = 0;
};

}
using namespace circular_buffer_details;
template <typename T, std::size_t Capacity = DYNAMIC_CAPACITY>
class CircularBuffer : private CapacityHelper<Capacity> {
 public:
  CircularBuffer() : head(0), tail(0), sz(0) {
    if constexpr (Capacity == DYNAMIC_CAPACITY) {
      static_assert(!std::is_same_v<T, T>,
                    "Dynamic capacity is forbidden in this constructor");
    }
  }

  explicit CircularBuffer(std::size_t capacity) : head(0), tail(0), sz(0) {
    if constexpr (Capacity == DYNAMIC_CAPACITY) {
      this->cp = capacity;
      arr = static_cast<T*>(operator new(capacity * sizeof(T)));
    } else {
      if (Capacity != capacity) {
        throw std::invalid_argument(
            "Capacity parameter must equal to template's parameter");
      }
    }
  }

  CircularBuffer(const CircularBuffer& other)
      : CircularBuffer(other.capacity()) {
    for (std::size_t i = 0; i < other.sz; ++i) {
      push_back(other.get_storage()[(other.head + i) % other.capacity()]);
    }
  }

  CircularBuffer& operator=(const CircularBuffer& other) {
    if (this != &other) {
      clear();

      if constexpr (Capacity == DYNAMIC_CAPACITY) {
        if (this->cp != other.cp) {
          operator delete(arr);
          this->cp = other.cp;
          arr = static_cast<T*>(operator new(this->cp * sizeof(T)));
        }
      }

      head = other.head;
      tail = other.tail;
      sz = other.sz;

      try {
        for (std::size_t i = 0; i < sz; ++i) {
          new (&get_storage()[(head + i) % capacity()])
              T(other.get_storage()[(other.head + i) % other.capacity()]);
        }
      } catch (...) {
        for (std::size_t i = 0; i < sz; ++i) {
          get_storage()[(head + i) % capacity()].~T();
        }
        sz = 0;
        throw;
      }
    }
    return *this;
  }

  void clear() {
    for (std::size_t i = 0; i < sz; ++i) {
      get_storage()[(head + i) % capacity()].~T();
    }
    sz = 0;
    head = 0;
    tail = 0;
  }

  ~CircularBuffer() {
    clear();
    if constexpr (Capacity == DYNAMIC_CAPACITY) {
      ::operator delete(arr);
    }
  }

  template <bool IsConst>
  class base_iterator {
   public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = std::conditional_t<IsConst, const T*, T*>;
    using reference = std::conditional_t<IsConst, const T&, T&>;
    using SpanType = std::conditional_t<
        Capacity == DYNAMIC_CAPACITY,
        std::conditional_t<IsConst, std::span<const T>, std::span<T>>,
        std::conditional_t<IsConst, std::span<const T, Capacity>,
                           std::span<T, Capacity>>>;

    base_iterator() = default;

    base_iterator(pointer data, std::size_t capacity, std::size_t index,
                  std::size_t count)
        : data_span_(data, capacity), index_(index), count_(count) {
      if constexpr (Capacity == DYNAMIC_CAPACITY) {
        data_span_ = SpanType(data, capacity);
      } else {
        data_span_ = SpanType(data, capacity);
      }
    }

    template <bool OtherIsConst,
              typename = std::enable_if_t<IsConst && !OtherIsConst>>
    base_iterator(const base_iterator<OtherIsConst>& other)
        : data_span_(other.data_span_),
          index_(other.index_),
          count_(other.count_) {}

    template <bool OtherIsConst,
              typename = std::enable_if_t<IsConst && !OtherIsConst>>
    base_iterator& operator=(const base_iterator<OtherIsConst>& other) {
      data_span_ = other.data_span_;
      index_ = other.index_;
      count_ = other.count_;
      return *this;
    }

    reference operator*() const { return data_span_[index_]; }
    pointer operator->() const { return &data_span_[index_]; }

    base_iterator& operator++() {
      index_ = (index_ + 1) % data_span_.size();
      ++count_;
      return *this;
    }

    base_iterator operator++(int) {
      base_iterator temp = *this;
      ++(*this);
      return temp;
    }

    base_iterator& operator--() {
      index_ = (index_ - 1 + data_span_.size()) % data_span_.size();
      --count_;
      return *this;
    }

    base_iterator operator--(int) {
      base_iterator temp = *this;
      --(*this);
      return temp;
    }

    base_iterator& operator+=(difference_type n) {
      index_ = (index_ + n) % data_span_.size();
      count_ += n;
      return *this;
    }

    base_iterator operator+(difference_type n) const {
      base_iterator tmp = *this;
      tmp += n;
      return tmp;
    }

    friend base_iterator operator+(difference_type n,
                                   const base_iterator& iter) {
      return iter + n;
    }

    base_iterator& operator-=(difference_type n) {
      index_ = (index_ - n + data_span_.size()) % data_span_.size();
      count_ -= n;
      return *this;
    }

    base_iterator operator-(difference_type n) const {
      base_iterator tmp = *this;
      tmp -= n;
      return tmp;
    }

    difference_type operator-(const base_iterator& other) const {
      return static_cast<difference_type>(count_) -
             static_cast<difference_type>(other.count_);
    }

    reference operator[](difference_type n) const noexcept {
      return *(*this + n);
    }
    std::partial_ordering operator<=>(const base_iterator& other) const {
      if (data_span_.data() != other.data_span_.data() ||
          data_span_.size() != other.data_span_.size()) {
        return std::partial_ordering::unordered;
      }
      return count_ <=> other.count_;
    }
    bool operator==(const base_iterator& other) const {
      return data_span_.data() == other.data_span_.data() &&
             data_span_.size() == other.data_span_.size() &&
             count_ == other.count_;
    }

    std::size_t get_index() const { return index_; }

   private:
    SpanType data_span_;
    std::size_t index_ = 0;
    std::size_t count_ = 0;

    template <bool OtherIsConst>
    friend class base_iterator;
  };

  using iterator = base_iterator<false>;
  using const_iterator = base_iterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() { return iterator(get_storage(), capacity(), head, 0); }

  const_iterator begin() const {
    return const_iterator(get_storage(), capacity(), head, 0);
  }

  const_iterator cbegin() const { return begin(); }

  iterator end() { return iterator(get_storage(), capacity(), tail, sz); }

  const_iterator end() const {
    return const_iterator(get_storage(), capacity(), tail, sz);
  }

  const_iterator cend() const { return end(); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return const_reverse_iterator(end());
  }

  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return const_reverse_iterator(begin());
  }

  void push_back(const T& elem) {
    if (full()) {
        get_storage()[head] = elem;
        head = (head + 1) % capacity();
    } else {
        try {
            new (&get_storage()[tail]) T(elem);
        } catch (...) {
            throw;
        }
        ++sz;
    }
    tail = (tail + 1) % capacity();
}

void push_front(const T& elem) {
  std::size_t new_head = (head - 1 + capacity()) % capacity();
  if (full()) {
      get_storage()[new_head] = elem;
      head = new_head;
      tail = (tail - 1 + capacity()) % capacity();
  } else {
      try {
          new (&get_storage()[new_head]) T(elem);
      } catch (...) {
          throw;
      }
      head = new_head;
      ++sz;
  }
}

  void pop_back() {
    if (empty()) return;

    get_storage()[(tail - 1 + capacity()) % capacity()].~T();
    tail = (tail - 1 + capacity()) % capacity();
    --sz;
  }

  void pop_front() {
    if (empty()) return;

    get_storage()[head].~T();
    head = (head + 1) % capacity();
    --sz;
  }

  void insert(iterator iter, const T& elem) {
    if (full() && iter == begin()) {
      return;
    }
    if (full()) {
      pop_front();
    }
    std::size_t insert_index = iter.get_index();
    for (std::size_t i = tail; i != insert_index;
         i = (i - 1 + capacity()) % capacity()) {
      std::size_t prev = (i - 1 + capacity()) % capacity();
      new (&get_storage()[i]) T(get_storage()[prev]);
      get_storage()[prev].~T();
    }
    new (&get_storage()[insert_index]) T(elem);
    tail = (tail + 1) % capacity();
    ++sz;
  }

  void erase(iterator pos) {
    if (empty()) return;

    std::size_t erase_pos = pos.get_index();
    for (std::size_t i = erase_pos; i != (tail - 1 + capacity()) % capacity();
         i = (i + 1) % capacity()) {
      std::size_t next = (i + 1) % capacity();
      get_storage()[i] = get_storage()[next];
    }
    get_storage()[(tail - 1 + capacity()) % capacity()].~T();
    tail = (tail - 1 + capacity()) % capacity();
    --sz;
  }

  T& at(std::size_t index) {
    if (index >= sz) throw std::out_of_range("Index out of range");
    return get_storage()[(head + index) % capacity()];
  }

  const T& at(std::size_t index) const {
    if (index >= sz) throw std::out_of_range("Index out of range");
    return get_storage()[(head + index) % capacity()];
  }

  T& operator[](std::size_t index) noexcept {
    return get_storage()[(head + index) % capacity()];
  }
  const T& operator[](std::size_t index) const noexcept {
    return get_storage()[(head + index) % capacity()];
  }

  bool full() const { return sz == capacity(); }
  bool empty() const { return sz == 0; }
  std::size_t capacity() const {
    if constexpr (Capacity == DYNAMIC_CAPACITY) {
      return this->cp;
    } else {
      return Capacity;
    }
  }
  std::size_t size() const { return sz; }

 private:
  T* get_storage() {
    if constexpr (Capacity == DYNAMIC_CAPACITY) {
      return arr;
    } else {
      return reinterpret_cast<T*>(arr.data());
    }
  }

  const T* get_storage() const {
    if constexpr (Capacity == DYNAMIC_CAPACITY) {
      return arr;
    } else {
      return reinterpret_cast<const T*>(arr.data());
    }
  }

  std::size_t head = 0;
  std::size_t tail = 0;
  std::size_t sz = 0;

  using BufferType =
      std::conditional_t<Capacity == DYNAMIC_CAPACITY, T*,
                         std::array<std::byte, sizeof(T) * Capacity>>;
  alignas(alignof(T) * 8) BufferType arr;                         
};