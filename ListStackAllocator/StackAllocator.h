#pragma once
#include <iterator>

template <std::size_t N> class StackStorage {
public:
  StackStorage(const StackStorage &) = delete;
  StackStorage() = default;

  void *allocate(std::size_t n, std::size_t alignment) {
    auto space = N - (ptr - arr);
    void *p = ptr;
    if (std::align(alignment, n, p, space)) {
      ptr = static_cast<char *>(p) + n;
      return p;
    }
    return nullptr;
  }

  void deallocate(void *, std::size_t) {}

private:
  char arr[N];
  char *ptr = arr;
};

template <typename T, std::size_t N> class StackAllocator {
public:
  using value_type = T;
  StackAllocator() = default;
  StackAllocator(StackStorage<N> &storage) : storage(&storage) {}

  template <typename U>
  StackAllocator(const StackAllocator<U, N> &other) : storage(other.storage) {}

  ~StackAllocator() = default;

  T *allocate(std::size_t n) {
    return static_cast<T *>(storage->allocate(n * sizeof(T), alignof(T)));
  }

  void deallocate(T *p, std::size_t n) {
    storage->deallocate(p, n * sizeof(T));
  }

  template <typename U> struct rebind { using other = StackAllocator<U, N>; };

  bool operator==(const StackAllocator &other) const {
    return storage == other.storage;
  }

  bool operator!=(const StackAllocator &other) const {
    return !(*this == other);
  }

private:
  StackStorage<N> *storage;

  template <typename U, std::size_t M> friend class StackAllocator;
};

template <typename T, typename Alloc = std::allocator<T>> class List {
private:
  std::size_t sz;

  struct Node;

  struct BaseNode {
    Node *prev;
    Node *next;
    BaseNode() : prev(nullptr), next(nullptr) {}
    BaseNode(Node *p, Node *n) : prev(p), next(n) {}
  };

  struct Node : BaseNode {
    T value;
    template <typename... Args>
    Node(Node *prev, Node *next, Args &&... args)
        : BaseNode(prev, next), value(std::forward<Args>(args)...) {}
  };

  using node_allocator =
      typename std::allocator_traits<Alloc>::template rebind_alloc<Node>;
  using alloc_traits = std::allocator_traits<node_allocator>;

  [[no_unique_address]] node_allocator allocator;
  BaseNode fake_node;

  Node *fake_to_node() { return reinterpret_cast<Node *>(&fake_node); }
  const Node *fake_to_node() const {
    return reinterpret_cast<const Node *>(&fake_node);
  }

  void initialize_empty() {
    fake_node.prev = fake_to_node();
    fake_node.next = fake_to_node();
    sz = 0;
  }

public:
  List() : sz(0) { initialize_empty(); }

  explicit List(const Alloc &alloc) : allocator(alloc) { initialize_empty(); }

  explicit List(std::size_t count, const T &value, const Alloc &alloc = Alloc())
      : allocator(alloc) {
    try {
      for (std::size_t i = 0; i < count; ++i) {
        push_back(value);
      }
    } catch (...) {
      clear();
      throw;
    }
  }

  explicit List(std::size_t count, const Alloc &alloc = Alloc())
      : allocator(alloc) {
    initialize_empty();
    try {
      for (std::size_t i = 0; i < count; ++i) {
        emplace_back();
      }
    } catch (...) {
      clear();
      throw;
    }
  }

  List(const List &other, const Alloc &alloc) : allocator(alloc) {
    initialize_empty();
    Fill(other);
  }

  List(List &&other) noexcept
      : allocator(std::move(other.allocator)), sz(other.sz) {
    if (sz > 0) {
      fake_node.next = other.fake_node.next;
      fake_node.prev = other.fake_node.prev;

      fake_node.next->prev = fake_to_node();
      fake_node.prev->next = fake_to_node();

      other.fake_node.next = other.fake_to_node();
      other.fake_node.prev = other.fake_to_node();
    } else {
      initialize_empty();
    }
  }

  List &operator=(List &&other) noexcept {
    if (this != &other) {
      clear();
      if (alloc_traits::propagate_on_container_move_assignment::value) {
        allocator = std::move(other.allocator);
        if (other.sz > 0) {
          fake_node = other.fake_node;
          fake_node.next->prev = fake_to_node();
          fake_node.prev->next = fake_to_node();
          sz = other.sz;
        } else {
          initialize_empty();
        }
        other.initialize_empty();
        other.sz = 0;
      } else {
        for (auto &elem : other) {
          push_back(std::move(elem));
        }
        other.clear();
      }
    }
    return *this;
  }
  ~List() { clear(); }

  List(const List &other)
      : allocator(alloc_traits::select_on_container_copy_construction(
            other.allocator)) {
    initialize_empty();
    Fill(other);
  }
  void Fill(const List &other) {
    List tmp(allocator);
    for (const auto &item : other) {
      tmp.push_back(item);
    }

    swap(tmp);
  }
  List &operator=(const List &other) {
    if (this == &other) {
      return *this;
    }

    Alloc new_alloc =
        alloc_traits::propagate_on_container_copy_assignment::value
            ? other.allocator
            : allocator;

    List temp(new_alloc);
    allocator = new_alloc;
    for (const auto &item : other) {
      temp.push_back(item);
    }
    swap(temp);
    return *this;
  }

  void swap(List &other) noexcept {
    using std::swap;
    swap(sz, other.sz);
    BaseNode temp_fake = fake_node;
    fake_node = other.fake_node;
    other.fake_node = temp_fake;

    if (sz > 0) {
      fake_to_node()->prev->next = fake_to_node();
      fake_to_node()->next->prev = fake_to_node();
    }
    if (other.sz > 0) {
      other.fake_to_node()->prev->next = other.fake_to_node();
      other.fake_to_node()->next->prev = other.fake_to_node();
    }
  }

  void clear() noexcept {
    while (sz > 0) {
      pop_back();
    }
  }

  template <bool IsConst> class base_iterator {
  public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = std::conditional_t<IsConst, const T *, T *>;
    using reference = std::conditional_t<IsConst, const T &, T &>;
    using node_type = std::conditional_t<IsConst, const Node, Node>;

    base_iterator() = default;
    explicit base_iterator(node_type *node) : node(const_cast<Node *>(node)) {}

    operator base_iterator<true>() const { return base_iterator<true>(node); }

    base_iterator<IsConst> &operator++() {
      node = node->next;
      return *this;
    }

    base_iterator<IsConst> operator++(int) {
      base_iterator<IsConst> cp = *this;
      node = node->next;
      return cp;
    }

    base_iterator &operator--() {
      node = node->prev;
      return *this;
    }

    base_iterator operator--(int) {
      base_iterator<IsConst> cp = *this;
      node = node->prev;
      return cp;
    }

    reference operator*() const {
      return static_cast<node_type *>(node)->value;
    }
    pointer operator->() const { return &node->value; }

    bool operator==(const base_iterator &other) const {
      return (node == other.node);
    }

    bool operator!=(const base_iterator &other) const {
      return (node != other.node);
    }

    Node *get_ptr() const { return node; }

  private:
    Node *node;
  };

  using iterator = base_iterator<false>;
  using const_iterator = base_iterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  std::size_t size() const noexcept { return sz; }
  bool empty() const noexcept { return (sz == 0); }
  node_allocator get_allocator() const noexcept { return allocator; }

  template <typename... Args> void emplace_back(Args &&... args) {
    Node *new_node = alloc_traits::allocate(allocator, 1);
    try {
      alloc_traits::construct(allocator, new_node, fake_node.prev,
                              fake_to_node(), std::forward<Args>(args)...);
    } catch (...) {
      alloc_traits::deallocate(allocator, new_node, 1);
      throw;
    }

    fake_node.prev->next = new_node;
    fake_node.prev = new_node;
    if (sz == 0) {
      fake_node.next = new_node;
    }
    ++sz;
  }

  void push_back(const T &elem) { emplace_back(elem); }
  void push_back(T &&elem) { emplace_back(std::move(elem)); }

  template <typename... Args> void emplace_front(Args &&... args) {
    Node *new_node = alloc_traits::allocate(allocator, 1);
    try {
      alloc_traits::construct(allocator, new_node, fake_to_node(),
                              fake_node.next, std::forward<Args>(args)...);
    } catch (...) {
      alloc_traits::deallocate(allocator, new_node, 1);
      throw;
    }

    fake_node.next->prev = new_node;
    fake_node.next = new_node;
    if (sz == 0) {
      fake_node.prev = new_node;
    }
    ++sz;
  }
  void push_front(const T &elem) { emplace_front(elem); }
  void push_front(T &&elem) { emplace_front(std::move(elem)); }
  template <typename... Args>
  iterator emplace(const_iterator iter, Args &&... args) {
    Node *pos = iter.get_ptr();
    Node *new_node = alloc_traits::allocate(allocator, 1);
    try {
      alloc_traits::construct(allocator, new_node, pos->prev, pos,
                              std::forward<Args>(args)...);
    } catch (...) {
      alloc_traits::deallocate(allocator, new_node, 1);
      throw;
    }

    pos->prev->next = new_node;
    pos->prev = new_node;
    ++sz;
    return iterator(new_node);
  }
  iterator insert(const_iterator iter, const T &elem) {
    return emplace(iter, elem);
  }
  iterator insert(const_iterator iter, T &&elem) {
    return emplace(iter, std::move(elem));
  }
  void pop_back() { erase(--end()); }
  void pop_front() { erase(begin()); }

  iterator erase(const_iterator iter) {
    Node *to_delete = iter.get_ptr();
    Node *next_node = to_delete->next;
    to_delete->prev->next = to_delete->next;
    to_delete->next->prev = to_delete->prev;
    alloc_traits::destroy(allocator, to_delete);
    alloc_traits::deallocate(allocator, to_delete, 1);
    --sz;
    return iterator(next_node);
  }

  iterator begin() noexcept { return iterator(fake_node.next); }
  iterator end() noexcept { return iterator(fake_to_node()); }
  const_iterator begin() const noexcept {
    return const_iterator(fake_node.next);
  }
  const_iterator end() const noexcept { return const_iterator(fake_to_node()); }
  const_iterator cbegin() const noexcept {
    return const_iterator(fake_node.next);
  }
  const_iterator cend() const noexcept {
    return const_iterator(fake_to_node());
  }
  reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
  reverse_iterator rend() noexcept { return reverse_iterator(begin()); }
  const_reverse_iterator rbegin() const noexcept {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator rend() const noexcept {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crbegin() const noexcept { return rbegin(); }
  const_reverse_iterator crend() const noexcept { return rend(); }
};
