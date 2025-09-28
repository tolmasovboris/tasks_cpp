#pragma once
#include <algorithm>
template<typename T>
class EnableSharedFromThis;
struct BaseControlBlock
{
  size_t shared_count;
  size_t weak_count;

  BaseControlBlock(size_t shared_count, size_t weak_count)
    : shared_count(shared_count)
    , weak_count(weak_count)
  {}
  BaseControlBlock()
    : BaseControlBlock(1, 0)
  {}
  virtual ~BaseControlBlock() = default;

  virtual void destroy_object(void* ptr) = 0;
  virtual void destroy_block() = 0;
  virtual void* get_ptr() const = 0;
};

template<typename T>
class WeakPtr;

template<typename T>
class SharedPtr
{
private:
  BaseControlBlock* control;

  template<typename U>
  friend class SharedPtr;
  template<typename U>
  friend class WeakPtr;

  template<typename Deleter, typename Alloc>
  struct RegularControlBlock : BaseControlBlock
  {
    T* ptr;
    [[no_unique_address]] Deleter deleter;
    [[no_unique_address]] Alloc allocator;

    RegularControlBlock(T* p, Deleter d = std::default_delete<T>(), Alloc a = std::allocator<T>())
      : BaseControlBlock(1, 0)
      , ptr(p)
      , deleter(std::move(d))
      , allocator(std::move(a))
    {}

    void* get_ptr() const override { return ptr; }

    void destroy_object(void* ptr) override { deleter(static_cast<T*>(ptr)); }

    void destroy_block() override
    {
      using ControlBlockAlloc = typename std::allocator_traits<
        Alloc>::template rebind_alloc<RegularControlBlock>;
      ControlBlockAlloc cb_alloc(allocator);
      this->~RegularControlBlock();
      std::allocator_traits<ControlBlockAlloc>::deallocate(cb_alloc, this, 1);
    }
  };

  template<typename Alloc>
  class MakeSharedControlBlock : public BaseControlBlock
  {
  public:
    mutable T object;
    [[no_unique_address]] Alloc alloc;
    using alloc_block = typename std::allocator_traits<
      Alloc>::template rebind_alloc<MakeSharedControlBlock>;
    template<typename... Args>
    explicit MakeSharedControlBlock(const Alloc& alloc, Args&&... args)
      : object(std::forward<Args>(args)...)
      , alloc(alloc)
    {}

    void* get_ptr() const override { return &object; }

    void destroy_object(void* ptr) override
    {
      std::allocator_traits<Alloc>::destroy(alloc, reinterpret_cast<T*>(ptr));
    }

    void destroy_block() override
    {
      alloc_block control_block = alloc;
      std::allocator_traits<alloc_block>::deallocate(control_block, this, 1);
      (&alloc)->~Alloc();
    }
  };

  explicit SharedPtr(BaseControlBlock* control)
    : control(control)
  {}

public:
  SharedPtr()
    : control(nullptr)
  {}
  template<typename U, typename Deleter = std::default_delete<U>, typename Allocator = std::allocator<U>,
         typename = std::enable_if_t<std::is_convertible_v<U*, T*>>>
explicit SharedPtr(U* ptr, Deleter d = Deleter(), Allocator a = Allocator())
    : control(nullptr)
{
    if constexpr (std::is_base_of_v<EnableSharedFromThis<T>, U>) {
        if (ptr) {
            ptr->weak_this = *this;
        }
    }

    using Block = RegularControlBlock<Deleter, Allocator>;
    using BlockAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<Block>;

    BlockAlloc block_alloc(a);
    Block* cb = block_alloc.allocate(1);
    try {
        new (cb) Block(ptr, std::move(d), std::move(a));
        control = cb;
    } catch (...) {
        block_alloc.deallocate(cb, 1);
        d(ptr);
        throw;
    }
}

  T* get() const
  {
    return control ? static_cast<T*>(control->get_ptr()) : nullptr;
  }
  SharedPtr(const SharedPtr& other)
    : control(other.control)
  {
    if (control) {
      ++control->shared_count;
    }
  }

  template<typename U>
  SharedPtr(const SharedPtr<U>& other)
    : control(other.control)
  {
    if (control) {
      ++control->shared_count;
    }
  }

  SharedPtr(SharedPtr&& other)
    : control(other.control)
  {
    other.control = nullptr;
  }

  template<typename U>
  SharedPtr(SharedPtr<U>&& other)
    : control(other.control)
  {
    other.control = nullptr;
  }

  ~SharedPtr()
  {
    if (control) {
      if (--control->shared_count == 0) {
        control->destroy_object(control->get_ptr());
        if (control->weak_count == 0) {
          control->destroy_block();
        }
      }
    }
  }

  SharedPtr& operator=(const SharedPtr& other)
  {
    SharedPtr(other).swap(*this);
    return *this;
  }

  SharedPtr& operator=(SharedPtr&& other)
  {
    SharedPtr(std::move(other)).swap(*this);
    return *this;
  }

  void swap(SharedPtr& other) { std::swap(control, other.control); }

  size_t use_count() const { return control ? control->shared_count : 0; }

  T& operator*() const { return *get(); }

  T* operator->() const { return get(); }

  void reset() { SharedPtr().swap(*this); }

  template<typename U,
           typename Deleter = std::default_delete<U>,
           typename Allocator = std::allocator<U>>
  void reset(U* new_ptr, Deleter d = Deleter(), Allocator a = Allocator())
  {
    SharedPtr(new_ptr, std::move(d), std::move(a)).swap(*this);
  }

  template<typename Alloc, typename... Args>
  static SharedPtr allocate_shared(const Alloc& alloc, Args&&... args)
  {
    using ControlBlock = MakeSharedControlBlock<Alloc>;
    using BlockAlloc = typename std::allocator_traits<
      Alloc>::template rebind_alloc<ControlBlock>;

    BlockAlloc block_alloc(alloc);
    ControlBlock* cb = block_alloc.allocate(1);
    try {
      std::allocator_traits<BlockAlloc>::construct(
        block_alloc, cb, block_alloc, std::forward<Args>(args)...);
      return SharedPtr(cb);
    } catch (...) {
      block_alloc.deallocate(cb, 1);
      throw;
    }
  }

  template<typename... Args>
  static SharedPtr make_shared(Args&&... args)
  {
    return allocate_shared(std::allocator<T>(), std::forward<Args>(args)...);
  }
};

template<typename T, typename... Args>
SharedPtr<T>
makeShared(Args&&... args)
{
  return SharedPtr<T>::make_shared(std::forward<Args>(args)...);
  ;
}

template<typename T, typename Alloc, typename... Args>
SharedPtr<T>
allocateShared(const Alloc& alloc, Args&&... args)
{
  return SharedPtr<T>::allocate_shared(alloc, std::forward<Args>(args)...);
}

template<typename T>
class WeakPtr
{
private:
  BaseControlBlock* control;
  template<typename U>
  friend class SharedPtr;
  template<typename U>
  friend class WeakPtr;

public:
  WeakPtr()
    : control(nullptr)
  {}

  WeakPtr(const WeakPtr& other)
    : control(other.control)
  {
    if (control) {
      ++control->weak_count;
    }
  }

  template<typename U>
  WeakPtr(const WeakPtr<U>& other)
    : control(other.control)
  {
    if (control) {
      ++control->weak_count;
    }
  }

  WeakPtr(WeakPtr&& other)
    : control(other.control)
  {
    other.control = nullptr;
  }

  template<typename U>
  WeakPtr(WeakPtr<U>&& other)
    : control(other.control)
  {
    other.control = nullptr;
  }

  WeakPtr(const SharedPtr<T>& other)
    : control(other.control)
  {
    if (control) {
      ++control->weak_count;
    }
  }

  template<typename U>
  WeakPtr(const SharedPtr<U>& other)
    : control(other.control)
  {
    if (control) {
      ++control->weak_count;
    }
  }

  ~WeakPtr()
  {
    if (control) {
      if (--control->weak_count == 0 && control->shared_count == 0) {
        control->destroy_block();
      }
    }
  }

  template<typename U>
  WeakPtr& operator=(const WeakPtr<U>& other)
  {
    WeakPtr tmp(other);
    swap(tmp);
    return *this;
  }

  WeakPtr& operator=(const WeakPtr& other)
  {
    WeakPtr tmp(other);
    swap(tmp);
    return *this;
  }

  template<typename U>
  WeakPtr& operator=(WeakPtr<U>&& other)
  {
    WeakPtr tmp(std::move(other));
    swap(tmp);
    return *this;
  }

  WeakPtr& operator=(WeakPtr&& other)
  {
    WeakPtr tmp(std::move(other));
    swap(tmp);
    return *this;
  }

  template<typename U>
  WeakPtr& operator=(const SharedPtr<U>& other)
  {
    WeakPtr tmp(other);
    swap(tmp);
    return *this;
  }

  void swap(WeakPtr& other) { std::swap(control, other.control); }

  T* get() const
  {
    return control ? static_cast<T*>(control->get_ptr()) : nullptr;
  }

  bool expired() const { return use_count() == 0; }

  size_t use_count() const { return control ? control->shared_count : 0; }

  SharedPtr<T> lock() const
  {
    if (!control || control->shared_count == 0) {
      return SharedPtr<T>();
    }
    ++control->shared_count;
    return SharedPtr<T>(control);
  }
};

template<typename T>
class EnableSharedFromThis
{
private:
  mutable WeakPtr<T> weak_this;
  template<typename U>
  friend class SharedPtr;
  friend SharedPtr<T> makeShared();

protected:
  EnableSharedFromThis() = default;
  EnableSharedFromThis(const EnableSharedFromThis&) = default;
  EnableSharedFromThis& operator=(const EnableSharedFromThis&) { return *this; }
  ~EnableSharedFromThis() = default;

public:
  SharedPtr<T> shared_from_this() { return weak_this.lock(); }

  SharedPtr<const T> shared_from_this() const { return weak_this.lock(); }

  WeakPtr<T> weak_from_this() noexcept { return weak_this; }

  WeakPtr<const T> weak_from_this() const noexcept { return weak_this; }
};