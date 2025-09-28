#pragma once

struct cat_tag{};

template<bool Repeat, typename T, typename... Ts>
struct contains_elem : std::false_type {};

template<typename T, typename Head, typename... Tail>
struct contains_elem<false, T, Head, Tail...> : 
std::conditional_t<std::is_same_v<T, Head>, contains_elem<true, T, Tail...>, contains_elem<false, T, Tail...>>{};

template<typename T, typename Head, typename... Tail>
struct contains_elem<true, T, Head, Tail...> : 
std::conditional_t<std::is_same_v<T, Head>, std::false_type, contains_elem<true, T, Tail...>>{};

template<typename T, typename Head>
struct contains_elem<false, T, Head> : 
std::conditional_t<std::is_same_v<T, Head>, std::true_type, std::false_type> {}; 

template<typename T, typename Head>
struct contains_elem<true, T, Head> : 
std::conditional_t<std::is_same_v<T, Head>, std::false_type, std::true_type> {};

template<typename T, typename... Ts>
bool const contains_elem_v = contains_elem<false, T, Ts...>::value;

template<size_t N, typename... Ts>
struct GetTypeByIndex;

template<typename Head, typename... Tail>
struct GetTypeByIndex<0, Head, Tail...> {
    using type = Head;
};

template<size_t N, typename Head, typename... Tail>
struct GetTypeByIndex<N, Head, Tail...> {
    using type = typename GetTypeByIndex<N-1, Tail...>::type;
};

template<typename T, typename = void>
struct is_cpinit : std::false_type {};

template<typename T>
struct is_cpinit<T, 
    std::void_t<decltype(T{} = {})>
> : std::true_type {};

template<typename T>
constexpr bool is_cpinit_v = is_cpinit<T>::value;

template<typename... Ts>
struct AllCpinit : std::false_type {};

template<typename Elem>
struct AllCpinit<Elem> : std::bool_constant<is_cpinit_v<Elem>> {};

template<typename Head, typename... Tail>
struct AllCpinit<Head, Tail...> : std::conditional_t<is_cpinit_v<Head>, AllCpinit<Tail...>, is_cpinit<Head>> {};


template<typename... Ts>
constexpr bool AllCpinit_v = AllCpinit<Ts...>::value;

template<typename... Ts>
class Tuple {};
template <typename...>
struct ConcatTypes;

template <typename T>
struct ConcatTypes<T> : std::type_identity<T> {};

template <typename... Args1, typename... Args2, typename... Tuples>
struct ConcatTypes<Tuple<Args1...>, Tuple<Args2...>, Tuples...> 
: std::type_identity<typename ConcatTypes<Tuple<Args1..., Args2...>, Tuples...>::type>{};

template<typename Head, typename... Tail>
class Tuple<Head, Tail...> {
    Head head;
    Tuple<Tail...> tail;

    template<size_t N, typename... Ts>
    friend typename GetTypeByIndex<N, Ts...>::type& get(Tuple<Ts...>& t);

    template<typename T, typename... Ts>
    friend T& get(Tuple<Ts...>& t);

    template<typename... UTypes>
    friend class Tuple;

    template <typename... Args>
    friend typename ConcatTypes<std::remove_reference_t<Args>...>::type tupleCat(Args&&... args);

    template <typename UHead, typename... Us,
    std::enable_if_t<std::is_same_v<decltype(std::declval<UHead>().tail), Tuple<>>, bool> = true>
    Tuple(cat_tag tag, UHead&& uhead, Us&&... tail)
    : head(std::forward<decltype(uhead.head)>(uhead.head)),
      tail(tag, std::forward<Us>(tail)...) {}

    template <typename UHead, typename... Rest,
    std::enable_if_t<!std::is_same_v<decltype(std::declval<UHead>().tail), Tuple<>>, bool> = true>  
    Tuple(cat_tag tag, UHead&& other, Rest&&... rest)
    : head(other.head),
      tail(tag, other.tail, std::forward<Rest>(rest)...) {}

    template <typename UHead>
    Tuple(cat_tag, UHead&& other)
        : head(other.head),
        tail(std::forward<decltype(other.tail)>(other.tail)) {
    }
public:
    template<typename NewHead = Head, 
    std::enable_if_t<std::conjunction_v<std::is_default_constructible<NewHead>, 
    std::is_default_constructible<Tail>...>, bool> = true>
    explicit(!AllCpinit_v<Head, Tail...>) Tuple(): head(), tail() {}
    
    template<typename NewHead = Head,
    typename std::enable_if_t<std::conjunction_v<std::is_copy_constructible<NewHead>,
    std::is_copy_constructible<Tail>...>, bool> = true> 
    explicit(!std::conjunction_v<std::is_convertible<const NewHead&, NewHead>,
    std::is_convertible<const Tail&, Tail>...>) 
    Tuple(const Head& head, const Tail&... tail): head(head), tail(tail...) {}

    template<typename Uhead, typename... Utail, typename NewHead = Head, 
    std::enable_if_t<sizeof...(Utail) == sizeof...(Tail) && std::conjunction_v<std::is_constructible<Head, Uhead>, 
    std::is_constructible<Tail, Utail>...>, bool> = true>
    explicit(!std::conjunction_v<std::is_convertible<NewHead, Uhead>, std::is_convertible<Tail, Utail>...>)
    Tuple(Uhead&& head, Utail&&... tail) : head(std::forward<Uhead>(head)), tail(std::forward<Utail>(tail)...) {}

    template<typename Uhead, typename... Utail,
    std::enable_if_t<sizeof...(Tail) == sizeof...(Utail) && 
    std::conjunction_v<std::is_constructible<Head, const Uhead&>, 
    std::is_constructible<Tail, const Utail&>...>
    && (sizeof...(Tail) > 0 || (!std::is_convertible_v<Uhead, Head> && 
    !std::is_constructible_v<Head, Uhead> && !std::is_same_v<Head, Uhead>)), bool> = true>
    explicit(!std::conjunction_v<std::is_convertible<Uhead, Head>, std::is_convertible<Utail, Tail>...>)
    Tuple(const Tuple<Uhead, Utail...>& other): head(other.head), tail(other.tail) {}
      
    template<typename Uhead, typename... Utail,
    std::enable_if_t<sizeof...(Tail) == sizeof...(Utail) && 
    std::conjunction_v<std::is_constructible<Head, Uhead>, 
    std::is_constructible<Tail, Utail>...>
    && sizeof...(Tail) == 0, bool> = true>
    explicit(!std::conjunction_v<std::is_convertible<Uhead, Head>, std::is_convertible<Utail, Tail>...>)
    Tuple(const Tuple<Uhead, Utail...>& other): head(other.head) {}

    template<typename Uhead, typename... Utail,
    std::enable_if_t<sizeof...(Tail) == sizeof...(Utail) && 
    std::conjunction_v<std::is_constructible<Head, Uhead>, 
    std::is_constructible<Tail, Utail>...>
    && (sizeof...(Tail) > 0 || (!std::is_convertible_v<Uhead, Head> && 
    !std::is_constructible_v<Head, Uhead> && !std::is_same_v<Head, Uhead>)), bool> = true>
    explicit(!std::conjunction_v<std::is_convertible<Uhead, Head>, std::is_convertible<Utail, Tail>...>)
    Tuple(Tuple<Uhead, Utail...>&& other): head(std::forward<decltype(other.head)>(other.head)), 
    tail(std::forward<decltype(other.tail)>(other.tail)) {}
      
    template<typename Uhead, typename... Utail,
    std::enable_if_t<sizeof...(Tail) == sizeof...(Utail) && 
    std::conjunction_v<std::is_constructible<Head, Uhead>, 
    std::is_constructible<Tail, Utail>...>
    && sizeof...(Tail) == 0, bool> = true>
    explicit(!std::conjunction_v<std::is_convertible<Uhead, Head>, std::is_convertible<Utail, Tail>...>)
    Tuple(Tuple<Uhead, Utail...>&& other): head(std::forward<decltype(other.head) >(other.head)) {}

    Tuple(const Tuple& other)
    requires((std::is_copy_constructible_v<Head> && ... && 
        std::is_copy_constructible_v<Tail>)) 
    :head(other.head), tail(other.tail) {}

    Tuple& operator=(const Tuple& other)
    requires(std::conjunction_v<
        std::is_copy_assignable<Head>,
        std::is_copy_assignable<Tail>...
    >)  
    {
        head = other.head;
        tail = other.tail;
        return *this;
    }

    Tuple& operator=(Tuple&& other)
    requires(std::conjunction_v<
        std::is_move_assignable<Head>,
        std::is_move_assignable<Tail>... 
    >)
    {  
        head = std::forward<Head>(other.head);
        tail = std::forward<Tuple<Tail...>>(other.tail);
        return *this;
    }
    
    template<typename Uhead, typename... Uother, 
    std::enable_if_t<std::conjunction_v<std::is_assignable<Head&, const Uhead&>, 
    std::is_assignable<Tail&, const Uother&>...> && sizeof...(Uother) == sizeof...(Tail), bool> = true>
    Tuple& operator=(const Tuple<Uhead, Uother...>& other) {
        head = other.head;
        tail = other.tail;
        return *this; 
    }

    template<typename Uhead, typename... Uother, 
    std::enable_if_t<std::conjunction_v<std::is_assignable<Head&, const Uhead&>, 
    std::is_assignable<Tail&, const Uother&>...> && sizeof...(Uother) == sizeof...(Tail), bool> = true>
    Tuple& operator=(Tuple<Uhead, Uother...>&& other) {
        head = std::forward<Uhead>(other.head);
        tail = std::forward<Tuple<Uother...>>(other.tail);
        return *this; 
    }

    template<typename T1, typename T2>
    Tuple(const std::pair<T1, T2>& other) : head(other.first), tail(other.second) {}

    template<typename T1, typename T2>
    Tuple(std::pair<T1, T2>&& other) : head(std::move(other.first)), tail(std::move(other.second)) {}
};

template<>
class Tuple<> {};

template<size_t N, typename... Ts>
typename GetTypeByIndex<N, Ts...>::type& get(Tuple<Ts...>& t)
{
    if constexpr(N == 0) {
        return t.head;
    }
    else
    {
        return get<N-1>(t.tail);
    }
}

template<typename T, typename... Ts>
const T& get(const Tuple<Ts...>& t)
{
    static_assert(contains_elem_v<T, Ts...>);
    if constexpr(std::is_same_v<T, decltype(t.head)>) {
        return t.head;
    }
    else
    {
        return get<T>(t.tail);
    }
}
template<typename T, typename... Ts>
T& get(Tuple<Ts...>& t)
{
    static_assert(contains_elem_v<T, Ts...>);
    if constexpr(std::is_same_v<T, decltype(t.head)>) {
        return t.head;
    }
    else
    {
        return get<T>(t.tail);
    }
}

template<typename T, typename... Ts>
T&& get(Tuple<Ts...>&& t)
{
    static_assert(contains_elem_v<T, Ts...>);
    if constexpr(std::is_same_v<T, decltype(t.head)>) {
        return t.head;
    }
    else
    {
        return get<T>(std::move(t.tail));
    }
}

template<typename T, typename... Ts>
const T&& get(const Tuple<Ts...>&& t)
{
    static_assert(contains_elem_v<T, Ts...>);
    if constexpr(std::is_same_v<T, decltype(t.head)>) {
        return t.head;
    }
    else
    {
        return get<T>(std::move(t.tail));
    }
}

template<typename T1, typename T2>
Tuple(const std::pair<T1, T2>&) -> Tuple<T1, T2>;

template<typename T1, typename T2>
Tuple(std::pair<T1, T2>&&) -> Tuple<T1, T2>;

template <typename... Ts>
Tuple<Ts...> makeTuple(Ts&&... ts) noexcept { return Tuple<Ts...>(std::forward<Ts>(ts)...); }

template <typename... Ts>
Tuple<Ts&...> tie(Ts&... ts) noexcept { return {ts...}; }

template <typename... Ts>
Tuple<Ts&&...> forwardAsTuple(Ts&&... ts) noexcept { return Tuple<Ts&&...>(std::move(ts)...); };

template<typename... Tuples>
auto tupleCat(Tuples&&... tuples)
-> typename ConcatTypes<std::remove_reference_t<Tuples>...>::type 
{
    return typename ConcatTypes<std::remove_reference_t<Tuples>...>::type (cat_tag(), tuples...);
}