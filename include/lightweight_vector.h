#ifndef LIGHTWEIGHTVECTOR_H
#define LIGHTWEIGHTVECTOR_H


//#include <std/shared_ptr>
#include <vector>

template <typename T, typename Allocator = std::allocator<T>>
struct lightweight_vector
{
  typedef T           value_type;
  typedef std::size_t size_type;

  typedef std::vector<T,Allocator> content_t;
private:
  typedef std::shared_ptr<content_t> storage_t;
  storage_t m_content;

public:
  typedef typename content_t::reference       reference;
  typedef typename content_t::const_reference const_reference;
  typedef typename content_t::iterator        iterator;
  typedef typename content_t::const_iterator  const_iterator;

  lightweight_vector(std::size_t size)
    : m_content(new content_t(size))
  { }

  lightweight_vector()
    : m_content(new content_t())
  { }

  lightweight_vector(const_iterator begin, const_iterator end)
    : m_content(new content_t(begin,end))
  { }


  size_t size(void) const
  {
    return m_content->size();
  }

  //////////////////////////////////////////////////////////////////////
  /// @brief Return an const_iterator to the beginning of
  /// the sequence of elements this container stores.
  //////////////////////////////////////////////////////////////////////
  const_iterator begin(void) const
  {
    return m_content->begin();
  }

  //////////////////////////////////////////////////////////////////////
  /// @brief Return an iterator to the beginning of
  /// the sequence of elements this container stores.
  //////////////////////////////////////////////////////////////////////
  iterator begin(void)
  {
    return m_content->begin();
  }

  //////////////////////////////////////////////////////////////////////
  /// @brief Return an const_iterator to the end (one past the last element) of
  /// the sequence of elements this container stores.
  //////////////////////////////////////////////////////////////////////
  const_iterator end(void) const
  {
    return m_content->end();
  }

  //////////////////////////////////////////////////////////////////////
  /// @brief Return an iterator to the end (one past the last element) of
  /// the sequence of elements this container stores.
  //////////////////////////////////////////////////////////////////////
  iterator end(void)
  {
    return m_content->end();
  }

  //////////////////////////////////////////////////////////////////////
  /// @brief Return reference to element in buffer at index @p idx.
  //////////////////////////////////////////////////////////////////////
  reference operator[](std::size_t idx)
  {
    return m_content->operator[](idx);
  }

  //////////////////////////////////////////////////////////////////////
  /// @brief Return reference to element in buffer at index @p idx.
  //////////////////////////////////////////////////////////////////////
  const_reference operator[](std::size_t idx) const
  {
    return m_content->operator[](idx);
  }

  bool operator==(lightweight_vector<T> const& other) const
  {
    return (m_content==other.m_content);
  }

  void push_back(T const& t)
  {
    m_content->push_back(t);
  }

  template <class... Args>
  void emplace_back(Args&&... args) 
  {
    m_content->emplace_back(args...);
  }

  void reserve(std::size_t size)
  {
    m_content->reserve(size);
  }

  void resize(std::size_t size)
  {
    m_content->resize(size);
  }

  bool empty() const
  {
    return m_content->empty();
  }

  const_reference front(void) const
  {
    return m_content->front();
  }

  reference front(void)
  {
    return m_content->front();
  }

  reference back(void) 
  {
    return m_content->back();
  }

  const_reference back(void) const 
  {
    return m_content->back();
  }

  void shrink_to_fit() 
  {
    m_content->shrink_to_fit();
  }

  template <class InputIterator>
  void insert (iterator position, InputIterator first, InputIterator last)
  {
    m_content->insert(position,first,last);
  }

  void insert (iterator position, size_type n, const value_type& val)
  {
    m_content->insert(position,n,val);
  }

  iterator insert (iterator position, const value_type& val)
  {
    return m_content->insert(position,val);
  }

  content_t* get()
  {
    return m_content.get();
  }

  void clear()
  {
    m_content->clear();
  }
  
  template <class InputIterator>
  void assign (InputIterator first, InputIterator last)
  {
    m_content->assign(first,last);
  }

  void assign (size_type n, const value_type& val)
  {
    m_content->assign(n,val);
  }
};

#endif