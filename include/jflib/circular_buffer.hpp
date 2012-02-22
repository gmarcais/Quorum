#ifndef _JFLIB_circular_buffer_H_
#define _JFLIB_circular_buffer_H_

#include <jflib/compare_and_swap.hpp>
#include <iostream>
#include <limits.h>

namespace jflib {
  template<typename T, unsigned int n, T g = ((T)1 << n) - 1>
  class basic_circular_buffer {
    static const unsigned int m = 8*sizeof(T)-n;
    struct splitT {
      T id:m;
      T val:n;
    };
    union elt {
      T      binary;
      splitT split;
    };
    size_t  _size;
    T      *_buffer;
    size_t  _head;
    size_t  _tail;
    bool    _closed;

  public:
    typedef T value_type;
    static const T guard = g;

    basic_circular_buffer(size_t size) :
      _size(size+1), _buffer(new T[_size]),  _head(0), _tail(0), _closed(false)
    { 
      elt init;
      init.split.id  = 0;
      init.split.val = guard;
      for(size_t i = 0; i < _size; ++i)
        _buffer[i] = init.binary;
    }
    virtual ~basic_circular_buffer() { }

    /** Enqueue an element.
     * @return false if the FIFO is full.
     */
    bool enqueue(const T &v);
    /** Dequeue an element.
     * @return 0 if the FIFO is empty.
     */
    T dequeue();
    bool is_closed() const { return a_load(_closed); }
    void close() { a_store(_closed, true); }
  };

#define BSIZEOF(x) (sizeof(T) * CHAR_BIT)
  template<typename T, T g = (T)-1>
  class circular_buffer : public basic_circular_buffer<uint64_t, BSIZEOF(T), g> {
  public:
    circular_buffer(size_t size) :
      basic_circular_buffer<uint64_t, BSIZEOF(T), g>(size) { }
    virtual ~circular_buffer() { }

    bool enqueue(const T &v) { 
      return basic_circular_buffer<uint64_t, BSIZEOF(T), g>::enqueue((uint64_t)v);
    }
    T dequeue() {
      return basic_circular_buffer<uint64_t, BSIZEOF(T), g>::dequeue();
    }
  };
}

template<typename T, unsigned int n, T guard>
bool jflib::basic_circular_buffer<T,n,guard>::enqueue(const T &v) {
  bool done = false;
 
  size_t chead = a_load(_head);
  while(!done) {
    size_t nhead = chead + 1;
    if(nhead % _size == a_load(_tail) % _size)
      return false;
    
    elt celt;
    celt.binary = a_load(_buffer[chead % _size]);
    size_t achead = a_load(_head);
    if(achead != chead) {
      chead = achead;
      continue;
    }
    if(celt.split.val == guard) {
      // entry is empty
      elt nelt;
      nelt.split.id  = celt.split.id + 1;
      nelt.split.val = v;
      done = cas(&_buffer[chead % _size], celt.binary, nelt.binary);
      // done == true <=> sucessfully written entry
    }
    cas(&_head, chead, nhead, &chead);
  }
 
  return true;
}

template<typename T, unsigned int n, T guard>
T jflib::basic_circular_buffer<T,n,guard>::dequeue() {
  bool done = false;
  elt res;

  size_t ctail = a_load(_tail);
  while(!done) {
    bool dequeued = false;
    do {
      if(ctail % _size == a_load(_head) % _size)
        return guard;
      size_t ntail = ctail + 1;
      dequeued = cas(&_tail, ctail, ntail, &ctail);
    } while(!dequeued);

    res.binary = a_load(_buffer[ctail % _size]);
    elt nres;
    nres.split.val = guard;
    while(true) {
      nres.split.id = res.split.id + 1;
      if(res.split.val == guard) {
        if(cas(&_buffer[ctail % _size], res.binary, nres.binary, &res.binary))
          break;
      } else {
        done = cas(&_buffer[ctail % _size], res.binary, nres.binary);
        break;
      }
    }
  }

  return res.split.val;
}
#endif /* _JFLIB_circular_buffer_H_ */
