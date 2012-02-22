#ifndef _JFLIB_MULTIPLEXED_IO_H_
#define _JFLIB_MULTIPLEXED_IO_H_

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ostream>
#include <stdexcept>
#include <jflib/pool.hpp>

/** Multiplexed io stream. Multiple threads write to a pool of buffers
 * which are then written by a single thread to an ostream. Similar to
 * asynchronous IO, potentially easier to use.
 *
 * This works in 2 steps: first create the o_multiplexer based on an
 * std::ostream. This will setup the pool. Then, each thread creates
 * an omstream based on the o_multiplexer. Each of this omstream
 * obtain and send buffers to the ostream using the pool. To ensure
 * that the output is not interleaved with other threads output, there
 * is a end_record() (or << endr) method which marks the end of a
 * record that can be send to the underlying ostream and written as
 * one atomic operation. A sync() operation (or << std::flush) implies
 * an end of record as well.
 */

namespace jflib {
  class o_multiplexer {
    class                buffer;
    class                multiplexedbuf;
    typedef pool<buffer> buffer_pool;
    struct writer_info {
      std::ostream *os;
      buffer_pool  *pool;
      writer_info(std::ostream *_os, buffer_pool *_pool) : os(_os), pool(_pool) { }
    };
  public:
    /** Constructor based on an ostream os and nb_buffers of size
     * buffer_size. For reliable performance, there should be more
     * nb_buffers than threads that will create omstream, otherwise
     * there will always be threads waiting on getting a
     * buffer. Probably, 2 to 3 times the number of threads.
     *
     * The size of a buffer should be larger than the record
     * length. If not, the buffer get enlarged (size doubled) when a
     * record does not fit into a buffer.
     */
    o_multiplexer(std::ostream *os, size_t nb_buffers, size_t buffer_size)
      : _os(os), _pool(nb_buffers), _info(_os, &_pool), _writer_started(false)
    {
      for(buffer_pool::iterator it = _pool.begin(); it != _pool.end();
          ++it) {
        if(!it->resize(buffer_size))
          throw std::runtime_error("Failed to allocate memory for buffers");
      }
      int res = pthread_create(&_writer_id, 0, writer_loop, (void*)&_info);
      if(res)
        throw std::runtime_error("Failed to create the writer thread");
      _writer_started = true;
      
    }
    /** Destructor. Sync is called on the underlying ostream, but it
     * is otherwise not changed (not closed or destroyed) by this
     * destructor.
     */
    virtual ~o_multiplexer() { close(); }
    /** Tell the writer thread to stop after the pool is empty. To be
     * safe, this should be called after every threads with an
     * omstream() associated have are done writing to the
     * stream. Unspecified behavior if a thread writes after close has
     * been called.
     *
     * TODO: should we implement a safe version that waits for every
     * omstream object to be destroyed?
     */
    void close() {
      _pool.close_A_to_B();
      if(_writer_started) {
        pthread_join(_writer_id, 0);
        _writer_started = false;
      }
      _os->flush();
    }

  private:
    static void *writer_loop(void *_d) {
      writer_info *d = (writer_info*)_d;

      while(true) {
        buffer_pool::elt e(d->pool->get_B());
        if(e.is_empty())
          break;
        if(e->write_len() > 0)
          d->os->write(e->base(), e->write_len());
        if(e->do_sync())
          d->os->flush();
      }
      d->os->flush();
      return 0;
    }

    class buffer {
    public:
      buffer() : _write_len(0), _sync(false), _data(0), _end(0) { }
      buffer(size_t size) 
        : _write_len(0), _sync(false), _data((char*)malloc(size)), _end(_data + size)
      {
        if(!_data)
          throw std::runtime_error("Failed to allocate multiplexer::buffer");
      }
      ~buffer() { free(_data); }
      char *base() const { return _data; }
      char *end() const { return _end; }
      std::streamsize write_len() const { return _write_len; }
      void set_write_len(std::streamsize s) { _write_len = s; }
      bool do_sync() const { return _sync; }
      void set_do_sync(bool sync) { _sync = sync; }
      char *resize(size_t size = 0) {
        size_t csize = _end - _data;
        if(!size)
          size = 2 * csize;
        if(size <= csize)
          return _data;
        char *ndata = (char*)realloc(_data, size);
        if(!ndata)
          return 0;
        _data = ndata;
        _end  = _data + size;
        return _data;
      }
    private:
      std::streamsize  _write_len;
      bool             _sync;
      char            *_data;
      char            *_end;
    };

    class basebuf : public std::streambuf {
    protected:
      basebuf() { }
    public:
      virtual ~basebuf() { }
      virtual void end_record() { }
    };

    class closedbuf : public basebuf {
    public:
      closedbuf() { }
      virtual ~closedbuf() { }
    };

    class multiplexedbuf : public basebuf {
    public:
      multiplexedbuf(buffer_pool *pool) : _pool(pool), _closed(false)
      { 
        update_buffer();
      }
      virtual ~multiplexedbuf() { }

      /** Send a sync signal. A sync implies the end of a record.
       */
      virtual int sync() {
        if(_closed)
          return -1;
        size_t s = pptr() - pbase();
        if(s == 0)
          return 0;

        _elt->set_write_len(s);
        _elt->set_do_sync(true);
    
        // Release and get new buffer
        if(!update_buffer())
          return -1;
        return 0;
      }

      void end_record() {
        if(_closed)
          return;
        size_t clen = pptr() - pbase();
        if(clen == 0)
          return;
        _elt->set_write_len(pptr() - pbase());
        ++_nbr;
        if(2 * clen >= (size_t)(epptr() - pptr()) * _nbr) {
          _elt->set_do_sync(false);
          update_buffer();
        }
      }

      virtual int overflow(int c) {
        if(_closed)
          return EOF;
        if(_nbr == 0) { // No record in its entirety
          size_t old_size = _elt->end() - _elt->base();
          if(!_elt->resize())
            return EOF;
          setp(_elt->base(), _elt->end());
          pbump(old_size);
        } else { // Save part not in a record and release rest.
          ptrdiff_t in_buf = pptr() - pbase();
          if(_elt->write_len() < in_buf) {
            size_t extra = in_buf - _elt->write_len();
            _overflow_buffer.resize(extra);
            memcpy(_overflow_buffer.base(), pptr() - extra, extra);
            _overflow_buffer.set_write_len(extra);
          } else
            _overflow_buffer.set_write_len(0);
          
          _elt->set_do_sync(false);
          if(!update_buffer())
            return EOF;
          if(_overflow_buffer.write_len() > 0) {
            memcpy(_elt->base(), _overflow_buffer.base(),
                   _overflow_buffer.write_len());
            pbump(_overflow_buffer.write_len());
          }
        }

        if(c != EOF) {
          *pptr() = c;
          pbump(1);
        }
        return !EOF;
      }

    private:
      // Release the current buffer and update it with a new one,
      // unless the stream is closed.
      bool update_buffer() {
        if(!_closed)
          _closed = _pool->is_closed_B_to_A();

        if(_closed) {
          _elt.release();
          setp(0, 0);
        } else {
          _elt = _pool->get_A();
          _elt->set_write_len(0);
          setp(_elt->base(), _elt->end());
        }
        _nbr = 0;
        return !_closed;
      }

      buffer_pool      *_pool;
      bool              _closed;
      unsigned int      _nbr; // number of records
      buffer_pool::elt  _elt;
      buffer            _overflow_buffer;
    };

    friend class omstream;
    multiplexedbuf *new_multiplexedbuf() { return new multiplexedbuf(&_pool); }

    std::ostream *_os;
    buffer_pool   _pool;
    writer_info   _info;
    pthread_t     _writer_id;
    bool          _writer_started;
  };

  class omstream : public std::ostream {
  public:
    omstream(o_multiplexer &om) : std::ostream(om.new_multiplexedbuf()) { }
    virtual ~omstream() {
      rdbuf()->pubsync();
      delete rdbuf();
    }
    /** Mark the end of a record
     */
    void end_record() { ((o_multiplexer::basebuf*)rdbuf())->end_record(); }
    void close() { 
      rdbuf()->pubsync();
      delete rdbuf(new o_multiplexer::closedbuf());
    }
  private:
  };
  struct endrT { };
  static const endrT endr = { };
}

jflib::omstream &operator<<(jflib::omstream &os, const jflib::endrT &rhs) {
  os.end_record();
  return os;
}


#endif /* _JFLIB_MULTIPLEXED_IO_H_ */
