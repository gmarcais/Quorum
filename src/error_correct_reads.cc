/*  This file is part of k_unitig.

    k_unitig is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    k_unitig is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with k_unitig.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <memory>

#ifndef typeof
#define typeof __typeof__
#endif

//#define DEBUG 1
#include <jellyfish/dbg.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/mer_counting.hpp>

#include <jflib/multiplexed_io.hpp>
#include <src/error_correct_reads.hpp>
#include <gzip_stream.hpp>

typedef uint64_t hkey_t;
typedef uint64_t hval_t;
//typedef jellyfish::invertible_hash::array<uint64_t,atomic::gcc,allocators::mmap> inv_hash_storage_t;
typedef std::vector<inv_hash_storage_t> hashes_t;

std::ostream &operator<<(std::ostream &os, const forward_counter &c) {
  return os << c._c;
}
std::ostream &operator<<(std::ostream &os, const backward_counter &c) {
  return os << c._c;
}

// Base class which gives all possible continuation given a k-mer and
// Jellyfish databases. One version uses many Jellyfish databases
// while the other uses a combined database.
//
// The best alternatives version gets all the possible continuation at
// the highest possible level. The best level is recorded. The
// get_alternatives version gets all the possible continuation at the
// currently recorded level.
class alternative_finder { 
public:
  virtual ~alternative_finder() { }

  // Get count at highest level for a k-mer. If the k-mer is not
  // found, 0 is returned
  virtual hval_t get_val(uint64_t m) = 0;
  virtual int get_best_alternatives(forward_mer& m, uint64_t counts[], uint64_t &ucode, int& level) = 0;
  virtual int get_best_alternatives(backward_mer& m, uint64_t counts[], uint64_t &ucode, int& level) = 0;
};

class alternative_multiple_dbs : public alternative_finder {
  const hashes_t*          hashes_;
  int                      min_count_;

public:
  alternative_multiple_dbs(const hashes_t* hashes, int min_count) :
    hashes_(hashes), min_count_(min_count) { }

  virtual ~alternative_multiple_dbs() { std::cerr << __PRETTY_FUNCTION__ << "\n"; }

  virtual hval_t get_val(uint64_t mer) {
    hashes_t::const_iterator chash = hashes_->begin();
    hval_t                   res   = 0;

    if(!chash->get_val(mer, res, true))
      return 0;
    return res;
  }
  virtual int get_best_alternatives(forward_mer& m, uint64_t counts[], uint64_t &ucode, int& level) {
    return get_best_alternatives__(m, counts, ucode, level);
  }
  virtual int get_best_alternatives(backward_mer& m, uint64_t counts[], uint64_t &ucode, int& level) {
    return get_best_alternatives__(m, counts, ucode, level);
  }

private:  
  template<typename dir_mer>
  int get_best_alternatives__(const dir_mer &mer, uint64_t counts[], uint64_t &ucode, int& level) {
    hashes_t::const_iterator chash = hashes_->begin();
    int                      count = 0;
    while(count == 0) {
      count = get_alternatives__(chash, mer, counts, ucode);
      if(count == 0 && ++chash == hashes_->end())
        return 0;
    }
    level = (hashes_->end() - chash) - 1 ;
    return count;
  }

  template <typename dir_mer>
  int get_alternatives__(hashes_t::const_iterator& chash, const dir_mer &mer, uint64_t counts[], uint64_t &ucode) {
    dir_mer  nmer(mer);
    int      count = 0;
    DBG << V(mer);
    for(uint64_t i = 0; i < (uint64_t)4; ++i) {
      nmer.replace(0, i);
      hval_t val = 0;
      if(chash->get_val(nmer.canonical(), val, true)) {
        counts[i] = val;
        if(val >= (uint64_t)min_count_) {
          count++;
          ucode = i;
        }
      } else {
        counts[i] = 0;
      }
    }
    return count;
  }
};

class alternative_combined_dbs : public alternative_finder {
  hashes_t::const_iterator hash_;
  const int                nb_levels_;
  int                      min_count_;

public:
  alternative_combined_dbs(const hashes_t* hashes, int levels,  int min_count) :
    hash_(hashes->begin()), nb_levels_(levels), min_count_(min_count) { }

  virtual ~alternative_combined_dbs() { }
  virtual hval_t get_val(uint64_t mer) {
    hval_t res = 0;
    bool found = hash_->get_val(mer, res, true, true);
    DBG << V(forward_mer(mer)) << V(found) << V(res) << V(res % nb_levels_) << V(res / nb_levels_);
    if(!found)
      return 0;
    if(res % nb_levels_ != (hval_t)(nb_levels_ - 1))
      return 0;
    return res / nb_levels_;
  }
  virtual int get_best_alternatives(forward_mer& m, uint64_t counts[], uint64_t& ucode, int& level) {
    return get_best_alternatives__(m, counts, ucode, level);
  }
  virtual int get_best_alternatives(backward_mer& m, uint64_t counts[], uint64_t& ucode, int& level) {
    return get_best_alternatives__(m, counts, ucode, level);
  }

private:
  template<typename dir_mer>
  int get_best_alternatives__(dir_mer& mer, uint64_t counts[], uint64_t& ucode, int& level) {
    int      nlevel;
    uint64_t val;
    dir_mer  nmer(mer);
    int      count = 0;
    level = 0;

    for(uint64_t i = 0; i < (uint64_t)4; ++i) {
      nmer.replace(0, i);
      if(!hash_->get_val(nmer.canonical(), val, true, true))
        val = 0;
      nlevel = val % nb_levels_;
      val    = val / nb_levels_;
      if(val == 0 || nlevel < level) {
        counts[i] = 0;
      } else {
        if(val >= (uint64_t)min_count_) {
          if(nlevel > level) {
            for(uint64_t j = 0; j < (uint64_t)i; ++j)
              counts[j] = 0;
            count  = 0;
            level = nlevel;
          }
          counts[i] = val;
          ucode     = i;
          ++count;
        }
      }
    }
    return count;
  }
};

// Contaminant database. If a Jellyfish database is given, return true
// iff the k-mer is in the database. With no database, it always
// returns false.
class contaminant_check {
public:
  virtual ~contaminant_check() { }
  
  virtual bool is_contaminant(uint64_t m) = 0;
  virtual void debug(const char*msg) = 0;
};

class contaminant_no_database : public contaminant_check {
public:
  virtual ~contaminant_no_database() { }
  virtual bool is_contaminant(uint64_t m) { return false; }
  virtual void debug(const char* msg) { std::cerr << msg << " no database" << std::endl; }
};

class contaminant_database : public contaminant_check {
  inv_hash_storage_t ary_;
public:
  contaminant_database(inv_hash_storage_t ary) : ary_(ary) { }
  virtual ~contaminant_database() { }
  virtual void debug(const char* msg) { std::cerr << msg << " block_len " << ary_.get_block_len() << std::endl; }
  virtual bool is_contaminant(uint64_t m) {
    hval_t val = 0;
    return ary_.get_val(m, val, false);
  }
};

template<class instance_t>
class error_correct_t : public thread_exec {
  jellyfish::parse_read* _parser;
  hashes_t*              _hashes;
  alternative_finder*    _af;
  int                    _mer_len;
  int                    _skip;
  int                    _good;
  int                    _anchor;
  std::string            _prefix;
  int                    _min_count;
  int                    _window;
  int                    _error;
  bool                   _gzip;
  int                    _combined;
  contaminant_check*     _contaminant;
  bool                   _trim_contaminant;
  jflib::o_multiplexer * _output;
  jflib::o_multiplexer * _log;
public:
  error_correct_t(jellyfish::parse_read* parser, hashes_t *hashes) :
    _parser(parser), _hashes(hashes),
    _mer_len(_hashes->begin()->get_key_len() / 2),
    _skip(0), _good(1), _min_count(1), _window(0), _error(0), _gzip(false),
    _combined(0), _contaminant(0), _trim_contaminant(false) { }

private:
  std::ostream *open_file(const char *suffix) {
    std::ostream *res;
    std::string file(_prefix);
    file += suffix;
    if(_gzip) {
      file += ".gz";
      res = new gzipstream(file.c_str());
    } else {
      res = new std::ofstream(file.c_str());
    }
      
    if(!res->good())
      eraise(std::runtime_error)
        << "Failed to open file '" << file << "'" << err::no;
    res->exceptions(std::ios::eofbit|std::ios::failbit|std::ios::badbit);
    return res;
  }

public:
  void do_it(int nb_threads) {
    std::auto_ptr<std::ostream> details(open_file(".log"));
    std::auto_ptr<std::ostream> output(open_file(".fa"));

    std::auto_ptr<jflib::o_multiplexer> 
      log_m(new jflib::o_multiplexer(details.get(), 3 * nb_threads, 1024));
    std::auto_ptr<jflib::o_multiplexer>
      output_m(new jflib::o_multiplexer(output.get(), 3 * nb_threads, 1024));
    _log    = log_m.get();
    _output = output_m.get();

    exec_join(nb_threads);
  }
  
  virtual void start(int id) {
    instance_t(this, id).start();
  }

  error_correct_t & skip(int s) { _skip = s; return *this; }
  error_correct_t & good(int g) { _good = g; return *this; }
  error_correct_t & anchor(int a) { _anchor = a; return *this; }
  error_correct_t & prefix(const char *s) { _prefix = s; return *this; }
  error_correct_t & prefix(const std::string s) { _prefix = s; return *this; }
  error_correct_t & mer_len(int l) { _mer_len = l; return *this; }
  error_correct_t & min_count(int m) { _min_count = m; return *this; }
  //  error_correct_t & advance(int a) { _advance = a; return *this; }
  error_correct_t & window(int w) { _window = w; return *this; }
  error_correct_t & error(int e) { _error = e; return *this; }
  error_correct_t & gzip(bool g) { _gzip = g; return *this; }
  error_correct_t & combined(int c) { _combined = c; return *this; }
  error_correct_t & contaminant(contaminant_check* c) { _contaminant = c; return *this; }
  error_correct_t & trim_contaminant(bool t) { _trim_contaminant = t; return *this; }

  jellyfish::parse_read* parser() const { return _parser; }
  int skip() const { return _skip; }
  int good() const { return _good; }
  int anchor() const { return _anchor; }
  const std::string & prefix() const { return _prefix; }
  int mer_len() const { return _mer_len; }
  int min_count() const { return _min_count; }
  //  int advance() const { return _advance; }
  int window() const { return _window ? _window : _mer_len; }
  int error() const { return _error ? _error : _mer_len / 2; }
  bool gzip() const { return _gzip; }
  int combined() const { return _combined; }
  contaminant_check* contaminant() const { return _contaminant; }
  bool trim_contaminant() const { return _trim_contaminant; }

  jflib::o_multiplexer &output() { return *_output; }
  jflib::o_multiplexer &log() { return *_log; }

  alternative_finder* new_af() {
    if(_combined == 0)
      return new alternative_multiple_dbs(_hashes, _min_count);
    else
      return new alternative_combined_dbs(_hashes, _combined, _min_count);
  }
};

class error_correct_instance {
public:
  typedef error_correct_t<error_correct_instance> ec_t ;

private:
  ec_t                     *_ec;
  int                       _id;
  size_t                    _buff_size;
  char                     *_buffer;
  alternative_finder*       _af;

  static const char* error_contaminant;
  static const char* error_no_starting_mer;
  
public:
  error_correct_instance(ec_t *ec, int id) :
    _ec(ec), _id(id), _buff_size(0), _buffer(0) { }

  void start() {
    jellyfish::parse_read::thread parser = _ec->parser()->new_thread();
    _af = _ec->new_af();
    
    const jellyfish::read_parser::read_t *read;

    jflib::omstream output(_ec->output());
    jflib::omstream details(_ec->log());

    uint64_t nb_reads = 0;
    bool parity = true;
    while((read = parser.next_read())) {
      parity = !parity;
      nb_reads++;
      insure_length_buffer(read->seq_e - read->seq_s);
      
      const char* error = "";
      kmer_t      mer;
      const char *input = read->seq_s + _ec->skip();
      char       *out   = _buffer + _ec->skip();
      //      DBG << V(_ec->skip()) << V((void*)read->seq_s) << V((void*)input);
      //Prime system. Find and write starting k-mer
      if(!find_starting_mer(mer, input, read->seq_e, out, &error)) {
        details << "Skipped " << substr(read->header, read->hlen) 
                << ": " << error << "\n";
        details << jflib::endr;
        output << jflib::endr;
        continue;
      }
      //      DBG << V((void*)read->seq_s) << V((void*)input) << V(kmer_t::k());
      DBG << V(std::string(read->header, 15));
      // Extend forward and backward
      forward_log fwd_log(_ec->window(), _ec->error());
      char *end_out = 
        extend(forward_mer(mer), forward_ptr<const char>(input),
               forward_counter(input - read->seq_s),
               forward_ptr<const char>(read->seq_e),
               forward_ptr<char>(out), fwd_log,
               &error);
      if(!end_out) {
        details << "Skipped " << substr(read->header, read->hlen) 
                << ": " << error << "\n";
        details << jflib::endr;
        output << jflib::endr;
        continue;
      }
      DBG << V((void*)end_out) << V((void*)read->seq_e);
      assert(input > read->seq_s + kmer_t::k());
      assert(out > _buffer + kmer_t::k());
      assert(input - read->seq_s == out - _buffer);
      backward_log bwd_log(_ec->window(), _ec->error());
      char *start_out =
        extend(backward_mer(mer), 
               backward_ptr<const char>(input - kmer_t::k() - 1),
               backward_counter(input - kmer_t::k() - read->seq_s - 1),
               backward_ptr<const char>(read->seq_s - 1),
               backward_ptr<char>(out - kmer_t::k() - 1), bwd_log,
               &error);
      if(!start_out) {
        details << "Skipped " << substr(read->header, read->hlen) 
                << ": " << error << "\n";
        details << jflib::endr;
        output << jflib::endr;
        continue;
      }
      DBG << V((void*)start_out) << V((void*)read->seq_s);
      start_out++;
      assert(start_out >= _buffer);
      assert(_buffer + _buff_size >= end_out);

      output << ">" << substr(read->header, read->hlen) 
             << " " << fwd_log << " " << bwd_log << "\n"
             << substr(start_out, end_out) << "\n";
      if(parity)
        output << jflib::endr;
    }
    details.close();
    output.close();
  }

private:

  // Extend and correct read. Copy from input to out. mer should be
  // represent a "good" starting k-mer at the input position.
  // out point to the next character to be written.
  template<typename dir_mer, typename in_dir_ptr, typename out_dir_ptr,
           typename counter, typename elog>
  char * extend(dir_mer mer, in_dir_ptr input, 
                counter pos, in_dir_ptr end,
                out_dir_ptr out, elog &log, const char** error) {
    counter cpos = pos;
    DBG << V((void*)input.ptr()) << V((void*)end.ptr()) << V(cpos);
    for( ; input < end; ++input) {
      char     base        = *input;
      DBG << V((void*)input.ptr()) << V((void*)end.ptr()) << V(base);
      if(base == '\n')
        continue;
      cpos = pos;
      ++pos;

      uint64_t ori_code;
      if(!mer.shift(base)) {
        ori_code = 5; // Invalid base
        mer.shift((uint64_t)0);
      } else {
        if(_ec->contaminant()->is_contaminant(mer.canonical())) {
          if(_ec->trim_contaminant()) {
            log.truncation(cpos);
            goto done;
          }
          *error = error_contaminant;
          return 0;
        }
        ori_code = mer.code(0);
      }
      uint64_t counts[4];
      uint64_t ucode = 0;
      int      count;
      int      level;

      count = _af->get_best_alternatives(mer, counts, ucode, level);
      DBG << V(*cpos) << V(mer) << V(count) << V(level) << V(counts[0]) << V(counts[1]) << V(counts[2]) << V(counts[3]);

      if(count == 0) {
        log.truncation(cpos);
        goto done;
      }
      if(count == 1) { // One continuation. Is it an error?
        if(ucode != ori_code) {
          mer.replace(0, ucode);
          if(_ec->contaminant()->is_contaminant(mer.canonical())) {
            if(_ec->trim_contaminant()) {
              log.truncation(cpos);
              goto done;
            }
            *error = error_contaminant;
            return 0;
          }
          if(log.substitution(cpos, base, mer.base(0)))
            goto truncate;
        }
        *out++ = mer.base(0);
        continue;
      }

      // Check that there is at least one more base in the
      // sequence. If not, leave it along
      if(input >= end) {
        log.truncation(cpos);
        goto done;
      }

      // Select the replacement base to try. Find the one with highest
      // coverage and check that it is at least 3 times as big as
      // current base count and that it can continue one more base. In
      // that case, make the switch.
      uint64_t check_code = ori_code;
      uint64_t ori_count  = ori_code < 4 ? counts[ori_code] : 0;
      if(ori_count < (uint64_t)_ec->min_count())
        ori_count = 0;

      uint64_t max_count = 10000000000;
      if(ori_count <= (uint64_t)_ec->anchor())
        max_count  = 3 * ori_count;

      for(int i = 0; i < 4; i++) {
        if(counts[i] < (uint64_t)_ec->min_count())
            continue;
        if(counts[i] > max_count && counts[i]-ori_count<200) {
          check_code = i;
          max_count  = counts[i];
        }
      }
      if(check_code == ori_code) {
        // Don't need to check that check_code == 5 as an alternative
        // would have been found by now.
        *out++ = base;
        continue;
      }

      // Check that it continues at least one more base with that quality
      dir_mer    nmer   = mer;
      //      in_dir_ptr ninput = input;
      nmer.replace(0, check_code);
      // char       nbase  = *ninput++;
      // Does not matter what we shift, check all alternative anyway.
      nmer.shift((uint64_t)0);

      uint64_t   ncounts[4];
      int        ncount;
      uint64_t   nucode = 0;
      int        nlevel;
      ncount = _af->get_best_alternatives(nmer, ncounts, nucode, nlevel);
      DBG << V(*cpos) << V(ncount) << V(nlevel) << V(level) << V(ncounts[0]) << V(ncounts[1]) << V(ncounts[2]) << V(ncounts[3]);
      if(ncount > 0 && nlevel >= level) { // TODO: Shouldn't we break if this test is false?
        mer.replace(0, check_code);
        if(_ec->contaminant()->is_contaminant(mer.canonical())) {
          if(_ec->trim_contaminant()) {
            log.truncation(cpos);
            goto done;
          }
          *error = error_contaminant;
          return 0;
        }
        *out++ = mer.base(0);
        if(check_code != ori_code)
          if(log.substitution(cpos, base, mer.base(0)))
            goto truncate;
        // if(ncount == 1) { // While we are at it, there is a uniq continuation
        //   mer    = nmer;
        //   input  = ninput;
        //   ++pos;
        //   if(nucode != mer.code(0)) {
        //     mer.replace(0, nucode);
        //     if(log.substitution(cpos, nbase, mer.base(0)))
        //       goto truncate;
        //   }
        //   *out++ = mer.base(0);
        // }
      }
    }
    
  done:
    return out.ptr();

  truncate:
    int diff = log.remove_last_window();
    out = out - diff;
    DBG << V(*cpos) << V(diff) << V(*(cpos - diff));
    log.truncation(cpos - diff);
    goto done;
  }

  void insure_length_buffer(size_t len) {
    if(len > _buff_size) {
      _buff_size = len > 2 * _buff_size ? len + 100 : 2 * _buff_size;
      _buffer = (char *)realloc(_buffer, _buff_size);

      if(!_buffer)
        eraise(std::runtime_error)
          << "Buffer allocation failed, size " << _buffer << err::no;
    }
  }

  bool find_starting_mer(kmer_t &mer, const char * &input, const char *end, char * &out,
                         const char** error) {
    while(input < end) {
      for(int i = 0; input < end && i < _ec->mer_len(); ++i) {
        char base = *input++;
        *out++ = base;
        DBG << V(base) << V(mer);
        if(!mer.shift_left(base))
          i = -1;        // If an N, skip to next k-mer
      }
      int found = 0;
      while(input < end) {
        bool contaminated = _ec->contaminant()->is_contaminant(mer.canonical());
        if(contaminated && !_ec->trim_contaminant()) {
          *error = error_contaminant;
          return false;
        }

        if(!contaminated) {
          hval_t val = _af->get_val(mer.canonical());
          
          found = (int)val >= _ec->anchor() ? found + 1 : 0;
          DBG << V(val) << V(mer) << V(_ec->anchor()) << V(*input) << V(found);
          if(found >= _ec->good())
            return true;
        }

        char base = *input++;
        *out++ = base;
        if(!mer.shift_left(base))
          break;
      }
    }

    *error = error_no_starting_mer;
    return false;
  }
};

const char* error_correct_instance::error_contaminant     = "Contaminated read";
const char* error_correct_instance::error_no_starting_mer = "No high quality mer";


int main(int argc, char *argv[])
{
  args_t args(argc, argv);

  if(args.combined_given && args.db_arg.size() > 1)
    die << "Only one Jellyfish database when using combined mode";

  // Open Jellyfish databases
  hashes_t hashes;
  unsigned int key_len = 0;
  for(auto it = args.db_arg.begin(); it != args.db_arg.end(); ++it) {
    mapped_file dbf(*it);
    dbf.random().will_need().load();
    hashes.push_back(*raw_inv_hash_query_t(dbf).get_ary());
    
    if(key_len == 0)
      key_len = hashes.front().get_key_len();
    else if(key_len != hashes.back().get_key_len())
      die << "Different key length (" << hashes.back().get_key_len() 
          << " != " << key_len
          << ") for hash '" << *it << "'";
  }

  // Open contaminant database
  contaminant_check* contaminant = 0;
  if(args.contaminant_given) {
    mapped_file dbf(args.contaminant_arg);
    dbf.random().will_need().load();
    inv_hash_storage_t ary = *raw_inv_hash_query_t(dbf).get_ary();
    if(ary.get_key_len() != key_len)
      die << "Contaminant hash must have same key length as other hashes ("
          << ary.get_key_len() << " != " << key_len << ")";
    contaminant = new contaminant_database(ary);
  } else {
    contaminant = new contaminant_no_database();
  }
  jellyfish::parse_read parser(args.file_arg.begin(), args.file_arg.end(), 100);

  kmer_t::k(key_len / 2);
  error_correct_instance::ec_t correct(&parser, &hashes);
  correct.skip(args.skip_arg).good(args.good_arg)
    .anchor(args.anchor_count_given ? args.anchor_count_arg : args.min_count_arg)
    .prefix(args.output_arg).min_count(args.min_count_arg)
    .window(args.window_given ? args.window_arg : kmer_t::k())
    .error(args.error_given ? args.error_arg : kmer_t::k() / 2)
    .gzip(args.gzip_flag)
    .combined(args.combined_arg)
    .contaminant(contaminant)
    .trim_contaminant(args.trim_contaminant_flag);
  correct.do_it(args.thread_arg);

  return 0;
}
