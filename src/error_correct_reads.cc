/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <vector>
#include <memory>
#include <limits>
#include <cmath>

#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_dna.hpp>

#include <jflib/multiplexed_io.hpp>
#include <gzip_stream.hpp>

#include <src/mer_database.hpp>
#include <src/error_correct_reads.hpp>
#include <src/error_correct_reads_cmdline.hpp>

using jellyfish::mer_dna;
typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> read_parser;


typedef uint64_t hkey_t;
typedef uint64_t hval_t;


// Contaminant database. If a Jellyfish database is given, return true
// iff the k-mer is in the database. With no database, it always
// returns false.
class contaminant_check {
public:
  virtual ~contaminant_check() { }

  virtual bool is_contaminant(const mer_dna& m) = 0;
  virtual void debug(const char*msg) = 0;
};

class contaminant_no_database : public contaminant_check {
public:
  virtual ~contaminant_no_database() { }
  virtual bool is_contaminant(const mer_dna& m) { return false; }
  virtual void debug(const char* msg) { std::cerr << msg << " no database" << std::endl; }
};

// class contaminant_database : public contaminant_check {
//   inv_hash_storage_t ary_;
// public:
//   contaminant_database(inv_hash_storage_t ary) : ary_(ary) { }
//   virtual ~contaminant_database() { }
//   virtual void debug(const char* msg) { std::cerr << msg << " block_len " << ary_.get_block_len() << std::endl; }
//   virtual bool is_contaminant(uint64_t m) {
//     hval_t val = 0;
//     return ary_.get_val(m, val, false);
//   }
// };

template<class instance_t>
class error_correct_t : public jellyfish::thread_exec {
  read_parser            _parser;
  int                    _mer_len;
  int                    _skip;
  int                    _good;
  int                    _anchor;
  std::string            _prefix;
  int                    _min_count;
  int			 _cutoff;
  int                    _window;
  int                    _error;
  bool                   _gzip;
  database_query*        _mer_database;
  contaminant_check*     _contaminant;
  bool                   _trim_contaminant;
  int                    _homo_trim;
  jflib::o_multiplexer * _output;
  jflib::o_multiplexer * _log;
public:
  error_correct_t(int nb_threads, stream_manager& streams) :
    _parser(4 * nb_threads, 100, 1, streams),
    _mer_len(mer_dna::k()),
    _skip(0), _good(1), _min_count(1), _cutoff(4), _window(0), _error(0), _gzip(false),
    _mer_database(0), _contaminant(0), _trim_contaminant(false),
    _homo_trim(std::numeric_limits<int>::min()) { }

private:
  // Open the data (error corrected reads) and log files. Default to
  // STDOUT and STDERR if none specified.
  std::ostream* open_file(const std::string prefix, const char* suffix,
                          const std::string def) {
    std::ostream* res;
    std::string file;

    if(prefix.empty())
      file = def;
    else {
      file = prefix;
      file += suffix;
    }
    if(_gzip) {
      if(!prefix.empty())
        file += ".gz";
      res = new gzipstream(file.c_str());
    } else
      res = new std::ofstream(file.c_str());
    if(!res->good())
      eraise(std::runtime_error)
        << "Failed to open file '" << file << "'" << jellyfish::err::no;
    res->exceptions(std::ios::eofbit|std::ios::failbit|std::ios::badbit);
    return res;
  }

public:
  void do_it(int nb_threads) {
    // Make sure they are deleted when done
    std::auto_ptr<std::ostream> details(open_file(_prefix, ".log", "/dev/fd/2"));
    std::auto_ptr<std::ostream> output(open_file(_prefix, ".fa", "/dev/fd/1"));
    // Multiplexers, same thing
    std::auto_ptr<jflib::o_multiplexer>
      log_m(new jflib::o_multiplexer(details.get(), 3 * nb_threads, 1024));
    std::auto_ptr<jflib::o_multiplexer>
      output_m(new jflib::o_multiplexer(output.get(), 3 * nb_threads, 1024));
    _log    = log_m.get();
    _output = output_m.get();

    exec_join(nb_threads);
  }

  virtual void start(int id) {
    instance_t(*this, id).start();
  }

  error_correct_t& skip(int s) { _skip = s; return *this; }
  error_correct_t& good(int g) { _good = g; return *this; }
  error_correct_t& anchor(int a) { _anchor = a; return *this; }
  error_correct_t& prefix(const char *s) { _prefix = s; return *this; }
  error_correct_t& prefix(const std::string s) { _prefix = s; return *this; }
  error_correct_t& mer_len(int l) { _mer_len = l; return *this; }
  error_correct_t& min_count(int m) { _min_count = m; return *this; }
  error_correct_t& cutoff(int p) { _cutoff = p; return *this; }
  //  error_corret_t & advance(int a) { _advance = a; return *this; }
  error_correct_t& window(int w) { _window = w; return *this; }
  error_correct_t& error(int e) { _error = e; return *this; }
  error_correct_t& gzip(bool g) { _gzip = g; return *this; }
  error_correct_t& mer_database(database_query* q) { _mer_database = q; return *this; }
  error_correct_t& contaminant(contaminant_check* c) { _contaminant = c; return *this; }
  error_correct_t& trim_contaminant(bool t) { _trim_contaminant = t; return *this; }
  error_correct_t& homo_trim(int t) { _homo_trim = t; return *this; }

  read_parser& parser() { return _parser; }
  int skip() const { return _skip; }
  int good() const { return _good; }
  int anchor() const { return _anchor; }
  const std::string & prefix() const { return _prefix; }
  int mer_len() const { return _mer_len; }
  int min_count() const { return _min_count; }
  int cutoff() const { return _cutoff; }
  //  int advance() const { return _advance; }
  int window() const { return _window ? _window : _mer_len; }
  int error() const { return _error ? _error : _mer_len / 2; }
  bool gzip() const { return _gzip; }
  database_query* mer_database() const { return _mer_database; }
  contaminant_check* contaminant() const { return _contaminant; }
  bool trim_contaminant() const { return _trim_contaminant; }
  bool do_homo_trim() const { return _homo_trim != std::numeric_limits<int>::min(); }
  int homo_trim() const { return _homo_trim; }

  jflib::o_multiplexer& output() { return *_output; }
  jflib::o_multiplexer& log() { return *_log; }
};

class error_correct_instance {
public:
  typedef error_correct_t<error_correct_instance> ec_t ;

private:
  ec_t&  _ec;
  int    _id;
  size_t _buff_size;
  char*  _buffer;

  static const char* error_contaminant;
  static const char* error_no_starting_mer;
  static const char* error_homopolymer;

public:
  error_correct_instance(ec_t& ec, int id) :
    _ec(ec), _id(id), _buff_size(0), _buffer(0) { }
  ~error_correct_instance() {
    free(_buffer);
  }

  void start() {
    jflib::omstream output(_ec.output());
    jflib::omstream details(_ec.log());
    kmer_t          mer, tmer;

    uint64_t nb_reads = 0;
    bool     parity   = true;
    while(true) {
      read_parser::job job(_ec.parser());
      if(job.is_empty()) break;
      for(size_t i = 0; i < job->nb_filled; ++i) {
        const std::string& header = job->data[i].header;
        const std::string& sequence = job->data[i].seq;
        const char* const seq_s = sequence.c_str();
        const char* const seq_e = seq_s + sequence.size();

        parity = !parity;
        nb_reads++;
        insure_length_buffer(sequence.size());

        const char* error = "";
        const char *input = seq_s + _ec.skip();
        char       *out   = _buffer + _ec.skip();
        //Prime system. Find and write starting k-mer
        if(!find_starting_mer(mer, input, seq_e, out, &error)) {
          details << "Skipped " << header
                  << ": " << error << "\n";
          details << jflib::endr;
          output << jflib::endr;
          continue;
        }
        // Extend forward and backward
        tmer = mer;
        forward_log fwd_log(_ec.window(), _ec.error());
        char *end_out =
          extend(forward_mer(tmer), forward_ptr<const char>(input),
                 forward_counter(input - seq_s),
                 forward_ptr<const char>(seq_e),
                 forward_ptr<char>(out), fwd_log,
                 &error);
        if(!end_out) {
          details << "Skipped " << header
                  << ": " << error << "\n";
          details << jflib::endr;
          output << jflib::endr;
          continue;
        }
        assert(input > seq_s + mer_dna::k());
        assert(out > _buffer + mer_dna::k());
        assert(input - seq_s == out - _buffer);
        tmer = mer;
        backward_log bwd_log(_ec.window(), _ec.error());
        char *start_out =
          extend(backward_mer(tmer),
                 backward_ptr<const char>(input - mer_dna::k() - 1),
                 backward_counter(input - mer_dna::k() - seq_s - 1),
                 backward_ptr<const char>(seq_s - 1),
                 backward_ptr<char>(out - mer_dna::k() - 1), bwd_log,
                 &error);
        if(!start_out) {
          details << "Skipped " << header
                  << ": " << error << "\n";
          details << jflib::endr;
          output << jflib::endr;
          continue;
        }
        start_out++;
        assert(start_out >= _buffer);
        assert(_buffer + _buff_size >= end_out);

        if(_ec.do_homo_trim()) {
          end_out = homo_trim(_buffer, start_out, end_out, fwd_log, bwd_log, &error);
          if(!end_out) {
            details << "Skipped " << header
                    << ": " << error << "\n";
            details << jflib::endr;
            output << jflib::endr;
            continue;
          }
        }
        assert(end_out >= _buffer);
        assert(_buffer + _buff_size >= end_out);

        output << ">" << header
               << " " << fwd_log << " " << bwd_log << "\n"
               << substr(start_out, end_out) << "\n";
        if(parity)
          output << jflib::endr;
      } // for(size_t i...  Loop over reads in job
    } // while(true)... loop over all jobs
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
    counter  cpos       = pos;
    uint32_t prev_count = _ec.mer_database()->get_val(mer.canonical());

    for( ; input < end; ++input) {
      char     base        = *input;

      if(base == '\n')
        continue;
      cpos = pos;
      ++pos;

      int ori_code;
      if(!mer.shift(base)) {
        ori_code = -1; // Invalid base
        mer.shift(0);
      } else {
        if(_ec.contaminant()->is_contaminant(mer.canonical())) {
          if(_ec.trim_contaminant()) {
            log.truncation(cpos);
            goto done;
          }
          *error = error_contaminant;
          return 0;
        }
        ori_code = mer.code(0);
      }
      uint64_t counts[4];
      int      ucode = 0;
      int      count;
      int      level;

      count = _ec.mer_database()->get_best_alternatives(mer, counts, ucode, level);

      // No coninuation whatsoever, trim.
      if(count == 0) {
        log.truncation(cpos);
        goto done;
      }

      if(count == 1) { // One continuation. Is it an error?
        prev_count = counts[ucode];
        if(ucode != ori_code) {
          mer.replace(0, ucode);
          if(_ec.contaminant()->is_contaminant(mer.canonical())) {
            if(_ec.trim_contaminant()) {
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
      // We get here if there is more than one alternative base to try
      // at some level. Note that if the current base is low quality
      // and all alternatives are higher quality, then the current
      // base will have a count of zero. If the current base is non N
      // and has a high count or previous count is low enough that
      // continuity does not apply, output current base. But if the current
      // base has count of zero, all alternatives are low quality and prev_count is low
      // then trim
      if(ori_code >= 0){ //if the current base is valid base (non N)
	if(counts[ori_code] > (uint64_t)_ec.min_count()) {
          if(counts[ori_code]>=(uint32_t)_ec.cutoff()) {
            *out++ = mer.base(0);
            continue;
          }
          // Now we ask for a probability of getting
          // counts[ori_code] errors with p=1/300 in sum_counts
          // trials.  If this probability is < 10e-6, do not correct
          double n = 0;
          for(int i = 0; i < 4; ++i)
            n += (double)counts[i];

          const double k = counts[ori_code];
          const double p = n / 300.;
          const double prob = pow(p / k, k) * exp(-p + k) / sqrt(2 * 3.1415927 * k);
          if(prob < 1e-6) {
            *out++ = mer.base(0);
            continue;
          }
	} else if(level == 0  && counts[ori_code] == 0) {
          // definitely an error and all alternatives are low quality
          log.truncation(cpos);
          goto done;
	}
      } else if(level == 0) { //if all alternatives are low quality
	log.truncation(cpos);
	goto done;
      }

      // We get here if there are multiple possible substitutions, the
      // original count is low enough and the previous count is high (good) or
      // the current base is an N
      // We find out all alternative bases
      // that have a continuation at the same or better level.  Then
      // if the current base is N, pick the one with the highest
      // count that is the most similar to the prev_count,
      // otherwise pick the one with the most similar count.
      // If no alternative continues, leave the base alone.
      int          check_code               = ori_code;
      bool         success                  = false;
      uint64_t     cont_counts[4]; //here we record the counts for the continuations
      bool         continue_with_correct_base[4];
      int          read_nbase_code          = -1;
      bool         candidate_continuations[4];
      unsigned int ncandidate_continuations = 0;

      //here we determine what the next base in the read is
      if(input + 1 < end)
        read_nbase_code = mer_dna::code(*(input + 1));

      for(int i = 0; i < 4; ++i) {
        cont_counts[i]                = 0;
        continue_with_correct_base[i] = false;
        if(counts[i] <= (uint64_t)_ec.min_count())
          continue;
        check_code = i;
        // Check that it continues at least one more base with that quality
        dir_mer    nmer   = mer;
        nmer.replace(0, check_code);
        // Does not matter what we shift, check all alternative anyway.
        nmer.shift(0);

        uint64_t   ncounts[4];
        int        ncount;
        int        nucode = 0;
        int        nlevel;
        ncount = _ec.mer_database()->get_best_alternatives(nmer, ncounts, nucode, nlevel);
        if(ncount > 0 && nlevel >= level) {
          continue_with_correct_base[i] = read_nbase_code >= 0 && ncounts[read_nbase_code] > 0;
          success                       = true;
          cont_counts[i]                = counts[i];
        }
      }

      if(success) {
        // We found at least one alternative base that continues now
        // we look for all alernatives that have a count closest to
        // the previous count first we determine the count that is the
        // closest to the current count but in the special case of
        // prev_count == 1 we simply pick the largest count
        check_code           = -1;
        uint32_t _prev_count = prev_count<=(uint64_t)_ec.min_count() ? std::numeric_limits<uint32_t>::max() : prev_count;
        int      min_diff    = std::numeric_limits<int>::max();
        for(int  i = 0; i < 4; ++i) {
          candidate_continuations[i] = false;
          if(cont_counts[i] > 0)
            min_diff = std::min(min_diff, abs(cont_counts[i] - _prev_count));
        }

        //we now know the count that is the closest, now we determine how many alternatives have this count
        for(uint32_t  i = 0; i < 4; i++) {
          if(abs(cont_counts[i] - _prev_count) == min_diff){
            candidate_continuations[i] = true;
            ++ncandidate_continuations;
            check_code=i;
          }
        }

        //do we have more than one good candidate? if we do then check which one continues with the correct base
        if(ncandidate_continuations>1 && read_nbase_code >= 0)
          for(int  i = 0; i < 4; ++i){
            if(candidate_continuations[i]){
              if(!continue_with_correct_base[i])
                --ncandidate_continuations;
              else
                check_code = i;
            }
          }

        //fail if we still have more than one candidate
        if(ncandidate_continuations != 1)
          check_code = -1;

        if(check_code >= 0){
          mer.replace(0, check_code);
          if(_ec.contaminant()->is_contaminant(mer.canonical())) {
            if(_ec.trim_contaminant()) {
              log.truncation(cpos);
              goto done;
            }
            *error = error_contaminant;
            return 0;
          }
          if(log.substitution(cpos, base, mer.base(0)))
            goto truncate;
        }
      }
      if(ori_code < 0 && check_code < 0) {// if invalid base and no good sub found
        log.truncation(cpos);
        goto done;
      }
      *out++ = mer.base(0);
    }

  done:
    return out.ptr();

  truncate:
    int diff = log.remove_last_window();
    out = out - diff;
    log.truncation(cpos - diff);
    goto done;
  }

  char* homo_trim(const char* start, char* out_start, char* out_end,
		  forward_log& fwd_log, backward_log& bwd_log, const char** error) {
    int   max_homo_score = std::numeric_limits<int>::min();
    char* max_pos        = 0;
    int   homo_score     = 0;
    char* ptr            = out_end - 1;
    char  pbase          = mer_dna::code(*ptr);

    for(--ptr; ptr >= out_start; --ptr) {
      char cbase = mer_dna::code(*ptr);
      homo_score += ((pbase == cbase) << 1) - 1; // Add 1 if same as last, -1 if not
      pbase       = cbase;
      if(homo_score > max_homo_score) {
	max_homo_score = homo_score;
	max_pos        = ptr;
      }
    }

    if(max_homo_score < _ec.homo_trim())
      return out_end; // Not a high score -> return without truncation
    assert(max_pos >= out_start);
    assert(max_pos >= start);
    if(max_pos == out_start) {
      *error = error_homopolymer;
      return 0;
    }
    fwd_log.force_truncate(forward_counter(max_pos - start));
    bwd_log.force_truncate(backward_counter(max_pos - start));
    fwd_log.truncation(forward_counter(max_pos - start));
    return max_pos;
  }

  void insure_length_buffer(size_t len) {
    if(len > _buff_size) {
      _buff_size = len > 2 * _buff_size ? len + 100 : 2 * _buff_size;
      _buffer    = (char *)realloc(_buffer, _buff_size);

      if(!_buffer)
	eraise(std::runtime_error)
	  << "Buffer allocation failed, size " << _buffer << jellyfish::err::no;
    }
  }

  bool find_starting_mer(kmer_t &mer, const char * &input, const char *end, char * &out,
			 const char** error) {
    while(input < end) {
      for(int i = 0; input < end && i < _ec.mer_len(); ++i) {
	char base = *input++;
	*out++ = base;
	if(!mer.shift_left(base))
	  i = -1;        // If an N, skip to next k-mer
      }
      int found = 0;
      while(input < end) {
	bool contaminated = _ec.contaminant()->is_contaminant(mer.canonical());
	if(contaminated && !_ec.trim_contaminant()) {
	  *error = error_contaminant;
	  return false;
	}

	if(!contaminated) {
	  hval_t val = _ec.mer_database()->get_val(mer.canonical());

	  found = (int)val >= _ec.anchor() ? found + 1 : 0;
	  if(found >= _ec.good())
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
const char* error_correct_instance::error_homopolymer     = "Entire read is an homopolymer";


int main(int argc, char *argv[])
{
  args_t args(argc, argv);

  database_query mer_database(args.db_arg);

  // Open contaminant database. Skipped for now. No contaminant.
  std::unique_ptr<contaminant_check> contaminant;
  contaminant.reset(new contaminant_no_database());
  // if(args.contaminant_given) {
  //   mapped_file dbf(args.contaminant_arg);
  //   dbf.random().will_need().load();
  //   inv_hash_storage_t ary = *raw_inv_hash_query_t(dbf).get_ary();
  //   if(ary.get_key_len() != key_len)
  //     die << "Contaminant hash must have same key length as other hashes ("
  //         << ary.get_key_len() << " != " << key_len << ")";
  //   contaminant.reset(new contaminant_database(ary));
  // } else {
  //   contaminant.reset(new contaminant_no_database());
  // }
  stream_manager streams(args.sequence_arg.cbegin(), args.sequence_arg.cend(), 1);

  error_correct_instance::ec_t correct(args.thread_arg, streams);
  correct.skip(args.skip_arg).good(args.good_arg)
    .anchor(args.anchor_count_arg)
    .prefix(args.output_given ? (std::string)args.output_arg : "")
    .min_count(args.min_count_arg)
    .cutoff(args.cutoff_arg)
    .window(args.window_given ? args.window_arg : mer_dna::k())
    .error(args.error_given ? args.error_arg : mer_dna::k() / 2)
    .gzip(args.gzip_flag)
    .mer_database(&mer_database)
    .contaminant(contaminant.get())
    .trim_contaminant(args.trim_contaminant_flag)
    .homo_trim(args.homo_trim_given ? args.homo_trim_arg : std::numeric_limits<int>::min());
  correct.do_it(args.thread_arg);

  return 0;
 }

