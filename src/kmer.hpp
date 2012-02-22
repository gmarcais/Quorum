#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <stdint.h>
#include <assert.h>
#include <iostream>
#include <string>

#ifndef bsizeof
#define bsizeof(x) (sizeof(x) * 8)
#endif
static uint64_t reverse_complement(uint64_t v, int length) {
  v = ((v >> 2)  & 0x3333333333333333UL) | ((v & 0x3333333333333333UL) << 2);
  v = ((v >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((v & 0x0F0F0F0F0F0F0F0FUL) << 4);
  v = ((v >> 8)  & 0x00FF00FF00FF00FFUL) | ((v & 0x00FF00FF00FF00FFUL) << 8);
  v = ((v >> 16) & 0x0000FFFF0000FFFFUL) | ((v & 0x0000FFFF0000FFFFUL) << 16);
  v = ( v >> 32                        ) | ( v                         << 32);
  return (((uint64_t)-1) - v) >> (bsizeof(v) - (length << 1));
}

class kmer_t {
  uint64_t _fmer, _rmer;

  static const uint64_t c3 = (uint64_t)0x3;
  static int            _k;
  static int            _lshift;
  static uint64_t       _mask;

public:
  static const uint64_t codes[256];
  static const char     rev_codes[4];

  kmer_t() : 
    _fmer(0), _rmer(0) {}

  kmer_t(uint64_t mer) :
    _fmer(mer), _rmer(reverse_complement(_fmer, _k)) {}

  static void k(int new_k) { 
    _k = new_k;
    _mask = ((uint64_t)1 << (2*_k)) - 1;
    _lshift = 2 * (_k - 1);
  }
  static int k() { return _k; }

  bool shift_left(char c) { 
    uint64_t x = codes[(int)c];
    if(x == (uint64_t)-1)
      return false;
    shift_left(x);
    return true;
  }

  void shift_left(uint64_t c) {
    c &= c3;
    _fmer = ((_fmer << 2) & _mask) | c;
    _rmer = (_rmer >> 2) | ((c3 - c) << _lshift);
  }

  bool shift_right(char c) {
    uint64_t x = codes[(int)c];
    if(x == (uint64_t)-1)
      return false;
    shift_right(x);
    return true;
  }

  void shift_right(uint64_t c) {
    c &= c3;
    _fmer = (_fmer >> 2) | (c << _lshift);
    _rmer = ((_rmer << 2) & _mask) | (c3 - c);
  }

  uint64_t canonical() const { return _fmer < _rmer ? _fmer : _rmer; }
  uint64_t fmer() const { return _fmer; }
  uint64_t rmer() const { return _rmer; }

  void replace(int i, uint64_t c) {
    assert(i >= 0 && i < _k);
    assert(c >= 0 && c < 4);
    c &= c3;
    i <<= 1;
    int j = (_k << 1) - i - 2;
    _fmer = (_fmer & ~(c3 << i)) | (c << i);
    _rmer = (_rmer & ~(c3 << j)) | ((c3 - c) << j);
    
    assert(code(i>>1) == c && rcode((_k-1)-(i>>1)) == (c3 - c));
  }

  char base(int i) const { assert(i >= 0 && i < _k); return rev_codes[(_fmer >> (2*i)) & c3]; }
  uint64_t code(int i) const { assert(i >= 0 && i < _k); return (_fmer >> (2*i)) & c3; }
  uint64_t rcode(int i) const { assert(i >= 0 && i < _k); return (_rmer >> (2*i)) & c3; }
  std::string str() const { return to_str(_fmer); }
  std::string rstr() const { return to_str(_rmer); }

  friend std::ostream &operator<<(std::ostream &os, const kmer_t &mer);
  friend class forward_mer;
  friend class backward_mer;
private:
  std::string to_str(uint64_t x) const;
};

std::ostream &operator<<(std::ostream &os, const kmer_t &mer);

class forward_mer {
  kmer_t _m;
public:
  forward_mer(kmer_t m) : _m(m) {}
  forward_mer(uint64_t m) : _m(m) {}
  bool shift(char c) { return _m.shift_left(c); }
  void shift(uint64_t x) { _m.shift_left(x); }
  bool rev_shift(char c) { return _m.shift_right(c); }
  void rev_shift(uint64_t x) { _m.shift_right(x); }
  uint64_t canonical() const { return _m.canonical(); }
  char base(int i) const { return _m.base(i); }
  uint64_t code(int i) const { return _m.code(i); }
  void replace(int i, uint64_t c) { _m.replace(i, c); }
  uint64_t rmer() const { return _m.rmer(); }
  friend std::ostream &operator<<(std::ostream &os, const forward_mer &mer);
};
std::ostream &operator<<(std::ostream &os, const forward_mer &mer);

class backward_mer {
  kmer_t _m;
public:
  backward_mer(kmer_t m) : _m(m) {}
  backward_mer(uint64_t m) : _m(m) {}
  bool shift(char c) { return _m.shift_right(c); }
  void shift(uint64_t x) { _m.shift_right(x); }
  bool rev_shift(char c) { return _m.shift_left(c); }
  void rev_shift(uint64_t x) { _m.shift_left(x); }
  uint64_t canonical() const { return _m.canonical(); }
  char base(int i) const { return _m.base(kmer_t::_k - i - 1); }
  uint64_t code(int i) const { return _m.code(kmer_t::_k - i - 1); }
  void replace(int i, uint64_t c) { _m.replace(kmer_t::_k - i - 1, c); }
  friend std::ostream &operator<<(std::ostream &os, const backward_mer &mer);
};
std::ostream &operator<<(std::ostream &os, const backward_mer &mer);

    // static uint64_t mer_string_to_binary(const char *in, uint_t klen) {
    //   uint64_t res = 0;
    //   for(uint_t i = 0; i < klen; i++) {
    //     const uint_t c = fasta_parser::codes[(uint_t)*in++];
    //     if(c & CODE_NOT_DNA)
    //       return 0;
    //     res = (res << 2) | c;
    //   }
    //   return res;
    // }
    // static void mer_binary_to_string(uint64_t mer, uint_t klen, char *out) {
    //   static const char table[4] = { 'A', 'C', 'G', 'T' };
      
    //   for(unsigned int i = 0 ; i < klen; i++) {
    //     out[klen-1-i] = table[mer & UINT64_C(0x3)];
    //     mer >>= 2;
    //   }
    //   out[klen] = '\0';
    // }

#endif
