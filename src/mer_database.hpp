/* Quorum
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

#ifndef __QUORUM_MER_DATABASE_HPP__
#define __QUORUM_MER_DATABASE_HPP__

#include <fstream>

#include <jellyfish/file_header.hpp>
#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/atomic_bits_array.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>

using jellyfish::mer_dna;
typedef jellyfish::large_hash::array<mer_dna> mer_array;
typedef jellyfish::large_hash::array_raw<mer_dna> mer_array_raw;
typedef jellyfish::atomic_bits_array<uint64_t> val_array;
typedef jellyfish::atomic_bits_array_raw<uint64_t> val_array_raw;

class database_header : public jellyfish::file_header {
public:
  database_header() : jellyfish::file_header() { }
  database_header(std::istream& is) : jellyfish::file_header(is) { }

  void bits(uint32_t b) { root_["bits"] = (Json::UInt)b; }
  uint32_t bits() const { return root_["bits"].asUInt(); }

  size_t value_bytes() const { return root_["value_bytes"].asLargestUInt(); }
  void value_bytes(size_t bytes) { root_["value_bytes"] = (Json::UInt64)bytes; }

  size_t key_bytes() const { return root_["key_bytes"].asLargestUInt(); }
  void key_bytes(size_t bytes) { root_["key_bytes"] = (Json::UInt64)bytes; }

  void set_format() {
    this->format("binary/quorum_db");
  }
  bool check_format() const {
    return "binary/quorum_db" == this->format();
  }
};

class hash_with_quality {
  mer_array      keys_;
  val_array      vals_;
  const uint64_t max_val_;

public:
  hash_with_quality(size_t size, uint16_t key_len, int bits, uint16_t reprobe_limit = 126) :
    keys_(size, key_len, 0, reprobe_limit),
    vals_(bits + 1, keys_.size()),
    max_val_((uint64_t)-1 >> (sizeof(uint64_t) * 8 - bits))
  { }

  bool add(const mer_dna& key, unsigned int quality) {
    bool is_new;
    size_t id;
    if(!keys_.set(key, &is_new, &id))
      return false;

    auto     entry = vals_[id];
    uint64_t oval;
    uint64_t nval = entry.get();;
    do {
      oval = nval;
      if((oval & 1) < quality)
        nval = 3;
      else if((oval >> 1) == max_val_ || (oval & 1) > quality)
        return true;
      else
        nval = oval + 2;
    } while(!entry.set(nval));
    return true;
  }

  void write(std::ostream& os, database_header* header = 0) const {
    if(header) {
      header->set_format();
      header->update_from_ary(keys_);
      header->bits(vals_.bits() - 1);
      header->key_bytes(keys_.size_bytes());
      header->value_bytes(vals_.size_bytes());
      header->write(os);
    }
    keys_.write(os);
    vals_.write(os);
  }

  mer_array& keys() { return keys_; }
  val_array& vals() { return vals_; }
};

class database_query {
  const database_header        header_;
  const jellyfish::mapped_file file_;
  const mer_array_raw          keys_;
  const val_array_raw          vals_;

  static database_header parse_header(const char* filename) {
    std::ifstream file(filename);
    if(!file.good())
      eraise(std::runtime_error) << "Can't open '" << filename << "' for reading";
    database_header res;
    if(!res.read(file))
      eraise(std::runtime_error) << "Can't parse header of file '" << filename << "'";
    if(!res.check_format())
      eraise(std::runtime_error) << "Wrong type '" << res.format() << "' for file '" << filename << "'";
    return res;
  }

public:
  database_query(const char* filename) :
  header_(parse_header(filename)),
  file_(filename),
  keys_(file_.base() + header_.offset(), header_.key_bytes(),
        header_.size(), header_.key_len(), header_.val_len(),
        header_.max_reprobe(), header_.matrix()),
  vals_(file_.base() + header_.offset() + header_.key_bytes(), header_.value_bytes(),
        header_.bits() + 1, header_.size())
  { }

  const database_header& header() const { return header_; }
  const mer_array_raw& keys() const { return keys_; }
  const val_array_raw& vals() const { return vals_; }

  std::pair<uint64_t, int> operator[](const mer_dna& m) const {
    std::pair<uint64_t, int> res = { 0, 0 };
    size_t                   id  = 0;
    if(keys_.get_key_id(m, &id)) {
      uint64_t v = vals_[id];
      res.first  = v >> 1;
      res.second = v & 0x1;
    }
    return res;
  }

  // Get value of m in the high quality database
  uint64_t get_val(const mer_dna& m) const {
    auto v = operator[](m);
    return v.second ? v.first : 0;
  }

  // Get all alternatives at the best level
  template<typename mer_type>
  int get_best_alternatives(mer_type& m, uint64_t counts[4], int& ucode, int& level) const {
    mer_dna tmp_mer;
    int     count = 0;
    memset(counts, '\0', sizeof(uint64_t) * 4);
    level = 0;
    int ori_code = m.code(0);

    for(int i = 0; i < 4; ++i) {
      m.replace(0, i);
      auto v = operator[](m.canonical());
      if(v.first > 0) {
        if(v.second >= level) {
          if(v.second > level && count > 0) {
            for(int j = 0; j < i; ++j)
              counts[j] = 0;
            count = 0;
          }
          counts[i] = v.first;
          ucode     = i;
          level     = v.second;
          ++count;
        }
      }
    }
    m.replace(0, ori_code); // Reset m to original value
    return count;
  }
};

#endif /* __QUORUM_MER_DATABASE_HPP__ */
