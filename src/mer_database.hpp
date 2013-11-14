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

#include <jellyfish/file_header.hpp>
#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/atomic_bits_array.hpp>
#include <jellyfish/mer_dna.hpp>

using jellyfish::mer_dna;
typedef jellyfish::large_hash::array<mer_dna> mer_array;
typedef jellyfish::atomic_bits_array<uint64_t> val_array;


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

class database_query {
public:
  database_query(const char* filename) { }
};

#endif /* __QUORUM_MER_DATABASE_HPP__ */
