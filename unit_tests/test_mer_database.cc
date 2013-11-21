#include <fstream>

#include <gtest/gtest.h>

#include <unit_tests/test_misc.hpp>
#include <src/mer_database.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/misc.hpp>

namespace {
std::string generate_sequence(const size_t s) {
  std::string res(s, 'A');

  for(size_t i = 0; i < s; ++i)
    res[i] = mer_dna::rev_code(jellyfish::random_bits(2));
  return res;
}

void insert_sequence(hash_with_quality& hash, const std::string& seq, const unsigned int quality) {
  mer_dna m;
  for(size_t i = 0; i <= seq.size() - mer_dna::k(); ++i) {
    m = seq.substr(i, mer_dna::k());
    hash.add(m, quality);
  }
}

void test_sequence(const database_query& hash, const std::string& seq, const uint64_t count, const int quality) {
  mer_dna m;
  for(size_t i = 0; i <= seq.size() - mer_dna::k(); ++i) {
    //  for(size_t i = 0; i <= 10; ++i) {
    m = seq.substr(i, mer_dna::k());
    auto res = hash[m];
    EXPECT_EQ(count, res.first);
    EXPECT_EQ(quality, res.second);
  }
}

TEST(MerDatabase, WriteRead) {
  file_unlink database_file("mer_database");
  database_file.do_unlink = false;

  // Insert in the database the following data sets:
  // hq2: high quality, inserted twice
  // hq1: high quality, inserted once
  // lqhq: inserted first as low quality, then as high quality
  // hqlq: inserted first as high quality, then as low quality
  // lq1
  // lq2
  static const size_t       sequence_len = 10000;
  static const unsigned int bits         = 4;
  std::string hq2  = generate_sequence(sequence_len);
  std::string hq1  = generate_sequence(sequence_len);
  std::string lqhq = generate_sequence(sequence_len);
  std::string hqlq = generate_sequence(sequence_len);
  std::string lq1  = generate_sequence(sequence_len);
  std::string lq2  = generate_sequence(sequence_len);

  mer_dna::k(20);

  size_t key_start, val_start, total_len;
  {
    hash_with_quality database(10 * sequence_len, mer_dna::k() * 2, bits);
    insert_sequence(database, hq2, 1);
    insert_sequence(database, hq2, 1);
    insert_sequence(database, hq1, 1);
    insert_sequence(database, lqhq, 0);
    insert_sequence(database, lqhq, 1);
    insert_sequence(database, hqlq, 1);
    insert_sequence(database, hqlq, 0);
    insert_sequence(database, lq1, 0);
    insert_sequence(database, lq2, 0);
    insert_sequence(database, lq2, 0);
    std::ofstream os(database_file.path.c_str());
    ASSERT_TRUE(os.good());
    database_header header;
    database.write(os, &header);
    EXPECT_TRUE(os.good());
    EXPECT_EQ((size_t)0, header.offset() % sizeof(uint64_t));
    key_start = header.offset();
    val_start = key_start + header.key_bytes();
    total_len = val_start + header.value_bytes();
    EXPECT_EQ((std::ofstream::pos_type)total_len, os.tellp());
  }

  database_query database(database_file.path.c_str());
  EXPECT_EQ(bits, database.header().bits());
  EXPECT_EQ(mer_dna::k() * 2, database.header().key_len());
  test_sequence(database, hq2, 2, 1);
  test_sequence(database, hq1, 1, 1);
  test_sequence(database, lqhq, 1, 1);
  test_sequence(database, hqlq, 1, 1);
  test_sequence(database, lq1, 1, 0);
  test_sequence(database, lq2, 2, 0);
}
}
