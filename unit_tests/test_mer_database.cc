#include <fstream>
#include <stdexcept>
#include <thread>

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

void insert_sequence(hash_with_quality* hash, const std::string& seq, const unsigned int quality) {
  mer_dna m;
  for(size_t i = 0; i <= seq.size() - mer_dna::k(); ++i) {
    m = seq.substr(i, mer_dna::k());
    if(!hash->add(m, quality))
      throw std::runtime_error("Hash is full");
  }
  hash->done();
}

void test_sequence(const database_query& hash, const std::string& seq, const uint64_t count, const int quality,
                   const char* msg) {
  SCOPED_TRACE(msg);
  mer_dna m;
  for(size_t i = 0; i <= seq.size() - mer_dna::k(); ++i) {
    m = seq.substr(i, mer_dna::k());
    SCOPED_TRACE(::testing::Message() << "i:" << i << " m:" << m);
    auto res = hash[m];
    EXPECT_EQ(count, res.first);
    EXPECT_EQ(quality, res.second);
  }
}

class MerDatabase : public ::testing::TestWithParam<int> { };

TEST_P(MerDatabase, WriteRead) {
  file_unlink database_file("mer_database");

  // Insert in the database the following data sets:
  // hq2: high quality, inserted twice
  // hq1: high quality, inserted once
  // lqhq: inserted first as low quality, then as high quality
  // hqlq: inserted first as high quality, then as low quality
  // lq1
  // lq2
  static const size_t       sequence_len = 100000;
  static const unsigned int bits         = 4;
  std::string hq2  = generate_sequence(sequence_len);
  std::string hq1  = generate_sequence(sequence_len);
  std::string lqhq = generate_sequence(sequence_len);
  std::string hqlq = generate_sequence(sequence_len);
  std::string lq1  = generate_sequence(sequence_len);
  std::string lq2  = generate_sequence(sequence_len);

  mer_dna::k(33);

  size_t key_start, val_start, total_len;
  {
    hash_with_quality database(GetParam() * sequence_len, mer_dna::k() * 2, bits, 10);
    std::thread th_hq2_1(insert_sequence, &database, hq2, 1);
    std::thread th_hq2_2(insert_sequence, &database, hq2, 1);
    std::thread th_hq1_1(insert_sequence, &database, hq1, 1);
    std::thread th_lqhq_0(insert_sequence, &database, lqhq, 0);
    std::thread th_lqhq_1(insert_sequence, &database, lqhq, 1);
    std::thread th_lqhq_2(insert_sequence, &database, hqlq, 1);
    std::thread th_lqhq_3(insert_sequence, &database, hqlq, 0);
    std::thread th_lq1_1(insert_sequence, &database, lq1, 0);
    std::thread th_lq2_1(insert_sequence, &database, lq2, 0);
    std::thread th_lq2_2(insert_sequence, &database, lq2, 0);

    th_hq2_1.join();
    th_hq2_2.join();
    th_hq1_1.join();
    th_lqhq_0.join();
    th_lqhq_1.join();
    th_lqhq_2.join();
    th_lqhq_3.join();
    th_lq1_1.join();
    th_lq2_1.join();
    th_lq2_2.join();

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
  test_sequence(database, hq2, 2, 1, "hq2");
  test_sequence(database, hq1, 1, 1, "hq1");
  test_sequence(database, lqhq, 1, 1, "lqhq");
  test_sequence(database, hqlq, 1, 1, "hqlq");
  test_sequence(database, lq1, 1, 0, "lq1");
  test_sequence(database, lq2, 2, 0, "lq2");
}

// Instantiate test for different size of mer databases
INSTANTIATE_TEST_CASE_P(MerDatabaseTest, MerDatabase, ::testing::Values(1, 10, 20, 40));
}
