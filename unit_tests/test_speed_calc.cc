#include <gtest/gtest.h>

#include <jellyfish/misc.hpp>
#include <src/mer_database.hpp>

namespace {
TEST(SpeedCalc, Simple) {
  uint64_t oids[4];

  for(int i = 0; i < 100; ++i) {
    const int r = jellyfish::random_bits(5) + 20;
    const int c = 2 * (jellyfish::random_bits(7) + 1);
    mer_dna::k(c >> 1);
    const uint64_t size_mask = ((uint64_t)1 << r) - 1;
    RectangularBinaryMatrix matrix(r, c);
    matrix.randomize(jellyfish::random_bits);
    oid_speed_calc osc(matrix, size_mask);

    mer_dna m;
    for(int j = 0; j < 100; ++j) {
      m.randomize();
      osc.calc(m, oids);
      for(int b = 0; b < 4; ++b) {
        m.base(0) = b;
        EXPECT_EQ(matrix.times(m) & size_mask, oids[b]);
      }
    }
  }
}
}
