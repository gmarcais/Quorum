#include <algorithm>

#include <jellyfish/err.hpp>
#include <src/mer_database.hpp>

int main(int argc, char* argv[])
{
  if(argc != 2)
    die << "Usage: " << argv[0] << " db";

  static const size_t hlen = 1001;
  uint64_t histos[hlen][2];
  memset(histos, '\0', sizeof(histos));
  const database_query mer_database(argv[1]);
  for(auto it = mer_database.vals().begin(); it != mer_database.vals().end(); ++it) {
    const uint64_t val = std::min(*it >> 1, hlen - 1);
    ++histos[val][*it & 0x1];
  }

  for(size_t i = 0; i < hlen; ++i) {
    if(histos[i][0] || histos[i][1])
      std::cout << i << " " << histos[i][0] << " " << histos[i][1] << "\n";
  }

  return 0;
}
