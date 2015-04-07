#include <armadillo>
#include <ctime>
#include <limits.h>
#include <algorithm>
#include <stdint.h>
#include <bitset>
#include "bitwise_input/bitwise_input_header.hpp"
#include "util/util_header.hpp"
#include "bitwise_kpca/bitwise_kpca_header.hpp"
#include "optimise/optimise_header.hpp"
using namespace std;

int main(int argc, char* argv[]){
  Mat<unit> dataset;
  char filename[512];
  strcpy(filename,"Example/OSRMatrix_Complete.txt.geno");
  long long N,L;
  read_snps(dataset,filename,N,L);
  cout<<getWindow(dataset.unsafe_col(1),43,200,L);
}
