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
#include "local_ancestry/local_ancestry_header.hpp"
using namespace std;

int main(int argc, char* argv[]){
  Mat<unit> dataset;
  char filename[512];
  strcpy(filename,"testfile.geno");
  long long N,L;
  arma::mat rDataset;
  vec evals;
  arma::mat evec;
  int newDim=3;
  read_snps(dataset,filename,N,L);

  Apply(dataset,rDataset,evals,evec,newDim);

  cout<<evec<<"\n";
  cout<<rDataset<<"\n";
  //arma::mat W=evec*(getDatasetWindow(dataset,0,15,L)).t();
  cout<<rDataset;
  arma::mat Q(3,3);
  Q(0,0)=0.97; Q(0,1)=0.97; Q(0,2)=0.97;
  Q(1,0)=0.97; Q(1,1)=0.03; Q(1,2)=0.0;
  Q(2,0)=0.0; Q(2,1)=0.97; Q(2,2)=0.97;
  cout<<applyWindow(getDatasetWindow(dataset,0,15,L),evec,evals,Q);

}
