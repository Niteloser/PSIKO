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
  //strcpy(filename,"Example/OSRMatrix_Complete.txt.geno");
  strcpy(filename,"sim1K3Actual.geno");
  long long N,L;
  arma::mat rDataset;
  vec evals;
  arma::mat evec;
  int newDim=3;
  read_snps(dataset,filename,N,L);

  Apply(dataset,rDataset,evals,evec,newDim);

  //cout<<evec<<"\n";
  //cout<<rDataset<<"\n";
  //arma::mat W=evec*(getDatasetWindow(dataset,0,15,L)).t();
  cout<<rDataset;
  arma::mat means(newDim-1,newDim,arma::fill::zeros);
  int ind=0; 
  initialise(rDataset,means);
  my_optimize(means,rDataset); 

  arma::mat A= arma::ones(1,means.n_cols); 
  A.insert_rows(1,means);
  A=arma::inv(A);
  arma::mat Q;
  Q=arma::zeros(rDataset.n_cols,means.n_cols);
  for(int i=0;i<rDataset.n_cols;i++)
   {
	  arma::vec ancestry=inferAncestry(rDataset.unsafe_col(i),A);
      for(int j=0;j<ancestry.n_elem;j++)
          {
                  Q(i,j)=ancestry(j);
          }
   }
  cout<<"Q=\n";
  cout<<Q;
  applyWindow(getDatasetWindow(dataset,0,15,L),evec,evals,Q);

}
