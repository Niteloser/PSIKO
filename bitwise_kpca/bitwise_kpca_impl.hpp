/*  This file is part of PSIKO.

    PSIKO Copyright (C) 2014  Popescu Andrei-Alin

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/.
    */

inline double shortEval(const Col<unit>& a, const Col<unit>& b)//, int ii, int jj)
  {
    double ret=0.0;
    unit prod;
    short ct=0;  

    for(int i=0;i<a.n_elem/2;i++)
      {
       /*a1,b1: first genotype copies of indiv a and b (e.g. mother)
	 a2,b2: second genotype copies of indiv a and b (e.g. father)
         Then the sum __builtin_popcount(a1&b1)+__builtin_popcount(a1&b2)+__builtin_popcount(a2&b1)+__builtin_popcount(a2&b2) corresponds  to the dot product a*b. This can be computed bit chunk by bit chunk and added up to the grand total ret to efficiently obtain the dot product we seek. 
       This property can be easily verified for a single-locus dataset. It generalises to arbitrarily many loci since loci are independent when bit-wise and-ing. Also provable by induction for the curios reader. 
       */
       unit a1=a(i),a2=a(i+a.n_elem/2),
       b1=b(i),b2=b(i+a.n_elem/2);       
        ret+=__builtin_popcountll(a1&b1)+__builtin_popcountll(a1&b2)+__builtin_popcountll(a2&b1)+__builtin_popcountll(a2&b2);
      }
   return ret;
  }

void GetKernelMatrix(const arma::Mat<unit>& data,
                                            arma::mat& kernelMatrix)//based off MLPACK, see http://www.mlpack.org/
{
  kernelMatrix.set_size(data.n_cols, data.n_cols);

  for (size_t i = 0; i < data.n_cols; ++i)
  {
    for (size_t j = i; j < data.n_cols; ++j)
    {
      kernelMatrix(i, j) = shortEval(data.unsafe_col(i),
                                           data.unsafe_col(j));
    }
  }

  for (size_t i = 1; i < data.n_cols; ++i)
    for (size_t j = 0; j < i; ++j)
      kernelMatrix(i, j) = kernelMatrix(j, i);
}


void Apply(const arma::Mat<unit>& data,
                                  arma::mat& transformedData,
                                  arma::vec& eigval,
                                  arma::mat& eigvec,int& nComp)//based off MLPACK, see http://www.mlpack.org/
{
  // Construct the kernel matrix.
  
  clock_t st_time,ed_time;
  arma::mat kernelMatrix;
  double tdiff;
  st_time=clock();
  GetKernelMatrix(data, kernelMatrix);
  ed_time=clock();
  tdiff=(ed_time-st_time)/(double) CLOCKS_PER_SEC;
  cout<<"get covar completed in "<<tdiff<<" \n";
  arma::rowvec rowMean = arma::sum(kernelMatrix, 0) / kernelMatrix.n_cols;
  arma::colvec coll=arma::sum(kernelMatrix, 1) / kernelMatrix.n_cols;
  kernelMatrix.each_row() -= rowMean;
  //arma::colvec coll=arma::sum(kernelMatrix, 1) / kernelMatrix.n_cols;
  kernelMatrix.each_col() -= coll;
  kernelMatrix += arma::sum(rowMean) / kernelMatrix.n_cols;  
  st_time=clock();
  arma::eig_sym(eigval, eigvec, kernelMatrix);
  ed_time=clock();
  tdiff=(ed_time-st_time)/(double) CLOCKS_PER_SEC;
  cout<<"Eigendecomposition completed in "<<tdiff<<" \n";
  
  for (size_t i = 0; i < floor(eigval.n_elem / 2.0); ++i)
    eigval.swap_rows(i, (eigval.n_elem - 1) - i);
  /*cout<<"Kernel Matrix\n";
  kernelMatrix.print();*/
  eigvec = arma::fliplr(eigvec);
  cout<<"eigenvectors\n";
  eigvec.print();
  //cout<<"eigenvalues\n";
  arma::rowvec lambda=sqrt(eigval.t());
  //lambda.print();
  transformedData=eigvec;
  transformedData.each_row() %= lambda;
  //cout<<"after division\n";
  //eigvec.print();
  transformedData = transformedData.t();//eigvec.t() * kernelMatrix;

    if(nComp<0)
     {
      nComp=findScree(eigval.subvec(0,eigval.n_elem-1));
      nComp++;
     }    
    
    cout<<"K= "<<nComp<<"\n";
    if(nComp<=transformedData.n_rows){
                                      transformedData.shed_rows(nComp-1,transformedData.n_rows-1);
                                      eigvec.shed_cols(nComp-1,eigvec.n_cols-1);
                                      eigval.shed_rows(nComp-1,eigval.n_rows-1);
                                     }
    //transformedData.shed_row(0);    

    scale(transformedData);
}
