/*  PSIKO Copyright (C) 2014  Popescu Andrei-Alin

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


#include <armadillo>
#include <ctime>
#include <limits.h>
#include <algorithm>
#include <stdint.h>
#include "bitwise_input/bitwise_input_header.hpp"
#include "util/util_header.hpp"
#include "bitwise_kpca/bitwise_kpca_header.hpp"
#include "optimise/optimise_header.hpp"


int main(int argc, char* argv[]){

  char qFile[512];
  char filename[512];
  char red[512];
  strcpy(red,"reduced_data.csv");
  int newDim=-1;
  char c;
  int iFlag=0;
  int lFlag=0;

  while ((c = getopt (argc, argv, "i:K:q:r:l")) != -1)
         switch (c)
           {
           case 'i':
             strcpy(filename,optarg);
             strcpy(qFile,filename);
             strcat(qFile,".PSI.Q");  
             iFlag=1;
             break;
           case 'K':
             newDim = atoi(optarg);
             break;
           case 'r':
              strcpy(red,optarg);
              break;
           case 'q':
             strcpy(qFile,optarg);
             break;
           case 'l':
              lFlag=1;
              break;
           default:
             abort ();
           }
  
  if(lFlag)
    {
      cout<<
      "PSIKO Copyright (C) 2014  Popescu Andrei-Alin\n"<<

      "This program is free software: you can redistribute it and/or modify\n"<<
      "it under the terms of the GNU General Public License as published by\n"<<
      "the Free Software Foundation, either version 3 of the License, or\n"<<
      "(at your option) any later version.\n"<<

      "This program is distributed in the hope that it will be useful,\n"<<
      "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"<<
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"<<
      "GNU General Public License for more details.\n"<<

      "You should have received a copy of the GNU General Public License\n"<<
      "along with this program.  If not, see http://www.gnu.org/licenses/.\n";
      exit(0);
    }
  
  cout<<"PSIKO  Copyright (C) 2014  Popescu Andrei-Alin\n";
  cout<<"This program comes with ABSOLUTELY NO WARRANTY; for details type \'./PSIKO -l\'.\n";
  cout<<"This is free software, and you are welcome to redistribute it\n"<<
    "under certain conditions; type \'./PSIKO -l\' for details.\n";

  if(!iFlag)
   {
    cout<<"ERROR: You need to specify input file\n";
    exit(1);
   }
  
  if(newDim<0)
   {
    cout<<" K unspecified, using TW test\n";
   }

  double m=25;
  Mat<unit> dataset;
  arma::mat rDataset;
  vec evals;
  arma::mat evec;
  clock_t st_time,ed_time;
  cout<<"reading "<<filename<<"\n";
  read_snps(dataset,filename);
  cout<<"Input done \n";
  
  st_time=clock();
  Apply(dataset,rDataset,evals,evec,newDim);
  ofstream out("means.csv");
  arma::mat means(newDim-1,newDim,arma::fill::zeros);
  int ind=0; 
  //cout<<"reduced data\n";
  //rDataset.print();
  initialise(rDataset,means);
  my_optimize(means,rDataset);
  for(int i=0;i<means.n_cols;i++)
   {
    for(int j=0;j<means.n_rows;j++)
      out<<means(j,i)<<" ";
    out<<"\n";
   }  
  ed_time=clock();
  double tdiff=(ed_time-st_time)/(double) CLOCKS_PER_SEC;
  cout<<"PCA completed in "<<tdiff<<" \n"; 
  ofstream qout(qFile);

  arma::mat A= arma::ones(1,means.n_cols); 
  A.insert_rows(1,means);
  A=arma::inv(A);

  for(int i=0;i<rDataset.n_cols;i++)
   {
	  arma::vec ancestry=inferAncestry(rDataset.unsafe_col(i),A);
        for(int j=0;j<ancestry.n_elem;j++)
          {
		qout<<ancestry(j)<<" ";
          }
	  qout<<"\n";
   }
  qout.close();
  out.close();
  ofstream rout(red);
  for(int i=0;i<rDataset.n_cols;i++)
  {
    for(int j=0;j<rDataset.n_rows;j++)
     {
       rout<<rDataset(j,i)<<(j<rDataset.n_rows-1?",":"");
     }
    rout<<"\n";
  }
  rout.close();
}
