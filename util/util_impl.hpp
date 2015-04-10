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


inline arma::vec inferAncestry(const arma::vec& x,const arma::mat& A)
{
  arma::mat od=arma::ones(x.n_elem+1,1);
  od(0)=1.0;
  for(int i=0;i<x.n_elem;i++)
   {
    od(i+1)=x(i);
   }
  
  arma::vec ret=A*od; 
  for(int i=0;i<ret.n_elem;i++)
   {
    if(ret(i)<0)ret(i)=0;
   }

  double norm=arma::norm(ret,1);
  
  for(int i=0;i<ret.n_elem;i++)
   {
    ret(i)=ret(i)/norm;
   }
  return ret;
}

unsigned int popcount64(unsigned long long x)
{
    x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    return (x * 0x0101010101010101ULL) >> 56;
}

int findScree(arma::vec evals)
{
  /*arma::vec l1(evals.n_elem-2,arma::fill::zeros);
  int ind=0;
  for(int i=1;i<evals.n_elem-1;i++)
    {
  double prev=0,next=0;
        prev=evals(i-1);
        next=evals(i+1);
        l1(ind++)=next-2*evals[i]+prev;
    }
  int max_ind=-1;
  double max=-1.0;
  for(int i=0;i<l1.n_elem;i++)
    {
  if(l1(i)>max)
         {
     max=l1(i);
     max_ind=i;
         }
    }
  return max_ind+1;*/
  int m=evals.n_elem;
  int ret=0;
  //evals.print();
  for(int j=0;j<m-1;j++)
    {
	//cout<<"j= "<<j<<"\n";
	double l1=0.0,l2=0.0;
	//cout<<"evals\n";
	//cout<<"len of evals "<<evals.n_elem<<"\n";
	for(int i=j;i<m;i++)
  	  {
	   //cout<<evals(i)<<" ";
	   l1+=evals(i);
	   l2+=evals(i)*evals(i);
	  }
	//cout<<"\n";
	//cout<<"l1= "<<l1<<"\n";
	//cout<<"l2= "<<l2<<"\n";
	double lambda=evals(j)*(m-(j+1))/l1;
        //cout<<"lambda= "<<lambda<<"\n";
	double nhat=l1*l1/l2;
	//cout<<"nhat= "<<nhat<<"\n";
        double mu=(sqrt(nhat-1)+sqrt(m-1))*(sqrt(nhat-1)+sqrt(m-1))/nhat;
        //cout<<"mu= "<<mu<<"\n";
	double sigma = (sqrt(nhat-1)+sqrt(m-1))/nhat*pow(1/sqrt(nhat-1)+1/sqrt(m-1),1.0/3);
        //cout<<"sigma= "<<sigma<<"\n";
	double twstat=(lambda-mu)/sigma;
	//cout<<"Tracy-Widom "<<twstat<<"\n";
	if( twstat <0.9794){return ret;}//tw 95% confidence interval no lower tail
        ret++;  
    }
  return ret;
}

inline int Factorial(int x) {
  return (x == 1 ? x : x * Factorial(x - 1));
}


arma::mat& getDatasetWindow(Mat<unit>& dataset,long long start, long long end, long long L)
 {
  end=min(end,L);
  arma::mat* ret=new arma::mat(end-start,dataset.n_cols);
  for(int i=0;i<dataset.n_cols;i++)
     {
      ret->unsafe_col(i)=getWindow(dataset.unsafe_col(i),start,end,L);
     }
  return *ret;
 }

arma::vec getWindow(const Col<unit>& col, long long start, long long end,long long L)
 {
   //return window of individual ind starting at start
   end=min(end,L);
   int wSize=end-start;
   arma::vec ret=arma::zeros<arma::vec>(wSize);
   int startChunk=start/sz;    
   int endChunk=end/sz;
   int n=endChunk-startChunk+1;
   int nr=col.n_elem;
   Col<unit> masks(n);
   for(int i=0;i<n;i++)
    {
      masks(i)=static_cast<unit>(-1);//all ones 
    }
   masks(0)&=~((static_cast<unit>(1)<<(start%sz))-1);
   masks(n-1)&=((static_cast<unit>(1)<<(end%sz))-1);
   //std::cout<<"masks"<<std::endl;
   //print_bits(masks);
   int id=0;
   int chunkStart=start-(start%sz);
   for(int i=0;id<wSize;i++)//loop over all chunks used
    {    
     unit mask=masks(i/sz);
     unit bit=static_cast<unit>(1)<<(i%sz);
     if((bit&mask)==0)continue;//outside window
       int iChunk=chunkStart+i;
       unit chunk=col(iChunk/sz);
       unit chunk2=col(iChunk/sz+nr/2);
       /*cout<<id<<" "<<i<<" here\n";
       cout<<"first chunk\n";
       print_bits(chunk&bit);
       cout<<endl;
       cout<<"second chunk\n";
       print_bits(chunk2&bit);
       cout<<endl;*/
       if((chunk&bit)!=0 && (chunk2&bit)!=0){ret(id)=2;}
       else if((chunk&bit)!=0){ret(id)=1;}
       else{ret(id)=0;}
       id++;
    }
   return ret;
 }

void print_bits(unit u)
  {
    std::bitset<sz> set(u);
    std::cout<<set;
  }

void print_bits(const Col<unit>& v)
 {
  for(int i=0;i<v.n_elem;i++){
     print_bits(v(i)); 
     std::cout<<std::endl;
   }
  std::cout<<std::endl;
 }

void scale(arma::mat& mt){
  arma::colvec mtMean = arma::mean(mt, 1);
  arma::colvec mtStd = arma::stddev(mt, 1, 1);
  mt.each_col() -= mtMean;
  mt.each_col() /=mtStd;
}
