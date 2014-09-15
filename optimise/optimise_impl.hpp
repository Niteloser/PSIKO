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


double cdf(double x, double m, double s)
{
   return 0.5 * erfc(-(x-m)/(s*sqrt(2)));
}

void initialise(const arma::mat& observations,
                  arma::mat& means)
{
  arma::Col<size_t> assignments;
  arma::vec d(observations.n_cols*(observations.n_cols-1)/2); 
  int k=0,n=observations.n_cols,el_max=-1;
  double mx=-1.0;
  for(int i=0;i<observations.n_cols-1;i++)
   for(int j=i+1;j<observations.n_cols;j++)
    {
    arma::vec dif=observations.unsafe_col(i)-observations.unsafe_col(j);
          d(k)=arma::norm(dif,2);
    if(d(k)>mx)
            {
    mx=d(k);
      el_max=i;
      }
    k++;
  }
  double mn=arma::mean(d);
  double s=arma::stddev(d);

  double cond=0.95;
  arma::vec klist(means.n_cols);
  klist(0)=el_max;
  k=1;
  while(k<means.n_cols)
    {
     double mx_sm=-1;
     int mx_el=-1;
     for(int i=0;i<n;i++)
       {
        arma::vec obs=observations.unsafe_col(i);
        bool sw=false;
  double sm=0;
  for(int j=0;j<k;j++)
   {
          arma::vec k_el=observations.unsafe_col(klist(j));
    double dist=arma::norm(obs-k_el,2);
    if(cdf(dist,mn,s)<cond)
             {
              sw=true;
        break;
             }
    sm+=dist;
   }
         if(sw)continue;
         if(mx_sm<sm)
     {
       mx_el=i;
             mx_sm=sm;
           }
        }
      if(mx_el==-1)
     {
      cond-=0.05;
      continue;
     }
      klist(k)=mx_el;
      k++;
      cond=0.95;
    }
  for(int i=0;i<klist.n_elem;i++)
    {
      means.unsafe_col(i)=observations.col(klist(i));
    }
}

void my_optimize(arma::mat& means, arma::mat& X,double alpha)
  {
      double prev=0.0;
      double cur=0.0;
      int iter=0;
      arma::mat Q(means.n_cols,X.n_cols,arma::fill::zeros);
      do
       {  //1: estimate Q-matrix
          arma::mat A= arma::ones(1,means.n_cols); 
          A.insert_rows(1,means);
          A=arma::inv(A);
          for(int i=0;i<Q.n_cols;i++)
          {
            arma::vec col=Q.unsafe_col(i);
            arma::vec c=X.unsafe_col(i);
            col=inferAncestry(c,A);
          }
          //2: estimate means
          arma::mat I=0.5*alpha*arma::eye(means.n_cols,means.n_cols);
          means=arma::inv(Q*Q.t()+I.t()*I)*Q*X.t();
          means=means.t();
          arma::mat content= arma::ones(1,means.n_cols);
          content.insert_rows(1,means);
          prev=cur;
          double penalty=1.0/(Factorial(means.n_rows)) * fabs(arma::det(content));
          cur=pow(arma::norm((X-means*Q),"fro"),2);
          //if(iter%20==0)std::cout<<"objective function "<<cur<<"\n";
          iter++;
       }while(abs(prev-cur)>1e-5 && iter<1000);
  }
