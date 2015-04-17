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

    arma::cube findStds(arma::mat projected, arma::mat Q)
      {
        arma::cube ret(projected.n_rows,projected.n_rows,Q.n_cols);
        for(int j=0;j<Q.n_cols;j++)
        {
          arma::mat localData(projected.n_rows,1);
          for(int i=0;i<Q.n_rows;i++)
             {
                if(Q(i,j)<0.9)continue;
                localData.insert_cols(0,projected.unsafe_col(i));
             }
          localData.shed_col(localData.n_cols-1);
          //cout<<j<<" subset\n";
          //cout<<localData;
          arma::mat cv=arma::cov(localData.t(),1);
          //cout<<"cv\n";
          //cout<<cv;
          ret.slice(j)=cv;
          
        }
        return ret;
      }

    arma::mat findMeans(arma::mat projected, arma::mat Q)
      {
        arma::mat ret(projected.n_rows,Q.n_cols);
        for(int j=0;j<Q.n_cols;j++)
        {
          arma::mat localData(projected.n_rows,1);
          for(int i=0;i<Q.n_rows;i++)
             {
                if(Q(i,j)<0.9)continue;
                localData.insert_cols(0,projected.unsafe_col(i));
             }
          localData.shed_col(localData.n_cols-1);
          /*cout<<j<<" subset\n";
          cout<<localData;*/
          ret.unsafe_col(j)=arma::mean(localData, 1);
        }
        return ret;
      }

    arma::vec applyWindow(const arma::mat& dataWindow, arma::mat& evec, arma::mat& evals,arma::mat Q){
        //compute principal components
        arma::vec ret(dataWindow.n_cols);
        arma::mat W=dataWindow*evec;
        W.each_row() /= sqrt(evals.t());
        //project data window onto PCs
        arma::mat projected=W.t()*dataWindow;
        //scale(projected);
        arma::mat means=findMeans(projected,Q);
        arma::cube stds=findStds(projected,Q);
        /*cout<<"means\n";
        cout<<means;
        cout<<"stds\n";
        cout<<stds;
        cout<<"projected\n";
        for(int i=0;i<10;i++)
          cout<<projected.col(i).t();*/

        for(int i=0;i<dataWindow.n_cols;i++)
          {
            double mxProb=0;
            int mxAncestor=-1;
            for(int j=0;j<means.n_cols;j++)
              {
                double prob=mvn(projected.col(i),means.col(j),stds.slice(j));
                if(prob>mxProb)
                 {
                    mxProb=prob;
                    mxAncestor=j;
                    //if(i==21){cout<<"current max "<<mxAncestor<<"with mxProb "<<mxProb<<" and prob "<<prob<<"\n";}
                 }
               //if(i==21){cout<<"prob("<<i<<","<<j<<")"<<prob<<endl;}                 
              }
            //if(i==21){cout<<"max ancestor= "<<mxAncestor<<"\n";}
            ret(i)=mxAncestor;
          }

        delete &dataWindow;

        return ret;
    }