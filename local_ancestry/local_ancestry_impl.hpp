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

    arma::mat applyWindow(const arma::mat& dataWindow, arma::mat& evec, arma::mat& evals,arma::mat Q){
        //compute principal components
        arma::mat W=dataWindow*evec;
        W.each_row() /= sqrt(evals.t());
        //project data window onto PCs
        arma::mat projected=W.t()*dataWindow;
        scale(projected);

        arma::mat means=findMeans(projected,Q);
        arma::mat stds=findStds(projected,Q);

        delete &dataWindow;

        return projected;
    }