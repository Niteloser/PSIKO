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

#include <armadillo>
#include <ctime>
#include <limits.h>
#include <algorithm>
#include <stdint.h>

using namespace std;
using namespace arma;

typedef  uint64_t unit;
extern const unit sz=sizeof(unit)*CHAR_BIT;

int get_cols(char* fname);

int get_rows(char* fname,int nc);

bool read_snps(Mat<unit>& data, char* fname);

#include "bitwise_input_impl.hpp"