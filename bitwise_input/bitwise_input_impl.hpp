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

int get_cols(char* fname)
{
  FILE *f = fopen(fname,"r");
  int cols=0;
  char c=fgetc(f);
  while ((c!=EOF) && (c!=10)) {
	  cols++;
	  c = fgetc(f);
	}
  return cols;
}


int get_rows(char* fname,int nc)
{
  FILE *f = fopen(fname,"r");
  int rows=0;
  char *bf=(char*) calloc(nc*3, sizeof(char));
  char *l=fgets(bf,nc*3,f);
  while (!feof(f) && l) {
	  rows++;
	  l = fgets(bf,nc*3,f);
	}
  free(bf);
  return rows;
}


bool read_snps(Mat<unit>& data, char* fname)
{     
     int nc=get_cols(fname);
     int nr=get_rows(fname,nc);
     FILE *f = fopen(fname,"r");
     cout<<nr<<" rows X "<<nc<<" columns\n";
     data.set_size((nr/sz+1)*2,nc);
     data.fill(0);
     char *bf=(char*) calloc(nc*3, sizeof(char));
     char *l=fgets(bf,nc*3,f);
     int i=0;
     while (i<nr && !feof(f)) {
          int j=0;
	        unit line=i/sz;//actual line
          unit msk=((unit)1)<<(i%sz);//bit in 'unit' stored at actual line
          for(char* it = l; j<nc; ++it) {
              char c=*it;
              
              if((short)(c-'0')>0 && c!='9'){//leave missing data set to zero
	         data.at(line,j)|=msk;//set for individual j if nonzero
                   }
              if(c=='2')
                {
                 data.at(line+(nr/sz+1),j)|=msk;
                }
              
              j++;
                  }

          i++;
          l = fgets(bf,nc*3,f);
	} 
     free(bf);
     return true;
}