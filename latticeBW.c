/***********************************INSTRUCTIONS********************************
                           
Author: Tiago Mendes-Santos (ICTP) Date: 04/2018
Code to generate the ALPS lattice.xml file for the 2D BW-EH Heisenberg model 

Parameters
"EH_or_not"
(1) EH_or_not = 0 - homogeneous lattice
(2) EH_or_not = 1 - BW-EH cylinder
(3) EH_or_not = 2 - BW-EH torus

"mm"
Cutoff for the SSE expansion. 
The choice of mm must be consistent with the value of $\beta_{EH}$

Instructions to run the ALPS code: 
http://alps.comp-phys.org/mediawiki/index.php/Documentation:qwl
http://alps.comp-phys.org/mediawiki/index.php/ALPS_2_Tutorials:MC-06_QWL

Basic commands
./parameter2xml "d2dlx*ly*EH*"
./qwl "d2dlx*ly*EH*.in.xml"
./qwl_evaluate --T_MIN 0.263 --T_MAX 10 --DELTA_T  0.1 "d2lx*ly*EH*.task*.out.xml"


*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>


void main()
{
 int s,x1,x2,y1,y2,limx,limy,s_aux,nn,lx,ly,bcond,EH_or_not,nb,mm,dimension=2;
 double coupling;
 char out1[100],out2[100];


 lx=6;
 ly=6;
 EH_or_not=1;
 
 bcond=1;
 mm= 1000;
 
 FILE *arq=NULL;
 FILE *arq2=NULL;
 
 sprintf(out1,"graphd%dlx%dly%dEH%d.xml",dimension,lx,ly,EH_or_not);
 sprintf(out2,"d%dlx%dly%dEH%d",dimension,lx,ly,EH_or_not);

 arq = fopen(out1,"w"); 
 arq2 = fopen(out2,"w");

 nn=lx*ly;
 
 if(EH_or_not == 1 || EH_or_not == 2){
       nb = 2*nn-ly;
       limx=2;
       limy=1;
 }
 if(EH_or_not == 0){
       if(bcond == 1){
         nb = 2*nn - ly;
         limx=2;
         limy=1;
       }
       if(bcond == 0){
         nb = 2*nn;
         limx=1;
         limy=1;
       }
 }
 
 fprintf(arq2,"LATTICE_LIBRARY=\"%s\"\n",out1);
 fprintf(arq2,"MODEL_LIBRARY=\"model.xml\"\n"); 
 fprintf(arq2,"GRAPH=\"graph\"\n"); 
 fprintf(arq2,"MODEL=\"spin\"\n");
 fprintf(arq2,"USE_ZHOU_BHATT_METHOD = 1\n");
 fprintf(arq2,"MEASURE_MAGNETIC_PROPERTIES = 0\n");
 fprintf(arq2,"local_S = 1/2\n");
 fprintf(arq2,"L = \%d\n",ly);
 fprintf(arq2,"CUTOFF = %d\n",mm);
 
   
 fprintf(arq,"<LATTICES>\n");
 fprintf(arq,"<GRAPH name=\"graph\" dimension=\"%d\" vertices=\"%d\" edges=\"%d\">\n",dimension,nn,nb);
 
 for(y1=1;y1<=nn;y1++){
          x2=(y1-1)%lx;
          y2=(y1-1)/lx;
          fprintf(arq," <VERTEX id=\"%d\" type=\"0\"><COORDINATE>%d %d</COORDINATE></VERTEX>\n",y1,x2,y2);
 }   
 for(y1=0;y1<=ly-1;y1++){
	  for(x1=0;x1<=lx-limx;x1++){
	    s=1+x1+y1*(lx - limx + 1);
            s_aux = 1+x1+y1*lx;
	    x2=(x1+1)%lx;
	    y2=y1;
            if (EH_or_not == 1) coupling=x1+1.0;
            if (EH_or_not == 2) coupling=(((double)x1+1.0)*((double)lx - ((double)x1+1.0)))/(double)lx;
            fprintf(arq," <EDGE source=\"%d\" target=\"%d\" id=\"%d\" type=\"%d\" vector=\"1 0\"/>\n",s_aux,1+x2+y2*lx,s-1,s-1);
            fprintf(arq2,"J%d=%lf\n",s-1,coupling);
          }
  }

  for(y1=0;y1<=ly-limy;y1++){
	  for(x1=0;x1<=lx-1;x1++){
	    s = 1+x1+y1*lx;
	    x2=x1;
	    y2=(y1+1)%ly;
            if (EH_or_not == 1) coupling = (x1 + 1.0 - 0.5);
            if (EH_or_not == 2) coupling=(((double)x1+1.0-0.5)*((double)lx - ((double)x1+1.0-0.5)))/(double)lx;
            fprintf(arq," <EDGE source=\"%d\" target=\"%d\" id=\"%d\" type=\"%d\" vector=\"0 1\"/>\n",s,1+x2+y2*lx,s+nn-ly*(limx-1)-1,s+nn-ly*(limx-1)-1);
            fprintf(arq2,"J%d=%lf\n",s+nn-ly*(limx-1)-1,coupling);
	   }
  }
 
  fprintf(arq,"</GRAPH>\n");
  fprintf(arq,"</LATTICES>");
  
 
}
