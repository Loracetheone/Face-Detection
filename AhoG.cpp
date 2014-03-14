#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include <vector>
//#include "stdafx.h"
#include <stdio.h>
//#include <iostream>
//#include "matrix.h"

using namespace std;   
//FILE *fpr; 

void HoG(double *pixels, double *params, int *img_size, double *dth_des)
{    
	double block_norm, pi = 3.1415926536;    
    int nb_bins = 9, bloques = 11, block_size = (int) params[2], orient = (int) params[3], img_width  = 50, img_height = 50, des_indx = 0, indx=0;    
	double Gsx[50][50]={0},Gsy[50][50]={0}, Gx[50][50]={0},Gy[50][50]={0},Orient=0, Magni[50][50]={0},vgsx,vgsy,vmagni;
          
    //Calculate gradients (zero padding)	
	   for( int y=0; y<img_width; y++)
	   {    //mascara en Y:[-1 0 1]'        
			for( int x=0; x<img_height; x++)
			{     //mascara en X: -1 0 1       
				if(x==0) Gx[y][x] = pixels[y +(x+1)*img_height];
				else
				{
					if (x==img_width-1) Gx[y][x] = -pixels[y + (x-1)*img_height];
					else Gx[y][x] = pixels[y+(x+1)*img_height] - pixels[y + (x-1)*img_height];
				}
				if(y==0) Gy[y][x] = -pixels[y+1+x*img_height];
				else
				{
					if (y==img_height-1) Gy[y][x] = pixels[y-1+x*img_height];
					else Gy[y][x] = -pixels[y+1+x*img_height] + pixels[y-1+x*img_height];
				}				
				Magni[y][x] = sqrt(pow(Gx[y][x],2)+pow(Gy[y][x],2));			
				Gsx[y][x]= pow(Gx[y][x],2)-pow(Gy[y][x],2);		
				Gsy[y][x]= 2*Gy[y][x]*Gx[y][x];				
			}
		}
	     
	   int i0, i1, j0, j1;
	   for(int x=0; x<img_width-bloques; x+=4)
	   {  
			for(int y=0; y<img_height-bloques; y+=4)
			{				
				//Dividimos en bloques a la imagen 11x11
				double bin[9]={0};

				for(int i=0; i<bloques; i++)
				{
					for(int j=0; j<bloques; j++)
					{
						//Calculamos las sumatorias con los vecinos locales 3x3 de cada pixel
						if ((x+i-1) >= 0) i0 = x+i-1;
						else   i0 = 0;

						if ((x+i+1) < img_width)  i1 = x+i+1;
						else i1 = img_width-1;

						if ((y+j-1) >= 0)  j0 = y+j-1;
						else j0= 0;

						if ((y+j+1) < img_height)  j1 = y+j+1;
						else j1 = img_height-1;						

						vgsx=0; vgsy=0; vmagni=0; 

						//Dividimos en bloques a la imagen 3x3
						for(int k=i0; k<=i1; k++)
						{
							for(int w=j0; w<=j1; w++)
							{
								vgsx+= Gsx[k][w];
								vgsy+= Gsy[k][w];
								vmagni+= pow(Magni[k][w],2);
							}
						}

						//Calculamos la orientacion promedio y Convertir de radianes a grados
						Orient= atan2(vgsy,vgsx)*180/pi;	
			
						//quantizamos la orientacion promedio son 180 grados máximo en 9 bins
						indx= (int)floor((Orient+180)/40)-1;
						bin[indx]+= vmagni;
					}
				}
			//Hacemos la normalizacion del histograma
			   block_norm=0;
			   
			   for(int i=0;i<nb_bins;i++)
				   block_norm+= bin[i]*bin[i];
			   
			   block_norm= sqrt(block_norm);
			   if (block_norm>0)
			   {
				   for(int i=0;i<nb_bins;i++)
				   {
						bin[i]= bin[i]/block_norm;
						if (bin[i]>0.2)
							bin[i]=0.2;
				   }
			   }
				//Concatenamos a los histogramas
				for (int i=0; i<nb_bins; i++)
				{
					dth_des[des_indx]= bin[i];
					des_indx++;
				}				
			}
	   }	 
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *pixels, *dth_des, *params;
    int nb_bins, block_size;
    int img_size[2];
        
    if (nlhs>1)  mexErrMsgTxt("Too many output arguments");
    if (nrhs==0) mexErrMsgTxt("No Image -> No HoG");
    
    if (mxGetClassID(prhs[0])!=6) mexErrMsgTxt("Matrix is not of type double");
    
    pixels     = mxGetPr(prhs[0]);    
    
    img_size[0] = mxGetM(prhs[0]);
    img_size[1]  = mxGetN(prhs[0]);
    
    if (nrhs>1){
        params     = mxGetPr(prhs[1]);
        if (params[0]<=0) mexErrMsgTxt("Number of orientation bins must be positive");
        if (params[1]<=0) mexErrMsgTxt("Cell size must be positive");
        if (params[2]<=0) mexErrMsgTxt("Block size must be positive");
    }
    else {
        params = new double[5];
        params[0]=8;
        params[1]=11;
        params[2]=2;
        params[3]=1;
        params[4]=0.2;
    }
    
    nb_bins       = (int) params[0];    
    block_size    = (int) params[2];     

	//int tam= pow(round((img_size[0]-params[1])/4),2)*params[0];
  
    plhs[0] = mxCreateDoubleMatrix(900, 1, mxREAL);
    dth_des = mxGetPr(plhs[0]);
    
    HoG(pixels, params, img_size, dth_des);
    if (nrhs==1) delete[] params;
}
