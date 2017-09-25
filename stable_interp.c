/*  Based on
Federico Simmross-Wattenberg, Marcos Martín-Fernández, Pablo Casaseca-de-la-Higuera, Carlos Alberola-López. Fast calculation of alpha-stable density functions based on off-line precomputations. Application to ML parameter estimation. Digital Signal Processing, 38, 1-12, 2015.
 
 Federico Simmross-Wattenberg, Marcos Martín-Fernández, Pablo Casaseca-de-la Higuera, Carlos Alberola-López. Fast calculation of stable density functions based on off-line precomputations. Sample implementation,
 http://es.mathworks.com/matlabcentral/fileexchange/44576-fast-calculation-of-stable-density-functions-based-on-off-line-precomputations
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#include "stable.h"
#include "stable_interp.h"


#include "tri_precalc.h"
#include "grid_ab_precalc.h"


#include <gsl/gsl_spline.h>



void create_global_grid_tri (tri_grid * GRID){
    
    int i;
    for(i = 0; i < DT_NUM_TRI; i++){
        GRID[i].v1.x = grid_alpha[DelaunayTri[i][0]-1];
        GRID[i].v1.y = grid_beta[DelaunayTri[i][0]-1];
        GRID[i].v2.x = grid_alpha[DelaunayTri[i][1]-1];
        GRID[i].v2.y = grid_beta[DelaunayTri[i][1]-1];
        GRID[i].v3.x = grid_alpha[DelaunayTri[i][2]-1];
        GRID[i].v3.y = grid_beta[DelaunayTri[i][2]-1];
        GRID[i].invDet = 1/((GRID[i].v2.y - GRID[i].v3.y)*(GRID[i].v1.x - GRID[i].v3.x) + (GRID[i].v3.x - GRID[i].v2.x)*(GRID[i].v1.y - GRID[i].v3.y));
        
    }
}

void create_global_grid_pdf (pdf_grid * GRID, char * filename_pdf){
    
    FILE * f_grid_pdf;
    double buff;
    int count_pdf;
    int offset = 0;
    
    if ((f_grid_pdf = fopen(filename_pdf,"rb")) == NULL)
        perror("Error al abrir el fichero de pdfs");
    
    
    for (count_pdf = 0; count_pdf < NUM_POINTS; count_pdf++)
    {
        
        fread(&buff, sizeof(double),1,f_grid_pdf);
        offset = buff * sizeof(double);
        GRID[count_pdf].tam = (int)buff;
        fread (GRID[count_pdf].xtable, offset, 1, f_grid_pdf);
        fread (GRID[count_pdf].ytable, offset, 1, f_grid_pdf);
    }

    fclose(f_grid_pdf);
}


double get_ms_time()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    
    return (double) t.tv_sec * 1000 + (double) t.tv_usec / 1000;
}


int bary_coord (Point v1, Point v2, Point v3, double x, double y,int i) {

    double lambda1;
    double lambda2;
    double lambda3;
    
    int tri_found = 0;

    
    lambda1 = ((v2.y - v3.y) * (x - v3.x) + (v3.x - v2.x) * (y - v3.y)) * GRID_TRI[i].invDet;
    if (lambda1 >= 0 &&  lambda1 <= 1){
        lambda2 = ((v3.y - v1.y) * (x - v3.x) + (v1.x - v3.x) * (y - v3.y)) * GRID_TRI[i].invDet;
        if (lambda2 >= 0 &&  lambda2 <= 1){
            
            lambda3 = 1.0 - lambda1 - lambda2;
            
            if (lambda3 >= 0 &&  lambda3 <= 1){

            resultTri.tri_idx = i + 1;
            resultTri.bar_coords[0] = lambda1;
            resultTri.bar_coords[1] = lambda2;
            resultTri.bar_coords[2] = lambda3;
            tri_found = i + 1;

            }
        }

    }
    return tri_found;
}


void tri_location (double alpha, double absbeta)

{
    int i, idtri;
    
   #pragma omp parallel for private (idtri)
    for (i = 0; i < DT_NUM_TRI; i++){
        idtri = bary_coord (GRID_TRI[i].v1, GRID_TRI[i].v2, GRID_TRI[i].v3, alpha, absbeta,i);

    }
    
}


void unnormalize_stored_pdf(double sigma, pdf_grid * pdf)
{
    int i;
    
    for (i = 0; i < pdf->tam; i++){
        pdf->xtable[i]= pdf->xtable[i] * sigma;
    }
}

    
double * cubic_interp(pdf_grid * pdf,double * x, const int length)
{
    
    double * yi_cubic;
    int i;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, pdf->tam);
    gsl_spline_init (spline, pdf->xtable, pdf->ytable, pdf->tam);
    
    yi_cubic = malloc(length * sizeof(double));
    
    for (i = 0; i < length; i++)
    {
        double xi = *(x + i);
        yi_cubic[i] = gsl_spline_eval(spline, xi, acc);
        
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return yi_cubic;
    
}


void stable_pdf_interp(StableDist *dist, const double data[], const int length, double *pdf){
    
    
    double * x;
    pdf_grid pdf1, pdf2, pdf3;
    int i;
    double * pdf1_interp;
    double * pdf2_interp;
    double * pdf3_interp;
    double sigma = dist->sigma;
    
    
    resultTri.tri_idx = 0;
    tri_location (dist->alpha,fabs(dist->beta));
    
    
    if(resultTri.tri_idx==0){
        printf("stable_pdf_interp: Requested point is outside the grid (alfa=%1.10f, beta=%1.10f)\n",dist->alpha,dist->beta);
        for (i = 0; i < length; i++){
            pdf[i]= 0;
        }
        return;
    }

    x = malloc(length*sizeof(double));
    memcpy (x, data, length*sizeof(double));
    
    for (i = 0; i < length; i++){
        *(x+i) = *(x + i) - dist-> mu_0;
        if(dist->beta < 0)
        {
            *(x + i) = - (*(x + i));
        }
        
    }
    
    
    #pragma omp parallel sections
    {
       #pragma omp section
       {
           pdf1 = GRID_PDF[DelaunayTri[resultTri.tri_idx-1][0]-1];
           unnormalize_stored_pdf(sigma, &pdf1);
           pdf1_interp = cubic_interp(&pdf1,x,length);
       }
     
      #pragma omp section
      {
     
          pdf2 = GRID_PDF[DelaunayTri[resultTri.tri_idx-1][1]-1];
          unnormalize_stored_pdf(sigma, &pdf2);
          pdf2_interp = cubic_interp(&pdf2,x,length);
      }
     
     #pragma omp section
     {
    
         pdf3 = GRID_PDF[DelaunayTri[resultTri.tri_idx-1][2]-1];
         unnormalize_stored_pdf(sigma, &pdf3);
         pdf3_interp = cubic_interp(&pdf3,x,length);
     }
    }
    
 
    
    
    for (i = 0; i < length; i++)
    {
        pdf [i] = (resultTri.bar_coords[0]*pdf1_interp[i] + resultTri.bar_coords[1]*pdf2_interp[i] + resultTri.bar_coords[2]*pdf3_interp[i])/sigma;

    }
    
    free(x);
    free(pdf1_interp);
    free(pdf2_interp);
    free(pdf3_interp);
    
}


void vector_step(double **x, double min, double max, double step, int * n)

{
    int i,m;
    double aux;
    
    aux=(max-min)/step;
    
    if (aux<0 || isnan(aux) || isinf(aux))
    {
        *n=0;
        (*x)=NULL;
        printf("Warning: Empty vector");
        return;
    }
    
    m=(int)(aux)+1;
    (*x)=(double*)malloc(m*sizeof(double));
    
    if ((*x) == NULL)
    {
        perror("Error while creating x array");
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<m;i++)
    {
        (*x)[i]=min+i*step;
    }
    
    *n=m;
    
    return;
}
