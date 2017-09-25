
#ifndef STABLE_INTERP_H
#define STABLE_INTERP_H

#include "stable.h"

#define NUM_POINTS 6076
#define DT_NUM_TRI 11709
#define MAX_PDF 210

typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    Point v1;
    Point v2;
    Point v3;
    double invDet;
} tri_grid;


typedef struct {
    int tam;
    double xtable [MAX_PDF];
    double ytable [MAX_PDF];
} pdf_grid;


typedef struct {
    double bar_coords[3];
    int tri_idx;
} TriFound;


tri_grid GRID_TRI[DT_NUM_TRI];
pdf_grid GRID_PDF[NUM_POINTS];

TriFound resultTri;


double get_ms_time();
void create_global_grid_tri (tri_grid GRID[]);
void create_global_grid_pdf (pdf_grid GRID[], char * filename_pdf);

int bary_coord (Point v1, Point v2, Point v3, double x, double y,int i);
void tri_location (double alfa, double absbeta);

void unnormalize_stored_pdf(double sigma, pdf_grid * pdf);
double * cubic_interp(pdf_grid * pdf, double * x, const int length);
void stable_pdf_interp(StableDist *dist, const double data[], const int length, double * pdf);

void vector_step(double **x, double min, double max, double step, int * n);

#endif
