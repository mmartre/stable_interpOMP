#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>

#include "stable.h"
#include "stable_interp.h"



int main(int argc, char *argv[])
{
	double alpha, beta, sigma, mu_0;
	double *data;
	int i;
	double start, end, ms_duration = 0;

	char * data_fname = "data_10e3.txt";
	int N = 1000;
    int Nexp = 100;
	
	StableDist *dist = NULL;
	
    char name_gridPDF[25]="grid_pdf.bin";

   	create_global_grid_tri (GRID_TRI);
    create_global_grid_pdf (GRID_PDF, name_gridPDF);

	alpha = 1.5;
	beta = 0.5;
	sigma = 1;
	mu_0 = 0;


	if ((dist = stable_create(alpha, beta, sigma, mu_0, 0)) == NULL) {
		printf("Error when creating the distribution");
		exit(1);
	}


    data = malloc(N* sizeof(double));
    data = load_rand_data(data_fname, N);
    
    for (i=0; i< Nexp; i++){
	
		start = get_ms_time();
		stable_fit_init(dist, data, N, NULL, NULL);
		stable_fit_mle(dist, data, N);
		end = get_ms_time();

		ms_duration += end - start;
	

    }
    ms_duration = ms_duration/Nexp;

    printf(" %d %lf %lf %lf %lf %lf\n", N, ms_duration, dist->alpha, dist->beta, dist->sigma, dist->mu_0);

	

free(data);

stable_free(dist);

return 0;
}
