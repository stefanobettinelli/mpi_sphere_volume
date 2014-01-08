/*
 TODO: 1) incapsulare centro e raggio in qualche modo provare a stimare il volume delle sfere con montecarlo
 DONE:calcolare il Vbb con i dati presi dal file
 DONE:riesco a leggere i dati dal file
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

double max(double a, double b){ return (a) > (b) ? (a) : (b); }
double min(double a,double b){ return (a) < (b) ? (a) : (b); }

/* mi serve un intervallo semiapert [min, max) */
double rand_point(double min, double max)
{
	int base_random = rand(); /* in [0, RAND_MAX] */
	if (RAND_MAX == base_random) return rand_point(min, max);
												
	double range = max - min;
	double remainder = RAND_MAX % (int)range;
	double bucket = RAND_MAX / range;
	
	if (base_random < RAND_MAX - remainder)
	{
		return min + base_random/bucket;
	}
	else
	{
		return rand_point(min, max);
	}
}

int main( int argc, char *argv[])
{
	if(argc != 3)
	{
		fprintf(stderr,"usage:\n\tmpirun -n NN sfere <filename> <pointsnumber>\n");
		exit(1);
	}
    
	FILE* sfere_file = fopen(argv[1],"r");
	if( sfere_file == NULL )
	{
		fprintf(stderr,"Error opening sfere file\n");
		exit(1);
	}
    
	int i = 0;
    
	double x_max, y_max, z_max, x_min, y_min, z_min;
	x_max = y_max = z_max = DBL_MIN;
	x_min = y_min = z_min = DBL_MAX;
    
	double Vbb = 0;
	double passo = 0;
	int sfere_n = 0;
	double cx, cy, cz = 0;
	double raggio = 0;
	double r_points[sfere_n][3];
	
	/*leggo il contenuto del file il numero di sfere e le coordinate dei centri+lunghezza_raggio */
	fscanf(sfere_file,"%d", &sfere_n);
	while(fscanf(sfere_file,"%lf %lf %lf %lf",&cx, &cy, &cz, &raggio) != EOF)
	{
		x_max = max(x_max,cx + raggio);
		x_min = min(x_min,cx - raggio);
		y_max = max(y_max,cy + raggio);
		y_min = min(y_min,cy - raggio);
		z_max = max(z_max,cz + raggio);
		z_min = min(z_min,cz - raggio);
		i++;
	}
	
   /*calcolo il volume del bounding box, lo divido per il numero di punti per ottenere il passo*/
	Vbb = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	passo = Vbb / sfere_n;
	
	/*genero i punti casuali all'interno del bouding box usando il passo*/
	srand((unsigned int)time(NULL));
	for (i=0; i<sfere_n; i++)
	{
		r_points[i][0] = rand_point(passo*i, passo*(i+1));
		r_points[i][1] = rand_point(passo*i, passo*(i+1));
		r_points[i][2] = rand_point(passo*i, passo*(i+1));
		printf("[%lf | %lf | %lf]\n", r_points[i][0], r_points[i][1], r_points[i][2]);
	}
    
	printf("(x_max %f;x_min %f), (y_max %f;y_min %f), (z_max %f;z_min %f)\n",x_max,x_min,y_max,y_min,z_max,z_min);
	printf("Vbb %f\n", Vbb);
	return 0;
}
