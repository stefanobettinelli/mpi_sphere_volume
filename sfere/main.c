#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#define MASTER 0
#define COORDS 3

double max(double a, double b){ return (a) > (b) ? (a) : (b); }
double min(double a,double b){ return (a) < (b) ? (a) : (b); }

/* mi serve un intervallo semiapert [min, max) */
double rand_point(double min, double max)
{
	double f = (double)rand() / RAND_MAX;
	double dart = min + f * (max - min);
	if (dart == max) {
		rand_point(min, max);
	}
	return dart;
}

int hit(double x, double y, double z, double cx, double cy, double cz, double r)
{
	int hit = 0;
	if( (pow(x-cx,2)+pow(y-cy,2)+pow(z-cz,2)) <= pow(r,2) ) hit=1;
	return hit;
}

int main( int argc, char *argv[])
{
	if(argc != 3)
	{
		fprintf(stderr,"usage:\n\tmpirun -n NN sfere <filename> <pointsnumber>\n");
		return -1;
	}
    
	FILE* sfere_file = fopen(argv[1],"r");
	if( sfere_file == NULL )
	{
		fprintf(stderr,"Error opening sfere file\n");
		return -1;
	}
	
	/* Obtain number of tasks and task ID */
	int numtasks, taskid;
	int rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	
	if ((atoi(argv[2]) % numtasks) != 0) {
		fprintf(stderr, "Attenzione in numero punti deve essere un multiplo del numero di processi\n");
		MPI_Finalize();
		return -1;
	}
    
	int i = 0;
	int j = 0;
	double x_max, y_max, z_max, x_min, y_min, z_min;
	x_max = y_max = z_max = DBL_MIN;
	x_min = y_min = z_min = DBL_MAX;
	double Vbb = 0;
	int hits = 0;
	double passo_x = 0;
	double passo_y = 0;
	double passo_z = 0;
	int sfere_n = 0;
	double cx, cy, cz = 0;
	double raggio = 0;
	int points = atoi(argv[2]);
	double r_points[points][3];
	
	/*leggo il contenuto del file il numero di sfere e le coordinate dei centri+lunghezza_raggio */
	fscanf(sfere_file,"%d", &sfere_n);
	/*
	 calcolo la radice cubica del numero di punti per ottenere
	 il numero di intervalli in cui dividere gli spigoli del bounding box
	*/
	int intervals = ceil(cbrt(points));
	if( intervals == 1 ){
		fprintf(stderr, "È stato fornito solo 1 punto...esco\n");
		return -1;
	}
	
	/*matrice di sfere cx,cy,cz,raggio NB: tutti i task hanno a disposizione la matrice*/
	double sfere[sfere_n][4];
	/*nel ciclo vengono calcolate le coordinare degli angoli estermi 
	 del parallelepipedo minimo contenete le sfere*/
	while(fscanf(sfere_file,"%lf %lf %lf %lf",&cx, &cy, &cz, &raggio) != EOF)
	{
		x_max = max(x_max,cx + raggio);
		x_min = min(x_min,cx - raggio);
		y_max = max(y_max,cy + raggio);
		y_min = min(y_min,cy - raggio);
		z_max = max(z_max,cz + raggio);
		z_min = min(z_min,cz - raggio);
		sfere[i][0] = cx;
		sfere[i][1] = cy;
		sfere[i][2] = cz;
		sfere[i][3] = raggio;
		i++;
	}
	
	/*Vbb: volume del bounding box. 
	 Lo divido per il numero di punti per ottenere il passo,
	 il passo è semplicemente la lunghezza degli intervalli in cui divido gli spigoli del BB
	 */
	Vbb = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	passo_x = (x_max-x_min) / intervals;
	passo_y = (y_max-y_min) / intervals;
	passo_z = (z_max-z_min) / intervals;
	
	/*genero i punti casuali all'interno del bouding box usando il passo*/
	if(taskid == MASTER)
	{
		srand((unsigned int)time(NULL));
		/*
		for (i=0; i<points; i++) {
			if( i % (intervals-1) == 0 ){
				r_points[i][0] = rand_point(x_min, x_min + passo_x);
				r_points[i][1] = rand_point(y_min, y_min + passo_y);
				r_points[i][2] = rand_point(z_min, z_min + passo_z);
			}
			else{
				r_points[i][0] = rand_point(x_min+(passo_x*(i%intervals)), x_min+(passo_x*(i%intervals)) + passo_x);
				r_points[i][1] = rand_point(y_min+(passo_y*(i%intervals)), y_min+(passo_y*(i%intervals)) + passo_y);
				r_points[i][2] = rand_point(z_min+(passo_z*(i%intervals)), z_min+(passo_z*(i%intervals)) + passo_z);
			}
		}*/
		/*
		int temp_intervals = intervals;
		int _min = x_min;
		int _max = x_min + passo_x;
		for (i=0; i<points; i++) {
			r_points[i][0] = rand_point(_min,_max);
			printf("[%d...%d] ==> %f\n",_min,_max,r_points[i][0]);
			temp_intervals--;
			if (temp_intervals == 0) {
				_min = x_min;
				_max = x_min + passo_x;
				temp_intervals = intervals;
			}
			else{
				_min = _max;
				_max = _min + passo_x;
			}
			
		}
		_min = y_min;
		_max = y_min + passo_y;
		temp_intervals = intervals;
		for (i=0; i<points; i++) {
			r_points[i][1] = rand_point(_min,_max);
			temp_intervals--;
			if (temp_intervals == 0) {
				_min = y_min;
				_max = y_min + passo_y;
				temp_intervals = intervals;
			}
			else{
				_min = _max;
				_max = _min + passo_y;
			}
		}
		_min = z_min;
		_max = z_min + passo_z;
		temp_intervals = intervals;
		i=0;
		for (i=0; i<points; i++) {
			r_points[i][2] = rand_point(_min,_max);
			temp_intervals--;
			if (temp_intervals == 0) {
				_min = z_min;
				_max = z_min + passo_z;
				temp_intervals = intervals;
			}
			else{
				_min = _max;
				_max = _min + passo_z;
			}
			
		}*/
		for (i=0; i<points; i++) {
			r_points[i][0] = rand_point(x_min,x_max);
			r_points[i][1] = rand_point(y_min,y_max);
			r_points[i][2] = rand_point(z_min,z_max);
		}
		
	}
	
	int rows = (int)(points / numtasks);
	double task_points[rows][3];
	MPI_Scatter(r_points,rows*COORDS,MPI_DOUBLE,task_points,rows*COORDS,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
	
	for (i = 0; i<rows; i++) {
		//printf("%d [%lf | %lf | %lf]\n", taskid, task_points[i][0], task_points[i][1], task_points[i][2]);
		for (j=0; j<sfere_n; j++) {
			if( hit(task_points[i][0], task_points[i][1], task_points[i][2], sfere[j][0], sfere[j][1], sfere[j][2], sfere[j][3]) == 1 ){
				hits++;
				break;
			}
		}
	}
	printf("hit %d from task %d\n",hits,taskid);
	
	int* hits_master = NULL;
	//if(taskid == MASTER){
		hits_master = (int*) malloc(numtasks*sizeof(int));
	//}
	MPI_Gather(&hits, 1, MPI_INT, hits_master, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	i=0;
	int total_hits = 0;
	if (taskid == MASTER) {
		for (i=0; i<numtasks; i++){
			total_hits += hits_master[i];
			printf("hits_array[%d] = %d\n",i ,hits_master[i]);
		}
		printf("Vbb %lf\ntotal_hits %d\npoints %d\nEST VOLUME %lf from task %d\n",Vbb,total_hits,points,Vbb*(1.0*total_hits/points), taskid);
	}
	MPI_Finalize();
	free(hits_master);
	return 0;
}
