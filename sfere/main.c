#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#define MASTER 0

double max(double a, double b){ return (a) > (b) ? (a) : (b); }
double min(double a,double b){ return (a) < (b) ? (a) : (b); }

/*genera un reale nell'intervallo [min, max] */
double rand_point(double min, double max)
{
	double f = (double)rand() / RAND_MAX;
	double dart = min + f * (max - min);
	return dart;
}
/*controllo se il punto generato (x,y,z) è all'interno di una sfera (cx,cy,cz,r)*/
int hit(double x, double y, double z, double cx, double cy, double cz, double r)
{
	int hit = 0;
	if( (pow(x-cx,2)+pow(y-cy,2)+pow(z-cz,2)) <= pow(r,2) ) hit=1;
	return hit;
}

int main( int argc, char *argv[])
{
	/*qualche controllo preliminare...*/
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
	
	/* inizio programma mpi */
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
	int sfere_n = 0;
	double cx, cy, cz = 0;
	double raggio = 0;
	int points = atoi(argv[2]);
	double r_points[points][3];
	
	/*leggo il contenuto del file il numero di sfere e le coordinate dei centri+lunghezza_raggio */
	fscanf(sfere_file,"%d", &sfere_n);
	
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
	
	/*Vbb: volume del bounding box*/
	Vbb = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
	/*genero i punti casuali all'interno del bouding box*/
	if(taskid == MASTER)
	{
		srand((unsigned int)time(NULL));
		for (i=0; i<points; i++) {
			r_points[i][0] = rand_point(x_min,x_max);
			r_points[i][1] = rand_point(y_min,y_max);
			r_points[i][2] = rand_point(z_min,z_max);
		}
		
	}
	/*rows sono le righe della matrice sottomesse ad ogni task*/
	int rows = (int)(points / numtasks);
	double task_points[rows][3];
	MPI_Scatter(r_points,rows*3,MPI_DOUBLE,task_points,rows*3,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
	/*ogni task controlla se la porzione di punti a lui assegnata ricade all'interno di almeno una sfera*/
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
	if(taskid == MASTER){
		hits_master = (int*) malloc(numtasks*sizeof(int));
	}
	/*roccolata degli hits nel MASTER TASK*/
	MPI_Gather(&hits, 1, MPI_INT, hits_master, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	int total_hits = 0;
	if (taskid == MASTER) {
		for (i=0; i<numtasks; i++){
			total_hits += hits_master[i];
		}
		/*Stampo il valore del volume stimato con il metodo monte carlo*/
		printf("Vbb %lf\ntotal_hits %d\npoints %d\nESTIMATED VOLUME %lf from task %d\n",Vbb,total_hits,points,Vbb*(1.0*total_hits/points), taskid);
	}
	MPI_Finalize();
	free(hits_master);
	return 0;
}
