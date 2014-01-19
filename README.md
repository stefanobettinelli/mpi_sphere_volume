% **Relazione di Algoritmi Avanzati**
% Stefano Bettinelli
% Gennaio 2014

Abstract
--------

In questa breve relazione viene spiegata l'implementazione che ho fornito al problema della stima del volume di *n* sfere attraverso l'uso di **open-mpi**.

----------


Descrizione ad alto livello
---------------------------

Vengono ora elencati in modo ordinato le fasi attraversate dall'algoritmo per il calcolo del volume delle sfere:

- Inizializzazione mpi con un determinato numero di tasks (ogni task viene eseguito logicamente da un singolo processore)
- Ogni processore (task) effettua la lettura da file dei dati delle sfere (coordiante centro e raggio)
- Calcolo del volume del **Bounding Box**
- Il task 0 genera una matrice $n*3$ di n punti casuali nello spazio a 3 dimensioni
- Le righe della matrice generata nel punto precedente vengono divise per il numero di processori, ottenendo $n/p$ sottomatrici
- Ad ogni processore (compreso lo zeresimo) viene sottomessa una sottomatrice di dimensione $n/p$, per ogni punto contenuto nella sottomatrice il processore controlla se esso ricade in almeno una sfera all'interno del Bounding Box
- Ogni processore produce un intero che indica quanti punti di quelli ricevuti in input ricadono all'interno di almeno una sfera
- Vengono aggregati gli output di ogni processore e nel task 0 si calcola il valore della stima del volume: $EstVol = Vbb * c/S$ in cui: *Vbb* è il volume del Bounding Box, *S* è il numero totale dei punti generati e *c* sono i punti che ricadono all'interno delle sfere

----------

Dettagli implementativi
-----------------------

Vediamo ora quali sono gli aspetti implementativi più rilevanti, attraverso l'analisi del codice prodotto.
Successivamente all'inizializzazione **mpi** tutti processori leggono i dati dal file di sfere passato come paramentro e costruiscono una matrice di sfere $NSfere*4$, ogni riga della matrice contiene le coordinate del centro e la lunghezza del raggio.  

```
fscanf(sfere_file,"%d", &sfere_n);
double sfere[sfere_n][4];
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
```

Nel ciclo `while si calcolano anche le coordinate degli estremi del *Bounding Box* (*BB*)

Il volume del *BB* è calcolato in questo modo:

```
Vbb = (x_max-x_min)*(y_max-y_min)*(z_max-z_min);
```
Ora il task 0 detto anche *MASTER* genera la matrice di punti random all'interno del *BB*.


```
if(taskid == MASTER)
	{
		srand((unsigned int)time(NULL));
		for (i=0; i<points; i++) {
			r_points[i][0] = rand_point(x_min,x_max);
			r_points[i][1] = rand_point(y_min,y_max);
			r_points[i][2] = rand_point(z_min,z_max);
		}	
	}
```

I punti vengono generati negli intervalli delle 3 diminsioni che corrispondono alla lunghezza degli spigoli del *BB*.

Quindi si divide in porzioni uguali la matrice di punti random e ad ogni task vengono sottomesse le porzioni ottenute.

```
int rows = (int)(points / numtasks);
	double task_points[rows][3];
	MPI_Scatter(r_points,rows*3,MPI_DOUBLE
				,task_points,rows*3,
				MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
	for (i = 0; i<rows; i++) {
		for (j=0; j<sfere_n; j++) {
			if( hit(task_points[i][0], task_points[i][1], task_points[i][2], 
				sfere[j][0], sfere[j][1], sfere[j][2], sfere[j][3]) == 1 ){
				hits++;
				break;
			}
		}
	}
```

I due *for* annidati consentono: di verificare per ogni punto della sottomatrice se esso ricade in almeno una delle sfere nel *BB*, se si la funzione *hit* restituisce 1 altrimenti 0.

Infine l'ultima fase dell'algoritmo si occupa di raccogliere il numero di 'hits' provenienti da ogni task, sommarli e usare il risultato per stimare il volume.

```
MPI_Gather(&hits, 1, MPI_INT, hits_master, 1, 
			MPI_INT, MASTER, MPI_COMM_WORLD);
int total_hits = 0;
if (taskid == MASTER) {
	double estVol = 0;
	for (i=0; i<numtasks; i++){
		total_hits += hits_master[i];
	}
	estVol = Vbb*(1.0*total_hits/points);
}
MPI_Finalize();
```

------------

Analisi prestazioni
-------------------

La strategia di parallelizzazione, seppur semplice permette comunque di avere un incremento di prestazioni rispetto alla versione seriale, con esecuzione su singolo processore.
In particolare quando l'algoritmo si trova nella penultima fase, deve effettuare un loop annidato verificare quale dei punti random generati ricade in almeno una sfera, e nel caso pessimo con singolo processore si ottiene una complessità $O(n*m)$ dove *n* è il numero di punti e *m* è il numero di sfere.
Mentre in presenza di *p* processori, dividendo la matrice di punti in *n* porzioni ogni processore esegue il calcolo con una complessità pari a $O(n/p * m)$.
Si ottiene uno **Speed Up** pari al numero di processori:

$$
\frac{n*m}{n/p * m} = p
$$
