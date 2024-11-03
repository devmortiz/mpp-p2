#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#include "../include/mh.h"

#define MUTATION_RATE 0.15
#define PRINT 1
#define MAX_NODES 2
#define NGM_DEFAULT 1
#define NEM_DEFAULT 1
#define PROCESO_MAESTRO 0

void crear_tipo_datos(int m, MPI_Datatype *individuo_type)
{
	int blocklen[2] = {m, 1};
	MPI_Datatype dtype[2] = { MPI_INT, MPI_DOUBLE };
	
	MPI_Aint disp[2];
	disp[0] = offsetof(Individuo, array_int);
	disp[1] = offsetof(Individuo, fitness);
	
	MPI_Type_create_struct(2, blocklen, disp, dtype, individuo_type); 
	MPI_Type_commit(individuo_type);
}

void crear_tipo_datos2(int m, MPI_Datatype *individuo_type)
{
	int blocklen[1] = {1};
	MPI_Datatype dtype[1] = { MPI_DOUBLE };

	
	MPI_Aint disp[1];
	//Individuo2 ind;
    //MPI_Aint base_address;
    //MPI_Get_address(&ind, &base_address);
    //MPI_Get_address(&ind.fitness, &disp[0]);
	//disp[0] = MPI_Aint_diff(disp[0], base_address);
	disp[0] = offsetof(Individuo2, fitness);
	
	MPI_Type_create_struct(1, blocklen, disp, dtype, individuo_type); 
	MPI_Type_commit(individuo_type);
}

int aleatorio(int n) {
	return (rand() % n);  // genera un numero aleatorio entre 0 y n-1
}

int find_element(int *array, int end, int element)
{
	int i=0;
	int found=0;
	
	// comprueba que un elemento no está incluido en el individuo (el cual no admite enteros repetidos)
	while((i < end) && ! found) {
		if(array[i] == element) {
			found = 1;
		}
		i++;
	}
	return found;
}

void crear_individuo(int n, int m, Individuo individuo)
{
	int i=0, value;
	
	// inicializa array de elementos
	memset(individuo.array_int, -1, m * sizeof(int));
	
	while(i < m) {
		value = aleatorio(n);
		// si el nuevo elemento no está en el array...
		if(!find_element(individuo.array_int, i, value)) {
			individuo.array_int[i] = value;  // lo incluimos
			i++;
		}
	}
}

int comp_array_int(const void *a, const void *b) {
	return (*(int *)a - *(int *)b);
}

int comp_fitness(const void *a, const void *b) {
	/* qsort pasa un puntero al elemento que está ordenando */
	return (*(Individuo *)b).fitness - (*(Individuo *)a).fitness;
}

void evolucionar(int g, const double *d, int n, int m, int n_gen, int tam_pob, Individuo *poblacion) {
	int i, mutation_start;
	// los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
	for(i = 0; i < (tam_pob/2) - 1; i += 2) {
		cruzar(poblacion[i], poblacion[i+1], poblacion[tam_pob/2 + i], poblacion[tam_pob/2 + i + 1], n, m);
	}
	
	// inicia la mutacion a partir de 1/4 de la poblacion
	mutation_start = tam_pob/4;
	
	// muta 3/4 partes de la poblacion
	for(i = mutation_start; i < tam_pob; i++) {
		mutar(poblacion[i], n, m, MUTATION_RATE);
	}
	
	// recalcula el fitness del individuo
	for(i = 0; i < tam_pob; i++) {
		fitness(d, poblacion[i], n, m);
	}
	
	// ordena individuos segun la funcion de bondad (mayor "fitness" --> mas aptos)
	qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);


	if (PRINT) {
		printf("Generacion %d - ", g);
		printf("Fitness = %.0lf - ", (poblacion[0].fitness));
	}
}

double aplicar_mh(const double *d, int n, int m, int n_gen, int tam_pob, int *sol, int argc, char** argv)
{
  	// para reiniciar la secuencia pseudoaleatoria en cada ejecución
  	srand(time(NULL) + getpid());
	int g, i;
	
	// crea poblacion inicial (array de individuos)
	Individuo *poblacion = (Individuo *) malloc(tam_pob * sizeof(Individuo));
	assert(poblacion);
	
	// crea cada individuo (array de enteros aleatorios)
	for(i = 0; i < tam_pob; i++) {
		/*
		poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
		memset(poblacion[i], 0, sizeof(Individuo));
    	poblacion[i]->array_int = crear_individuo(n, m);
		
		// calcula el fitness del individuo
		fitness(d, poblacion[i], n, m);
		*/
		Individuo individuo;
		crear_individuo(n, m, individuo);
		fitness(d, individuo, n, m);
		poblacion[i] = individuo;
	}
	// printf("array[0]: %d\n", poblacion[0]);
	// printf("array[1]: %d\n", poblacion[1]);
	
	// ordena individuos segun la funcion de bondad (mayor "fitness" --> mas aptos)
	qsort(&poblacion, tam_pob, sizeof(Individuo *), comp_fitness);

	int np, rank;
	int ngm = NGM_DEFAULT;
	int nem = NEM_DEFAULT;

	MPI_Status status;
	//MPI_Request request, *p_requests;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//p_requests = (MPI_Request *) malloc(np*sizeof(MPI_Request));

	// dividir los individuos entre los procesos
	Individuo *subpoblacion = (Individuo *) malloc(tam_pob * sizeof(Individuo));
	assert(subpoblacion);
	// crea cada individuo (array de enteros aleatorios)
	for(i = 0; i < tam_pob; i++) {
		Individuo individuo;
		crear_individuo(n, m, individuo);
		fitness(d, individuo, n, m);
		subpoblacion[i] = individuo;
	}

	/*
	int *distribucion = (int *) malloc(np * sizeof(int));
	int *desplazamientos = (int *) malloc(np * sizeof(int));
	int count = 0;
	double aux = tam_pob*1.0/np*1.0;
	for (int i = 0; i < np; i++) {
		distribucion[i] = (aux*(i+1)-count);
		desplazamientos[i] = count;
		count = aux*(i+1);
	}
	*/

	Individuo *nuevos_individuos = (Individuo *) malloc(nem * sizeof(Individuo));
	assert(nuevos_individuos);

	Individuo2 individuo_prueba;
	Individuo2 nuevo_individuo_prueba;

	individuo_prueba.fitness = 2.0;

	MPI_Datatype MPI_individuo_type;
	MPI_Datatype MPI_individuo2_type;
	crear_tipo_datos(m, &MPI_individuo_type);
	crear_tipo_datos2(m, &MPI_individuo2_type);
	//MPI_Scatterv(poblacion, distribucion, desplazamientos, MPI_individuo_type, subpoblacion, distribucion, MPI_individuo_type, 0, MPI_COMM_WORLD);
	if(rank == 0) {
		MPI_Send(&individuo_prueba, 1, MPI_individuo2_type, 1, 0, MPI_COMM_WORLD);
		//MPI_Send(poblacion, nem, MPI_individuo_type, 1, 0, MPI_COMM_WORLD);
		printf("El proceso 1 ha enviado este individuo:\n"); // imprime bien
		printf("%f\n", individuo_prueba.fitness);
		
	} else {
		MPI_Recv(&nuevo_individuo_prueba, 1, MPI_individuo2_type, PROCESO_MAESTRO, 0, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_individuo2_type, &count);
		printf("El proceso 2 ha recibido %d datos, el primer individuo es:\n", count); // imprime nem
		printf("%f\n", nuevo_individuo_prueba.fitness); // imprime 0.000000
		/*
		MPI_Recv(nuevos_individuos, nem, MPI_individuo_type, PROCESO_MAESTRO, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_individuo_type, &count);
		printf("El proceso 2 ha recibido %d datos, el primer individuo es:\n", count); // imprime nem
		printf("%f\n", nuevos_individuos[0].fitness); // imprime 0.000000
		*/
	}
	
	printf("hola5\n");

  	// evoluciona la poblacion durante un numero de generaciones
	for(g = 0; g < n_gen; g++) {
    	evolucionar(g, d, n, m, n_gen, tam_pob, subpoblacion);
		if (g != 0 && g%ngm == 0) {
			if (rank > 0) {
				//MPI_Send(subpoblacion, nem, MPI_individuo_type, PROCESO_MAESTRO, MPI_ANY_TAG, MPI_COMM_WORLD);
				//Individuo* nuevos_individuos = (Individuo *) malloc(nem * sizeof(Individuo));
				//MPI_Recv(nuevos_individuos, nem, MPI_individuo_type, PROCESO_MAESTRO, 0, MPI_COMM_WORLD, NULL);
			} else {
				// nodo maestro
				// Individuo** nuevaPoblacion = (Individuo **) malloc(nem * np * sizeof(Individuo *));
				// añadir individuos propios
				for (int p = 1; p < np; p++) {
					// MPI_RECV(poblacion)
					// nuevaPoblacion.add(poblacion)
				}
				// mezclar(nuevaPoblacion)
				// subpoblaciones = dividir(nuevaPoblacion, np)
				for (int p = 1; p < np; p++) {
					// MPI_SEND(poblacion)
				}
			}
		}
  	}
	printf("hola6\n");
	// Recibir las poblaciones finales de los demás procesos
	if (rank > 0) {
		// MPI_SEND(poblacion)
	} else {
		for (int p = 1; p < np; p++) {
			// MPI_RECV(poblacion)
			// poblacion.add(poblacion)
		}
	}
	MPI_Finalize();
	printf("hola7\n");
	// ordena el array solucion
	qsort(poblacion[0].array_int, m, sizeof(int), comp_array_int);

    // y lo mueve a sol para escribirlo
	memmove(sol, poblacion[0].array_int, m*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (poblacion[0].fitness);
	
   	// Liberamos la memoria asociada a cada individuo
   	for (int i = 0; i < tam_pob; i++)
   	{
     	// Liberamos la memoria asociada al array dentro de cada individuo
     	//free(poblacion[i].array_int);
   	}
	// se libera la memoria reservada
	// Hay que liberar el array de enteros de cada individuo, luego cada individuo y luego la población
	free(poblacion);
	
	// devuelve el valor obtenido para el fitness
	printf("hola8\n");
	return value;
}

void factibilizar(Individuo individuo, int n, int m) {
	
	int * elementosPresentes = (int *) malloc(n * sizeof(int));
	memset(elementosPresentes, -1, n * sizeof(int));
	for(int i = 0; i < m; i++) {
		while (elementosPresentes[individuo.array_int[i]] == 1) {
			individuo.array_int[i] = aleatorio(n);
		}
		elementosPresentes[individuo.array_int[i]] = 1;
	}
	free(elementosPresentes);
}

	// CRUZAR: Para cada par de individuos, selecciona un punto aleatorio a partir
	// del que se produce la mezcla y los nuevos descendientes se crean intercambiando
	// los genes de los padres entre sí.
void cruzar(Individuo padre1, Individuo padre2, Individuo hijo1, Individuo hijo2, int n, int m)
{
	// Elegir un "punto" de corte aleatorio a partir del que se realiza el intercambio de los genes
	int pc = aleatorio(m-1); // El -1 es porque si sale n no se realizaría ningún corte
	
	// Los primeros genes del padre1 van al hijo1. Idem para el padre2 e hijo2.
	int i = 0;
	for (; i < pc; i++) {
		hijo1.array_int[i] = padre1.array_int[i];
		hijo2.array_int[i] = padre2.array_int[i];
	}
	
	// Y los restantes son del otro padre, respectivamente.
	for (; i < m; i++) { // Intercambio
		hijo1.array_int[i] = padre2.array_int[i];
		hijo2.array_int[i] = padre1.array_int[i];
	}
	
	// Factibilizar: eliminar posibles repetidos de ambos hijos
	// Si encuentro alguno repetido en el hijo1, lo cambio por otro que no este en el conjunto
	factibilizar(hijo1, n, m);
	factibilizar(hijo2, n, m);
}

void mutar(Individuo actual, int n, int m, int m_rate)
{
	// Decidir cuantos elementos mutar: Si el valor es demasiado pequeño la 
	// convergencia es muy pequeña y si es demasiado alto diverge
	
	// Cambia el valor de algunos elementos de array_int de forma aleatoria
	// teniendo en cuenta que no puede haber elementos repetidos:
	// una posibilidad podría ser usar una variable, m_rate, para establecer la intensidad de la mutación 
  	// (un bucle for con un número de iteraciones que dependa, por ejemplo, de m_rate*m)
	for(int i = 0; i < m; i++) {
		double rndm = aleatorio(100)/100;
  		if (rndm <= m_rate) {
	  		int elemento = aleatorio(n);
	  		while (find_element(actual.array_int, m, elemento))
		  		elemento = aleatorio(n);
		  	actual.array_int[i] = aleatorio(n);
    	}
  	}
}

int getIndiceMatrizCompacta(int i, int j, int n) {
	int k = 0;
	k = (n*n-n)/2;
	k = k-(((n-i)*(n-i))-(n-i))/2;
	k = k + j - i + 1;
	return k;
}

double distancia_ij(const double *d, int i, int j, int n)
{
	// Devuelve la distancia entre dos elementos i, j de la matriz 'd'
	if (i > j) {
		int aux = i;
		i = j;
		j = aux;
	}
	int k = getIndiceMatrizCompacta(i, j, n);
	//printf("%lf\n", d[k]);
	return d[k];
}


void fitness(const double *d, Individuo individuo, int n, int m)
{
	// Determina la calidad del individuo calculando la suma de la distancia entre cada par de enteros
  	double fitness = 0;
	for (int i = 0; i < m; i++) {
		for (int j = i; j < m; j++) {
			fitness += distancia_ij(d, individuo.array_int[i], individuo.array_int[j], n);
		}
	}
	individuo.fitness = fitness;
}

