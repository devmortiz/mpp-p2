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
#define NGM_DEFAULT 2
#define NEM_DEFAULT 2
#define PROCESO_MAESTRO 0
#define TAG 0

void crear_tipo_datos(int m, MPI_Datatype *individuo_type)
{
	int blocklen[2] = {1000, 1};
	MPI_Datatype dtype[2] = { MPI_INT, MPI_DOUBLE };
	
	MPI_Aint disp[2];
	disp[0] = offsetof(Individuo, array_int);
	disp[1] = offsetof(Individuo, fitness);
	
	MPI_Type_create_struct(2, blocklen, disp, dtype, individuo_type); 
	MPI_Type_commit(individuo_type);
}

void imprimirPoblacion(Individuo* poblacion, int n, int m) {
	for(int i = 0; i < n; i++) {
		printf("Individuo %d, array: ", i);
		for(int j = 0; j < m; j++) {
			printf("%d ", poblacion[i].array_int[j]);
		}
		printf("\nFitness: %.2f\n", poblacion[i].fitness);
	}
}

/*
void empaquetar(Individuo *individuo, int m, char *buffer) {
	int posicion = 0;
	MPI_Pack(individuo->array_int, m, MPI_INT, buffer, 4*m+8, &posicion, MPI_COMM_WORLD);
	MPI_Pack(&individuo->fitness, 1, MPI_DOUBLE, buffer, 4*m+8, &posicion, MPI_COMM_WORLD);
	
}

void desempaquetar(char* buffer, int m, Individuo *individuo) {
	int posicion = 0;
	MPI_Unpack(buffer, 4*m+8, &posicion, individuo->array_int, m, MPI_INT, MPI_COMM_WORLD);
	MPI_Unpack(buffer, 4*m+8, &posicion, &individuo->fitness, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
}
*/

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

void crear_individuo(int n, int m, Individuo *individuo)
{
	int i=0, value;
	
	// inicializa array de elementos
	memset(individuo->array_int, -1, m * sizeof(int));
	
	while(i < m) {
		value = aleatorio(n);
		// si el nuevo elemento no está en el array...
		if(!find_element(individuo->array_int, i, value)) {
			individuo->array_int[i] = value;  // lo incluimos
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

double aplicar_mh(const double *d, int n, int m, int n_gen, int tam_pob, int *sol, int argc, char** argv, int np, int rank, int MODO)
{
	int ngm = NGM_DEFAULT;
	int nem = NEM_DEFAULT;

	MPI_Status status;
	MPI_Request p_requests[np];

	// para reiniciar la secuencia pseudoaleatoria en cada ejecución
  	srand(time(NULL) + getpid());
	int g, i;
	
	// crea poblacion inicial (array de individuos)
	Individuo poblacion[tam_pob];
	
	// crea cada individuo (array de enteros aleatorios)
	for(i = 0; i < tam_pob; i++) {
		Individuo individuo;
		crear_individuo(n, m, &individuo);
		fitness(d, &individuo, n, m);
		poblacion[i] = individuo;
	}
	
	// ordena individuos segun la funcion de bondad (mayor "fitness" --> mas aptos)
	qsort(poblacion, tam_pob, sizeof(Individuo), comp_fitness);

	// dividir los individuos entre los procesos
	
	MPI_Datatype MPI_individuo_type;
	crear_tipo_datos(m, &MPI_individuo_type);
	
	int distribucion[np];
	int desplazamientos[np];
	int count = 0;
	int cont = 0; // para comprobar el correcto funcionamiento
	double aux = tam_pob*1.0/np*1.0;
	for (int i = 0; i < np; i++) {
		distribucion[i] = (aux*(i+1)-count);
		desplazamientos[i] = count;
		count = aux*(i+1);
		cont += distribucion[i];
	}
	// comprobación
	assert(cont == tam_pob);
	assert(desplazamientos[np-1] == count-distribucion[np-1]);

	int tam = distribucion[rank];
	Individuo subpoblacion[tam];
	
	if (MODO == 2) {
		MPI_Scatterv(&poblacion, distribucion, desplazamientos, MPI_individuo_type, &subpoblacion, tam, MPI_individuo_type, PROCESO_MAESTRO, MPI_COMM_WORLD);
	}
	else {
		if(rank == PROCESO_MAESTRO) {
			for (int i = 0; i < tam; i++) {
				subpoblacion[i] = poblacion[i]; // Guardar su parte
			}
			for (int p = 1; p < np; p++) {
				/* Paquetes
				char buffer[4*m*distribucion[p]+8];
				empaquetar(&poblacion[desplazamientos[p]], distribucion[p], m, buffer);
				MPI_Send(buffer, 4*m*distribucion[p]+8, MPI_PACKED, p, TAG, MPI_COMM_WORLD);
				*/
				
				if (MODO == 0)
					MPI_Send(&poblacion[desplazamientos[p]], distribucion[p], MPI_individuo_type, p, TAG, MPI_COMM_WORLD); // Síncrona
				else {
					MPI_Isend(&poblacion[desplazamientos[p]], distribucion[p], MPI_individuo_type, p, TAG, MPI_COMM_WORLD, &p_requests[rank]);
					MPI_Wait(&p_requests[rank], &status);
				}
			}

		} else {
			/* Paquetes
			char buffer[4*m*distribucion[rank]+8];
			MPI_Recv(buffer, 4*m*distribucion[rank]+8, MPI_PACKED, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD, &status);
			desempaquetar(buffer, m, subpoblacion);
			*/
			if (MODO == 0)
				MPI_Recv(&subpoblacion, tam, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD, &status); // Síncrona
			else {
				MPI_Irecv(&subpoblacion, tam, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD, &p_requests[rank]);
				MPI_Wait(&p_requests[rank], &status);
			}
		}
	}
	

  	// evoluciona la poblacion durante un numero de generaciones
	for(g = 0; g < n_gen;) {
		//evolucionar
		int mutation_start;
		// los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
		for(int i = 0; i < (tam/2) - 1; i += 2) {
			cruzar(&subpoblacion[i], &subpoblacion[i+1], &subpoblacion[tam/2 + i], &subpoblacion[tam/2 + i + 1], n, m);
		}
		
		mutation_start = tam/4; // inicia la mutacion a partir de 1/4 de la poblacion
		if (mutation_start == 0)
			mutation_start = 1; // Para el caso en que el proceso tenga 3 individuos o menos
		
		// muta 3/4 partes de la poblacion
		for(int i = mutation_start; i < tam; i++) {
			mutar(&subpoblacion[i], n, m, MUTATION_RATE);
		}
		
		// recalcula el fitness del individuo
		for(int i = 0; i < tam; i++) {
			fitness(d, &subpoblacion[i], n, m);
		}
		
		// ordena individuos segun la funcion de bondad (mayor "fitness" --> mas aptos)
		qsort(subpoblacion, tam, sizeof(Individuo), comp_fitness);

		if (PRINT && rank == PROCESO_MAESTRO) {
			printf("Generacion %d - ", g);
			printf("Fitness = %.2lf - \n", (subpoblacion[0].fitness));
		}

		g++;
		// Migrar individuos
		if (g%ngm == 0) {
			if (MODO == 2) {
				Individuo nuevaPoblacion[np*nem];
				MPI_Gatherv(&subpoblacion, nem, MPI_individuo_type, &nuevaPoblacion, distribucion, desplazamientos, MPI_individuo_type, PROCESO_MAESTRO, MPI_COMM_WORLD);
				qsort(nuevaPoblacion, np * nem, sizeof(Individuo), comp_fitness);

				if (rank == PROCESO_MAESTRO) {
					for (int i = 0; i < nem; i++)
						subpoblacion[i] = nuevaPoblacion[i];
				}
				MPI_Bcast(subpoblacion, nem, MPI_individuo_type, PROCESO_MAESTRO, MPI_COMM_WORLD);
			}
			else {
				if (rank != PROCESO_MAESTRO) {
					// nodos hijos
					if (MODO == 0) {
						MPI_Send(&subpoblacion, nem, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD); // Síncrona
					}
					else {
						MPI_Isend(&subpoblacion, nem, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD, &p_requests[rank]);
						MPI_Wait(&p_requests[rank], &status);
					}
					
					if (MODO == 0)
						MPI_Recv(&subpoblacion, nem, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD, &status); // Síncrona
					else {
						MPI_Irecv(&subpoblacion, nem, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD, &p_requests[rank]);
						MPI_Wait(&p_requests[rank], &status);
					}

				} else {
					// nodo maestro
					Individuo nuevaPoblacion[np*nem];

					// añadir individuos propios
					for (int i = 0; i < nem; i++) {
						nuevaPoblacion[i] = subpoblacion[i];
					}
					// añadir el resto de individuos
					for (int p = 1; p < np; p++) {
						if (MODO == 0)
							MPI_Recv(&nuevaPoblacion[p*nem], nem, MPI_individuo_type, p, TAG, MPI_COMM_WORLD, &status); // Síncrona
						else {
							MPI_Irecv(&nuevaPoblacion[p*nem], nem, MPI_individuo_type, p, TAG, MPI_COMM_WORLD, &p_requests[rank]);
							MPI_Wait(&p_requests[rank], &status);
						}
					}
					qsort(nuevaPoblacion, np * nem, sizeof(Individuo), comp_fitness);

					for (int i = 0; i < nem; i++) {
						subpoblacion[i] = nuevaPoblacion[i];
					}

					for (int p = 1; p < np; p++) {
						if (MODO == 0)
							MPI_Send(&nuevaPoblacion, nem, MPI_individuo_type, p, TAG, MPI_COMM_WORLD); // Síncrona
						else {
							MPI_Isend(&nuevaPoblacion, nem, MPI_individuo_type, p, TAG, MPI_COMM_WORLD, &p_requests[rank]);
							MPI_Wait(&p_requests[rank], &status);
						}
					}
				}
			}
		}
  	}
	// Recibir las poblaciones finales de los demás procesos
	if (MODO == 2) {
		MPI_Gatherv(&subpoblacion, tam, MPI_individuo_type, &poblacion, distribucion, desplazamientos, MPI_individuo_type, PROCESO_MAESTRO, MPI_COMM_WORLD);
		qsort(poblacion, tam_pob, sizeof(Individuo), comp_fitness);
	}
	else {
		if (rank != PROCESO_MAESTRO) {
			if (MODO == 0)
				MPI_Send(&subpoblacion, tam, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD); // Síncrona
			else {
				MPI_Isend(&subpoblacion, tam, MPI_individuo_type, PROCESO_MAESTRO, TAG, MPI_COMM_WORLD, &p_requests[rank]);
				MPI_Wait(&p_requests[rank], &status);
			}
		} else {
			// añadir individuos propios
			for (int i = 0; i < tam; i++) {
				poblacion[i] = subpoblacion[i];
			}
			// añadir el resto de individuos
			for (int p = 1; p < np; p++) {
				if (MODO == 0)
					MPI_Recv(&poblacion[desplazamientos[p]], distribucion[p], MPI_individuo_type, p, TAG, MPI_COMM_WORLD, &status); // Síncrona
				else {
					MPI_Irecv(&poblacion[desplazamientos[p]], distribucion[p], MPI_individuo_type, p, TAG, MPI_COMM_WORLD, &p_requests[rank]);
					MPI_Wait(&p_requests[rank], &status);
				}
			}
			
			qsort(poblacion, tam_pob, sizeof(Individuo), comp_fitness);
		}
	}

	MPI_Finalize();
	if (rank != PROCESO_MAESTRO)
		exit(EXIT_SUCCESS);
	// ordena el array solucion
	qsort(poblacion[0].array_int, m, sizeof(int), comp_array_int);

    // y lo mueve a sol para escribirlo
	memmove(sol, poblacion[0].array_int, m*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (poblacion[0].fitness);
	
	// devuelve el valor obtenido para el fitness
	return value;
}

void factibilizar(Individuo *individuo, int n, int m) {
	
	int elementosPresentes[n];
	memset(elementosPresentes, -1, n * sizeof(int));
	for(int i = 0; i < m; i++) {
		while (elementosPresentes[individuo->array_int[i]] == 1) {
			individuo->array_int[i] = aleatorio(n);
		}
		elementosPresentes[individuo->array_int[i]] = 1;
	}
}

	// CRUZAR: Para cada par de individuos, selecciona un punto aleatorio a partir
	// del que se produce la mezcla y los nuevos descendientes se crean intercambiando
	// los genes de los padres entre sí.
void cruzar(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n, int m)
{
	// Elegir un "punto" de corte aleatorio a partir del que se realiza el intercambio de los genes
	int pc = aleatorio(m-1); // El -1 es porque si sale n no se realizaría ningún corte
	
	// Los primeros genes del padre1 van al hijo1. Idem para el padre2 e hijo2.
	int i = 0;
	for (; i < pc; i++) {
		hijo1->array_int[i] = padre1->array_int[i];
		hijo2->array_int[i] = padre2->array_int[i];
	}
	
	// Y los restantes son del otro padre, respectivamente.
	for (; i < m; i++) { // Intercambio
		hijo1->array_int[i] = padre2->array_int[i];
		hijo2->array_int[i] = padre1->array_int[i];
	}
	
	// Factibilizar: eliminar posibles repetidos de ambos hijos
	// Si encuentro alguno repetido en el hijo1, lo cambio por otro que no este en el conjunto
	factibilizar(hijo1, n, m);
	factibilizar(hijo2, n, m);
}

void mutar(Individuo *actual, int n, int m, double m_rate)
{
	// Decidir cuantos elementos mutar: Si el valor es demasiado pequeño la 
	// convergencia es muy pequeña y si es demasiado alto diverge
	
	// Cambia el valor de algunos elementos de array_int de forma aleatoria
	// teniendo en cuenta que no puede haber elementos repetidos:
	// una posibilidad podría ser usar una variable, m_rate, para establecer la intensidad de la mutación 
  	// (un bucle for con un número de iteraciones que dependa, por ejemplo, de m_rate*m)
	for(int i = 0; i < m; i++) {
		double rndm = aleatorio(100)/100.0;
  		if (rndm <= m_rate) {
	  		int elemento = aleatorio(n);
	  		while (find_element(actual->array_int, m, elemento))
		  		elemento = aleatorio(n);
		  	actual->array_int[i] = elemento;
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


void fitness(const double *d, Individuo *individuo, int n, int m)
{
	// Determina la calidad del individuo calculando la suma de la distancia entre cada par de enteros
  	double fitness = 0;
	for (int i = 0; i < m; i++) {
		for (int j = i+1; j < m; j++) {
			assert(individuo->array_int[i] != individuo->array_int[j]); // Debug
			int res = distancia_ij(d, individuo->array_int[i], individuo->array_int[j], n);
			assert(res >= 0);
			fitness += res;
		}
	}
	individuo->fitness = fitness;
}
