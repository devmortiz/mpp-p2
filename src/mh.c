#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>

#include "../include/mh.h"

#define MUTATION_RATE 0.15
#define PRINT 1

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

int *crear_individuo(int n, int m)
{
	int i=0, value;
	int *individuo = (int *) malloc(m * sizeof(int));
	
	// inicializa array de elementos
	memset(individuo, -1, m * sizeof(int));
	
	while(i < m) {
		value = aleatorio(n);
		// si el nuevo elemento no está en el array...
		if(!find_element(individuo, i, value)) {
			individuo[i] = value;  // lo incluimos
			i++;
		}
	}
	return individuo;
}

int comp_array_int(const void *a, const void *b) {
	return (*(int *)a - *(int *)b);
}

int comp_fitness(const void *a, const void *b) {
	/* qsort pasa un puntero al elemento que está ordenando */
	return (*(Individuo **)b)->fitness - (*(Individuo **)a)->fitness;
}

double aplicar_mh(const double *d, int n, int m, int n_gen, int tam_pob, int *sol)
{

  

  	// para reiniciar la secuencia pseudoaleatoria en cada ejecución
  	srand(time(NULL) + getpid());
	int i, g, mutation_start;
	
	// crea poblacion inicial (array de individuos)
	Individuo **poblacion = (Individuo **) malloc(tam_pob * sizeof(Individuo *));
	assert(poblacion);
	
	// crea cada individuo (array de enteros aleatorios)
	for(i = 0; i < tam_pob; i++) {
		poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
		memset(poblacion[i], 0, sizeof(Individuo));
    	poblacion[i]->array_int = crear_individuo(n, m);
		
		// calcula el fitness del individuo
		fitness(d, poblacion[i], n, m);
	}
	
	// ordena individuos segun la funcion de bondad (mayor "fitness" --> mas aptos)
	qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);

  
	double best_fitness = 0;  // El fitness de las últimas 5 iteraciones.
  	int last_change = 0;
  	// evoluciona la poblacion durante un numero de generaciones
	for(g = 0; g < n_gen; g++)
	{
    
		
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
			printf("Fitness = %.0lf - ", (poblacion[0]->fitness));

    }
  }
	// ordena el array solucion
	qsort(poblacion[0]->array_int, m, sizeof(int), comp_array_int);

    // y lo mueve a sol para escribirlo
	memmove(sol, poblacion[0]->array_int, m*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (poblacion[0]->fitness);
	
   	// Liberamos la memoria asociada a cada individuo
   	for (int i = 0; i < tam_pob; i++)
   	{
     	// Liberamos la memoria asociada al array dentro de cada individuo
     	free(poblacion[i]->array_int);

     	free(poblacion[i]);
   	}
	// se libera la memoria reservada
	// Hay que liberar el array de enteros de cada individuo, luego cada individuo y luego la población
	free(poblacion);
	
	// devuelve el valor obtenido para el fitness
	return value;
}

void factibilizar(Individuo *individuo, int n, int m) {
	
	int * elementosPresentes = (int *) malloc(n * sizeof(int));
	memset(elementosPresentes, -1, n * sizeof(int));
	for(int i = 0; i < m; i++) {
		while (elementosPresentes[individuo->array_int[i]] == 1) {
			individuo->array_int[i] = aleatorio(n);
		}
		elementosPresentes[individuo->array_int[i]] = 1;
	}
	free(elementosPresentes);
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

void mutar(Individuo *actual, int n, int m, int m_rate)
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
	  		while (find_element(actual->array_int, m, elemento))
		  		elemento = aleatorio(n);
		  	actual->array_int[i] = aleatorio(n);
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
		for (int j = i; j < m; j++) {
			fitness += distancia_ij(d, individuo->array_int[i], individuo->array_int[j], n);
		}
	}
	individuo->fitness = fitness;
}
