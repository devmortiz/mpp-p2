#ifndef _MH
#define _MH
	
	typedef struct {
		int array_int[1000];
		double fitness;
	} Individuo;
	
	void cruzar(Individuo, Individuo, Individuo, Individuo, int, int);
	void mutar(Individuo, int, int, int);
	void fitness(const double *, Individuo, int, int);
#endif
