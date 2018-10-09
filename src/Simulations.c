/*
    'MolTi' and 'molti-console' detects communities from multiplex networks / 'bonf' computes q-values of annotations enrichment of communities / 'test' simulates random multiplex to test community detection approaches
    Copyright (C) 2015  Gilles DIDIER

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "Utils.h"
#include "Graph.h"
#include "ERMG.h"
#include "RandomGraph.h"
#include "Partition.h"
//#include "PartitionConsensus.h"
#include "Louvain.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define NAME_OUTPUT "out"
#define EXT_OUTPUT "csv"
#define EXT_OUTPUT_TREE "_tree.newick"
#define EXT_OUTPUT_FOSSIL "_fossil.csv"
#define INTMIN 0.6
#define INTMAX 0.7
#define EXTMIN 0.4
#define EXTMAX 0.5
#define NTABLE 3
#define NCLASS 2
#define NITER 4
#define NTRIALS 2
#define NP_MAX 20


#define INFTY 1E99
#define TRIAL 10
#define MAX_ITER 1000


//./testthbis -i 500 -t 10 -c 20 -p 0.5 0.2 -g 1000
//./testthbis -i 500 -t 10 -c 20 -p 0.1 0.01 -g 1000

#define HELPMESSAGE "\nsimul first simulates random multiplex networks with g vertices and from 1 to t layers with a balanced community structure of c communities (i simulations for each number of layers). It next detects communities by using aggregation and multiplex-modularity approaches and computes the adjusted Rand index between the reference community structure and the detected one. It writes a '.csv' file containing the means and the standard deviations of the adjusted Rand indexes for each method, with the following format:\n Column 1: number of layers\n Column 2 and 3: mean and standard deviation obtained with the multiplex-modularity approach\n Column 4 and 5: mean and standard deviation obtained with the sum-aggregation approach\n Column 6 and 7: mean and standard deviation obtained with the intersection\n Column 6 and 7: mean and standard deviation obtained with the union aggregation approach. \n\nusage: simul <options> <output file name>\n\nOptions:\n	-g <number>	set the number of vertices of the random graphs (50 as default)\n	-t <number>	set the max number of random graphs (10 as default)\n	-c <number>	set the number of classes (3 as default)\n	-i <number>	set the number of iterations (4 as default)\n	-p <prob intra> <prob inter>	add a new pair of probas\n	-a  <number>	set the modularity parameter (1 as default)\n	-m <number>	set the number of simultaneous threads\n	-h	display help message\n"



typedef struct {
	int niter, iter;
	TypePartition part1;
	TypeLouvainInfo info;
	double *randIndX, meanX, stdX;
	double *randInd, mean, std;
	double *randIndS, meanS, stdS;
	double *randIndI, meanI, stdI;
	double *randIndU, meanU, stdU;
	double *randIndC, meanC, stdC;
	TypeIterFunc iterFunct;
} TypeThreadArgument;


typedef struct
{
	int number, maxThreads;
	TypeThreadArgument arg;
	pthread_mutex_t mutex_rand;
	pthread_mutex_t mutex_number;
	pthread_mutex_t mutex_arg;
	pthread_cond_t cond_number;
} TypeThreadData;

 
TypeThreadData set =
{
	.number = 0,
	.maxThreads = 40,
	.mutex_rand = PTHREAD_MUTEX_INITIALIZER,
	.mutex_number = PTHREAD_MUTEX_INITIALIZER,
	.mutex_arg = PTHREAD_MUTEX_INITIALIZER,
	.cond_number = PTHREAD_COND_INITIALIZER
};


static void *simulThread(void *data) {
	double  cmp, cmpS, cmpI, cmpU, cmpC, cmpX, max;
	TypeMultiGraph *graph;
	TypePartition part2;
	int i;
	
	pthread_mutex_lock (&set.mutex_arg);
	pthread_mutex_unlock (&set.mutex_arg);

	graph = (TypeMultiGraph*) data;
	
	part2 = getPartition(graph, LouvainType, iterateGlobal, (void*) &(set.arg.info));
	pthread_mutex_lock (&set.mutex_arg);
		cmp = comparePart(correctedRandIndex, &set.arg.part1, &part2);
	pthread_mutex_unlock (&set.mutex_arg);
	free((void*) part2.atom);
	
	part2 = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
	pthread_mutex_lock (&set.mutex_arg);
		cmpS = comparePart(correctedRandIndex, &set.arg.part1, &part2);
	pthread_mutex_unlock (&set.mutex_arg);
	free((void*) part2.atom);

	part2 = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
	max = getModularityMulti(graph, &part2);
	for(i=1; i<2; i++) {
		TypePartition partTmp;
		double modu;
		partTmp = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
		modu = getModularityMulti(graph, &partTmp);
		if(modu>max) {
			max = modu;
			free((void*) part2.atom);
			part2 = partTmp;
		} else
			free((void*) partTmp.atom);
	}
	pthread_mutex_lock (&set.mutex_arg);
		cmpI = comparePart(correctedRandIndex, &set.arg.part1, &part2);;
	pthread_mutex_unlock (&set.mutex_arg);
	free((void*) part2.atom);
	
	part2 = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
	max = getModularityMulti(graph, &part2);
	for(i=1; i<3; i++) {
		TypePartition partTmp;
		double modu;
		partTmp = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
		modu = getModularityMulti(graph, &partTmp);
		if(modu>max) {
			max = modu;
			free((void*) part2.atom);
			part2 = partTmp;
		} else
			free((void*) partTmp.atom);
	}
	pthread_mutex_lock (&set.mutex_arg);
		cmpU = comparePart(correctedRandIndex, &set.arg.part1, &part2);
	pthread_mutex_unlock (&set.mutex_arg);
	free((void*) part2.atom);
	
	part2 = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
	max = getModularityMulti(graph, &part2);
	for(i=1; i<4; i++) {
		TypePartition partTmp;
		double modu;
		partTmp = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
		modu = getModularityMulti(graph, &partTmp);
		if(modu>max) {
			max = modu;
			free((void*) part2.atom);
			part2 = partTmp;
		} else
			free((void*) partTmp.atom);
	}
	pthread_mutex_lock (&set.mutex_arg);
		cmpC = comparePart(correctedRandIndex, &set.arg.part1, &part2);
	pthread_mutex_unlock (&set.mutex_arg);
	free((void*) part2.atom);
	
	part2 = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
	max = getModularityMulti(graph, &part2);
	for(i=1; i<5; i++) {
		TypePartition partTmp;
		double modu;
		partTmp = getPartition(graph, LouvainType, iterateGlobalRandomized, (void*) &(set.arg.info));
		modu = getModularityMulti(graph, &partTmp);
		if(modu>max) {
			max = modu;
			free((void*) part2.atom);
			part2 = partTmp;
		} else
			free((void*) partTmp.atom);
	}
	pthread_mutex_lock (&set.mutex_arg);
		cmpX = comparePart(correctedRandIndex, &set.arg.part1, &part2);
	pthread_mutex_unlock (&set.mutex_arg);
	free((void*) part2.atom);
	
	freeMultiGraph(graph);
	pthread_mutex_lock (&set.mutex_arg);
		if(set.arg.iter<set.arg.niter) {
			set.arg.randInd[set.arg.iter] = cmp;
			set.arg.mean += cmp;
			set.arg.randIndS[set.arg.iter] = cmpS;
			set.arg.meanS += cmpS;
			set.arg.randIndI[set.arg.iter] = cmpI;
			set.arg.meanI += cmpI;
			set.arg.randIndU[set.arg.iter] = cmpU;
			set.arg.meanU += cmpU;
			set.arg.randIndC[set.arg.iter] = cmpC;
			set.arg.meanC += cmpC;
			set.arg.randIndX[set.arg.iter] = cmpX;
			set.arg.meanX += cmpX;
			set.arg.iter++;
		}
	pthread_mutex_unlock (&set.mutex_arg);
	pthread_mutex_lock (&set.mutex_number);
	set.number--;
	pthread_cond_signal (&set.cond_number);
	pthread_mutex_unlock (&set.mutex_number);
	return NULL;
}


int main(int argc, char **argv) {		
	char outputFileName[STRING_SIZE], option[256], kind = 'S';
	FILE *fo;
	int i, j, nTable, maxTable = 10, size = 50, niter = NITER, sizeAtom = 3, np = 0;
	double *alpha,  pin[NP_MAX], pex[NP_MAX], alp = 1., thre = 0., pres = 1.;
	TypeIterFunc iterFunct = iterateGlobal;
//	srand(time(NULL));
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['g']) {
			option['g'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &size) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -g");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxTable) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &sizeAtom) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &thre) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -f");
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['p']) {
			option['p'] = 0;
			if(np >= NP_MAX) {
				fprintf(stderr, "No more than %d probas\n", NP_MAX-1);
				exit(1);
			}
			if((i+1)<argc && sscanf(argv[i+1], "%lf", pin+np) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", pex+np) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -p");
			np++;
		}
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &alp) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -a");
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &pres) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -x");
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &(set.maxThreads)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -m");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(np == 0) {
		pin[0] = 0.2;
		pex[0] = 0.1;
		np = 1;
	}

	if(!(i<argc && sscanf(argv[i++], "%s", outputFileName) == 1)) {
		sprintf(outputFileName, "M_%c_%s_%d_%d_%d_%d.%s", kind, NAME_OUTPUT, sizeAtom, size, maxTable, niter, EXT_OUTPUT);
	}
	alpha = (double*) malloc(sizeAtom*sizeof(double));
	for(i=0;i<sizeAtom; i++)
		alpha[i] = 1./((double)sizeAtom);

	if((fo = fopen(outputFileName, "w"))) {
		double std, stdS, stdI, stdU, stdC, stdX;
		set.arg.iterFunct = iterFunct;
		set.arg.randInd = (double*) malloc(niter*sizeof(double));
		set.arg.randIndS = (double*) malloc(niter*sizeof(double));
		set.arg.randIndI = (double*) malloc(niter*sizeof(double));
		set.arg.randIndU = (double*) malloc(niter*sizeof(double));
		set.arg.randIndC = (double*) malloc(niter*sizeof(double));
		set.arg.randIndX = (double*) malloc(niter*sizeof(double));
		setBalancedPartition(sizeAtom, size, &set.arg.part1);
		set.arg.info.weight = (double*) malloc(maxTable*sizeof(double));
		for(nTable=0; nTable<maxTable; nTable++)
			set.arg.info.weight[nTable] = 1.;
		set.arg.info.alpha = alp;
		for(nTable=1; nTable<maxTable; nTable++) {
			int cont = 1;
			i = 0;
			set.arg.niter = niter;
			set.arg.iter = 0;
			set.arg.mean = 0.;
			set.arg.meanS = 0.;
			set.arg.meanI = 0.;
			set.arg.meanU = 0.;
			set.arg.meanX = 0.;
			printf("simulating %d graph\n", nTable);
			while(cont) {
				pthread_mutex_lock(&set.mutex_number);
				while(i < niter && set.number < set.maxThreads) {
					TypeMultiERMG *model;
					pthread_t thread;	
					int ret = 0;
					model = getBasicMultiERMGX(sizeAtom, nTable, alpha, pin, pex, np);
//					if((ret = pthread_create(&thread, NULL, simulThread, (void*) getRandomMultiGraph(model, size, &set.arg.part1))) == 0) {
					if((ret = pthread_create(&thread, NULL, simulThread, (void*) getRandomMultiGraphMissing(model, size, &set.arg.part1, pres))) == 0) {
						int err;
						if((err = pthread_detach(thread)) == 0) {
							set.number++;
							i++;
						} else {
							fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror (err));
//							pthread_kill(thread, 0);
						}
					} else
						fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror (ret));
					freeMultiERMG(model);
fprintf (stderr, "i %d threads %d\n", i, set.number);
				}
				cont = (set.number > 0);
				if(cont)
					pthread_cond_wait (& set.cond_number, & set.mutex_number);
				pthread_mutex_unlock (& set.mutex_number);
			}
	pthread_mutex_lock (&set.mutex_arg);
			set.arg.mean /= (double) set.arg.iter;
			std = 0.;
			for(i=0;i<set.arg.iter; i++)
				std += pow(set.arg.randInd[i]-set.arg.mean, 2.);
			std /= (double) set.arg.iter -1.;
			std = sqrt(std);
			
			set.arg.meanS /= (double) set.arg.iter;
			stdS = 0.;
			for(i=0;i<set.arg.iter; i++)
				stdS += pow(set.arg.randIndS[i]-set.arg.meanS, 2.);
			stdS /= (double) set.arg.iter -1.;
			stdS = sqrt(stdS);
			
			set.arg.meanI /= (double) set.arg.iter;
			stdI = 0.;
			for(i=0;i<set.arg.iter; i++)
				stdI += pow(set.arg.randIndI[i]-set.arg.meanI, 2.);
			stdI /= (double) set.arg.iter -1.;
			stdI = sqrt(stdI);
			
			set.arg.meanU /= (double) set.arg.iter;
			stdU = 0.;
			for(i=0;i<set.arg.iter; i++)
				stdU += pow(set.arg.randIndU[i]-set.arg.meanU, 2.);
			stdU /= (double) set.arg.iter -1.;
			stdU = sqrt(stdU);

			set.arg.meanC /= (double) set.arg.iter;
			stdC = 0.;
			for(i=0;i<set.arg.iter; i++)
				stdC += pow(set.arg.randIndC[i]-set.arg.meanC, 2.);
			stdC /= (double) set.arg.iter -1.;
			stdC = sqrt(stdC);
			
			set.arg.meanX /= (double) set.arg.iter;
			stdX = 0.;
			for(i=0;i<set.arg.iter; i++)
				stdX += pow(set.arg.randIndX[i]-set.arg.meanX, 2.);
			stdX /= (double) set.arg.iter -1.;
			stdX = sqrt(stdX);

			fprintf(fo, "%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE", nTable, set.arg.mean, std, set.arg.meanS, stdS, set.arg.meanI, stdI, set.arg.meanU, stdU, set.arg.meanC, stdC, set.arg.meanX, stdX);
			fprintf(fo, "\n");
			fprintf(stdout, "%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", nTable, set.arg.mean, set.arg.meanS, set.arg.meanI, set.arg.meanU, set.arg.meanC, set.arg.meanX);
	pthread_mutex_unlock (&set.mutex_arg);
		}
		free((void*)set.arg.info.weight);
		free((void*) set.arg.part1.atom);
		free((void*)set.arg.randInd);
		free((void*)set.arg.randIndS);
		free((void*)set.arg.randIndI);
		free((void*)set.arg.randIndU);
		free((void*)set.arg.randIndC);
		fclose(fo);
	} else
		exitProg(ErrorWriting, "Can't read the file");
	exit(0);
	return 0;
}
