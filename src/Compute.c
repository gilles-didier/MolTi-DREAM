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




#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Utils.h"
#include "Graph.h"
#include "Partition.h"
#include "Louvain.h"
#include "EdgesComposition.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define NAME_OUTPUT "outFile"
#define EXT_OUTPUT "cls"
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
//./molti-console -g 5 ../../../../DREAM/DREAM11_subchallenge2_initial_networks/Raw/1_ppi_anonym_aligned_v2.txt -g 10 ../../../../DREAM/DREAM11_subchallenge2_initial_networks/Raw/2_ppi_anonym_aligned_v2.txt ../../../../DREAM/DREAM11_subchallenge2_initial_networks/Raw/3_signal_anonym_aligned_directed_v3.txt
#define HELPMESSAGE "\n 'molti-console' takes as input several graph files in the simple \"ncol\" format: one interaction per line, the 2 interactors being separated by a tab, an additional number separated by a tab is taken as the weight of the interaction (if not present the weight is 1). Communities are detected by maximizing their Multiplex-modularity. One can assign a weight to a graph by adding it with option \"-g\" (see below). By default graphs have weight 1. It returns a file containing the community partition. Option '-p' sets the resolution parameter and Option '-s' provides the community partitions of all the individual graphs and that of their sum and union graphs (for comparison purpose).\n\nusage: molti-console [options] <input file1> <input file2> ...\n\noptions are\n	-o <file name>      		set the output file name prefix\n	-p <real number>    		set the Newman modularity resolution parameter (set to 1 as default)\n	-g <weight> <input file> 	add the graph from <input file> with weight <weight>\n	-s                  		compute partition on each graph individually and on the sum graph\n	-h                  		display help\n\n\n"

void fillMultiOne(TypeGraph *g, TypeMultiGraph *m) {
	m->name = g->name;
	m->edge[0] = g->edge;
	m->sizeGraph = g->sizeGraph;
}


int main(int argc, char **argv) {		
	char **inputNameTable, **name, outputFileName[STRING_SIZE], effectifFileName[STRING_SIZE], 
	outputPrefix[STRING_SIZE], outputFileNameRand[STRING_SIZE], *tmp, option[256], format = 'c';
	FILE *fi, *fo;
	int i, j, t, sizeTable = 0, singleF = 0, nRandom = 5;
	TypeMultiGraph *graph;
	TypeGraph **tableGraph;
	TypePartition  part, *tablePart;
	double alpha = 1., *weight;
	TypeIterFunc iterFunct = iterateGlobal;
	TypeLouvainInfo info;
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	sprintf(outputFileName, "%s.%s", NAME_OUTPUT, EXT_OUTPUT);

	tableGraph = (TypeGraph**) malloc((argc+3)*sizeof(TypeGraph*));
	weight = (double*) malloc((argc+3)*sizeof(double));
	inputNameTable = (char**) malloc((argc+3)*sizeof(char*));
	sizeTable = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", outputFileName) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -f");
		}
		if(option['s']) {
			option['s'] = 0;
			singleF = 1;
		}
		if(option['p']) {
			option['p'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &alpha) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -p");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a real number is required after option -p");
		}
		if(option['r']) {
			option['r'] = 0;
			iterFunct = iterateGlobalRandomized;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nRandom) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -r");
		}
/*		if(option['m']) {
			option['m'] = 0;
			if(!(sscanf(argv[i+1], "%s", inputFileName) == 1))
				exitProg(ErrorArgument, "wrong file name");
			i++;
			inputNameTable[sizeTable] = (char*) malloc((strlen(inputFileName)+1)*sizeof(char));
			strcpy(inputNameTable[sizeTable], inputFileName);
			printf("Reading file %s\n", inputNameTable[sizeTable]);
			if((fi = fopen(inputNameTable[sizeTable], "r"))) {
				tableGraph[sizeTable++] = readMatrixGraph(fi);
				fclose(fi);
			} else
				exitProg(ErrorReading, inputNameTable[sizeTable]);
		}
*/		if(option['g']) {
			option['g'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(weight[sizeTable])) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -g");
			if((i+1)<argc) {
				inputNameTable[sizeTable] = strdpl(argv[i+1]);
				i++;
			} else
				exitProg(ErrorArgument, "a file name is required after option -g");
			sizeTable++;
		}

		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	for(j=i; j<argc; j++) {
		inputNameTable[sizeTable] = strdpl(argv[j]);
		weight[sizeTable] = 1.;
		sizeTable++;
	}
	for(j=0; j<sizeTable; j++) {
printf("Reading file %s\n", inputNameTable[j]);
		if((fi = fopen(inputNameTable[j], "r"))) {
			tableGraph[j] = readGraph(fi);
			fclose(fi);
/***************************************************/
/*ligne à commenter pour avoir des graphes pondérés*/
/***************************************************/
//			fixEdgeGraph(tableGraph[j]);
		} else
			exitProg(ErrorReading, inputNameTable[j]);
printf("%d nodes\n", tableGraph[j]->sizeGraph);
	}
	if(sizeTable <= 0)
		exitProg(ErrorArgument, "at least one graph is required.");
	strcpy(outputPrefix, outputFileName);
	if((tmp = strrchr(outputPrefix, '.')) != NULL)
		tmp[0] = '\0';
	info.alpha = alpha;
	info.weight = weight;
	graph = toMultiGraph(tableGraph, sizeTable);
	if(iterFunct == iterateGlobalRandomized) {
		int i;
		double max;
		part = getPartition(graph, LouvainType, iterFunct, &info);
		max = getModularityMultiWeighted(graph, &part, weight, alpha);
		for(i=1; i<nRandom; i++) {
			TypePartition partTmp;
			double modu;
			partTmp = getPartition(graph, LouvainType, iterFunct, &info);
			modu = getModularityMultiWeighted(graph, &partTmp, weight, alpha);
			if(modu>max) {
				max = modu;
				free((void*) part.atom);
				part = partTmp;
			} else
				free((void*) partTmp.atom);
		}
	} else
		part = getPartition(graph, LouvainType, iterFunct, &info);
	tableGraph[sizeTable] = sumMultiGraph(graph);
	inputNameTable[sizeTable] = (char*) malloc((strlen("Sum")+1)*sizeof(char));
	strcpy(inputNameTable[sizeTable], "Sum");
	sizeTable++;
	tableGraph[sizeTable] = unionMultiGraph(graph);
	inputNameTable[sizeTable] = (char*) malloc((strlen("Union")+1)*sizeof(char));
	strcpy(inputNameTable[sizeTable], "Union");
	sizeTable++;
	tableGraph[sizeTable] = interMultiGraph(graph);
	inputNameTable[sizeTable] = (char*) malloc((strlen("Intersection")+1)*sizeof(char));
	strcpy(inputNameTable[sizeTable], "Intersection");
	sizeTable++;
	name = (char**) malloc(sizeTable*sizeof(char*));
	for(t=0; t<sizeTable; t++) {
		char *tmp;
		if((tmp = strrchr(inputNameTable[t], '.')) != NULL)
			tmp[0] = '\0';
		if((tmp=strrchr(inputNameTable[t], '/')) == NULL)
			tmp = inputNameTable[t];
		else
			tmp++;
		name[t] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(name[t], tmp);
//		printf("name[%d]\t%s\n", t, name[t]);
	}
	tablePart = (TypePartition*) malloc(sizeTable*sizeof(TypePartition));
	if(singleF) {
		TypeMultiGraph *gtmp;
		gtmp = (TypeMultiGraph*) malloc(sizeof(TypeMultiGraph));
		gtmp->sizeTable = 1;
		gtmp->edge = (TypeEdgeG***) malloc(sizeof(TypeEdgeG**));
		gtmp->present = (int**) malloc(sizeof(int*));
		gtmp->present[0] = (int*) malloc(graph->sizeGraph*sizeof(int));
		for(i=0; i<graph->sizeGraph; i++)
			gtmp->present[0][i] = 1;
		for(t=0; t<sizeTable; t++) {
			char output[STRING_SIZE];
			fillMultiOne(tableGraph[t], gtmp);
			sprintf(output, "%s_%s.%s", outputPrefix, name[t], EXT_OUTPUT);
printf("Computing %s\n", output);
			tablePart[t] = getPartition(gtmp, LouvainType, iterFunct, &info);
			if((fo = fopen(output, "w"))) {
				fprintPartitionClustNSee(fo, &(tablePart[t]), gtmp->name);
				fclose(fo);
			} else
				exitProg(ErrorWriting, output);
		}
		free((void*)gtmp->present[0]);
		free((void*)gtmp->present);
		free((void*)gtmp->edge);
		free((void*)gtmp);
	}
	if((fo = fopen(outputFileName, "w"))) {
		switch(format) {
			case 's':
				fprintPartitionStandard(fo, &part, graph->name);
				break;
			case 'l':
				fprintPartitionList(fo, &part, graph->name);
				break;
			case 'd':
				fprintPartitionDream(fo, &part, graph->name);
				break;
			case 'c':
			default:
				fprintPartitionClustNSee(fo, &part, graph->name);
		}
		fclose(fo);
	} else
		exitProg(ErrorWriting, outputFileName);
	sprintf(effectifFileName, "%s_effectif.csv", outputPrefix);
	if((fo = fopen(effectifFileName, "w"))) {
		int **eff, c, t, *cs;
		eff = getEdgesNumbers(&part, graph);
		cs = getClassSize(&part);
		fprintf(fo, "\tSize");
		for(t=0; t<graph->sizeTable; t++)
			fprintf(fo, "\t%s", name[t]);
		fprintf(fo, "\n");
		for(c=0; c<part.sizeAtom; c++) {
			fprintf(fo, "ClusterID:%d", c+1);
			fprintf(fo, "\t%d", cs[c]);
			for(t=0; t<graph->sizeTable; t++) {
				fprintf(fo, "\t%d", eff[c][t]);
			}
			fprintf(fo, "\n");
		}
		fclose(fo);
		for(c=0; c<part.sizeAtom; c++)
			free((void*)eff[c]);
		free((void*)eff);
		free((void*)cs);
	} else
		exitProg(ErrorWriting, effectifFileName);
		
	if(singleF) {
		for(t=sizeTable-3; t<sizeTable; t++) {
			sprintf(effectifFileName, "%s_%s_effectif.csv", outputPrefix, name[t]);
			if((fo = fopen(effectifFileName, "w"))) {
				int **eff, c, t, *cs;
				eff = getEdgesNumbers(&(tablePart[sizeTable-1]), graph);
				cs = getClassSize(&(tablePart[sizeTable-1]));
				fprintf(fo, "\tSize");
				for(t=0; t<graph->sizeTable; t++)
					fprintf(fo, "\t%s", name[t]);
				fprintf(fo, "\n");
				for(c=0; c<tablePart[sizeTable-1].sizeAtom; c++) {
					fprintf(fo, "ClusterID:%d", c+1);
					fprintf(fo, "\t%d", cs[c]);
					for(t=0; t<graph->sizeTable; t++) {
						fprintf(fo, "\t%d", eff[c][t]);
					}
					fprintf(fo, "\n");
				}
				fclose(fo);
				for(c=0; c<tablePart[sizeTable-1].sizeAtom; c++)
					free((void*)eff[c]);
				free((void*)eff);
				free((void*)cs);
			} else
				exitProg(ErrorWriting, effectifFileName);
			sprintf(effectifFileName, "%s_%s_graph.csv", outputPrefix, name[t]);
			if((fo = fopen(effectifFileName, "w"))) {
				fprintGraph(fo, tableGraph[t]);
			} else
				exitProg(ErrorWriting, effectifFileName);
		}
	}
	if(singleF) {
		int tl, tc;
		sprintf(outputFileNameRand, "%s_Rand.csv", outputPrefix);
		if((fo = fopen(outputFileNameRand, "w"))) {
			fprintf(fo,"All\t%lf\n", comparePartDiff(correctedRandIndex, &(part), graph->name, &(part), graph->name));
			for(tl=0; tl<sizeTable; tl++) {
				fprintf(fo, "%s\t%lf", name[tl], comparePartDiff(correctedRandIndex, &(part), graph->name, &(tablePart[tl]), tableGraph[tl]->name));
				for(tc=0; tc<=tl; tc++)
					fprintf(fo, "\t%lf", comparePartDiff(correctedRandIndex, &(tablePart[tc]),  tableGraph[tc]->name, &(tablePart[tl]), tableGraph[tl]->name));
				fprintf(fo, "\n");
			}
			fprintf(fo, "\tAll");
			for(tl=0; tl<sizeTable; tl++) {
				fprintf(fo, "\t%s", name[tl]);
			}
			fprintf(fo, "\n");
			fclose(fo);
		} else
			exitProg(ErrorWriting, outputFileNameRand);
	}
	for(t=0; t<sizeTable; t++) {
		freeGraph(tableGraph[t]);
		free((void*) name[t]);
		free((void*) inputNameTable[t]);
	}
	free((void*) tableGraph);
	free((void*) name);
	free((void*) inputNameTable);
	exit(0);
	return 0;
}
