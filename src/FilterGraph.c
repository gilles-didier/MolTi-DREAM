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
#include <math.h>
#include <time.h>
#include "Utils.h"
#include "Graph.h"
#include "Partition.h"


#define OUTPUT "output"

#define MIN_TIME 0.001

#define SIZE_BUFFER_CHAR 300
#define HELPMESSAGE "\n 'filter_graph' discards all the vertices with a signle edge.\n\nusage: filter_graph <input file> [output file]\n"

static TypeGraph *newGraphSize(TypeSizeG sizeGraph);

/*simulate a random graph following model*/
TypeGraph *filterGraph(TypeGraph *g) {
	TypeSizeG n, m, *present;
	TypeGraph *r;
	present = (TypeSizeG*) malloc(g->sizeGraph*sizeof(TypeSizeG));
	r = (TypeGraph*) malloc(sizeof(TypeGraph));
	r->sizeGraph = 0;
	for(n=0; n<g->sizeGraph; n++) {
		TypeSizeG d = 0;
		for(m=0; m<g->sizeGraph; m++)
			if(g->edge[n][m]>0)
				d++;
		if(d>1)
			present[n] = r->sizeGraph++;
		else
			present[n] = NO_VERTEX;
	}
	r->name = (char**) malloc(r->sizeGraph*sizeof(char*));
	r->edge = (TypeEdgeG**) malloc(r->sizeGraph*sizeof(TypeEdgeG*));
	for(n=0; n<g->sizeGraph; n++) {
		if(present[n] != NO_VERTEX) {
			r->name[present[n]] = strdpl(g->name[n]);
			r->edge[present[n]] = (TypeEdgeG*) malloc(r->sizeGraph*sizeof(TypeEdgeG));
			r->edge[present[n]][present[n]] = 0.;
			for(m=0; m<n; m++)
				if(present[m] != NO_VERTEX) {
					r->edge[present[n]][present[m]] = g->edge[n][m];
					r->edge[present[m]][present[n]] = g->edge[m][n];
				}
		}
	}
	free((void*)present);
	return r;
}

int main(int argc, char **argv) {		
	char inputFileNameGraph[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], option[256], *ctmp;
	int i;
	FILE *fi, *fo;
	TypeGraph *g, *h;
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}

	if (i>=argc || sscanf(argv[i++], "%s", inputFileNameGraph) != 1) exitProg(ErrorArgument, "input file tree");
	if((fi = fopen(inputFileNameGraph, "r"))) {
		g = readGraph(fi);
	} else
		exitProg(ErrorReading, "input file graph missing");
	strcpy(outputFileName, inputFileNameGraph);
	ctmp = strrchr(outputFileName, '.');
	if(ctmp != NULL)
		ctmp[0] = '\0';
	sprintf(outputFileName, "%s_filtered.txt", outputFileName);
	printf("%d vertices\n", g->sizeGraph);
	if((fo = fopen(outputFileName, "w"))) {
		h = filterGraph(g);
		while(h->sizeGraph < g->sizeGraph) {
			freeGraph(g);
			g = h;
	printf("%d vertices\n", g->sizeGraph);
			h = filterGraph(g);
		}
		printf("%d vertices\n", h->sizeGraph);
		fprintGraphSpecial(fo, h);
		fclose(fo);
	} else
		exitProg(ErrorReading, "Error");
	return 0;
}

