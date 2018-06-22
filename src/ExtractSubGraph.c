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


#define OUTPUT "output"

#define MIN_TIME 0.001

#define SIZE_BUFFER_CHAR 300
#define HELPMESSAGE "\n 'extract_graph' takes as input two files containing a graph and a list of vertices respectively and extract the corresponding subgraph.\n\nusage: extract_graph <graph file> <list of vertices file>\n"

static TypeGraph *newGraphSize(TypeSizeG sizeGraph);
static TypeLexiTree *readDict(FILE *f);


/*simulate a random graph following model*/
TypeGraph *newGraphSize(TypeSizeG sizeGraph) {
	TypeGraph *g;
	TypeSizeG n, m;
	g = (TypeGraph*) malloc(sizeof(TypeGraph));
	g->sizeGraph = sizeGraph;
	g->name = (char**) malloc(g->sizeGraph*sizeof(char*));
	g->edge = (TypeEdgeG**) malloc(g->sizeGraph*sizeof(TypeEdgeG*));
	for(n=0; n<g->sizeGraph; n++) {
		g->name[n] = NULL;
		g->edge[n] = (TypeEdgeG*) malloc(g->sizeGraph*sizeof(TypeEdgeG));
		g->edge[n][n] = 0.;
		for(m=0; m<n; m++) {
			g->edge[n][m] = 0.;
			g->edge[m][n] = g->edge[n][m];
		}
	}
	return g;
}


/*read partition in Line format*/
TypeLexiTree *readDict(FILE *f) {
    char c;
    TypeLexiTree *dict;
    dict = newLexiTree();

    for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
    while(c != EOF) {
        char tmp[SIZE_BUFFER_CHAR+1];
        int i;
		if(c == '\'' || c == '"') {
			c = fgetc(f);
			for(i=0; i<SIZE_BUFFER_CHAR && c != EOF && c != '\'' && c != '"'; i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
			if(c == '\'' || c == '"')
				c = fgetc(f);
			else
				exitProg(ErrorReading, "Missing closing \" or '...");
		} else {
			for(i=0; i<SIZE_BUFFER_CHAR && c !=EOF && !issepline(c); i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
		}
		if(i == SIZE_BUFFER_CHAR)
			exitProg(ErrorExec, "Name too much long...");
		if(i>0) {
			tmp[i++] = '\0';
			addWordLexiTree(tmp, dict);
		}
		for(; c != EOF && issepline(c); c = fgetc(f));
    }
    return dict;
}

int main(int argc, char **argv) {		
	char inputFileNameGraph[SIZE_BUFFER_CHAR], inputFileNameList[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], option[256];
	int i;
	FILE *fi1, *fi2, *fo;

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
	if (i>=argc || sscanf(argv[i++], "%s", inputFileNameList) != 1) exitProg(ErrorArgument, "input file tree");
	if (i>=argc || sscanf(argv[i++], "%s", outputFileName) != 1) {
		char *ctmp;
		strcpy(outputFileName, inputFileNameGraph);
		ctmp = strchr(outputFileName, '.');
		if(ctmp != NULL)
			ctmp[0] = '\0';
		sprintf(outputFileName, "%s_fixed.gra", outputFileName);
	}
	if((fi1 = fopen(inputFileNameGraph, "r")) && (fi2 = fopen(inputFileNameList, "r")) && (fo = fopen(outputFileName, "w"))) {
		TypeGraph *g, *r;
		TypeSizeG *index, size = 0, n, m;
		TypeLexiTree *dict;

		g = readGraph(fi1);
		fclose(fi1);
		dict = readDict(fi2);
		fclose(fi2);
		index = (TypeSizeG*) malloc(g->sizeGraph*sizeof(TypeSizeG));
		for(n=0; n<g->sizeGraph; n++)
			if(g->name != NULL  && g->name[n] != NULL && findWordLexi(g->name[n], dict)>=0)
				index[n] = size++;
			else
				index[n] = -1;
		r = newGraphSize(size);
		for(n=0; n<g->sizeGraph; n++)
			if(index[n] >= 0) {
				r->name[index[n]] = strdpl(g->name[n]);
				for(m=0; m<g->sizeGraph; m++)
					if(index[m] >= 0)
						r->edge[index[n]][index[m]] = g->edge[n][m];
			}
printf("write file %s\n",outputFileName);
		fprintWeightedGraph(fo, r);
		fclose(fo);
	} else
		exitProg(ErrorReading, "Error");
	return 0;
}

