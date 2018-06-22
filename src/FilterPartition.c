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
#include "Partition.h"


#define OUTPUT "output"

#define MIN_TIME 0.001

#define SIZE_BUFFER_CHAR 300
#define HELPMESSAGE "\n 'filter_size' filters classes of a partition according to their sizes.\n\nusage: filter_size [options] <input file> [output file]\n\noptions are\n	-f <format>      	set the partition format of <input file>\n	-m <number>    		set the minimum size of classes to <number>\n	-M <number>    		set the maximum size of classes to <number>\n"

int main(int argc, char **argv) {		
	char inputFileNamePartition[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], outputPrefix[SIZE_BUFFER_CHAR], option[256], format = 'c';
	int i, max, min;
	FILE *fi, *fo;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "l or c is required after option -f");
		}
		if(option['M']) {
			option['M'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &max) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -m");
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &min) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -m");
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}

	if (i>=argc || sscanf(argv[i++], "%s", inputFileNamePartition) != 1) exitProg(ErrorArgument, "input file tree");
	if (i>=argc || sscanf(argv[i++], "%s", outputPrefix) != 1) {
		char *ctmp;
		strcpy(outputPrefix, inputFileNamePartition);
		ctmp = strrchr(outputPrefix, '.');
		ctmp[0] = '\0';
	}
	sprintf(outputFileName, "%s_filtered_size.cls", outputPrefix);
	if((fi = fopen(inputFileNamePartition, "r")) && (fo = fopen(outputFileName, "w"))) {
		TypePartition p;
		TypePartitionList *pl;
		char **name;
		int c, i;
		if(format == 'l')
			p = readPartitionLine(fi, &name);
		else
			p = readPartition(fi, &name);
		fclose(fi);
		pl = getPartitionList(&p);
		fprintf(fo, "#ClustnSee analysis export\n");
		for(c=0; c<pl->sizeAtom; c++) {
			int size = 0;
			for(i=pl->start[c]; i>=0; i=pl->next[i])
				size++;
			if(size<=max && size>=min) {
				fprintf(fo, "ClusterID:%d||\n", c+1);
				for(i=pl->start[c]; i>=0; i=pl->next[i])
					fprintf(fo, "%s\n", name[i]);
				fprintf(fo, "\n");
			}
		}
		fclose(fo);
		freePartitionList(pl);	
	} else
		exitProg(ErrorReading, "Error");
	return 0;
}

