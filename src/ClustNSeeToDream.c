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
#define HELPMESSAGE "\n 'cNsToDream' translates partition from Clust'N See to Dream formats.\n\nusage: cNsToDream <input file> [output file]\n"

int main(int argc, char **argv) {		
	char inputFileNamePartition[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], option[256];
	int i;
	FILE *fi, *fo;

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

	if (i>=argc || sscanf(argv[i++], "%s", inputFileNamePartition) != 1) exitProg(ErrorArgument, "input file tree");
	if (i>=argc || sscanf(argv[i++], "%s", outputFileName) != 1) {
		char *ctmp;
		strcpy(outputFileName, inputFileNamePartition);
		ctmp = strchr(outputFileName, '.');
		ctmp[0] = '\0';
		sprintf(outputFileName, "%s_dream.cls", outputFileName);
	}
	if((fi = fopen(inputFileNamePartition, "r")) && (fo = fopen(outputFileName, "w"))) {
		TypePartition p;
		char **name;
		p = readPartition(fi, &name);
		fclose(fi);
		fprintPartitionDream(fo, &p, name);
		fclose(fo);
	} else
		exitProg(ErrorReading, "Error");
	return 0;
}

