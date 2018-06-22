# list of executable files to produce
MOLTI = molti-console
EXTRACT=extract_graph
PARTITION = extract_part
STDTOCNS = stdTocNs
CNSTODREAM = cNsToDream
FILTER_SIZE = filter_size
FILTER_GRAPH = filter_graph

ALL = $(MOLTI) $(EXTRACT) $(PARTITION) $(FILTER_SIZE) $(FILTER_GRAPH) $(STDTOCNS) $(CNSTODREAM)
# .o files necessary to build the executables

OBJ_MOLTI = EdgesComposition.o Louvain.o Partition.o Graph.o  Utils.o Compute.o 
OBJ_EXTRACT = Partition.o Graph.o  Utils.o ExtractSubGraph.o 
OBJ_PARTITION = Partition.o  Utils.o ExtractPartition.o 
OBJ_STDTOCNS = Partition.o  Utils.o StandardToClustNSee.o 
OBJ_CNSTODREAM = Partition.o  Utils.o ClustNSeeToDream.o 
OBJ_FILTER_SIZE = Partition.o Utils.o FilterPartition.o 
OBJ_FILTER_GRAPH = Partition.o Graph.o Utils.o FilterGraph.o 

########### MODIFY ONLY TO CHANGE OPTIONS ############

# compiler and its options
CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE=200809L -std=c99 -Wall -pedantic -g -march=native



# linker and its options
LD = $(CC)
LDFLAGS = -g -L/usr/local/lib -L/usr/lib/ -lm -O3



############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(ALL)

# build the executable
$(MOLTI): $(OBJ_MOLTI)
	$(LD) $^ -o $@ $(LDFLAGS)
$(STDTOCNS): $(OBJ_STDTOCNS)
	$(LD) $^ -o $@ $(LDFLAGS)
$(CNSTODREAM): $(OBJ_CNSTODREAM)
	$(LD) $^ -o $@ $(LDFLAGS)
$(EXTRACT): $(OBJ_EXTRACT)
	$(LD) $^ -o $@ $(LDFLAGS)
$(FILTER_SIZE): $(OBJ_FILTER_SIZE)
	$(LD) $^ -o $@ $(LDFLAGS)
$(FILTER_GRAPH): $(OBJ_FILTER_GRAPH)
	$(LD) $^ -o $@ $(LDFLAGS)
$(PARTITION): $(OBJ_PARTITION)
	$(LD) $^ -o $@ $(LDFLAGS)

############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release

-include Makefile.d
