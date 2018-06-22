The MolTi software suite allows identifying communities from multiplex networks, and annotated the obtained clusters. MolTi performs both clustering and annotation enrichment tests from multiplex networks, and allow an easy exploration of the results.

--------------------

 'molti-console' takes as input several graph files in the simple "ncol" format: one interaction per line, the 2 interactors being separated by a tab, an additional number separated by a tab is taken as the weight of the interaction (if not present the weight is 1). Communities are detected by maximizing their Multiplex-modularity. One can assign a weight to a graph by adding it with option "-g" (see below). By default graphs have weight 1. It returns a file containing the community partition. Option '-p' sets the resolution parameter and Option '-s' provides the community partitions of all the individual graphs and that of their sum and union graphs (for comparison purpose).

usage: molti-console [options] <input file1> <input file2> ...

options are
	-o <file name>      		set the output file name prefix
	-p <real number>    		set the Newman modularity resolution parameter (set to 1 as default)
	-g <weight> <input file> 	add the graph from <input file> with weight <weight>
	-s                  		compute partition on each graph individually and on the sum graph
	-h                  		display help


--------------------

 'extract_part' takes as input a partition file and write the elements of each class greater than a given threshold in a different file (used in the recursion script).

usage: extract_part [options] <input file>

options are
	-f <format>      	set the partition format of <input file>
	-m <number>    		set threshold to <number>

--------------------

 'extract_graph' takes as input two files containing a graph and a list of vertices respectively and extract the corresponding subgraph.

usage: extract_graph <graph file> <list of vertices file>


--------------------

 'stdTocNs' translates partition from "flat" to Clust'N See formats.

usage: stdTocNs <input file> [output file]

--------------------

 'cNsToDream' translates partition from Clust'N See to Dream formats.

usage: cNsToDream <input file> [output file]

--------------------

 'filter_size' filters classes of a partition according to their sizes.

usage: filter_size [options] <input file> [output file]

options are
	-f <format>      	set the partition format of <input file>
	-m <number>    		set the minimum size of classes to <number>
	-M <number>    		set the maximum size of classes to <number>

--------------------

 'filter_graph' discards all the vertices with a signle edge.

usage: filter_graph <input file> [output file]

