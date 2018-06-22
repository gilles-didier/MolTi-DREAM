scriptStd.sh allows to apply molti-console in a recursive way. It needs the partition file (*.cls) and the folder with the network(s). For each cluster in the partition file that contains more than x proteins (to be set in the scriptStd.sh), it will create one or more subnetworks, extracting interaction from the network(s) folder. Then, molti-console (with the parameters set in the scriptStd.sh) will be re-applied on this or these networks.

>scriptStd.sh partition.cls network_folder/

scriptWeight.sh is the same script but for multiplex networks with different weights for the different graphs. It takes as input the partition, and then the different graphs with their corresponding weights.

>scriptWeight.sh partition.cls w1 network1.txt w2 network2.txt w3 network3.txt
w1, w2, w3 being the weights associated to each graph