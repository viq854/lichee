LICHeE: Fast and scalable inference of multi-sample cancer lineages
============

### About
LICHeE is a combinatorial method designed to reconstruct multi-sample cell lineage trees and infer the subclonal composition of the given samples based on variant allele frequencies (VAFs) of deep-sequencing somatic single nucleotide variants (SSNVs). The program accepts as input a list of SNVs with specified per-sample VAFs and outputs the inferred cell lineage tree(s) and the sample subclone decomposition. It provides an optional simple GUI to allow users to interact with the trees dynamically.

At a high level, LICHeE's execution can be broken down into the following steps: (1) SSNV calling across input samples, (2) SSNV clustering using VAFs (each group of SSNVs present in the same set of samples is clustered separately), (3) construction of the evolutionary constraint network (where the nodes are the clusters obtained in step (2) and the edges represent valid pairwise ancestry relationships), (4) search for lineage trees embedded in the network that satisfy all the phylogenetic constraints, and (5) output visualization.

For more information about the algorithm please see the following publication:  
Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S.  
*Fast and scalable inference of multi-sample cancer lineages*. Genome Biology 2015, 16:91.

### Program Parameters

For best results users are advised to explore the parameters exposed by the method and customize them to their specific datasets. The default values for several parameters are set fairly conservatively, assuming noisy real data, and relaxing these thresholds (especially when testing on simulated data) can produce more granular results. For example, lowering ```-maxClusterDist```, which controls the collapsing of nearby clusters, can order additional SSNVs by keeping them in separate clusters; similarly, lowering ```-minClusterSize``` to 1 will keep single-SSNV clusters in the network. More information on parameter tuning is provided below.

##### COMMANDS

```-build``` lineage tree reconstruction

##### INPUT/OUTPUT AND DISPLAY OPTIONS

```-i <arg>``` Input file path (*required*)  
```-o <arg>``` Output file path where the results should be written (default: input file name with the suffix .trees.txt)  
```-cp``` Input data represents cell prevalence (CP) values (as opposed to default VAF values)  
```-sampleProfile``` Input file contains the SSNV sample presence-absence profile (this will disable the default SSNV calling step)  
```-n,--normal <arg>``` Normal sample column id in the list of samples, 0-based (e.g. 0 is the first column) (*required*\*)  
```-clustersFile <arg>``` SSNV clusters file path  
```-s,--save <arg>``` Maximum number of output trees to save, if any (default: 1)  
```-showNetwork,--net``` Display the constraint network  
```-showTree,--tree <arg>``` Display the top ranking lineage tree(s) (default: 0)  
```-color``` Enable lineage tree visualization in color mode  
```-dot``` Enable DOT file export of the top-scoring tree for Graphviz visualization (saved by default to: input file with suffix .dot)  
```-dotFile <arg>``` DOT file path  

##### SSNV FILTERING AND CALLING

```-maxVAFAbsent,--absent <arg>``` Maximum VAF to consider an SSNV as robustly absent from a sample (*required*\*)  
```-minVAFPresent,--present <arg>``` Minimum VAF to consider an SSNV as robustly present in a sample (*required*\*)  
```-maxVAFValid <arg>``` Maximum allowed VAF in a sample (default: 0.6)  
```-minProfileSupport <arg>``` Minimum number of *robust*\*\* SSNVs required for a group presence-absence profile to be labeled robust during SNV calling: SNVs from non-robust groups can be re-assigned to existing robust groups (default: 2)

\* *these parameters are required unless the -sampleProfile option is specified*  
\*\* *robust SNVs have VAFs < maxVAFAbsent or > minVAFPresent across all the samples*

##### PHYLOGENETIC NETWORK CONSTRUCTION AND TREE SEARCH

```-minClusterSize <arg>``` Minimum number of SSNVs required per cluster (default: 2)  
```-minPrivateClusterSize <arg>``` Minimum number of SSNVs required for a private cluster (i.e. with SSNVs occurring only in one sample) (default: 1)  
```-minRobustNodeSupport <arg>``` Minimum number of robust SSNVs required for a node to be labeled robust during tree search: non-robust nodes can be removed from the network when no valid lineage trees are found (default: 2)  
```-maxClusterDist <arg>``` Maximum mean VAF difference on average per sample up to which two SSNV clusters can be collapsed (default: 0.2)  
```-c,--completeNetwork``` Add all possible edges to the constraint network, by default private nodes are connected only to closest level parents and only nodes with no other parents are descendants of root  
```-e <arg>``` VAF error margin (default: 0.1)  
```-nTreeQPCheck <arg>``` Number of top-ranking trees on which the QP consistency check is run, we have not seen this check to fail in practice (default: 0, for best performance)

##### OTHER
```-v,--verbose``` Verbose mode, prints more information about each step of the algorithm  
```-h,--help``` Print program usage information

### How to Run

From the /release directory:

```
./lichee -build -i <input_file_path> [-minVAFPresent <VAF1> -maxVAFAbsent <VAF2> -n <normal_sample_id>] [other options]
```
### Examples

From the /release directory (for other command-line settings used on the ccRCC and HGSC datasets see the README file in the data/ directory):

\#Show the top ranking tree <br>
```
./lichee -build -i ../data/ccRCC/RK26.txt -maxVAFAbsent 0.005 -minVAFPresent 0.005 -n 0 -showTree 1
```

\#Eliminate private clusters/nodes that have fewer than 2 SSNVs, show and save to file the top-ranking tree <br>
```
./lichee -build -i ../data/ccRCC/RMH008.txt -maxVAFAbsent 0.005 -minVAFPresent 0.005 -n 0 -minPrivateClusterSize 2 -showTree 1 -s 1
```

\#Reduce the VAF cluster centroid distance, which determines when the clusters are collapsed <br>
```
./lichee -build -i ../data/hgsc/case6.txt -maxVAFAbsent 0.005 -minVAFPresent 0.01 -n 0 -maxClusterDist 0.1 -showTree 1
```

### Input File Types

LICHeE accepts three different file format types. The main file format is composed of a list of SSNVs with their associated VAF or CP values per sample: one SSNV entry per line.  
The file contains the following header with fields separated by tabs:
```
#chr position description <sample names separated by tabs> 
```

For example (the following file contains 5 samples and 3 SSNVs):
```
#chr    position    description    Normal    S1    S2    S3    S4                     
17      123456      A/T DUSP19     0.0       0.1   0.2   0.25  0.15                   
11      341567      C/G MUC16      0.0       0.4   0.09  0.38  0.24                   
9       787834      A/C OR2A14     0.0       0.35  0.14  0.17  0.48                   
```

Users can also optionally provide pre-computed SSNV calls per sample, by adding one more column to the above format before the sample frequency information, which can specify the binary presence-absence pattern of this SSNV across samples. For example, for a file with 5 samples, a pattern of 01001 implies that the SSNV was called in the second and fifth sample (column id 1 and 4, since we start counting at 0). In order to use this file type (and disable the default calling mechanism), users should include the ```-sampleProfile``` flag.  
An example is shown below:
```
#chr    position    description    profile        Normal    S1    S2    S3    S4      
1       184306474   A/G HMCN1      01111          0.0       0.1   0.2   0.25  0.15    
1       18534005    C/A IGSF21     01111          0.0       0.1   0.25  0.2   0.1     
1       110456920   G/A UBL4B      01111          0.0       0.4   0.4   0.45  0.45    
10      26503064    C/G MYO3A      01001          0.0       0.4   0.0   0.0   0.24    
```

Finally, users can also specify pre-computed SSNV clustering information by providing an additional input file 
containing the clusters (with the corresponding centroid VAFs per sample and the member SSNVs): one cluster per line.
The file should contain the following fields separated by tabs (corresponding to the primary SSNV input file):

```
profile   <cluster VAFs per sample separated by tabs> <comma-separated list of SSNVs> 
```

For example (the following file contains 3 clusters for the SSNV example file shown above; the SSNVs are specified as line numbers in the SSNV input file ignoring the header line, starting from 1):

```
01111     0.0   0.1   0.23  0.23  0.13    1,2                                         
01111     0.0   0.4   0.4   0.45  0.45    3                                           
01001     0.0   0.4   0.0   0.0   0.24    4                                           
```

### Output Visualization

The resulting trees and sample decomposition information produced by LICHeE can be written to a text file (using the ```-s``` option that specifies up to how many top trees should be saved; it is recommended to evaluate all the trees that achieved the best score) and visualized via the interactive LICHeE Lineage Tree Viewer GUI (using the ```-showTree``` option that specifies how many trees should be displayed). It is also possible to export the best-scoring tree as a DOT file for Graphviz visualization (using the ```-dot``` or ```-dotFile``` options).  

The GUI allows users to dynamically remove nodes from the tree, collapse clusters of the same SSNV group, and view information about each node (e.g. SSNV composition of cluster nodes or the subclone decomposition of sample nodes). The Snapshot button can be used anytime to capture the current state of the tree as a vector graphic PDF file (please note that it takes a bit of time to write out the image to file).

We currently support two display modes: plain (default) and color (enabled with the ```-color``` flag). In the color mode, each cluster node is assigned a unique color and each sample node is decorated with the colors corresponding to the clusters of mutations present in the sample. The sample is decomposed by color according to the (approximate) prevalence of each cluster in the sample. The contribution of a cluster to each sample is highlighted (in purple) when the cluster node is selected.

A few useful tips for working with the GUI: one or multiple nodes can be selected and dragged to the desired position, the size (zoom) and position of the graph can be adjusted using the trackpad.

Example 1. Visualization for ccRCC patient RK26

```
./lichee -build -i ../data/ccRCC/RK26.txt -maxVAFAbsent 0.005 -minVAFPresent 0.005 -n 0 -showTree 1 -color -dot
```

Display using Graphviz (Graphviz must be installed separately):

```
dot -Tpdf ../data/ccRCC/RK26.txt.dot -O
```

<p align="center">
<img src="https://github.com/viq854/lichee/blob/master/img_demo/RK26.txt.dot.png" width="65%" height="65%" />
</p>

Example 2. Visualization for ccRCC patient RMH008

```
./lichee -build -i ../data/ccRCC/RMH008.txt -maxVAFAbsent 0.005 -minVAFPresent 0.005 -n 0 -minPrivateClusterSize 2 -showTree 1 -color -dot
```

```
dot -Tpdf ../data/ccRCC/RMH008.txt.dot -O
```

<p align="center">
<img src="https://github.com/viq854/lichee/blob/master/img_demo/RMH008.txt.color.dot.png" width="65%" height="65%" />
</p>

Plain mode simple look (withot ```-color``` flag):

<p align="center">
<img src="https://github.com/viq854/lichee/blob/master/img_demo/RMH008.txt.dot.png" width="65%" height="65%" />
</p>

GUI interaction examples:

Cluster node 10 is selected, sample constributions highlighted in purple.

<p align="center">
<img src="https://github.com/viq854/lichee/blob/master/img_demo/lichee_cluster_demo.png" width="65%" height="65%" />
</p>


Sample node R5 is selected, lineages highlighted in purple:

<p align="center">
<img src="https://github.com/viq854/lichee/blob/master/img_demo/lichee_sample_demo.png" width="65%" height="65%" />
</p>

### Parameter Tuning and Diagnostics

In some cases, LICHeE may not find a valid tumor lineage tree for an input dataset given a specific parameter setting. In some other cases, multiple alternative lineage trees might be valid under different parameter settings. Therefore, it is recommended to explore various parameters when analyzing a particular dataset. 

For instance, since LICHeE uses a heuristic method to call SSNVs that heavily relies on the values of the ```-maxVAFAbsent``` and ```-minVAFPresent``` parameters, adjusting these parameters to reflect the expected noise levels in the data, or supplying pre-computed calls can be very useful. Furthermore, it might be useful to adjust the criteria for incorporating clusters into the constraint network. For example, clusters that contain only a few SSNVs are more likely to represent mis-called presence patterns and can be filtered out by increasing the ```-minClusterSize``` and ```-minPrivateClusterSize``` parameters. The parameter ```-minRobustNodeSupport``` (which determines how many robustly-called SSNVs are required for a node to be non-removable) can be increased to iteratively remove nodes from the network while no valid trees are found automatically. For very noisy data, the ```-e``` parameter can be increased to relax the VAF constraint enforcement (although this should be done sparingly). On the other hand, adjusting these parameters in the opposite direction can result in more granular trees and is advisable on less noisy datasets in order to get the most informative results.

For diagnostics, LICHeE outputs a log detailing the execution of each step of the algorithm (with more information provided using the verbose, ```-v```, flag). This log can be very useful when trying to diagnose the performance of the program and view any of its intermediate results. In particular, it provides information about SSNV calling, clustering, the structure of the resulting constraint network, tree scoring, and any other operations controlled by various parameter settings. Using the log, the user can also trace why a particular SSNV was not included in the final output tree(s) (e.g. due to filtering based on the maximum VAF allowed or due to cluster size constraints) by searching the log for the "Filtered" keyword or for the unique descriptor of the SSNV (the log will output the SSNV entry line for each filtered SSNV as it appears in the input file). Furthermore, when no valid trees are found, examining the cluster centroid VAFs and the topology of the constraint network can be helpful to determine why at least one phylogenetic constraint is violated in each candidate embedded spanning trees. 


### System Requirements

Java Runtime Environment (JRE) 1.6 

### License

MIT License 

### Support

For help running the program or any questions/suggestions/bug reports, please contact viq@cs.stanford.edu
