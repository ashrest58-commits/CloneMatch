**CloneMatch**

**Overview**

A computational framework for identifying candidate genetic lineages from single-cell RNA sequencing (scRNA-seq) data by assessing concordance between clustering-based and phylogeny-based analyses of expressed single-nucleotide variants (SNVs). It addresses challenges in lineage inference arising from uncertain cluster assignments and poorly supported phylogenetic branching patterns and is demonstrated using acute myeloid leukemia scRNA-seq datasets.

CloneMatch requires user-provided cluster assignments and a cellular phylogeny and consists of three steps: (i) assessing the distribution of cells from each cluster on a given cellular phylogeny, (ii) merging overlapping clusters on the phylogeny, and (iii) filtering cells that are inconsistent with their assigned lineages.

#### ***Step 1: Evaluating cluster separation using phylogenetic distances***

CloneMatch tests whether each pair of clusters occupies a distinct clade on the cellular phylogeny, since cells from the same genetic lineage should be more closely related than cells from different lineages.

#### ***Step 2: Merging overlapping clusters***

Overlapping cluster pairs are connected in a graph (nodes \= clusters, edges \= overlapping pairs). Each connected component becomes one merged lineage..

#### ***Step 3: Filtering cells for phylogenetic consistency***

For each cell, CloneMatch compares its mean distance to the other cells of its own lineage with its mean distance to each other lineage (one-sided t-test).

***Outputs***

- pairwise inter-cluster distances, null-distribution summaries, (ΔD) values, overlap classifications, permutation Benjamini–Hochberg-adjusted *P*\-values  
- cell-level within-lineage and between-lineage distances, the original and merged cluster IDs, and the status of filtering in the final annotation  
- filtered phylogeny in Newick format and phylogenetic visualizations colored by original cluster, merged cluster and final lineage assignments.

**Processing data**

We provide a script to process the input files (.iqtree.tree files and \_cluster.tsv files). The pipeline matches the clusters to the tips (individual cells) of the tree, merges the cluster pairs that overlap based on the phylogenetic distances, filters out the cells that are farther than the cells of its own lineage as compared to the cells in other lineages, and generates the final tree with lineage annotation that is supported by both the approaches. The resulting sample output data is stored in the CloneMatch projects output directory (see directory tree structure below).

**Example input**

Input  
├── Initial  
│   └── test\_clusters[.tsv.gz](http://.tsv.gz)      	 (cluster assignment)  
├── test.iqtree.tree	 		 (cell phylogeny)

test\_clusters.tsv and test.iqtree.tree are outputs of clustering-based and cell-phylogeny based lineage annotation, which are used as inputs for CloneMatch. These files are generated from scRNAseq raw data using the scPhylogenomics pipeline, which uses SNPmanifold and IQ-TREE to generate the clusters and inferred tree, respectively ([https\://github.com/ewafula/scPhylogenomics](https://github.com/ewafula/scPhylogenomics)).

**Example usage:**

bash  
CloneMatch.R Tree=\<TREE FILE\> CloID=\<CLUSTER FILE\>

**Run one sample:**

bash  
CloneMatach.R Tree="Input/test.iqtree.tree" CloID="Input/Initial/test\_clusters[.tsv.gz](http://.tsv.gz)"

	

**Example output**

Output  
└── LE1  
    ├── test.iqtree.tree.png		          		(tree colored by original clusters)  
    ├── test.iqtree.treeConcordance1.txt       	(cluster-pair test)  
    ├── test.iqtree.treeMerge.png                  	(tree colored by merged clustering)  
    ├── test.iqtree.treeConcordance.txt         	(per-cell distances and filtering status)  
    ├── test.iqtree.treeFilter.png	           		(tree after filtering)  
    ├── test.iqtree.treeConcordance.nwk      		(tree after filtering, Newick)  
    ├── test.iqtree.treeCellCount.txt        		(cells per cluster before and after filtering)  
    └── test.iqtree.treeConcordanceCellAnno.txt 	(merged clone ID of each retained cell)

test.iqtree.tree.png: Visualization of an initial cell phylogeny-based tree. This gives an preliminary idea of whether the clones grouped together in a cluster-based result indeed match biological phylogeny or not. 

test.iqtree.treeConcordance1.txt: A text file with an information about the distance between two clone pairs, observed distance, null mean distance, difference between observed and null mean distance, and classification of the clone pairs (overlapped, weakly overlapped, and separated)

test.iqtree.treeMerge.png: Visualization of the merged tree based on the cluster-pair test. Each clone in the clusters are matched to the tips (cells) of the tree and the merged tree is generated in which the clones are usually separated into more than one lineages.

test.iqtree.treeConcordance.txt: A text file containing the information about distance of a cell from other cells within its own lineage vs other lineages, based on which the cells are kept or filtered out from the final tree. Only the cells that are closest to the cells in its own lineages are retained.

test.iqtree.treeFilter.png: The final filtered tree with only the cells that passed the per-cell distance filter.

test.iqtree.treeConcordance.nwk: Newick format of the final filtered tree

test.iqtree.treeCellCount.txt: A text file containing information about the original clone\_ID, Merged Clone\_ID, and original number of cells and retained number of cells in the final tree.

test.iqtree.treeConcordanceCellAnno.txt: A text file that contains information about the lineage annotation or clone\_ID for each cell.

## License
CloneMatch is distributed under the GNU GPL v3. For more information, see [license(LICENSE)
