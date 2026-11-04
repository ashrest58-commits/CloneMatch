**CloneMatch**

**Overview**

A computational framework for identifying candidate genetic lineages from single-cell RNA sequencing (scRNA-seq) data by assessing concordance between clustering-based and phylogeny-based analyses of expressed single-nucleotide variants (SNVs). It addresses challenges in lineage inference arising from uncertain cluster assignments and poorly supported phylogenetic branching patterns and is demonstrated using acute myeloid leukemia scRNA-seq datasets.

CloneMatch requires user-provided cluster assignments and a cellular phylogeny and consists of three steps: (i) assessing the distribution of cells from each cluster on a given cellular phylogeny, (ii) merging overlapping clusters on the phylogeny, and (iii) filtering cells that are inconsistent with their assigned lineages.

#### ***Step 1: Evaluating cluster separation using phylogenetic distances***

CloneMatch tests whether each pair of clusters occupies a distinct clade on the cellular phylogeny, since cells from the same genetic lineage should be more closely related than cells from different lineages.

For clusters A and B, the observed distance (*Dobs*,) is the mean topological distance (number of edges) between all cell pairs across the two clusters. Edges are used instead of branch lengths, which can be unreliable in phylogenies from sparse SNV data. Cluster labels are then permuted between the two clusters (preserving cluster sizes) to obtain the null distance (*Dperm*), and the effect size is:

 *d*(A, B) \= *Dobs*(A, B) – mean\[*Dperm*(A, B)\]

Cluster pairs are classified as overlapping (*d* \< 1), weakly overlapping (1 ≤ *d* ≤ 2), or separated (*d* \> 2\)

#### ***Step 2: Merging overlapping clusters***

Overlapping cluster pairs are connected in a graph (nodes \= clusters, edges \= overlapping pairs). Each connected component becomes one merged lineage. This allows indirect merging: clusters A and C are merged if both overlap with B, even if A and C overlap only weakly or not at all. If all clusters merge into one group, CloneMatch reports a single candidate lineage, meaning the phylogeny does not support subdivision under the specified criteria.

#### ***Step 3: Filtering cells for phylogenetic consistency***

For each cell, CloneMatch compares its mean distance to the other cells of its own lineage with its mean distance to each other lineage (one-sided t-test). A cell is retained only if it is significantly closer to its own lineage than to every other lineage (P \< 0.05); otherwise it is removed from the phylogeny. The remaining tree and cell assignments define the final lineages.

***Outputs***

- pairwise inter-cluster distances, null-distribution summaries, (ΔD) values, overlap classifications, permutation Benjamini–Hochberg-adjusted *P*\-values  
- cell-level within-lineage and between-lineage distances, the original and merged cluster IDs, and the status of filtering in the final annotation  
- filtered phylogeny in Newick format and phylogenetic visualizations colored by original cluster, merged cluster and final lineage assignments.

**Data Description**

**Study Data**

Acute Myeloid Leukemia (AML) dataset

Single-cell RNA sequencing datasets from the acute myeloid leukemia (AML) patients with tumor heterogeneity to assess how CloneMatch can be used to integrate two lineage annotation approaches \- clustering-based and phylogeny-based, each with certain limitations, to get a genetic lineages supported by both approaches.  

* LE1 sample from an individual with acute myeloid leukemia (AML)

**Processing data**

We provide a script to process the input files (.iqtree.tree files and \_cluster.tsv files) of sample LE1. The pipeline matches the clusters to the tips (individual cells) of the tree, merges the cluster pairs that overlap based on the phylogenetic distances, filters out the cells that are farther than the cells of its own lineage as compared to the cells in other lineages, and generates the final tree with lineage annotation that is supported by both the approaches. The resulting sample output data is stored in the CloneMatch projects output directory (see directory tree structure below).  
**Example usage:**

bash  
Rscript ConcorTest\_clean1.r Tree=\<TREE FILE\> CloID=\<CLUSTER FILE\>

**Example input**

Input  
├── Initial  
│   └── LE1\_clusters[.tsv.gz](http://.tsv.gz)      	  (cluster assignment)  
├── LE1.all-cell-types.iqtree.tree	                      (cell phylogeny)	

###write about input files and where they come form.also write that the input files can be generated using scPhylogeny pipeline (link it)

**Run one sample:**
bash  
Rscript ConcorTest\_clean1.r Tree="Input/LE1.all-cell-types.iqtree.tree" CloID="Input/Initial/LE1\_clusters[.tsv.gz](http://.tsv.gz)"
#change file name, it is to long

**Example output**
Output  
└── LE1  
    ├── LE1.all-cell-types.iqtree.tree.png		          (tree colored by original clusters)  
    ├── LE1.all-cell-types.iqtree.treeConcordance1.txt       (per-cell distances and filtering status)  
    ├── LE1.all-cell-types.iqtree.treeMerge.png                  (tree colored by merged clustering)  
    ├── LE1.all-cell-types.iqtree.treeConcordance.txt         (cluster-pair test)  
    ├── LE1.all-cell-types.iqtree.treeFilter.png	           (tree after filtering)  
    ├── LE1.all-cell-types.iqtree.treeConcordance.nwk       (tree after filtering, Newick)  
    ├── LE1.all-cell-types.iqtree.treeCellCount.txt        (cells per cluster before and after filtering)  
    └── LE1.all-cell-types.iqtree.treeConcordanceCellAnno.txt (merged clone ID of each retained cell) 
#change output name, this too long
#explain outcomes
