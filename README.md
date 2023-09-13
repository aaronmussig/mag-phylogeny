# MAG Phylogeny

This pipeline incorporates multiple programs to infer the novelty of 
**de-replicated** MAGs using the [GTDB](https://gtdb.ecogenomic.org/) taxonomy.

## 1. Installation

1a. Clone the repoisitory to a local directory, e.g. `/tmp/mag-phylogeny`

```bash
# Clone the repository to /tmp/mag-phylogeny
cd /tmp
git clone https://github.com/aaronmussig/mag-phylogeny.git
```

1b. Create a new conda environment and install the dependencies.

```bash
conda env create -f env.yaml
conda activate mag-phylogeny
```

1c. Some rules require the following conda environments to exist and 
be configured as per the respective documentation.

* `gtdbtk-2.3.0`
* `fasttree-2.1.11`
* `genometreetk-0.1.8`
* `phylorank-0.1.12`

These can be changed by modifying the rules.

## 2. Pipeline

The full pipeline can be accessed by running the magphylogeny module, e.g.:

```bash
conda activate mag-phylogeny
cd /tmp/mag-phylogeny
python -m magphylogeny --help
```

### 2a. Phylogeny

The GTDB-Tk _de novo_ pipeline has been run using individual steps to duplicating
the identify and align steps when running for different domains.

1. `gtdbtk identify`
2. `gtdbtk align`
    - Using the`--debug` flag to produce individual marker gene alignments.
3. `gtdbtk infer`
4. `gtdbtk root`
    - Archaeal outgroup: p__Undinarchaeota 
    - Bacterial outgroup: p__Chloroflexota
5. `gtdbtk decorate`
6. `workflow/scripts/bootstrap_msa.py`
    - This will subsample the MSA for 100 replicates per-domain.
7. `FastTree`
    - Used to infer the bootstrap replicates phylogeny.
8. `genometreetk bootstrap`
    - Computes the support values from the replicates and decorates the GTDB-Tk tree.
9. `workflow/scripts/phylorank_trusted_taxa.py`
    - This will generate a list of taxonomically stable groups 
     (as per the decorate output) to infer RED ranges.
10. `phylorank outliers`
    - Determine the RED ranges for stable groups, and RED scale the tree.

### 2b. Marker genes

The marker gene pipeline computes the individual marker gene trees for all taxa.
These are later used for identifying questionable taxa.

1. Download GTDB [representative marker genes](https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/)
2. Include marker genes from `gtdbtk align` for each marker alignment.
3. Remove columns with low parsimony and blank sequences.
4. Create marker gene trees for each alignment using FastTree.

### 2c. Summary

This pipeline generates plots, and analyses the novelty of each taxon.
See each script / output for more details.

1. `workflow/scripts/summary_novelty_of_genomes.py`
2. `workflow/scripts/plot_marker_summary.py`
3. `workflow/scripts/plot_novelty_at_red_ranks.py`

## 3. References

* [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)
* [GenomeTreeTk](https://github.com/donovan-h-parks/GenomeTreeTk)
* [Phylorank](https://github.com/donovan-h-parks/PhyloRank)
* [PhyloDM](https://github.com/aaronmussig/PhyloDM)
