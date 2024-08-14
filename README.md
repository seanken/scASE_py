# scASE_py

This is a basic python package to enable one to work with the output of our single cell allele specific expression pipeline (ASE). The pipeline itself is present at: https://github.com/seanken/ASE_pipeline

An R package to work with the pipeline, as well as other downstream analysis, is at: https://github.com/seanken/scAlleleExpression

## Install package

Need to install anndata and pandas.

Once that is done can install with `pip install git+https://github.com/seanken/scASE_py`

## Basic usage

The basic object around which this package is built is a meta data table implemented as a pandas data frame. Can be constructed with many tools, including scanpy. Each row represents a cell, each column a different piece of meta data (cell type, sample of origin, etc). We assume a basis meta data object has been created all ready.

The first step is to add ASE information to the meta data. In particular, let `meta` be the data frame with the meta data information, and `sampleCol` the column with sample of origin information (which 10X run/individual it comes from). Then one needs a dictionary, here called `dirs`, that maps from the sample names in sampleCol to the associated output directories from the ASE pipeline above. Also requires a `tag`, which is a character string added to the new column names being created as a prefix. Then one can run:

```
from scASE_py import *
meta=UpdateMeta(meta,dirs,sampleCol,tag)
```

to update the meta data with new columns pointing to the outputs of the ASE pipeline.

One can then use this package to load ASE data. To load SNP level ASE need a list `genes` of genes and a list `snps` of snps (must be same length, where the ith snp corresponds to the ith gene). One can then run:

```
dat=LoadSNPLevel(meta,snps,genes,tag)
```

To get a data frame, `dat`, with SNP level ASE (so for a particular snp/gene pair gets the allele specific expression of that gene, aligned so the haplotype with the reference value for the given SNP is the reference value for that gene). Can also get gene level ASE (same as SNP level except uses a default orientation for each haplotype) using:

```
dat=LoadGeneLevel(meta,genes,tag)
```

One can also load non-allele specific expression information. In particular to load pseudobulked expression information (one entry per sample) for a givem list of genes one can run:

```
dat=loadExpression(meta,genes,tag)
```

where dat is an AnnData object. There is also a method, `readStarSolo`, to load StarSolo outputs as AnnData objects. 

Finally there is also experimental support for running scDali (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02593-8) on the meta data we created above with the method `Run_scDali`, though this approach often crashes.




