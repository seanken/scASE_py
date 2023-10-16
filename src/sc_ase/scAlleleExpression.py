import pandas as pd
pd.options.mode.chained_assignment = None

##
##A python package to load the output of our ASE pipeline
##


##Modifies a given meta data pandas DataFrame to work with this package by adding new columns related to the pipeline output
def UpdateMeta(meta: pd.DataFrame,dirs: dict[str,str],sampleCol: str,tag: str) -> pd.DataFrame:

    ##Check sample column is present
    if not sampleCol in meta.columns:
        raise Exception("sampleCol needs to be the name of a column in the meta object")

    ##Add pipeline information to meta data
    meta[tag+"dirASE"]=[dirs.get(x,"NA") for x in meta[sampleCol]]
    meta[tag+"Allele"]=[x+"/output/AlleleCounts/counts.txt" for x in meta[tag+"dirASE"]]
    meta[tag+"SNP"]=[x+"/output/SNPLevelCounts/snps.bed" for x in meta[tag+"dirASE"]]
    meta[tag+"Raw"]=[x+"/output/STARSolo/output/resultsSolo.out/GeneFull/raw" for x in meta[tag+"dirASE"]]
    meta[tag+"Sample"]=meta[sampleCol]
    return(meta)

##Given a meta data table prepperd with UpdateMeta and a list of genes, get the ASE information for those genes
def LoadGeneLevel(meta: pd.DataFrame,genes: list[str],tag: str) -> pd.DataFrame:
    tab=AggPerSamp(meta,tag)

    ##Split meta data into pieces, one for each sample
    splitMetaBySamp=[meta[meta[tag+"Sample"]==x] for x in tab[tag+"Sample"]]

    ##Get the Gene level ASE and combined into on data frame
    aseList=[getGeneASE(metaSplit,genes,tag) for metaSplit in splitMetaBySamp]
    dat=pd.concat(aseList,0)

    return(dat)

##Given a meta data table prepperd with UpdateMeta and a list of genes and a list of snps, get the ASE information for those genes with ref vs alt alignment based on the SNP
def LoadSNPLevel(meta: pd.DataFrame,snps: list[str],genes: list[str],tag: str) -> pd.DataFrame:
    tab=AggPerSamp(meta,tag)

    ##Split meta data into pieces, one for each sample        
    splitMetaBySamp=[meta[meta[tag+"Sample"]==x] for x in tab[tag+"Sample"]]

    ##Get the SNP level ASE and combined into on data frame
    aseList=[getSNPASE(metaSplit,snps,genes,tag) for metaSplit in splitMetaBySamp]
    dat=pd.concat(aseList,0)
    return(dat)

##Gets a table with one entry per sample, with columns corresponding to sample name, SNP file, and Allele file from the pipeline
def AggPerSamp(meta: pd.DataFrame,tag: str) -> pd.DataFrame:
    ##Check for columns
    if not tag+"Sample" in meta.columns or not tag+"SNP" in meta.columns or not tag+"Allele" in meta.columns or not tag+"Raw" in meta.columns:
        raise Exception("Make sure tag is correct and the meta data was prepared with UpdateMeta")

    return(meta[[tag+"Sample",tag+"SNP",tag+"Allele",tag+"Raw"]].drop_duplicates())

##Given the meta data for one sample, gets the gene level ASE for that sample in the specified genes
def getGeneASE(meta: pd.DataFrame,genes: list[str],tag: str) -> pd.DataFrame:
    ##Checks for columns
    if not tag+"Sample" in meta.columns or not tag+"SNP" in meta.columns or not tag+"Allele" in meta.columns:
        raise Exception("Make sure tag is correct and the meta data was prepared with UpdateMeta")
    if not "CBC" in meta.columns:
        raise Exception("Need CBC column in meta data")

    ##Prep
    countFile=meta[tag+"Allele"].tolist()[0]
    samp=meta[tag+"Sample"].tolist()[0]
    cells=meta["CBC"].tolist()
    print(samp)

    ##Load gene level expression from pipeline
    cts=pd.read_csv(countFile,sep=" ",names=["CBC","Gene","Allele","Count"])
    cts=cts[cts["CBC"].isin(cells)]
    cts=cts[cts["Gene"].isin(genes)]
    cts.loc[:,"Sample"]=[samp for i in cts["CBC"]]
    cts=cts[cts["Allele"]!="Ambig"]
    print(cts.shape)
    return(cts)


##Given the meta data for one sample, gets the gene level ASE for that sample in the specified genes aligned with the specified SNPs
def getSNPASE(meta: pd.DataFrame, snps: list[str],genes: list[str],tag: str) -> pd.DataFrame:
    ##Get gene level counts
    cts_gene=getGeneASE(meta,genes,tag)

    ##Load SNP genotype information
    snpFile=meta[tag+"SNP"].tolist()[0]
    snpTable=pd.read_csv(snpFile,sep="\t",names=["chr","start","end","SNP","Phase"])
    snpTable=snpTable[snpTable["SNP"].isin(snps)]

    ##Remove genes and SNPs not found in this sample
    genes_new=[genes[i] for i in range(0,len(genes)) if genes[i] in cts_gene["Gene"].tolist() and snps[i] in snpTable["SNP"].tolist()]
    snps_new=[snps[i] for i in range(0,len(genes)) if genes[i] in cts_gene["Gene"].tolist() and snps[i] in snpTable["SNP"].tolist()]
    genes=genes_new
    snps=snps_new

    ##Reorder alleles for each gene based on associated SNP genotype
    cts_by_gene=[cts_gene[cts_gene["Gene"]==gene] for gene in genes]
    cts_by_snp=[toSNPCounts(cts_by_gene[i],snps[i],snpTable) for i in range(0,len(genes))]
    cts_snp=pd.concat(cts_by_snp)
    print(cts_snp.shape)
    print(cts_snp.head())
    return(cts_snp)


##For a gene ASE count table (cts) and a snp genotype (value in 0|1 or 1|0), realigns the Allele in the cts column to match the SNP
def toSNPCounts(cts: pd.DataFrame,snp: str,snpTable: pd.DataFrame) -> pd.DataFrame:
    ##Gets the phased of this snp (0|1 or 1|0)
    listSNP=snpTable[snpTable["SNP"]==snp]["Phase"].tolist()
    phaseSNP=listSNP[0]

    ##If needed flips Allele so matches SNP phase
    if phaseSNP=="1|0":
        swapAll=lambda a: "All1" if a=="All2" else "All2"
        cts.loc[:,"Allele"]=[swapAll(allele) for allele in cts["Allele"]]
    
    ##If not het SNP removes (shouldn't be an issue)
    if phaseSNP!="1|0" and phaseSNP!="0|1":
        cts["Count"]=[0 for i in cts["Count"]]
    cts.loc[:,"SNP"]=[snp for i in cts["Count"]]
    return(cts)