import pandas as pd
import anndata as ad
import scanpy as sc
from anndata import read_mtx 
from .scAlleleExpression import AggPerSamp
##Will load expression

def loadExpression(meta: pd.DataFrame,genes: list[str],tag:str) -> ad.AnnData:
    ##Checks correct format
    if not tag+"Sample" in meta.columns or not tag+"Raw" in meta.columns:
        raise Exception("Make sure tag is correct and the meta data was prepared with UpdateMeta")
    if not "CBC" in meta.columns:
        raise Exception("Need CBC column in meta data")

    ##Get sample level information
    tab=AggPerSamp(meta,tag)

    #split meta data and load anndata
    splitMetaBySamp=[meta[meta[tag+"Sample"]==x] for x in tab[tag+"Sample"]]
    listExpression=[getExpression(metaSamp,tag,genes) for metaSamp in splitMetaBySamp]
    dat = ad.concat(listExpression)

    return(dat)

def getExpression(meta: pd.DataFrame,tag: str,genes: list[str]) -> ad.AnnData:
    print("Load dataset!")
    ##prep
    cells=meta["CBC"].tolist()
    rawCountFile=meta[tag+"Raw"].tolist()[0]

    #load data
    dat=readStarSolo(rawCountFile)
    dat.var_names_make_unique()
    dat=dat[cells,:]
    dat.obs=meta
    if len(genes)>0:
        dat=dat[:,genes]

    return(dat)




def readStarSolo(dir: str) -> ad.AnnData:
    ##Load matrix
    adata=read_mtx(dir+"/matrix.mtx").T

    ##Get gene names
    genes=pd.read_csv(dir+"/features.tsv",sep="\t",names=["gene_id","gene_name","type"])
    adata.var_names = genes["gene_name"].tolist()

    ##Get barcodes
    barcodes= pd.read_csv(dir + '/barcodes.tsv',names=["CBC"],sep="\t")["CBC"].tolist()
    adata.obs_names=barcodes

    return(adata)


