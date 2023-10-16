import pandas as pd
import anndata
import scanpy as sc
##Will load expression

def loadExpression(meta: pd.DataFrame,genes: list[str],tag:str) -> anndata:
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

def getExpression(meta: pd.DataFrame,tag: str,genes: list[str]) -> anndata.AnnData:
    ##prep
    cells=meta["CBC"].tolist()
    rawCountFile=meta[tag+"Raw"].tolist()[0]

    #load data
    dat=sc.read_10x_mtx(rawCountFile)
    data.var_names_make_unique()
    dat=dat[:,cells]
    if len(genes)>0:
        dat=dat[genes,:]

    return(dat)