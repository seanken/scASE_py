from scdali import run_tests
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


def Run_scDali(meta: pd.DataFrame,dat: pd.DataFrame,tag: str,scDaliMode: str="scDALI-Het",cellType: str="CellType") -> pd.DataFrame:
    ##Prep data

    datTemp=dat.pivot_table(index=['CBC','Sample','Gene'],columns='Allele',values='Count',fill_value=0).reset_index()
    datTemp["Tot"]=datTemp["All1"]+datTemp["All2"]
    #rowInfo=datTemp.pivot_table(index=['CBC','Sample'],columns='Gene',values='All1',fill_value=0).reset_index()[["CBC","Sample"]]
    #rowInfo.columns=["CBC",tag+"Sample"]
    A=datTemp.pivot_table(index=['CBC','Sample'],columns='Gene',values='All1',fill_value=0)
    D=datTemp.pivot_table(index=['CBC','Sample'],columns='Gene',values='Tot',fill_value=0).to_numpy()
    meta=meta.set_index(["CBC",tag+"Sample"]).loc[A.index]
    genes=[i for i in A.columns]
    A=A.to_numpy()

    #meta[cellType]=meta[cellType].astype('category')
    #meta["cellState"]=meta[cellType].cat.codes
    cell_state=pd.get_dummies(meta[cellType],columns=[cellType],dtype=float).to_numpy()
    sum1=D.sum(axis=0)
    keep=[i for i in range(0, len(sum1)) if sum1[i]>100]
    A=A[:,keep]
    D=D[:,keep]
    genes=[genes[i] for i in keep]

    ##Run scdali
    pval = run_tests(A=A, D=D, model=scDaliMode, cell_state=cell_state, n_cores=1,base_rate=.5)['pvalues']
    mrk=pd.DataFrame(genes)
    mrk.columns=["Gene"]
    mrk["pval"]=pval.tolist()
    mrk=mrk.sort_values(['pval'])
    mrk["FDR"]=fdrcorrection(mrk["pval"])[1].tolist()
    return(mrk)