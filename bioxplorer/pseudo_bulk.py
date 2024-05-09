
import anndata
import numpy as np


def convert_to_pseudobulk(adata, patient_col='donor_id'):
    """
    This function is biased towards the scRNA-seq data structure of the adata object.

    Args:
        adata (anndata.AnnData): AnnData object with scRNA-seq data.
        patient_col (str): Column name of the patient ID.
    """
    # Group data by patient and aggregate the metadata
    pbulk_obs = adata.obs.groupby(patient_col).agg(lambda x: list(x.unique()))

    # Group data by patient and mean the expression matrices
    pbulk_x = adata.obs.groupby(patient_col).apply(lambda sub_df: adata[sub_df.index].X.mean(axis=0))

    # Convert the aggregated data into an AnnData object
    pdata = anndata.AnnData(
        X=np.vstack(pbulk_x.values),
        obs=pbulk_obs,
        var=adata.var.copy()
    )

    # to avoid anndata save error
    pdata.X = pdata.X.tolist()

    return pdata