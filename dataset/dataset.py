import GEOparse
import pandas as pd
import numpy as np

def get_formated_data_frames():
  gse = GEOparse.get_GEO(geo="GSE22552")
  dataframes = []
  for key, val in gse.gsms.items():
    table = val.table.set_index("ID_REF").rename(columns={"VALUE": key})
    dataframes.append(table)
    
  data = pd.concat(dataframes, axis=1)

  genes_comprehensible_data = gse.gpls["GPL570"].table.rename(columns={"ID": "ID_REF"}).set_index('ID_REF').loc[: , ('Gene Symbol', 'ENTREZ_GENE_ID')]

  def remove_trailing_slashes(col):
    arr = np.array(col)
    return [str(gene).split('/')[0] for gene in arr]
      
  refined_gpl = genes_comprehensible_data[['Gene Symbol', 'ENTREZ_GENE_ID']].apply(remove_trailing_slashes)

  refined_merged_dataframe = pd.merge(left=refined_gpl, right=data, left_index=True, right_index=True).reset_index().set_index(['ID_REF','Gene Symbol', 'ENTREZ_GENE_ID'])
  phenotypes = gse.phenotype_data
  removed_data_sets = list(phenotypes[phenotypes['title'].str.contains('Uns_')].index.values)
  dataset_pre_processed = refined_merged_dataframe.drop(removed_data_sets, axis=1)

  dataset_pre_processed = dataset_pre_processed[dataset_pre_processed[dataset_pre_processed < 6].count(axis=1) < 10]
  df_wiht_new_titles = dataset_pre_processed.rename(columns={ old: gse.gsms[old].metadata["title"][0].split('-')[0] for old in gse.gsms.keys() })

  return df_wiht_new_titles