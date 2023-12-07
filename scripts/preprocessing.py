import pandas as pd

def preprocess_disease(ppi_file, dga_file):

    # Load PPI data
    ppi = pd.read_csv(ppi_file)
    ppi = (ppi[['Symbol_A','Symbol_B']]
           .drop_duplicates()
           .dropna())

    # Load disease-gene association data
    dga = pd.read_csv(dga_file, sep='\t')

    # Lowercase disease names
    dga['diseaseName'] = dga['diseaseName'].str.lower()

    # Filter out disease types that are 'group' or 'phenotype'
    cleaned_dga = dga[dga.diseaseType != 'group']
    cleaned_dga = cleaned_dga[cleaned_dga.diseaseType != 'phenotype']
    cleaned_dga = cleaned_dga[['geneSymbol', 'diseaseName']].drop_duplicates()

    # Count the number of genes associated with each disease
    num_genes = (cleaned_dga.groupby('diseaseName')
                 .agg('count')
                 .reset_index()
                 .rename(columns={'geneSymbol':'count_genes'}))

    # Merge with the original data and filter diseases with at least 10 genes
    dga = cleaned_dga.merge(num_genes, on='diseaseName', how='inner')
    dga = dga[dga.count_genes >= 10]

    return dga, ppi

