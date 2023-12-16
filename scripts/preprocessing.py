import pandas as pd


def preprocess_disease(ppi_file, dga_file):
    # Load PPI data
    ppi = pd.read_csv(ppi_file)
    ppi = ppi[["Symbol_A", "Symbol_B"]].drop_duplicates().dropna()

    # Load disease-gene association data
    dga = pd.read_csv(dga_file, sep="\t")

    # Filter out disease types that are 'group' or 'phenotype'
    dga = dga[(dga.diseaseType != "group") & (dga.diseaseType != "phenotype")]

    # Lowercase disease names
    dga["diseaseName"] = dga["diseaseName"].str.lower()

    # Count the number of genes associated with each disease
    num_genes = (
        dga.groupby("diseaseName")
        .agg("count")
        .sort_values(by="geneSymbol")
        .reset_index()
        .rename(columns={"geneSymbol": "count_genes"})
    )

    filtered_dga = dga.merge(num_genes, on="diseaseName", how="inner")

    dga = (
        filtered_dga[filtered_dga.count_genes > 10][["geneSymbol", "diseaseName"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    return dga, ppi
