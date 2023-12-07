import seaborn as sns
import localization
import separation
import proximity
import preprocessing
import matplotlib.pylab as plt


def get_disease_module_info(dis_name, gda, ppi, gppi):
    
    d = gda[gda.diseaseName == dis_name]
    genes = list(d.geneSymbol.unique())
    genes_in_ppi = [n for n in gppi.nodes if n in genes]
    genes_in_lcc = localization.get_lcc(ppi, genes_in_ppi)
    
    print('Number of disease genes: ',len(genes))
    print('Number of disease genes in the PPI: ',len(genes_in_ppi))
    print('Number of disease genes in the LCC: ',len(genes_in_lcc))
    
    return genes, genes_in_ppi, genes_in_lcc

def plot_disease_module_info(ppi, genes):
    
    sig_dict = localization.get_random_comparison(ppi, genes[1], 1000)
    random_lcc = sig_dict['LCC_list']
    print('Full randomization')
    print('Mean: ',sig_dict['mean'])
    print('Std: ',sig_dict['std'])
    print('z-score: ',sig_dict['z_score'])
    print('p-value: ',sig_dict['p_value'])

    sig_dp_dict = localization.get_random_comparison(ppi, genes[1], 1000, degree_preserving=True)
    random_dp_lcc = sig_dp_dict['LCC_list']
    print('\nDegree preserving randomization')
    print('Mean: ',sig_dp_dict['mean'])
    print('Std: ',sig_dp_dict['std'])
    print('z-score: ',sig_dp_dict['z_score'])
    print('p-value: ',sig_dp_dict['p_value'], '\n')
    
    fig, axs = plt.subplots(1, 2, figsize=(7,3))

    sns.histplot(data=random_lcc,
                 bins=10, 
                 ax=axs[0])
    axs[0].axvline(len(genes[2]), color='r')
    axs[0].set_xlabel('LCC value')
    axs[0].set_ylabel('Count iterations')
    axs[0].set_title('Full randomization')

    sns.histplot(data=random_dp_lcc,
                 bins=10, 
                 ax=axs[1])
    axs[1].axvline(len(genes[2]), color='r')
    axs[1].set_xlabel('LCC value')
    axs[1].set_ylabel('Count iterations')
    axs[1].set_title('Degree preserving randomization')

    plt.tight_layout()
    plt.show()