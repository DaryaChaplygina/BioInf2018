import pandas as pd


def check_mut(s):
    if s == 'SYNONYMOUS_CODING':
        return 'syn'
    else:
        return 'nonsyn'


if __name__ == '__main__':
    table = pd.read_csv('ESP_SNP_data.tsv', sep='\t')
    table['MUT'] = table['EFF'].apply(check_mut)
    sum_freq = table.groupby(by=['GENE', 'MUT']).sum().AF.reset_index(level=['GENE', 'MUT'])
    dnds = pd.merge(sum_freq.where(sum_freq.MUT == 'syn').dropna(how='all'),
                    sum_freq.where(sum_freq.MUT == 'nonsyn').dropna(how='all'),
                    how='left', on=['GENE'])
    dnds['dNdS'] = dnds.AF_y / dnds.AF_x
    dnds.sort_values('dNdS', inplace=True)
    dnds[['GENE', 'dNdS']].head(200).to_csv('top_genes.csv', index=False)
    print(dnds.GENE.head(200))
