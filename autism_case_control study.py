#!/usr/bin/env python
# coding: utf-8

# # Steps to be followed prior 
#    
#    ```1.extract the common variants from all the vcf files.```
#    
#    `2.create a snp matrix of genotypes.`
#    
#    `3.create a pedigree file.`
#      
#    `4.create a  Map file.`  
#    
# 
#    
#   

# ### Format of SNP matrix
# ![Snp_matrix%20sample.png](attachment:Snp_matrix%20sample.png)

# ### Format of pedigree file
# ![Screenshot%202024-01-20%20114220.png](attachment:Screenshot%202024-01-20%20114220.png)

# ### Format of map file
# 
# ![Screenshot%202024-01-20%20114936.png](attachment:Screenshot%202024-01-20%20114936.png)

# ## From the available VCF files on autism.
# 
# ``` cases```   :- 2,9,7,17,11,13
# 
# ```Controls``` :- rest of cases

# In[1]:


# Here is the gene list of interest for Autism

gene_list = ["AMT", "NAGA", "PEX7", "SYNE1", "VPS13B", "PAH", "POMGNT1", "CYB5R3", "TCF20", "CHD8", "NF1", "FMR1", "TSC1",
            "TSC2", "UBE3A", "NIPBL", "RAD21", "SMC3", "HDAC8", "SMC1A", "MECP2", "PTPN11", "KRAS", "SOS1", "RIT1", "RAF1",
            "CHD7", "CEP290", "CPLANE1", "CC2D2A", "INPP5E", "KIAA0586", "MKS1", "RPGRIP1L", "TCTN2", "TMEM67", "TMEM216", 
            "CACNA1C", "ARID1B", "DYRK1A", "POGZ", "SHANK3", "SYNGAP1", "ASH1L", "CHD2", "KMT2C", "CELF4", "KMT2E", "KMT5B",
            "DEAF1", "EIF3G", "MKX", "ELAVL3", "NCOA1", "PAX5", "TBR1", "ZMYND8", "PHF2", "KDM6B", "ANK2", "AP2S1", "RFX3",
            "RORB", "SATB1", "SKI", "SMARCC2", "KCNMA1", "SHANK2", "CACNA2D3", "PTEN", "NRXN1", "DSCAM", "CORO1A", "GFAP",
            "PTK7", "DPYSL2", "MAP1A", "SPAST", "KIAA0232", "NUP155", "PPP5C", "TEK", "TM9SF4", "UBR1", "GNAI1", "HECTD4",
            "ADNP", "IRF2BPL", "SETD5", "ANKRD11", "MBD5", "SIN3A", "MED13L", "TBL1XR1", "ASXL3", "TCF4", "BCL11A", "NACC1",
            "TCF7L2", "NSD1", "CREBBP", "NR3C2", "TLK2", "CTNNB1", "PHF12", "TRAF7", "DNMT3A", "TRIP12", "FOXP1", "PPP2R5D",
            "VEZF1", "FOXP2", "RAI1", "WAC", "CACNA1E", "KCNQ3", "SLC6A1", "GABRB2", "LRRC4C", "STXBP1", "GABRB3", "PRR12", 
            "GRIN2B", "SCN2A", "DYNC1H1", "TAOK1", "SUPT16H", "FOXG1", "ALG14", "YY1", "TMLHE", "ADSL", "MAN1B1", "H4C5", 
            "CHRNA7", "NIPA1", "TUBG1", "NIPA2", "NAA20", "HDAC4", "EHMT1", "MKRN3", "HERC2", "NPAP1", "SNORD115-1", "PWAR1",
            "SNORD116-1", "IPW", "MAGEL2", "PWRN1", "GAMT", "DPYD", "AHDC1", "GNB1", "SATB2", "EXT2", "DHCR7", "SETD2",
            "SEMA3E", "KMT2A", "NHS", "REV3L", "PLXND1", "SNRPN", "NDN", "OCA2", "ALG13", "GRIA3", "MAOA", "FTSJ1", "FLCN",
            "STS", "ATRX", "COMT", "JMJD1C", "GP1BB", "RREB1", "TBX1", "HIRA", "SEC24C", "UFD1", "ARVCF", "SEC23B", "SDHB",
            "SDHD", "SDHC", "AKT1", "USF3", "PIK3CA", "KLLN", "CHD1", "ARCN1", "WFS1", "TBCK", "PDE4D", "KCNA2", "DNM1",
            "AP3B2", "ATP6V1A", "SLC38A3", "FZR1", "EEF1A2", "FBXO28", "GABRA5", "GRIN2D", "NUS1", "CLTC", "WWOX", "SZT2",
            "HCN1", "FGF13", "NECAP1", "PACS2", "AARS1", "SCN8A", "GABBR2", "GABRA2", "DALRD3", "CYFIP2", "NTRK2", "PARS2",
            "CDK19", "ARV1", "FGF12", "PPP3CA", "ATP1A3", "TRAK1", "GABRG2", "YWHAG", "UBA5", "SCN3A", "CACNA2D1", "CELF2", 
            "CACNA1B", "ACTL6B", "SLC1A2", "SYNJ1", "ATP1A2", "KCNB1", "SLC13A5", "CNKSR2", "CACNA1A", "DHDDS", "PIGL", 
            "ALDH5A1", "GATM", "BRD4", "DPAGT1", "FTCD", "LHX1", "HNF1B", "PTCHD1", "SPTBN1", "EP300", "CYP27A1", "TCF12",
            "IFNG", "STX1A", "CLIP2", "DNAJC30", "FKBP6", "RFC2", "GTF2I", "LIMK1", "GTF2IRD1", "VPS37D", "METTL27", "BUD23", 
            "EIF4H", "GTF2IRD2", "BCL7B", "TBL2", "NCF1", "TMEM270", "MLXIPL", "BAZ1B", "ELN", "MED12", "PRKCZ", "GABRD",
            "SPEN", "MMP23B", "LUZP1", "PDPN", "HSPG2", "UBE4B", "KCNAB2", "RERE", "PRDM16", "CASZ1", "RPL10", "CNTNAP2",
            "SLC35C1", "DICER1", "ZBTB20", "FLG", "SMAD4", "NONO", "IQSEC2", "NR2F1", "ALAD", "SLC9A6", "KDM5C", "NFIB",
            "BCKDK", "DOCK7", "NDP", "GLRA2", "IL1RAPL1", "NLGN3", "SH2B1", "CDK13", "GJA8", "GJA5", "MEIS2", "PAK2", "FGFR1",
            "PROKR2", "SOX2", "SOX3", "ARNT2", "HESX1", "OTX2"]


# In[ ]:





# In[2]:


import warnings
warnings.filterwarnings('ignore')


# In[3]:


# Import pandas and numpy libraries
import pandas as pd
import numpy as np

# Define a function that takes a file path and a sheet name as parameters
def get_set_from_file(file_path, sheet_name,gene_list):
    # Read the excel file into a dataframe
    df = pd.read_excel(file_path, sheet_name=sheet_name, skiprows=1)
  
    # Filter the columns that start with 'gene'
    gene_cols = [col for col in df.filter(like='ene')]
    
     # Convert all genes to a single column
    df['d_genes'] = df[gene_cols].apply(lambda x: ','.join(map(str, x)), axis=1)
    
    # Filter the rows that contain any of the genes in the gene list
    df2 = df[[any(x in gene_list for x in y) for y in df['d_genes'].str.split(',')]]

    # Create a new column 'C/p_R/A'
    df2['C/p_R/A'] = df2.iloc[:, [0, 1]].apply(lambda x: '_'.join(x.astype(str)), axis=1)
    # Get the set of variants
    vcf_variants = set(df2['C/p_R/A'])
    
    return df2,vcf_variants


# In[4]:


vcf1,vcf1_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP1 - KHAFAMGPCSP1-WES-Pool-101_S27_L004_R.xlsx',0,gene_list)


# In[5]:


vcf1.shape


# In[6]:


len(vcf1_var)


# In[7]:


vcf2,vcf2_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP2 - KHAFAMGPCSP2_S28_L004_R.xlsx',0,gene_list)
print(vcf2.shape)
print(vcf2_var)


# In[8]:


vcf3,vcf3_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP3 - KHAFAMGPCSP3_S29_L004_R.xlsx',0,gene_list)
print(vcf3.shape)
print(len(vcf3_var))


# In[9]:


vcf4,vcf4_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP4 - KHAFAMGPCSP4_S30_L004_R.xlsx',0,gene_list)
print(vcf4.shape)
print(len(vcf4_var))


# In[10]:


vcf5,vcf5_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP5 - KHAFAMGPCSP5_S31_L004_R.xlsx',0,gene_list)
print(vcf5.shape)
print(len(vcf5_var))


# In[11]:


vcf6,vcf6_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP6 - KHAFAMGPCSP6_S32_L004_R.xlsx',0,gene_list)
print(vcf6.shape)
print(len(vcf6_var))


# In[12]:


vcf7,vcf7_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP7 - KHAFAMGPCSP7_S33_L004_R.xlsx',0,gene_list)
print(vcf7.shape)
print(len(vcf7_var))


# In[13]:


vcf8,vcf8_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP8 - KHAFAMGPCSP8_S34_L004_R.xlsx',0,gene_list)
print(vcf8.shape)
print(len(vcf8_var))


# In[14]:


vcf9,vcf9_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP9 - KHAFAMGPCSP9_S35_L004_R.xlsx',0,gene_list)
print(vcf9.shape)
print(len(vcf9_var))


# In[15]:


vcf10,vcf10_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP10 - KHAFAMGPCSP10_S36_L004_R.xlsx',0,gene_list)
print(vcf10.shape)
print(len(vcf10_var))


# In[16]:


vcf11,vcf11_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP11 - KHAFAMGPCSP11_S37_L004_R.xlsx',0,gene_list)
print(vcf11.shape)
print(len(vcf11_var))


# In[17]:


vcf12,vcf12_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP12 - KHAFAMGPCSP12_S11_L003_R.xlsx',0,gene_list)
print(vcf12.shape)
print(len(vcf12_var))


# In[18]:


vcf13,vcf13_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP13 - KHAFAMGPCSP13_S12_L003_R.xlsx',0,gene_list)
print(vcf13.shape)
print(len(vcf13_var))


# In[19]:


vcf14,vcf14_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP14 - KHAFAMGPCSP14_S2_L004_R.xlsx',0,gene_list)
print(vcf14.shape)
print(len(vcf14_var))


# In[20]:


vcf15,vcf15_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP15 - KHAFAMGPCSP15_S3_L004_R.xlsx',0,gene_list)
print(vcf15.shape)
print(len(vcf15_var))


# In[21]:


vcf16,vcf16_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP16 - KHAFAMGPCSP16_S4_L004_R.xlsx',0,gene_list)
print(vcf16.shape)
print(len(vcf16_var))


# In[22]:


vcf17,vcf17_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP17 - KHAFAMGPCSP17_S5_L004_R.xlsx',0,gene_list)
print(vcf17.shape)
print(len(vcf17_var))


# In[23]:


vcf18,vcf18_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP18 - KHAFAMGPCSP18_S6_L004_R.xlsx',0,gene_list)
print(vcf18.shape)
print(len(vcf18_var))


# In[24]:


vcf19,vcf19_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP19 - KHAFAMGPCSP19_S7_L004_R.xlsx',0,gene_list)
print(vcf19.shape)
print(len(vcf19_var))


# In[25]:


vcf20,vcf20_var = get_set_from_file('/Users/Dell/Downloads/autism_vcf/Germline_KHAFAMGPCSP20 - KHAFAMGPCSP20_S42_L004_R.xlsx',0,gene_list)
print(vcf20.shape)
print(len(vcf20_var))


# In[26]:


common_variants = set.intersection(vcf1_var,vcf2_var,vcf3_var,vcf4_var,vcf5_var,vcf6_var,vcf7_var,vcf8_var,vcf9_var,vcf10_var,
                                  vcf11_var,vcf12_var,vcf13_var,vcf14_var,vcf15_var,vcf16_var,vcf17_var,vcf18_var,vcf19_var,vcf20_var)
cmv=list(common_variants)
len(cmv)


# In[27]:


def filter_vcf(vcf, cmv):
    com_vcf = vcf[vcf['C/p_R/A'].isin(cmv)]
    return com_vcf


# In[28]:


com_vcf1 = filter_vcf(vcf1,cmv)


# In[29]:


com_vcf2 = filter_vcf(vcf2,cmv)


# In[30]:


com_vcf3 = filter_vcf(vcf3,cmv)


# In[31]:


com_vcf4 = filter_vcf(vcf4,cmv)


# In[32]:


com_vcf5 = filter_vcf(vcf5,cmv)


# In[33]:


com_vcf6 = filter_vcf(vcf6,cmv)


# In[34]:


com_vcf7 = filter_vcf(vcf7,cmv)


# In[35]:


com_vcf8 = filter_vcf(vcf8,cmv)


# In[36]:


com_vcf9 = filter_vcf(vcf9,cmv)


# In[37]:


com_vcf10 = filter_vcf(vcf10,cmv)


# In[38]:


com_vcf11 = filter_vcf(vcf11,cmv)


# In[39]:


com_vcf12 = filter_vcf(vcf12,cmv)


# In[40]:


com_vcf13 = filter_vcf(vcf13,cmv)


# In[41]:


com_vcf14 = filter_vcf(vcf14,cmv)


# In[42]:


com_vcf15 = filter_vcf(vcf15,cmv)


# In[43]:


com_vcf16 = filter_vcf(vcf16,cmv)


# In[44]:


com_vcf17 = filter_vcf(vcf17,cmv)


# In[45]:


com_vcf18 = filter_vcf(vcf18,cmv)


# In[46]:


com_vcf19 = filter_vcf(vcf19,cmv)


# In[47]:


com_vcf20 = filter_vcf(vcf20,cmv)


# In[48]:


com_df = [com_vcf1,com_vcf2,com_vcf3,com_vcf4,com_vcf5,com_vcf6,com_vcf7,com_vcf8,com_vcf9,com_vcf10,
         com_vcf11,com_vcf12,com_vcf13,com_vcf14,com_vcf15,com_vcf16,com_vcf17,com_vcf18,com_vcf19,com_vcf20]


# In[49]:


vcf1.columns


# In[50]:


# Initialize an empty dataframe for the result
merged_df_GT= pd.DataFrame()

# Loop through the dataframes and extract the column
for i, df in enumerate(com_df):
    merged_df_GT['vcf'+str(i+1)] = df['0/1 Genotypes (GT)'].values


# In[51]:


merged_df_GT


# In[52]:


header = list(com_vcf1['Identifier.1'])
len(header)


# In[53]:


com_vcf1['Identifier.1'].isna().sum()


# In[54]:


df.drop(["Courses", "Fee"], axis = 1, inplace=True)


# In[55]:


GT_v=merged_df_GT.T
GT_v.columns = header
GT_v.drop(["rs3114914", "rs4819208"],axis=1,inplace=True) 


# In[56]:


GT_v


# In[60]:


def replace_values(df):
    df = df.applymap(lambda x: 1 if x == '1/1' else 2)
    return df


# In[61]:


snpmatrix = replace_values(GT_v)


# In[62]:


genotypes = snpmatrix


# In[63]:


genotypes_gi = GT_V2.filter(GI_ppi)
genotypes_gi.columns


# In[ ]:


diff_columns = set(GI_ppi) - set(genotypes_gi.columns)
diff_columns


# In[ ]:


genotypes['rs311']


# In[ ]:


vcf1.columns


# In[ ]:


GI_m=merged_df_GT.T
GI_m.columns = header


# In[ ]:


GI_m


# In[ ]:


genotypes_gi = GI_m.filter(GI_ppi)
genotypes_gi = replace_values(genotypes_gi)
genotypes_gi


# In[ ]:


Map = com_vcf1[['C/p_R/A','Identifier.2']]


# In[ ]:


Map


# In[ ]:


new_columns = Map['C/p_R/A'].str.split(r'[:_/]', expand=True)


# In[ ]:


new_columns


# In[ ]:


new_columns.columns = ['Chr', 'pos', 'Ref', 'Alt']


# In[ ]:


new_columns


# In[ ]:


PPi_genes = com_vcf1['Gene Names'].unique()


# In[ ]:


len(PPi_genes)


# In[ ]:


PPi_genes.to_text('ppi_genes.txt')


# In[ ]:


np.savetxt('ppi_genes.txt', PPi_genes, fmt='%s', delimiter=',')



# In[ ]:


map2 = pd.concat([new_columns, com_vcf1['Identifier.1']], axis=1)


# In[ ]:


map2['Identifier.1']


# In[ ]:


map2['cm'] = 0


# In[ ]:


map2


# In[ ]:


map2[['Chr','Identifier.1','cm','pos','Ref','Alt']]


# In[ ]:


# Define a dictionary mapping old column names to new ones
column_names = {
    'Chr': 'chr',
    'Identifier.1': 'snp.name',
    'pos': 'gpos',
    'Ref': 'allele.1',
    'Alt': 'allele.2'
}

# Use the .rename() method to change the column names
map2 = map2.rename(columns=column_names)


# In[ ]:


MAP=map2[['chr','snp.name','cm','gpos','allele.1','allele.2']]
MAP[MAP['snp.name'] == 'rs3114914']


# In[ ]:


map_gi = MAP[MAP['snp.name'].isin(GI_ppi)]
map_gi


# In[ ]:


diff_columns = set(GI_ppi) - set(map_gi['snp.name'])
diff_columns


# In[ ]:


S_G = pd.concat([MAP['snp.name'],com_vcf1[ 'Gene Names']],axis=1)


# In[ ]:


column_names = {
   'snp.name': 'snp',
    'Gene Names': 'gene'
}

# Use the .rename() method to change the column names
S_G = S_G.rename(columns=column_names)


# In[ ]:


S_G=S_G.loc[~S_G['snp'].isin(["rs3114914", "rs4819208"])]


# In[ ]:


S_G.to_csv('snp_gene.csv',index=False)


# In[ ]:


S_G


# In[ ]:


ppi = pd.read_table('string_interactions.tsv')


# In[ ]:


ppi= ppi[['#node1','node2']]


# In[ ]:


column_names = {
   '#node1': 'gene1',
    'node2': 'gene2'
}

# Use the .rename() method to change the column names
ppi = ppi.rename(columns=column_names)


# In[ ]:


ppi.to_csv('gene_gene.csv',index=False)


# In[ ]:


ppi.shape


# In[ ]:


unique_genes = pd.concat([ppi['gene1'], ppi['gene2']]).unique()
unique_genes


# In[ ]:


s_g = S_G[S_G['gene'].isin(unique_genes)]
GI_ppi=(s_g['snp'])


# In[ ]:


s_g.to_csv('/Users/Dell/Downloads/autism_vcf/snpmapping.csv',index=False)


# In[ ]:


GI_ppi


# In[ ]:


fam = pd.read_excel('autism_fam.xlsx')


# In[ ]:


fam


# In[ ]:


import pickle

# Assuming file1, file2, file3 are your data
data = [genotypes, map2, fam]

with open('data.pkl', 'wb') as f:
    pickle.dump(data, f)


# In[ ]:


MAP = MAP.loc[~MAP['snp.name'].isin(["rs3114914", "rs4819208"])]


# In[ ]:


genotypes.to_csv('/Users/Dell/Downloads/autism_vcf/genotypes.csv',index_label=0)
fam.to_csv('/Users/Dell/Downloads/autism_vcf/fam.csv')
MAP.to_csv('/Users/Dell/Downloads/autism_vcf/map.csv')
genotypes_gi.to_csv('/Users/Dell/Downloads/autism_vcf/genotypes_gi.csv',index_label=0)
map_gi.to_csv('/Users/Dell/Downloads/autism_vcf/map_gi.csv')


# In[ ]:


MAP[MAP['snp.name'].isin(["rs3114914", "rs4819208"])]


# In[ ]:


covars = genotypes.index.tolist()
covars = pd.DataFrame(covars)
covars.columns = ['sample' for _ in range(len(covars.columns))]

covars.to_csv('covars.csv',index=False)


# In[ ]:


MAP.loc[~MAP['snp.name'].isin(["rs3114914", "rs4819208"])]

