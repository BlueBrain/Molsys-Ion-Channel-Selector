#####
import pandas as pd
import numpy as np
import json
import mygene

# Resoures

# path_to_scRNAseq_data = './input/medians.csv'
path_to_scRNAseq_data = './input/Yao_et_al_trim_mean25.csv'
path_to_selected_channels = './input/Channels_genes_(correspondance_channels)_v2.csv'
path_to_BBP_m_type_list = './input/BBP_mtype_list.csv'

path_to_inh_map_L1 = './input/P(BBPmarker_metype)_L1_(Gouw+pseq_BBP)April_16_2021.csv'
path_to_inh_map_L26 = './input/P(BBPmarker_metype)_L26_(Gouw+pseq_BBP)April_16_2021.csv'


def unique_element(array):
    unique_ = []
    for x in array:
        if x not in unique_:
            unique_.append(x)
    return np.asarray(unique_)

def count_elements(array):
    """Return as and pandas DataFrame unique elements and the associated counts of an array."""
    
    unq = unique_element(array)
    
    return pd.DataFrame([len(array[[y==x for y in array]]) for x in unq], index=unq, columns=['counts'])

def sigmoid_f(x, a, mu):
    return 1 / (1 + np.exp(-a*(x - mu)))

def gaussian_f(x, a, mu, sigma):
    return a * np.exp(-((x - mu)**2 / (2*sigma**2)))

def generate_panda(data, ttype):
    """ generate a panda data frame for a given t-type from the json"""
    df = pd.DataFrame(data[ttype]['values'],
                      index=data[ttype]['index'],
                      columns=data[ttype]['columns'])
    return df

# def preprocess_X_df(X_df):
    
#     dict_class={}
#     for x in X_df.index:
#         if ('Lamp5' in x)|('Sncg' in x)|('Serpinf1' in x)|('Vip' in x)|('Sst' in x)|('Pvalb' in x):
#             dict_class[x]='GABAergic'
#         elif ('IT' in x)|('Car3' in x)|('ET' in x)|('NP' in x)|('CT' in x)|('L6b' in x)|('Ly6g6e' in x):
#             dict_class[x]='Glutamatergic'
#         elif ('CA' in x):
#             dict_class[x]='Hippocampus'
#         else:
#             dict_class[x]='Non-neuronal'

#     msk_0 = [(c == 'Glutamatergic' )|(c == 'GABAergic') for c in X_df.rename(index=dict_class).index]
#     msk_1 = [c == 'GABAergic' for c in X_df.rename(index=dict_class).index]
#     msk_2 = [c == 'Glutamatergic' for c in X_df.rename(index=dict_class).index]
#     msk_3 = [c == 'Non-neuronal' for c in X_df.rename(index=dict_class).index]

#     Xnrn_df = X_df[msk_0]
#     Xinh_df = X_df[msk_1]
#     Xexc_df = X_df[msk_2]
#     Xnon_nrn_df = X_df[msk_3]
    
#     return Xnrn_df, Xinh_df, Xexc_df, Xnon_nrn_df

def preprocess_X_df(X_df):
    
    dict_class={}
    for x in X_df.index:
        if ('Gaba' in x):
            dict_class[x]='GABAergic'
        elif ('Glut' in x):
            dict_class[x]='Glutamatergic'
        elif ('Chol' in x)|('Dopa' in x)|('Sero' in x)|('Gly' in x)|('Glyc' in x)|('Hist' in x):
            dict_class[x]='modulatory'
        else:
            dict_class[x]='Non-neuronal'

    msk_0 = [(c == 'Glutamatergic' )|(c == 'GABAergic') for c in X_df.rename(index=dict_class).index]
    msk_1 = [c == 'GABAergic' for c in X_df.rename(index=dict_class).index]
    msk_2 = [c == 'Glutamatergic' for c in X_df.rename(index=dict_class).index]
    msk_3 = [c == 'Non-neuronal' for c in X_df.rename(index=dict_class).index]

    Xnrn_df = X_df[msk_0]
    Xinh_df = X_df[msk_1]
    Xexc_df = X_df[msk_2]
    Xnon_nrn_df = X_df[msk_3]
    
    return Xnrn_df, Xinh_df, Xexc_df, Xnon_nrn_df


def make_binary_1Perc(X_df):
    Binary_1Perc = []
    
    for ttype in X_df.index:

        channels_expr_ = X_df.T[ttype]
        threshold = 0.01*np.sum(channels_expr_.values)

        Binary_1Perc.append(np.where(channels_expr_>threshold, 1, 0))
        
    Binary_1Perc_df = pd.DataFrame(np.asarray(Binary_1Perc),
                              index = X_df.index,
                              columns = X_df.columns)
    return Binary_1Perc_df

def make_binary_0thresh(X_df):
    Binary_0thresh = []
    
    for ttype in X_df.index:

        channels_expr_ = X_df.T[ttype]

        Binary_0thresh.append(np.where(channels_expr_ > 0., 1, 0))
        
    Binary_0thresh_df = pd.DataFrame(np.asarray(Binary_0thresh),
                                  index = X_df.index,
                                  columns = X_df.columns) 
    return Binary_0thresh_df

## TO DO: Implement gaussian thresh

def compute_default_distribution_file(Binary_df):
    Distribution_df = pd.DataFrame(
        np.asarray([['uniform']*len(Binary_df.columns)]*len(Binary_df.index)),
        columns = Binary_df.columns,
        index = Binary_df.index)
    return Distribution_df

def compute_IC_data(Binary_df, Distribution_df):
    
    IC_data = {}
    for ttype in Binary_df.index:
        df = pd.concat([Binary_df.T[ttype], 
                         Distribution_df.T[ttype], 
                         Distribution_df.T[ttype], 
                         Distribution_df.T[ttype]], 
                        axis=1)
        df.columns = ['presence', 'dendrites distribution', 'axon distribution', 'soma distribution']
#         df.rename(Channels_genes['Channel'].T, axis=1) # if you want to have protein names instead of channel symbols
        IC_data[ttype] = {'values': df.values.tolist(), 'index': df.index.tolist(), 'columns': df.columns.tolist()}
    
    return IC_data

def compute_default_exc_map(path_to_BBP_m_type_list, ttype_list):
    df = pd.read_csv(path_to_BBP_m_type_list, index_col=0)
    BBP_mtypes = df['m-type'].values
    msk_exc_bbp = np.asarray([('TPC' in x)|('BPC' in x)|('UPC' in x)|('IPC' in x)|('SSC' in x)|('HPC' in x)
                              for x in BBP_mtypes])
    
    lay_list_dict = {
        'L2' : ['L2/3'],
        'L23': ['L2/3'],
        'L3': ['L2/3'],
        'L4' : ['L4', 'L4/5'],
        'L5' : ['L4/5', 'L5', 'NP'],
        'L6' : ['L5/6', 'L6']
    }

    map_exc_ttype = {}

    for mtype in BBP_mtypes[msk_exc_bbp]:
        lay = mtype.split('_')[0]
        ttype_match = []
        for lay_el in lay_list_dict[lay]:
            msk_tmp = [lay_el in x for x in ttype_list]
            ttype_match += np.asarray(ttype_list)[msk_tmp].tolist()
        ttype_match = unique_element(np.asarray(ttype_match))
        map_exc_ttype [mtype] = ttype_match.tolist()
        
    return map_exc_ttype

def make_inh_map_binary(path_to_inh_map_L1, path_to_inh_map_L26):
    
    map_inh_ttype_L1 = pd.read_csv(path_to_inh_map_L1, index_col=0)
    map_inh_ttype_L26 = pd.read_csv(path_to_inh_map_L26, index_col=0)
    map_inh_ttype_L26 = map_inh_ttype_L26.reindex(['Vip', 'Lamp5', 'Pvalb', 'Sst', 'Sncg', 'Serpinf1'])

    map_inh_ttype = pd.concat([map_inh_ttype_L1, map_inh_ttype_L26], axis=1, sort=True).fillna(0)
    map_inh_ttype_binary = pd.DataFrame(np.where(map_inh_ttype>0., 1, 0),
                                         index = map_inh_ttype.index,
                                         columns = map_inh_ttype.columns)
    return map_inh_ttype_binary

if __name__ == "__main__":
    
    medians = pd.read_csv(path_to_scRNAseq_data, index_col=0)
    Channels_genes = pd.read_csv(path_to_selected_channels, index_col = 'gene_symbol')
    Channels_genes = Channels_genes.rename({'Scn2a1': 'Scn2a', 'Clcn4-2': 'Clcn4'}, axis=0).drop('Clca4c-ps', axis=0)
    msk_genes = [g in Channels_genes.index.tolist() for g in medians.index]
    X_df = medians.T[medians.index[msk_genes]]
    print(X_df)
    
    Xnrn_df = preprocess_X_df(X_df)[0]
    
    Binary_0thresh_df = make_binary_0thresh(Xnrn_df)
    
    Distribution_df = compute_default_distribution_file(Binary_0thresh_df)
    
    IC_data = compute_IC_data(Binary_0thresh_df, Distribution_df)
    

    with open('./output/t_types_IC_expression.json', 'w') as fp:
        json.dump(IC_data, fp)
        
    
    ttype_list = [x for x in IC_data.keys()]
    subclass_list = unique_element([x.split('_')[0] for x in IC_data.keys()])

    map_exc_ttype = compute_default_exc_map(path_to_BBP_m_type_list, ttype_list)
    
    with open('./output/map_exc_ttype.json', 'w') as fp:
        json.dump(map_exc_ttype, fp)
    
    compatible_profiles = {}
    for m_type in map_exc_ttype.keys():
        compatible_profiles[m_type] = {}
        for ttype in map_exc_ttype[m_type]:
            df = generate_panda(IC_data, ttype)
            compatible_profiles[m_type][ttype] = {'values': df.values.tolist(), 'index': df.index.tolist(), 'columns': df.columns.tolist()}
            
    with open('./output/exc_compatible_profiles.json', 'w') as fp:
        json.dump(compatible_profiles, fp)
    
    map_inh_ttype_binary = make_inh_map_binary(path_to_inh_map_L1, path_to_inh_map_L26)
    
    msk_clstr_Lamp5 = [('Lamp5' in c) for c in medians.columns]
    msk_clstr_Sncg = [('Sncg' in c) for c in medians.columns]
    msk_clstr_Serpinf = [('Serpinf' in c) for c in medians.columns]
    msk_clstr_Vip = [('Vip' in c) for c in medians.columns]
    msk_clstr_Sst = [('Sst' in c) for c in medians.columns]
    msk_clstr_Pvalb = [('Pvalb' in c) for c in medians.columns]

    dict_mkr_ttype = {
        'Lamp5': medians.columns[msk_clstr_Lamp5],
        'Pvalb': medians.columns[msk_clstr_Pvalb],
        'Serpinf1': medians.columns[msk_clstr_Serpinf],
        'Sncg': medians.columns[msk_clstr_Sncg],
        'Sst': medians.columns[msk_clstr_Sst],
        'Vip': medians.columns[msk_clstr_Vip],
    }
    
    compatible_inh_profiles = {}

    for me_type in map_inh_ttype_binary.columns:
        compatible_inh_profiles[me_type] = {}
        compatible_markers = map_inh_ttype_binary.index[np.where(map_inh_ttype_binary[me_type] == 1)[0]]

        compatible_ttype = []
        for mkr in compatible_markers:
            compatible_ttype += dict_mkr_ttype[mkr].tolist()
        compatible_ttype

        for ttype in compatible_ttype:
            df = generate_panda(IC_data, ttype)
            compatible_inh_profiles[me_type][ttype] = {'values': df.values.tolist(), 'index': df.index.tolist(), 'columns': df.columns.tolist()}
    
    with open('./output/inh_compatible_profiles.json', 'w') as fp:
        json.dump(compatible_inh_profiles, fp)