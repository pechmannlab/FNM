import sys, os
import pandas as pd
import numpy as np
import networkx as nx



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_gi_cdt(fileIN):
    """
    load and parse Toronto genetic interaction (GI) 
    data from cdt tree view file
    """
    
    data = pd.read_csv(fileIN, header=2, index_col=False, sep='\t', low_memory=False)
    data = data.drop([0,1,2], axis=0)
    
    #formatting of cdt file
    current_cols = data.columns
    current_sel = current_cols[[0,1,3,4,5]]

    data = data.drop(current_sel, axis=1)
    data = data.set_index(current_cols[2])
    data.index = data.index.rename('')

    return data


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def big_matrix(GImats):
    """
    combine data to full matrix
    scores are averages of all pairwise interactions
    """

    orf_list = []
    GIindmats = []
    GIdict = {}
    for dataset in range(len(GImats)):
        current_data = GImats[dataset]
        current_rows = list(current_data.index)
        current_cols = list(current_data.columns)
        current_data = np.array(current_data, dtype=float)

        current_rows_orf = [orf.split('.')[0] for orf in current_rows]
        current_cols_orf = [orf.split('.')[0] for orf in current_cols]
        GIindmats.append( (current_rows_orf, current_cols_orf) )
        
        for i in current_rows_orf:
            orf_list.append(i)
        for i in current_cols_orf:
            orf_list.append(i)

        current_orfs = list(set( current_rows_orf + current_cols_orf ))
        current_N = len(current_orfs)

        for i in range(current_N):
            orf_i = current_orfs[i]
            print(i, current_N, orf_i)
            for j in range(i, current_N):
                orf_j = current_orfs[j]
                orfs_sorted = sorted([orf_i, orf_j])
                current_key = (orfs_sorted[0], orfs_sorted[1])

                if orf_i in current_rows_orf and orf_j in current_cols_orf:
                    current_orf_i = np.where(np.array(current_rows_orf)==orf_i)[0]
                    current_orf_j = np.where(np.array(current_cols_orf)==orf_j)[0]
                    current_scores = current_data[current_orf_i,:]
                    current_scores = current_scores[:, current_orf_j]
                    current_scores = list( current_scores.flatten() )      
                    if current_key not in GIdict.keys():
                        GIdict[current_key] = current_scores
                    elif current_key in GIdict.keys():
                        current_list = GIdict[current_key]
                        current_list = current_list + current_scores
                        GIdict[current_key] = current_list

                if orf_i in current_cols_orf and orf_j in current_rows_orf:
                    current_orf_i = np.where(np.array(current_cols_orf)==orf_i)[0]
                    current_orf_j = np.where(np.array(current_rows_orf)==orf_j)[0]     
                    current_scores = current_data[current_orf_j,:]
                    current_scores = current_scores[:, current_orf_i]
                    current_scores = list( current_scores.flatten() ) 
                    if current_key not in GIdict.keys():
                        GIdict[current_key] = current_scores
                    elif current_key in GIdict.keys():
                        current_list = GIdict[current_key]
                        current_list = current_list + current_scores
                        GIdict[current_key] = current_list
  
    orfs = list(set(orf_list))
    N = len(orfs)
    print(N)
    GI = np.zeros(( N, N ))
 
    for i in range(N):
        orf_i = orfs[i]
        for j in range(i, N):
            orf_j = orfs[j]
            orfs_sorted = sorted([orf_i, orf_j])
            current_key = (orfs_sorted[0], orfs_sorted[1])

            if current_key in GIdict.keys():
                current_interactions = np.array( GIdict[current_key] ) 
                current_interactions = current_interactions[~np.isnan(current_interactions)]
                if len(current_interactions) > 0:
                    GI[i,j] = GI[j,i] = np.around(np.nanmean(np.array(current_interactions)), 3)
                else:
                    GI[i,j] = GI[j,i] = np.nan

    #np.savetxt("../data/processed/GI.txt", GI, fmt='%.3f' )
    GIdf = pd.DataFrame(data=GI, index=orfs, columns=orfs)
    GIdf.to_csv("../data/processed/GI.txt", header=True, index=True, sep='\t')

    return GIdf


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def sig_GInetwork(gi):
    """
    filter top/bottom 5% as signification interactions
    save as networkx file
    """

    #gi = pd.read_csv("../data/processed/GI.txt", header=0, index_col=0, sep='\t')
    #print(gi)

    data = np.array(gi)
    data = data[~np.isnan(data)]
    #print( np.sum(data > 0.051) / len(data), np.sum(data < -0.064) / len(data) )
    #print(np.percentile(data, 5), np.percentile(data, 95) )

    orfs = list(gi.index)
    N = len(orfs)
    gi_data = np.array(gi)
    Ggi = nx.Graph()
    for i in range(N):
        orf_i = orfs[i]
        for j in range(i+1,N):
            orf_j = orfs[j]

            score = gi_data[i,j] 
            if ~np.isnan(score):
                if score > 0.051:       # top 5%
                    Ggi.add_edge(orf_i, orf_j, weight=score)
                elif score < -0.064:    # bottom 5%
                    Ggi.add_edge(orf_i, orf_j, weight=score)

    nx.write_gpickle(Ggi, "../data/processed/GI.nx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_ppi(BIOGRID_file):
    """
    load BIOGRID yeast PPI
    write to smaller 2-columns file and nx.Graph
    """

    ppi = pd.read_csv(BIOGRID_file, header=0, index_col=False, sep='\t', low_memory=False)
    ppi = ppi[ppi['Experimental System Type'] == 'physical']
    column_selection = ppi.columns[[5,6]]
    PPI = ppi[column_selection]
    PPI.columns=['ORF1', 'ORF2'] 
    PPI.to_csv("../data/processed/yeast_ppi.txt", header=True, index=False, sep='\t')

    Gppi = nx.Graph()
    for i in range(len(PPI)):
        orf1, orf2 = list(PPI.loc[i]) 
        print(orf1, orf2)
        orf1 = PPI.loc[i]['ORF1']
        orf2 = PPI.loc[i]['ORF2']
        Gppi.add_edge(orf1, orf2)
    nx.write_gpickle(Gppi, '../data/processed/yeast_ppi.nx')








if __name__ == '__main__':

    # Generate smaller pre-processed input files from publicly available raw data

    # parse PPI network
    parse_ppi("../data/yeast/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.185.tab3.txt")


    # parse GI network
    gi_ExE = load_gi_cdt("../data/yeast/GI/SGA_ExE_clustered.cdt")
    gi_NxN = load_gi_cdt("../data/yeast//GI/SGA_NxN_clustered.cdt")
    gi_ExN = load_gi_cdt("../data/yeast/GI/SGA_ExN_NxE_clustered.cdt")
    gi_DAmP = load_gi_cdt("../data/yeast/GI/SGA_DAmP_clustered.cdt")
    GI = big_matrix([gi_ExE, gi_NxN, gi_ExN, gi_DAmP])
    sig_GInetwork(GI)