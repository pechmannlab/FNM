import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import itertools as it
import time
import multiprocessing as mp
import glob
import subprocess
import shutil
import random

from scipy.cluster.hierarchy import average, dendrogram, fcluster
from scipy.spatial.distance import squareform


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def orf2name():
    """
    Utility function to convert between ORFs and gene names
    """

    f_names = '../data/yeast/orf_coding_names.txt'
    yeast_names = {}
    yeast_names_rev = {}
    for line in open(f_names, 'r'):
        current_line = line.split()
        orf = current_line[0]
        if len(current_line) > 1:
            name = current_line[1]
            yeast_names[orf] = name
            yeast_names_rev[name] = orf
        else:
            yeast_names[orf] = orf
            yeast_names_rev[orf] = orf

    return yeast_names, yeast_names_rev
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter_degree_ppi(Gppi, max_degree):
    """
    Utility function to filter networks based on maximum degree in PPI
    removes nodes and associated edges with degree > max_degree
    """

    list_keep = []
    for entry in list(Gppi.degree()):
        if entry[1] <= max_degree:
            list_keep.append(entry[0])
    Gppi = nx.subgraph(Gppi, list_keep)
   
    return Gppi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_noclust(orflist):
    """
    Function to write motifs not part of clusters to dataframe
    """

    motifDF = pd.DataFrame(columns=['ORF'+str(i+1) for i in range(6)])

    for k in [3,4,5,6]:
        motif_data = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        current_columns = motif_data.columns
        current_columns = list(current_columns[1:-1])
      
        for i in range( len(motif_data) ):
            current_motif = list(motif_data.loc[i][current_columns])
            current_idx = np.zeros(( k )) * np.nan
            current_orflist1_filter = np.zeros(( k ), dtype=bool)
          
            for mx, current_orf in enumerate(current_motif) :
                if current_orf in orflist:
                    current_orflist1_filter[mx] = True

            if np.all( current_orflist1_filter ):
                current_noclust = current_motif
                while len(current_noclust) < 6:
                    current_noclust.append('none')
                motifDF.loc[len(motifDF)] = current_noclust
 
    motifDF.to_csv("../data/figures/figure3/motifs_noclust.txt", header=True, index=False, sep='\t')

    return motifDF




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_suppressors(category):
    """
    Check for enrichment of suppressor interactions as proxy 
    for 'functionally related', i.e. functional motifs
    """

    # load Boone/Andrews data on yeast suppressors
    supp = pd.read_csv("../data/yeast/yeast_suppressors.txt", header=0, index_col=False, sep='\t')
    supp = supp[['Query ORF', 'Suppressor ORF']]
    supp = supp.drop_duplicates()

    supp_nodes = sorted(list(set( list(supp['Query ORF']) + list(supp['Suppressor ORF']) )))
    supp_mat = np.zeros(( len(supp), len(supp_nodes) ))

    for i in range(len(supp)):
        current_query = supp.iloc[i]['Query ORF']
        current_supp  = supp.iloc[i]['Suppressor ORF']
        supp_mat[i, supp_nodes.index(current_query)] = 1
        supp_mat[i, supp_nodes.index(current_supp)] = 1

    # load Boone/Andrews data on bypass suppressors (too few to really analyze)
    bypa = pd.read_csv("../data/yeast/yeast_bypass.txt", header=0, index_col=False, sep='\t')
    bypa = bypa[['QueryORF', 'SuppressorORF']]
    bypa = bypa.drop_duplicates()
    
    bypa_nodes = list(bypa['QueryORF'])
    for orf in bypa['SuppressorORF']:
        if orf != 'unknown' and ';' not in orf:
            bypa_nodes.append(orf)
        elif ';' in orf:
            orf = orf.split(';')
            bypa_nodes.append(orf[0])
            bypa_nodes.append(orf[1])
    bypa_nodes = sorted(bypa_nodes)
    bypa_mat = np.zeros(( len(bypa), len(bypa_nodes) ))

    for i in range(len(bypa)):
        current_query = bypa.iloc[i]['QueryORF']
        current_bypa  = bypa.iloc[i]['SuppressorORF']
        if current_bypa != 'unknown' and ';' not in current_bypa:
            bypa_mat[i, bypa_nodes.index(current_query)] = 1
            bypa_mat[i, bypa_nodes.index(current_bypa)] = 1
        elif ';' in current_bypa:
            current_bypa = current_bypa.split(';')
            bypa_mat[i, bypa_nodes.index(current_query)] = 1
            bypa_mat[i, bypa_nodes.index(current_bypa[0])] = 1
            bypa_mat[i, bypa_nodes.index(current_bypa[1])] = 1

    idx = 0
    idx_supp = 0
    idx_bypa = 0
    for k in [3,4,5,6]:
        if category == "FNM":
            motif_data = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        elif category == "ALL":
            motif_data = pd.read_csv('../data/motif/result_d50_all_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')

        current_columns = motif_data.columns
        current_columns = list(current_columns[1:-1])    
        for i in range( len(motif_data) ):
            current_motif = list(motif_data.loc[i][current_columns])
            current_idx = []
            current_idx2 = []
            idx += 1          
            for mx, current_orf in enumerate(current_motif) :
                if current_orf in supp_nodes:
                    current_idx.append( supp_nodes.index(current_orf) )
                if current_orf in bypa_nodes:
                    current_idx2.append( bypa_nodes.index(current_orf) )
            current_idx = np.array(current_idx, dtype=int)
            current_idx2 = np.array(current_idx2, dtype=int)

            current_mat = supp_mat[:,current_idx]
            current_hits = current_mat[ np.sum(current_mat, 1) == 2, : ]
            if current_hits.shape[0] > 0:
                idx_supp += current_hits.shape[0]

            current_bypa = bypa_mat[:,current_idx2]
            sel_bypa = (np.sum(current_bypa, 1) == np.sum(bypa_mat,1)) * (np.sum(current_bypa,1) > 0)
            current_hits2 = current_bypa[ sel_bypa, : ]
            if current_hits2.shape[0] > 0:
                idx_bypa += current_hits2.shape[0]

        print(k, 'done')

    return idx, idx_supp, idx_bypa




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def expression_divergence():
    """
    Std of expression profile as proxy for regulation (not really informative)
    """

    gasch = pd.read_csv("../data/yeast/gasch_stress.txt", header=0, index_col=False, sep='\t')
   
    columns_gasch = gasch.columns
    orfs_gasch = list(gasch['UID'])
    data = np.array(gasch[columns_gasch[1:]])
    data_nr, data_nc = np.shape(data)

    reg_std = np.nanstd(data, 1)
    reg_abs = np.nansum(np.abs(data),1) / np.sum(np.isnan(data)==False, 1)

    result = pd.DataFrame({'ORF': orfs_gasch, 'SD': list( np.around(reg_std,3) ), 'ABS': list( np.around(reg_abs,3) )})
    result.to_csv("../data/figures/figure3/transcript_reg.txt", header=True, index=False, sep='\t')

    return result





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_gencat(category):
    """
    test for enrichment of essential genes and complex subunits in motifs
    """


    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')
    list_subunits = list(set(complexes['ORF']))

    
    yeast_essential = []
    for line in open("../data/yeast/yeast_essential.txt", 'r'):
        current_line = line.split()
        current_orf = current_line[0]
        yeast_essential.append(current_orf)

    idx = 0
    idx_ess = 0
    idx_sub = 0

    result_ess = np.zeros(( 7 ), dtype=int)
    result_sub = np.zeros(( 7 ), dtype=int)

    result_scores_ess = np.zeros(( 11000000 ))
    result_scores_sub = np.zeros(( 11000000 ))

    df_ess = pd.DataFrame(columns=['orf1', 'orf2', 'orf3', 'orf4', 'orf5', 'orf6'])
    df_sub = pd.DataFrame(columns=['orf1', 'orf2', 'orf3', 'orf4', 'orf5', 'orf6'])


    for k in [3,4,5,6]:
        if category == "FNM":
            motif_data = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        elif category == "ALL":
            motif_data = pd.read_csv('../data/motif/result_d50_all_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')

        current_columns = motif_data.columns
        current_columns = list(current_columns[1:-1])
      
        for i in range( len(motif_data) ):
            current_motif = list(motif_data.loc[i][current_columns])
            current_idx = []
            current_idx2 = []
            current_ess = 0
            current_sub = 0
      
            for mx, current_orf in enumerate(current_motif) :
                if current_orf in yeast_essential:
                    current_ess += 1
                if current_orf in list_subunits:
                    current_sub += 1

            if current_ess > 0:
                idx_ess += 1             
            if current_sub > 0:
                idx_sub += 1

            if current_ess == 6:
                df_ess.loc[len(df_ess)] = current_motif
            if current_sub == 6:
                df_sub.loc[len(df_sub)] = current_motif

            result_ess[current_ess] += 1
            result_sub[current_sub] += 1

            result_scores_ess[idx] = current_ess / len(current_motif)
            result_scores_sub[idx] = current_sub / len(current_motif)
            idx += 1

        print(k, 'done')

    result_scores_ess = result_scores_ess[:idx]
    result_scores_sub = result_scores_sub[:idx]

    np.savetxt("../data/figures/figure2/"+category+".scores.ess.txt", np.around(result_scores_ess,3), fmt='%.3f' )
    np.savetxt("../data/figures/figure2/"+category+".scores.sub.txt", np.around(result_scores_sub,3), fmt='%.3f' )

    df_ess.to_csv("../data/figures/figure2/"+category+".full6.ess.txt", header=True, index=False, sep='\t')
    df_sub.to_csv("../data/figures/figure2/"+category+".full6.sub.txt", header=True, index=False, sep='\t')


    return result_ess, result_sub




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def check_motif_noclust(inputDF):
    """
    check for essential genes and complex subunits in no-cluster motifs
    """

    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')
    list_subunits = list(set(complexes['ORF']))

    yeast_essential = []
    for line in open("../data/yeast/yeast_essential.txt", 'r'):
        current_line = line.split()
        current_orf = current_line[0]
        yeast_essential.append(current_orf)


    supp = pd.read_csv("../data/yeast/yeast_suppressors.txt", header=0, index_col=False, sep='\t')
    supp = supp[['Query ORF', 'Suppressor ORF']]
    supp = supp.drop_duplicates()

    supp_nodes = sorted(list(set( list(supp['Query ORF']) + list(supp['Suppressor ORF']) )))
    supp_mat = np.zeros(( len(supp), len(supp_nodes) ))

    for i in range(len(supp)):
        current_query = supp.iloc[i]['Query ORF']
        current_supp  = supp.iloc[i]['Suppressor ORF']
        supp_mat[i, supp_nodes.index(current_query)] = 1
        supp_mat[i, supp_nodes.index(current_supp)] = 1


    bypa = pd.read_csv("../data/yeast/yeast_bypass.txt", header=0, index_col=False, sep='\t')
    bypa = bypa[['QueryORF', 'SuppressorORF']]
    bypa = bypa.drop_duplicates()

    bypa_nodes = list(bypa['QueryORF'])
    for orf in bypa['SuppressorORF']:
        if orf != 'unknown' and ';' not in orf:
            bypa_nodes.append(orf)
        elif ';' in orf:
            orf = orf.split(';')
            bypa_nodes.append(orf[0])
            bypa_nodes.append(orf[1])
    bypa_nodes = sorted(bypa_nodes)
    bypa_mat = np.zeros(( len(bypa), len(bypa_nodes) ))

    for i in range(len(bypa)):
        current_query = bypa.iloc[i]['QueryORF']
        current_bypa  = bypa.iloc[i]['SuppressorORF']
        if current_bypa != 'unknown' and ';' not in current_bypa:
            bypa_mat[i, bypa_nodes.index(current_query)] = 1
            bypa_mat[i, bypa_nodes.index(current_bypa)] = 1
        elif ';' in current_bypa:
            current_bypa = current_bypa.split(';')
            bypa_mat[i, bypa_nodes.index(current_query)] = 1
            bypa_mat[i, bypa_nodes.index(current_bypa[0])] = 1
            bypa_mat[i, bypa_nodes.index(current_bypa[1])] = 1

    idx = 0
    idx_supp = 0
    idx_bypa = 0
    idx_ess = 0
    idx_sub = 0


    for i in range( len(inputDF) ):
        current_motif = np.array(inputDF.loc[i])
        current_motif = list( current_motif[current_motif!='none'] )
        current_idx = []
        current_idx2 = []
        current_ess = 0
        current_sub = 0
            
        for mx, current_orf in enumerate(current_motif) :
            
            if current_orf in yeast_essential:
                current_ess += 1

            if current_orf in list_subunits:
                current_sub += 1

            if current_orf in supp_nodes:
                current_idx.append( supp_nodes.index(current_orf) )

            if current_orf in bypa_nodes:
                current_idx2.append( bypa_nodes.index(current_orf) )

        if current_ess > 0:
            idx_ess += 1

        if current_sub > 0:
            idx_sub += 1

        current_idx = np.array(current_idx, dtype=int)
        current_idx2 = np.array(current_idx2, dtype=int)

        current_mat = supp_mat[:,current_idx]
        current_hits = current_mat[ np.sum(current_mat, 1) == 2, : ]
        if current_hits.shape[0] > 0:
            idx_supp += current_hits.shape[0]


        current_bypa = bypa_mat[:,current_idx2]
        sel_bypa = (np.sum(current_bypa, 1) == np.sum(bypa_mat,1)) * (np.sum(current_bypa,1) > 0)
        current_hits2 = current_bypa[ sel_bypa, : ]
        if current_hits2.shape[0] > 0:
            idx_bypa += current_hits2.shape[0]
     
        idx += 1

    f = open("../data/figures/figure2/noclust_gencat.txt", 'w')
    f.write("total\t"+str(idx)+"\n")
    f.write("suppressors\t"+str(idx_supp)+'\n')
    f.write("bypa\t"+str(idx_bypa)+"\n")
    f.write("essential\t"+str(idx_ess)+"\n")
    f.write("subunit\t"+str(idx_sub)+"\n")
    f.close()

    #print(idx, idx_supp, idx_bypa, idx_ess, idx_sub)

    return None


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def complex_network(orflist1, theta_subu):
    """
    Statistics on number of complexes and interactions as function of theta_subu
    """

    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')

    min_subunit = theta_subu
    complex_dict = {}
    for i in  list(set(list(complexes['ORF']))):
        current_entry = complexes[complexes['ORF']==i]
        if len(current_entry) == 1:
            current_complex = current_entry['Complex'].item()
        #elif len(current_entry) > 1:
        #    current_complex = 'multi'

            current_subu = complexes[complexes['Complex']==current_complex]
            if len(current_subu) >= min_subunit:
                complex_dict[i] = current_complex 

    N_complex = 0
    for i in list(set(list(complexes['Complex']))):
        current_entry = complexes[complexes['Complex']==i]
        if len(current_entry) >= min_subunit:
            N_complex += 1

    counter = 0
    counter_compcnx = 0
    for k in [3,4,5,6]:
        motif_data = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        current_columns = motif_data.columns
        current_columns = list(current_columns[1:-1])
      
        for i in range( len(motif_data) ):
            current_motif = list(motif_data.loc[i][current_columns])
            current_idx = np.zeros(( k )) * np.nan
            current_orflist1_filter = np.zeros(( k ), dtype=bool)

            counter += 1
            current_complex = []
            for mx, current_orf in enumerate(current_motif) :
                current_membership = complex_dict.get(current_orf, 'none')

                if current_orf in orflist1:
                    current_orflist1_filter[mx] = True

                if current_membership != 'none': # and current_membership != 'multi':
                    current_complex.append(current_membership)

            current_complex = list(set(current_complex))
            if np.all(current_orflist1_filter) and len(current_complex) > 1:
                counter_compcnx += 1
           
    print(theta_subu, N_complex, counter_compcnx)

    return theta_subu, N_complex, counter_compcnx






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_complex_net(inputDF, theta_complex, limit_subset=True):
    """
    Function to generate a protein complex - protein complex interaction
    network from motif data
    inputDF: dataframe with motifs
    theta_complex: minimum size of complex for consideration
    limit_subset: logical whether to filter based on theta_complex
    """

    yeast_names, yeast_names_rev = orf2name()
    
    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')
    G_ppi = nx.read_gpickle("../data/processed/yeast_ppi.nx")
    G = nx.Graph()

    transreg = expression_divergence()

    complex_names = {}
    names_short = pd.read_csv("../data/yeast/MIPS_shortnames.txt", header=0, index_col=False, sep=' ')
    for i in range(len(names_short)):
        current_df = names_short.iloc[i]
        current_long = str( current_df['Long'] )
        current_short = str( current_df['Short'] )
        complex_names[current_long] = current_short
   
    def get_complex_name(orf, cmplx, ynames):
        current_name = orf
        if orf in list(cmplx['ORF']):
            current_entry = cmplx[cmplx['ORF']==orf]
            if len(current_entry) == 1:
                current_name = current_entry['Complex'].item()
            elif len(current_entry) > 1:
            	current_name = list(current_entry['Complex'])[0]
        #else:
        #    current_name = ynames.get(orf, orf)

        return current_name

    for i in range( len(inputDF) ):
        current_motif = np.array(inputDF.loc[i])
        current_motif = list( current_motif[current_motif!='none'] )

       
        G_motif = nx.subgraph(G_ppi, current_motif)
        current_edges = G_motif.edges()
        for edge in current_edges:
            n1 = edge[0]
            n2 = edge[1]

            c1 = get_complex_name(n1, complexes, yeast_names)
            c2 = get_complex_name(n2, complexes, yeast_names)
            G.add_edge(c1, c2)


    if limit_subset:
        sel_complexes = np.zeros(( len(complexes) ), dtype=bool)
        for i in list(set(list(complexes['Complex']))):
            if len(complexes[complexes['Complex']==i]) >= theta_complex:
                sel_complexes[ complexes['Complex']==i] = True
        complexes2 = complexes.loc[sel_complexes]

        list_complexes_filtered = list(set(list(complexes2['Complex'])))
    
        list_subset = []
        for i in list(G.nodes()):
            if i in list(yeast_names.keys()):
                list_subset.append(complex_names.get(i, i))
            elif i in list_complexes_filtered:
                list_subset.append(i)

        G2 = nx.subgraph(G, list_subset)
        G0  = max(nx.connected_components(G2), key=len)
        G = nx.subgraph(G, G0)


    list_nodes = list(G.nodes())
    list_nodes = sorted(list_nodes)

    print(len(list_nodes))

    gmat = nx.to_numpy_array(G, nodelist=list_nodes)
    np.savetxt("../data/figures/figure4/network_"+str(theta_complex)+".txt", gmat, fmt='%i')


    list_anchors = []
    for i in list_nodes:
        if i in list(complexes['Complex']):
            list_anchors.append(i)

    #betw = nx.betweenness_centrality(G)				# all pathws
    betw = betweenness_anchornodes(G, list_anchors)		# only paths between complexes

    node_type = []
    node_betw = np.zeros(( len(list_nodes) ))
    node_names = []
    node_orfs = []
    node_transreg = np.zeros(( len(list_nodes) ))
    for ix, i in enumerate(list_nodes):
        node_betw[ix] = betw[ix]

        if i in list(complexes['Complex']):
            node_type.append('C')
            node_names.append( complex_names.get(i, i))
            node_orfs.append( complex_names.get(i, i))
            node_transreg[ix] = np.nan
        else:
            node_type.append('S')
            node_names.append(yeast_names.get(i, i))
            node_orfs.append(i)
            if i in list(transreg['ORF']):
                current_transreg = transreg[transreg['ORF']==i]['SD']
                node_transreg[ix] = float(current_transreg)
            else:
            	node_transreg[ix] = np.nan 

    graphDF = pd.DataFrame({'orf': list(node_orfs), 'name':node_names, 'type':list(node_type), 'betw': list(np.around(node_betw,5)), 'regulation': list(node_transreg)})
    graphDF.to_csv("../data/figures/figure4/network_annotation_"+str(theta_complex)+".txt", header=True, index=False, sep='\t', na_rep='NA')


    return None

    
    
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def betweenness_anchornodes(G, anchors):
    """
    computes betweenness centrality only for
    pathways between anchornodes
    anchors: input list of anchor nodes
    """

    L = len(anchors)

    list_nodes = list(G.nodes())

    sigma = np.zeros(( len(list_nodes) ))
    n_sigma = 0.

    for i in range(L):
        anchor_source = anchors[i]
        for j in range(i+1, L):
            anchor_target = anchors[j]
            current_shortest_paths = nx.all_shortest_paths(G, source=anchor_source, target=anchor_target)

            sigma_s_t = float(0)
            sigma_s_t_v = np.zeros(( len(list_nodes) ), dtype=int)

            for path in current_shortest_paths:
                sigma_s_t += 1
                if len(path) > 2:
                    current_middle_nodes = path[1:-1]
                    for node in current_middle_nodes:
                        sigma_s_t_v[ list_nodes.index(node) ] += 1

            current_sigma = ( sigma_s_t_v / sigma_s_t )
            if np.any(current_sigma > 1):
                print(current_sigma[current_sigma > 0])
            sigma += current_sigma
            n_sigma += 1

    betweenness = sigma / n_sigma

    return betweenness


  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def reference_complex_net(G_ppi, pre, theta_complex, limit_subset=True):
    """
    Complex interaction network based on PPI not FNMs
    """
   
    yeast_names, yeast_names_rev = orf2name()  
    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')
    
    G = nx.Graph()

    complex_names = {}
    names_short = pd.read_csv("../data/yeast/MIPS_shortnames.txt", header=0, index_col=False, sep=' ')
    for i in range(len(names_short)):
        current_df = names_short.iloc[i]
        current_long = str( current_df['Long'] )
        current_short = str( current_df['Short'] )
        complex_names[current_long] = current_short
    
    def get_complex_name(orf, cmplx, ynames):
        current_name = orf
        if orf in list(cmplx['ORF']):
            current_entry = cmplx[cmplx['ORF']==orf]
            if len(current_entry) == 1:
                current_name = current_entry['Complex'].item()
            elif len(current_entry) > 1:
            	current_name = list(current_entry['Complex'])[0]
        #else:
        #    current_name = ynames.get(orf, orf)
        return current_name

    for edge in list(G_ppi.edges()):
        n1 = edge[0]
        n2 = edge[1]
        c1 = get_complex_name(n1, complexes, yeast_names)
        c2 = get_complex_name(n2, complexes, yeast_names)       
        G.add_edge(c1, c2)

    if limit_subset:
        sel_complexes = np.zeros(( len(complexes) ), dtype=bool)
        for i in list(set(list(complexes['Complex']))):
            if len(complexes[complexes['Complex']==i]) >= theta_complex:
                sel_complexes[ complexes['Complex']==i] = True
        complexes2 = complexes.loc[sel_complexes]

        list_complexes_filtered = list(set(list(complexes2['Complex'])))   
        list_subset = []
        for i in list(G.nodes()):
            if i in list(yeast_names.keys()):
                list_subset.append(complex_names.get(i, i))
            elif i in list_complexes_filtered:
                list_subset.append(i)

        G2 = nx.subgraph(G, list_subset)
        G0  = max(nx.connected_components(G2), key=len)
        G = nx.subgraph(G, G0)

    list_nodes = list(G.nodes())
    list_nodes = sorted(list_nodes)
    gmat = nx.to_numpy_array(G, nodelist=list_nodes)
    np.savetxt("../data/figures/figure4/refnet_"+str(pre)+"_"+str(theta_complex)+".txt", gmat, fmt='%i')

    list_anchors = []
    for i in list_nodes:
        if i in list(complexes['Complex']):
            list_anchors.append(i)

    #betw = nx.betweenness_centrality(G)				# all pathws
    betw = betweenness_anchornodes(G, list_anchors)		# only paths between complexes

    node_type = []
    node_betw = np.zeros(( len(list_nodes) ))
    node_names = []
    node_orfs = []
    for ix, i in enumerate(list_nodes):
        node_betw[ix] = betw[ix]
        if i in list(complexes['Complex']):
            node_type.append('C')
            node_names.append( complex_names.get(i, i))
            node_orfs.append( complex_names.get(i, i))
        else:
            node_type.append('S')
            node_names.append(yeast_names.get(i, i))
            node_orfs.append(i)
  
    graphDF = pd.DataFrame({'orf': node_orfs, 'name':node_names, 'type':list(node_type), 'betw': list(np.around(node_betw,5))})
    graphDF.to_csv("../data/figures/figure4/refnet_annotation_"+str(pre)+"_"+str(theta_complex)+".txt", header=True, index=False, sep='\t')

    return graphDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def full6_gi(motifDF, G, GI):
    """
    Network of motifs that have complex subunits as nodes
    PPI graph: all edges from full6 motifs, edge weight is number of occurances
    GI graph: all GI edges between motif nodes, edge weight is average of per-complex interaction
    """

    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')
    shortnames = pd.read_csv("../data/yeast/MIPS_shortnames.txt", header=0, index_col=False, sep=' ')

    list_orf = []
    G_complex = nx.Graph()
    G_complex_gi = nx.Graph()
    G_complex_gi_pos = nx.Graph()
    G_complex_gi_neg = nx.Graph()
    current_gi_dict = {}

    for i in range(len(motifDF)):
        current_motif = list( motifDF.iloc[i] )

        current_g = nx.subgraph(G, current_motif)
        for edge in list( current_g.edges() ):
            complex1 = list(complexes[complexes['ORF']==edge[0]]['Complex'])
            complex2 = list(complexes[complexes['ORF']==edge[1]]['Complex'])
            if len(complex1) == 1 and len(complex2) == 1:
                node1 = complex1[0]
                node2 = complex2[0]
                if node1 != node2:
                    if G_complex.has_edge(node1, node2):
                        G_complex.edges[node1, node2]['weight'] += 1
                    else:
                        G_complex.add_edge( node1, node2, weight=1)

        current_gi = nx.subgraph(GI, current_motif)
        for gi_edge in list( current_gi.edges() ):
            complex1 = list(complexes[complexes['ORF']==gi_edge[0]]['Complex'])
            complex2 = list(complexes[complexes['ORF']==gi_edge[1]]['Complex'])
            if len(complex1) == 1 and len(complex2) == 1:
                current_sorted = sorted( [complex1[0], complex2[0]] )
                current_key = current_sorted[0]+"#"+current_sorted[1]
                current_weight = GI.edges[gi_edge]['weight']
                if current_key not in list(current_gi_dict.keys()):
                    current_gi_dict[current_key] = [current_weight]
                else:
                    current_weightlist = current_gi_dict[current_key]
                    current_weightlist.append(current_weight)
                    current_gi_dict[current_key] = current_weightlist

    G_complex = G_complex.subgraph( max(nx.connected_components(G_complex), key=len) )
    list_nodes = sorted(list(G_complex.nodes()))

    for key in list(current_gi_dict.keys()):
        current_nodes = key.split("#")
        current_weight = np.mean(np.array(current_gi_dict[key]))
        current_n_pos = np.sum(np.array(current_gi_dict[key]) > 0)
        current_n_neg = np.sum(np.array(current_gi_dict[key]) < 0)
        if current_nodes[0] != current_nodes[1]:
            G_complex_gi.add_edge( current_nodes[0], current_nodes[1], weight=current_weight)
            G_complex_gi_pos.add_edge( current_nodes[0], current_nodes[1], weight=current_n_pos )
            G_complex_gi_neg.add_edge( current_nodes[0], current_nodes[1], weight=current_n_neg )


    gmat_complex = nx.to_numpy_array(G_complex, nodelist=list_nodes)
    gmat_gi = nx.to_numpy_array(G_complex_gi, nodelist=list_nodes)
    gmat_gi_pos = nx.to_numpy_array(G_complex_gi_pos, nodelist=list_nodes)
    gmat_gi_neg = nx.to_numpy_array(G_complex_gi_neg, nodelist=list_nodes)
    np.savetxt("../data/figures/figure3/complex_network_ppi.txt", np.array(gmat_complex, dtype=int))
    np.savetxt("../data/figures/figure3/complex_network_gi.txt", gmat_gi, fmt='%.4f')
    np.savetxt("../data/figures/figure3/complex_network_gi_pos.txt", gmat_gi_pos, fmt='%.4f')
    np.savetxt("../data/figures/figure3/complex_network_gi_neg.txt", gmat_gi_neg, fmt='%.4f')


    #---- heuristic hierarchical cluster of nodes for the edge bundeling ---- 
    gmat2 = np.copy(gmat_complex)	# copy without edge weights
    gmat2[gmat2>0] = 1
    deg = np.sum(gmat2, 1)

    group = deg > 6	 	# arbitrary threshold used for clustering used in edge bundeling, just optics
    intab = np.zeros(( np.shape(gmat2)[1] )) * np.nan
    for ix, i in enumerate(np.where(group)[0]):
        intab[i] = ix
    
    while np.sum(np.isnan(intab)) > 0:
        for i in range(np.shape(gmat2)[1]):
            if i not in np.where(np.isnan(intab) == False)[0]: #group[i]:
                list_id = []
                list_score = []
                for j in np.where(np.isnan(intab)==False)[0]:
                    if gmat2[i,j] == 1:
                        current_overlap = (gmat2[i,:]==1)*(gmat2[j,:]==1)
                        list_id.append(j)
                        list_score.append( np.sum(current_overlap))
                if len(list_score) > 0:
                    best_j = np.argmax(np.array(list_score))
                    best_group = intab[list_id[best_j]]
                    intab[i] = best_group
   

    # write pseudo hierarchical clustering to data frame
    Rhierarchy = pd.DataFrame(columns=['from', 'to'])
    for ix, i in enumerate( list(set(intab)) ):
        Rhierarchy.loc[len(Rhierarchy)] = ("origin", "group_"+str(ix))
    for ix, i in enumerate( list(set(intab)) ):
        current_group = np.where(intab == i)[0]
        for j in current_group:
            Rhierarchy.loc[len(Rhierarchy)] = ("group_"+str(ix), "subgroup_"+str(j))
    Rhierarchy.to_csv("../data/figures/figure3/network_hierarchy.txt", header=True, index=False, sep='\t')


    # graphs sorted the same way at cost for another second of compute
    L = np.shape(gmat_complex)[1]
    connect = pd.DataFrame(columns=['from', 'to', 'weight'])
    conn_gi = pd.DataFrame(columns=['from', 'to', 'weight'])

    for i in range(L):
        for j in range(i+1, L):
            if gmat_complex[i,j] > 0:
                connect.loc[len(connect)] = ( "subgroup_"+str(i), "subgroup_"+str(j), gmat_complex[i,j] )
            if gmat_gi[i,j] != 0:
                conn_gi.loc[len(conn_gi)] = ( "subgroup_"+str(i), "subgroup_"+str(j), gmat_gi[i,j])

    connect.to_csv("../data/figures/figure3/network_connection.txt", header=True, index=False, sep='\t')
    conn_gi.to_csv("../data/figures/figure3/network_connection_gi.txt", header=True, index=False, sep='\t')

    # write vertex annotations to file
    deg = dict( G_complex.degree() )
    vertexDF = pd.DataFrame(columns=['name', 'full', 'short', 'degree'])
    vertexDF.loc[len(vertexDF)] = ("origin", "NA", "NA", "NA")
    for ix, i in enumerate( list(set(intab)) ):
        vertexDF.loc[len(vertexDF)] = ("group_"+str(ix), "NA", "NA", "NA")
    for ix, current_name in enumerate(list_nodes):
        current_deg = deg[current_name]
        current_short = list(shortnames[shortnames['Long']==current_name]['Short'])[0]
        vertexDF.loc[len(vertexDF)] = ("subgroup_"+str(ix), current_name, current_short, current_deg)
    vertexDF.to_csv("../data/figures/figure3/network_vertex.txt", header=True, index=False, sep='\t')

    
    return gmat_complex


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def score_full6(motifDF, G, GI):
    """
    check input motifs
    PPI graph: all edges from full6 motifs, edge weight is number of occurances
    GI graph: all GI edges between motif nodes, edge weight is average of per-complex interaction
    """

    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')
    shortnames = pd.read_csv("../data/yeast/MIPS_shortnames.txt", header=0, index_col=False, sep=' ')
    yeast_names, yeast_names_rev = orf2name()

    list_orf = []
    motif_columns = motifDF.columns
    for i in motif_columns:
        list_orf += list(motifDF[i])

    dict_orf = {}
    for i in list(set(list_orf)):
        dict_orf[i] = np.sum(np.array(list_orf) == i)

    best_idx = None
    best_score = 0
    for i in range(len(motifDF)):
        current_motif = list( motifDF.iloc[i] )
        current_score = 0
        current_gi = nx.subgraph(GI, current_motif)
        current_gimat = nx.to_numpy_array(current_gi)
        if np.any(current_gimat > 0):					# also pos GI
            for j in current_motif:
                current_score += dict_orf.get(j, 0)
            if current_score > best_score:
                best_score = current_score
                best_idx = i

    best_motif = sorted(list(motifDF.iloc[best_idx]))
    best_ppi = nx.subgraph(G, best_motif)
    best_gi  = nx.subgraph(GI, best_motif)

    gmat_ppi = nx.to_numpy_array(best_ppi, nodelist=best_motif)
    gmat_gi = nx.to_numpy_array(best_gi, nodelist=best_motif)
    np.savetxt("../data/figures/figure3/bestmotif_ppi.txt", gmat_ppi, fmt='%i')
    np.savetxt("../data/figures/figure3/bestmotif_gi.txt", gmat_gi, fmt='%.3f')

    networkDF = pd.DataFrame(columns=['from', 'to', 'weight', 'type', 'strength'])
    for i in range(len(best_motif)):
        for j in range(i+1, len(best_motif)):
            if gmat_ppi[i,j] == 1:
                networkDF.loc[len(networkDF)] = (best_motif[i], best_motif[j], 0.35, 'PPI', 0)
            if gmat_gi[i,j] > 0:
                networkDF.loc[len(networkDF)] = (best_motif[i], best_motif[j], gmat_gi[i,j], 'GIpos', 0.25)
            elif gmat_gi[i,j] < 0:
                networkDF.loc[len(networkDF)] = (best_motif[i], best_motif[j], gmat_gi[i,j], 'GIneg', 0.25)
  
    networkDF.to_csv("../data/figures/figure3/bestmotif_df.txt", header=True, index=False, sep='\t')

    # add name and complex_short
    nodesDF = pd.DataFrame(columns=['ORF', 'name', 'short'])

    for i in sorted(best_motif):
        current_complex = list(complexes[complexes['ORF']==i]['Complex'])[0]
        current_short = list(shortnames[shortnames['Long']==current_complex]['Short'])[0]
        current_name = yeast_names[i]
        nodesDF.loc[len(nodesDF)] = (i, current_name, current_short ) 
    nodesDF.to_csv("../data/figures/figure3/bestmotif_orfs.txt", header=True, index=False, sep='\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_yeastract():
    """
    Parse YEASTRACT transcription regulatory network (TRN)
    Data obtained by webserver query of all consensus TFs + targets
    Returns: TRN dataframe
    """

    yeast_names, yeast_names_rev = orf2name()

    yeastract = pd.DataFrame(columns=['TF', 'TG'])
    for line in open('../data/yeast/yeastract.txt', 'r'):
        current_line = line.split()

        current_tf = current_line[0]
        current_tg = current_line[1]

        current_tf_orf = yeast_names_rev.get(current_tf, 'none')
        current_tg_orf = yeast_names_rev.get(current_tg, 'none')

        if current_tf_orf != 'none' and current_tg_orf != 'none':
            yeastract.loc[len(yeastract)] = ( current_tf_orf, current_tg_orf )


    return yeastract


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def feedback_candidates(inputDF, G, GI):
    """
    Function to extract all motifs that link the yeast 
    transcription regulatory and metabolic networks
    inputDF: dataframe with list of ORFs for each query motif
    G: PPI network
    GI: GI network
    """

    yeast_names, yeast_names_rev = orf2name()

    # load YEASTRACT transcription regulatory network
    yeast_reg = parse_yeastract()
    list_tf = list(set(list(yeast_reg['TF']))) 

    # load yeast metabolic network
    metab = pd.read_csv("../data/yeast/metab/SysBioChalmers-yeast-GEM-c11e8dd/ModelFiles/txt/yeastGEM.txt", header=0, index_col=False, sep='\t')
    list_metab = []
    for i in list(metab['Gene-reaction association']):
        if len(str(i)) > 3:
            current_metab = i.replace('(', '').replace(')', '').replace('and','or').split('or')
            list_metab += current_metab
    list_metab = list(set(list_metab))
    
    netDF = pd.DataFrame(columns=['from', 'to', 'type', 'weight', 'strength'])
    idx_total = 0
    idx_count = 0
    visited = []
    for i in range( len(inputDF) ):
        current_motif = np.array(inputDF.loc[i])
        current_motif = list( current_motif[current_motif!='none'] )

        current_g = nx.subgraph(G, current_motif)
        current_gi = nx.subgraph(GI, current_motif) 
          
        current_tf = []
        current_metab = []

        for mx, current_orf in enumerate(current_motif) :
            if current_orf in list_tf and current_orf not in list_metab:
                current_tf.append(current_orf)
            elif current_orf not in list_tf and current_orf in list_metab:
                current_metab.append(current_orf)

        if len(current_tf) == 1 and len(current_metab) == 1 :       
            current_tf = current_tf[0]
            current_metab = current_metab[0]
            current_key = current_tf+"_"+current_metab

            if current_gi.has_edge(current_tf, current_metab) and current_key not in visited:
                visited.append(current_key)
                idx_count += 1

                for edge in list(current_g.edges()):
                    if edge[0] != edge[1]:
                        netDF.loc[len(netDF)] = (edge[0], edge[1], 'PPI', 1, 0)

                for edge in list(current_gi.edges()):
                    if edge[0] != edge[1]:
                        if current_gi[edge[0]][edge[1]]['weight'] > 0:
                            netDF.loc[len(netDF)] = (edge[0], edge[1], 'GIpos', current_gi[edge[0]][edge[1]]['weight'], 0.15)
                        else:
                            netDF.loc[len(netDF)] = (edge[0], edge[1], 'GIneg', current_gi[edge[0]][edge[1]]['weight'], 0.15)
        idx_total += 1  
    netDF.drop_duplicates(inplace=True)
    netDF.to_csv("../data/figures/figure4/tfmetab_network.txt", header=True, index=False, sep='\t')

    vertexDF = pd.DataFrame(columns=['name', 'gene', 'type', 'size'])
    for ix, i in enumerate(sorted(list(set( list(netDF['from']) + list(netDF['to']) )))):
        current_name = yeast_names.get(i, i)
        if i in list_tf:
            current_type = 'TF'
            current_size = 6
        elif i in list_metab:
            current_type = 'metab'
            current_size = 6
        else:
            current_type = 'other'
            current_size = 4
        vertexDF.loc[len(vertexDF)] = (i, current_name, current_type, current_size)
    vertexDF.to_csv("../data/figures/figure4/tfmetab_vertices.txt", header=True, index=False, sep='\t')














if __name__ == '__main__':

    ## Load network data
    G_ppi = nx.read_gpickle("../data/processed/yeast_ppi.nx")
    G_ppi_50 = filter_degree_ppi(G_ppi, 50)
    G_gi  = nx.read_gpickle("../data/processed/GI.nx")

    ## Extract no-cluster motifs
    list_orfs_noclust = pd.read_csv("../data/figures/figure2/list_orfs_noclust.txt", header=0, index_col=False)
    list_orfs_noclust = list(list_orfs_noclust['ORF'])  
    if os.path.exists("../data/figures/figure3/motifs_noclust.txt"):
        noclustDF = pd.read_csv("../data/figures/figure3/motifs_noclust.txt", header=0, index_col=False, sep='\t')
    else:
        noclustDF = motif_noclust(list_orfs_noclust)


    ## FIGURE 2 analyses
    fnm_total, fnm_supp, fnm_bypa = motif_suppressors("FNM")
    all_total, all_supp, all_bypa = motif_suppressors("ALL")
    suppDF = pd.DataFrame(data=np.array([[fnm_total, fnm_supp, fnm_bypa],[all_total, all_supp, all_bypa]]), columns=['all', 'supp', 'bypa'])
    suppDF.to_csv("../data/figures/figure2/yeast_suppressors.txt", header=True, index=False, sep='\t')

    fnm_ess, fnm_sub = motif_gencat("FNM")
    all_ess, all_sub = motif_gencat("ALL")
    gencatDF = pd.DataFrame({'FNM_ess': list(fnm_ess), 'FNM_sub': list(fnm_sub), "ALL_ess": list(all_ess), "ALL_sub": list(all_sub)})
    gencatDF.to_csv("../data/figures/figure2/gencatDF.txt", header=True, index=False, sep='\t')

    check_motif_noclust(noclustDF)


    ## FIGURE 3 analyses
#    df_ess = pd.read_csv("../data/figures/figure2/FNM.full6.ess.txt", header=0, index_col=False, sep='\t')
#    df_sub = pd.read_csv("../data/figures/figure2/FNM.full6.sub.txt", header=0, index_col=False, sep='\t')
#    gmat = full6_gi(df_sub, G_ppi_50, G_gi)
#    score_full6(df_sub, G_ppi_50, G_gi)


    ## FIGURE 4 analyses
    # complex network
#    parse_complex_net(noclustDF, 2)
#    parse_complex_net(noclustDF, 5)
#    parse_complex_net(noclustDF, 10)
#    reference_complex_net(G_ppi, 0, 2)
#    reference_complex_net(G_ppi, 0, 5)
#    reference_complex_net(G_ppi_50, 50, 2)
#    reference_complex_net(G_ppi_50, 50, 5)

    # TF - metabolism link
#    feedback_candidates(noclustDF, G_ppi_50, G_gi)

