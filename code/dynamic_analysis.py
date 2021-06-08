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

from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def orf2name():
    """
    utility function to convert between ORFs and gene names
    """

    #f_names = '../data/yeast/metab/SysBioChalmers-yeast-GEM-c11e8dd/ComplementaryData/databases/SGDgeneNames.tsv'
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
    filter networks based on maximum degree in PPI
    rationale: omit 'nonspecifically' interacting nodes
    removes nodes and associated edges with degree > max_degree
    """

    list_keep = []
    for entry in list(Gppi.degree()):
        if entry[1] <= max_degree:
            list_keep.append(entry[0])

    Gppi = nx.subgraph(Gppi, list_keep)
   

    return Gppi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_cossimi(packaged_input):
    """
    Function to test transcriptional co-regulation of motifs
    paralell cosine similarity (cossimi) calculations
    """

    expmat = packaged_input[0]
    colidx = packaged_input[1]

    NR, NC = np.shape(expmat)

    score = np.nan

    if colidx < NC:
        current_ref = np.ones(( NR ))
        current_score = expmat[:, colidx]           
        check_nan = np.isnan(current_score)
        if np.any(check_nan):
            current_score = current_score[check_nan == False]
            current_ref   = current_ref[check_nan == False]
          
        if len(current_ref) > 1 and len(current_score) > 1:
            current_score = current_score.reshape(1,-1)
            current_ref = current_ref.reshape(1,-1)
            score = cosine_similarity(current_ref, current_score)
            score = float(np.around(score, 3))
    
    return score


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def dynmat(motifDF):
    """
    Project Gasch expression data onto motifs, test for co-regulation
    orflist: limit ORFs to those on inputlist
    Returns: matrix of pairwise co-expression scores, dataframe of corresponding motifs
    """


    gasch = pd.read_csv("../data/yeast/gasch_stress.txt", header=0, index_col=False, sep='\t')

    # group time course data 
    columns_gasch = gasch.columns
    gasch_categories = [
        ['Heat Shock 05 minutes hs-1', 'Heat Shock 10 minutes hs-1', 'Heat Shock 15 minutes hs-1', 'Heat Shock 20 minutes hs-1', 'Heat Shock 30 minutes hs-1', 'Heat Shock 40 minutes hs-1', 'Heat Shock 60 minutes hs-1', 'Heat Shock 80 minutes hs-1'],
        ['Heat Shock 000 minutes hs-2', 'Heat Shock 000 minutes  hs-2', 'Heat Shock 000 minutes  hs-2.1', 'Heat Shock 005 minutes  hs-2', 'Heat Shock 015 minutes  hs-2', 'Heat Shock 030inutes  hs-2','Heat Shock 060 minutes  hs-2'],
        ['37C to 25C shock - 15 min', '37C to 25C shock - 30 min', '37C to 25C shock - 45 min', '37C to 25C shock - 60 min', '37C to 25C shock - 90 min'],
        ['heat shock 17 to 37, 20 minutes', 'heat shock 21 to 37, 20 minutes', 'heat shock 25 to 37, 20 minutes', 'heat shock 29 to 37, 20 minutes', 'heat shock 33 to 37, 20 minutes', '29C to 33C - 5 minutes', '29C to 33C - 15 minutes', '29C to 33C - 30 minutes', '33C vs. 30C - 90 minutes'],
        ['29C +1M sorbitol to 33C + 1M sorbitol - 5 minutes', '29C +1M sorbitol to 33C + 1M sorbitol - 15 minutes', '29C +1M sorbitol to 33C + 1M sorbitol - 30 minutes', '29C +1M sorbitol to 33C + *NO sorbitol - 5 minutes', '29C +1M sorbitol to 33C + *NO sorbitol - 15 minutes', '29C +1M sorbitol to 33C + *NO sorbitol - 30 minutes'],
        ['constant 0.32 mM H2O2 (10 min) redo', 'constant 0.32 mM H2O2 (20 min) redo', 'constant 0.32 mM H2O2 (30 min) redo', 'constant 0.32 mM H2O2 (40 min) rescan', 'constant 0.32 mM H2O2 (50 min) redo', 'constant 0.32 mM H2O2 (60 min) redo', 'constant 0.32 mM H2O2 (80 min) redo', 'constant 0.32 mM H2O2 (100 min) redo', 'constant 0.32 mM H2O2 (120 min) redo', 'constant 0.32 mM H2O2 (160 min) redo'],
        ['1 mM Menadione (10 min)redo', '1 mM Menadione (20 min) redo', '1 mM Menadione (30 min) redo', '1mM Menadione (40 min) redo', '1 mM Menadione (50 min)redo', '1 mM Menadione (80 min) redo', '1 mM Menadione (105 min) redo', '1 mM Menadione (120 min)redo', '1 mM Menadione (160 min) redo'],
        ['2.5mM DTT 005 min dtt-1', '2.5mM DTT 015 min dtt-1', '2.5mM DTT 030 min dtt-1', '2.5mM DTT 045 min dtt-1', '2.5mM DTT 060 min dtt-1', '2.5mM DTT 090 min dtt-1', '2.5mM DTT 120 min dtt-1', '2.5mM DTT 180 min dtt-1'],
        ['dtt 000 min  dtt-2', 'dtt 015 min dtt-2', 'dtt 030 min  dtt-2', 'dtt 060 min dtt-2', 'dtt 120 min dtt-2', 'dtt 240 min dtt-2', 'dtt 480 min dtt-2'],
        ['1.5 mM diamide (5 min)', '1.5 mM diamide (10 min)', '1.5 mM diamide (20 min)', '1.5 mM diamide (30 min)', '1.5 mM diamide (40 min)', '1.5 mM diamide (50 min)', '1.5 mM diamide (60 min)', '1.5 mM diamide (90 min)'],
        ['1M sorbitol - 5 min', '1M sorbitol - 15 min', '1M sorbitol - 30 min', '1M sorbitol - 45 min', '1M sorbitol - 60 min', '1M sorbitol - 90 min', '1M sorbitol - 120 min', 'steady-state 1M sorbitol'],
        ['Hypo-osmotic shock - 5 min', 'Hypo-osmotic shock - 15 min', 'Hypo-osmotic shock - 30 min', 'Hypo-osmotic shock - 45 min', 'Hypo-osmotic shock - 60 min'],
        ['aa starv 0.5 h', 'aa starv 1 h', 'aa starv 2 h', 'aa starv 4 h', 'aa starv 6 h'],
        ['Nitrogen Depletion 30 min.', 'Nitrogen Depletion 1 h', 'Nitrogen Depletion 2 h', 'Nitrogen Depletion 4 h', 'Nitrogen Depletion 8 h', 'Nitrogen Depletion 12 h', 'Nitrogen Depletion 1 d', 'Nitrogen Depletion 2 d', 'Nitrogen Depletion 3 d', 'Nitrogen Depletion 5 d'],
        ['Diauxic Shift Timecourse - 0 h', 'diauxic shift timecourse 9.5 h', 'diauxic shift timecourse11.5 ', 'diauxic shift timecourse 13.5 h', 'diauxic shift timecourse 15.5 h', 'diauxic shift timecourse 18.5 h', 'diauxic shift timecourse 20.5 h'],
        ['YPD 2 h ypd-2', 'YPD 4 h ypd-2', 'YPD 6 h ypd-2', 'YPD 8 h ypd-2', 'YPD 10 h  ypd-2', 'YPD 12 h ypd-2', 'YPD 1 d ypd-2', 'YPD 2 d ypd-2', 'YPD 3 d ypd-2', 'YPD 5 d ypd-2'],
        ['YPD stationary phase 2 h ypd-1', 'YPD stationary phase 4 h ypd-1', 'YPD stationary phase 8 h ypd-1', 'YPD stationary phase 12 h ypd-1', 'YPD stationary phase 1 d ypd-1', 'YPD stationary phase 2 d ypd-1', 'YPD stationary phase 3 d ypd-1', 'YPD stationary phase 5 d ypd-1', 'YPD stationary phase 7 d ypd-1', 'YPD stationary phase 13 d ypd-1', 'YPD stationary phase 22 d ypd-1', 'YPD stationary phase 28 d ypd-1'],
        ['DBY7286 37degree heat - 20 min', 'DBYmsn2-4- 37degree heat - 20 min', 'DBYmsn2/4 (real strain) + 37degrees (20 min)', 'DBYyap1- 37degree heat - 20 min (redo)', 'DBYyap1 + 37degree heat (repeat)', 'DBY7286 + 0.3 mM H2O2 (20 min)', 'DBYmsn2msn4 (good strain) + 0.32 mM H2O2', 'DBYmsn2/4 (real strain) + 0.32 mM H2O2 (20 min)', 'DBYyap1- + 0.3 mM H2O2 (20 min)', 'DBYyap1 + 0.32 mM H2O2 (20 min)'],
        ['ethanol vs. reference pool car-1', 'galactose vs. reference pool car-1', 'glucose vs. reference pool car-1', 'mannose vs. reference pool  car-1', 'raffinose vs. reference pool car-1', 'sucrose vs. reference pool car-1'],
        ['YP ethanol vs reference pool car-2', 'YP fructose vs reference pool car-2', 'YP galactose vs reference pool car-2', 'YP glucose vs reference pool car-2', 'YP mannose vs reference pool car-2', 'YP raffinose vs reference pool car-2', 'YP sucrose vs reference pool car-2']
    ]
    gasch_cat_names = ['Heat_shock_1', 'Heat_shock_2', 'Cold_shock', 'Temperature_div', 'HS+sorbitol', 'H2O2', 'Menadione', 'DTT_1', 'DTT_2', 'Diamide', 'Sorbitol', 'Hypo-osmotic_shock', 'Starvation', 'Nitrogen_depletion', 'Diauxic_shift', 'Ypd_2', 'Ypd_1', 'DBY', 'Carbon_1', 'Carbon_2']

    orfs_gasch = list(gasch['UID'])
    data = np.array(gasch[columns_gasch[1:]])
    data_nr, data_nc = np.shape(data)

    resmat = np.zeros(( len(motifDF), len(gasch_categories) )) * np.nan
    motifL = np.zeros(( len(motifDF) ))
    for i in range(len(motifDF)):
        current_motif = np.array( motifDF.iloc[i] )
        current_motif = list( current_motif[current_motif!='none'] )
        motifL[i] = len(current_motif)

        current_idx = np.zeros(( len(current_motif) )) * np.nan
        for mx, current_orf in enumerate(current_motif) :
            if current_orf in orfs_gasch:
                current_idx[mx] = int( orfs_gasch.index(current_orf) )
        current_idx = current_idx[np.isnan(current_idx)==False]
        current_idx = np.array(current_idx, dtype=int)

        for j in range(len(gasch_categories)):
            current_gasch = np.array(gasch[gasch_categories[j]])
            current_data = current_gasch[current_idx,:]  
            current_data = current_data 

            k = current_data.shape[0]
            current_combinations = list(it.combinations(range(k), 2))            
            current_output = np.zeros(( len(current_combinations) )) * np.nan

            for cx, comb in enumerate(current_combinations):
                vec1 = current_data[comb[0],:]
                vec2 = current_data[comb[1],:]
                sel_nonan = (np.isnan(vec1) == False) * (np.isnan(vec2) == False)

                if np.sum(sel_nonan) > 2:		# at least 2 non-nan data poits ... 
                    current_output[cx] = cosine_similarity( vec1[None, sel_nonan], vec2[None, sel_nonan] )

            current_output = current_output[np.isnan(current_output) == False]
            if len(current_output) > 0:
                resmat[i,j] = float(np.mean(current_output))
            else:
                resmat[i,j] = np.nan
   
    sel_nonan = np.sum(np.isnan(resmat),1) == 0
    resmat = resmat[sel_nonan,:]	   		# ignore rows with missing values, it's only for visualization
    motifDF = motifDF.iloc[sel_nonan]
    motifL = motifL[sel_nonan]

    dist = pdist(resmat)
    Z = linkage(dist, method="complete")    
    leaves_order = leaves_list(Z)
    resmat = resmat[leaves_order,:]
   
    resDF = pd.DataFrame(data=resmat, columns=gasch_cat_names)
    resDF.to_csv("../data/figures/figure5/cossimi.txt", header=True, index=False, sep='\t', na_rep='NA')
    np.savetxt("../data/figures/figure5/motifL.txt", motifL, fmt='%i')

    return resmat, motifDF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def analyze_dynmat(data, dataDF):
    """
    Extract possibly interesting clusters from heatmap
    """

    resDF = pd.read_csv("../data/figures/figure5/cossimi.txt", header=0, index_col=False, sep='\t')
    data = np.array(data)

    list_sel = []
    list_con = []

    for i in range(data.shape[1]):
        sel_pos = np.where(data[:,i] > 0.9)[0]
        list_sel += list(sel_pos)
        list_con += [resDF.columns[i],]*len(sel_pos)
    list_sel = np.array(list_sel, dtype=int)

    motifDF = dataDF.iloc[list_sel]
    motifDF['stress'] = list_con
    motifDF.to_csv("../data/figures/figure5/motifs_top_up.txt", header=True, index=False, sep='\t')

    return motifDF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def gasch_motifs(inputDF, G, GI):
    """
    Merge connected motifs from inputDF
    """

    yeast_names, yeast_names_rev = orf2name()
    
    list_stress = list(set(list(inputDF['stress'])))
    for current_stress in list_stress:
        dataDF = inputDF[inputDF['stress']==current_stress]
        dataDF = dataDF[dataDF.columns[:6]]

        G_new = nx.Graph()
        for motif_idx in range( len(dataDF) ):
            current_motif = np.array(dataDF.iloc[motif_idx])
            current_motif = list( current_motif[current_motif!='none'] )
            current_g = nx.subgraph(G, current_motif)

            for edge in list( current_g.edges() ):
                if edge[0] != edge[1]:
                    G_new.add_edge(edge[0], edge[1])

        conncomp = list( nx.connected_components(G_new) )
        for gx, g in enumerate(conncomp):
            G2 = nx.subgraph(G, list(g))
            GI2 = nx.subgraph(GI, list(G2.nodes()))

            if len( list(G2.nodes()) ) < 20:		# no super big clusters
                netDF = pd.DataFrame(columns=['from', 'to', 'type', 'weight', 'strength'])
                for edge in list(G2.edges()):
                    if edge[0] != edge[1]:
                        netDF.loc[len(netDF)] = (edge[0], edge[1], 'PPI', 1, 0)

                for edge in list(GI2.edges()):
                    if edge[0] != edge[1]:
                        if GI2[edge[0]][edge[1]]['weight'] > 0:
                            netDF.loc[len(netDF)] = (edge[0], edge[1], 'GIpos', GI2[edge[0]][edge[1]]['weight'], 0.15)
                        else:
                            netDF.loc[len(netDF)] = (edge[0], edge[1], 'GIneg', GI2[edge[0]][edge[1]]['weight'], 0.15)
                netDF.drop_duplicates(inplace=True)
                netDF.to_csv("../data/figures/figure5/motifs/"+current_stress+"_motif"+str(gx)+"_network.txt", header=True, index=False, sep='\t')

                vertexDF = pd.DataFrame(columns=['name', 'gene', 'type', 'size'])
                for ix, i in enumerate(sorted(list(set( list(netDF['from']) + list(netDF['to']) )))):
                    current_name = yeast_names.get(i, i)
                    current_type = 'other'
                    current_size = 4
                    vertexDF.loc[len(vertexDF)] = (i, current_name, current_type, current_size)
                vertexDF.to_csv("../data/figures/figure5/motifs/"+current_stress+"_motif"+str(gx)+"_vertices.txt", header=True, index=False, sep='\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ubr_motif(G, GI):
    """
    Explore one pre-defined motif (TF-metab-ubr motif)
    """

    yeast_names, yeast_names_rev = orf2name()
    orfs = ["YGR184C", "YGL058W", "YLR024C", "YBL041W", "YOR236W", "YDL020C"]

    current_g = nx.subgraph(G, orfs)
    current_gi = nx.subgraph(GI, orfs)

    # parse motif network
    netDF = pd.DataFrame(columns=['from', 'to', 'type', 'weight', 'strength'])
    for edge in list(current_g.edges()):
        if edge[0] != edge[1]:
            netDF.loc[len(netDF)] = (edge[0], edge[1], 'PPI', 1, 0)

    for edge in list(current_gi.edges()):
        if edge[0] != edge[1]:
            if current_gi[edge[0]][edge[1]]['weight'] > 0:
                netDF.loc[len(netDF)] = (edge[0], edge[1], 'GIpos', current_gi[edge[0]][edge[1]]['weight'], 0.15)
            else:
                netDF.loc[len(netDF)] = (edge[0], edge[1], 'GIneg', current_gi[edge[0]][edge[1]]['weight'], 0.15)

    netDF.drop_duplicates(inplace=True)
    netDF.to_csv("../data/figures/figure5/ubr_motif_network.txt", header=True, index=False, sep='\t')

    vertexDF = pd.DataFrame(columns=['name', 'gene', 'type', 'size'])
    for ix, i in enumerate(sorted(list(set( list(netDF['from']) + list(netDF['to']) )))):
        current_name = yeast_names.get(i, i)
        current_type = 'other'
        current_size = 4
        vertexDF.loc[len(vertexDF)] = (i, current_name, current_type, current_size)
    vertexDF.to_csv("../data/figures/figure5/ubr_motif_vertices.txt", header=True, index=False, sep='\t')

    # extract corresponding expression data
    gasch = pd.read_csv("../data/yeast/gasch_stress.txt", header=0, index_col=False, sep='\t')
    gasch_hs1 = ['Heat Shock 05 minutes hs-1', 'Heat Shock 10 minutes hs-1', 'Heat Shock 15 minutes hs-1', 'Heat Shock 20 minutes hs-1', 'Heat Shock 30 minutes hs-1', 'Heat Shock 40 minutes hs-1', 'Heat Shock 60 minutes hs-1', 'Heat Shock 80 minutes hs-1']
     
    orfs_gasch = list(gasch['UID'])
    current_idx = np.zeros(( len(orfs) )) * np.nan
    list_gene = []
    for mx, current_orf in enumerate(orfs) :
        list_gene.append(yeast_names.get(current_orf, current_orf))
        if current_orf in orfs_gasch:
            current_idx[mx] = int( orfs_gasch.index(current_orf) )
    current_idx = current_idx[np.isnan(current_idx)==False]
    current_idx = np.array(current_idx, dtype=int)

    data_ubr = gasch.iloc[current_idx][gasch_hs1]
    data_ubr.columns = ['hs5min', 'hs10min', 'hs15min', 'hs20min', 'hs30min', 'hs40min', 'hs60min', 'hs80min']
    data_ubr['ORF'] = orfs
    data_ubr['gene'] = list_gene
    data_ubr.to_csv("../data/figures/figure5/ubr_df.txt", header=True, index=False, sep='\t', na_rep='NA')







if __name__ == '__main__':

    G_ppi = nx.read_gpickle("../data/processed/yeast_ppi.nx")
    G_ppi_50 = filter_degree_ppi(G_ppi, 50)
    G_gi  = nx.read_gpickle("../data/processed/GI.nx")

    list_orfs_noclust = pd.read_csv("../data/figures/figure2/list_orfs_noclust.txt", header=0, index_col=False)
    list_orfs_noclust = list(list_orfs_noclust['ORF'])
    noclustDF = pd.read_csv("../data/figures/figure3/motifs_noclust.txt", header=0, index_col=False, sep='\t')

    # Co-regulation of motifs
    motif_gasch, motif_gaschDF = dynmat(noclustDF)

    # Top co-regulated motifs
    motif_topupDF = analyze_dynmat(motif_gasch, motif_gaschDF)
    gasch_motifs(motif_topupDF, G_ppi_50, G_gi)

    # TF-Ub-Metab motif analysis
    ubr_motif(G_ppi_50, G_gi)

