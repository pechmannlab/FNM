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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_summary():
    """
    This function pulls together all data on enumerated motifs, 
    merges the counts of isomorphic graphlets, and generates the
    data for Figures 1B and 1C. 
    """

    resultDF = pd.DataFrame(columns=['topo', 'k', 'count', 'count_rand', 'count_rand_sd', 'count_gi', 'count_gi_rand', 'count_gi_rand_sd'])

    for k in [3,4,5,6]:

        files_counts = glob.glob('../data/motif/rand/*_result_rand_d50_k'+str(k)+'_counts.txt')
        files_topo = glob.glob('../data/motif/rand/*result_rand_d50_k'+str(k)+'_topo.txt')
        files_topogi = glob.glob('../data/motif/rand/*result_randgi_d50_k'+str(k)+'.txt')

        result_counts = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        result_topo = pd.read_csv('../data/motif/result_d50_k'+str(k)+'_topo.txt', header=0, index_col=False, sep='\t')
    
        list_df = []
        list_topo = list(result_counts['G6'])
        for ix, i in enumerate(files_topo):
            data = pd.read_csv(i, header=0, index_col=False, sep='\t')
            list_topo = list(data['topo'])
            list_df.append(data)
        topo = sorted(list(set(list_topo)))

        iso_dict = {}
        iso_list = []

        for i in topo:              # loop through graphlet topology
            if i not in iso_list:
                g_i = nx.from_graph6_bytes( bytes(i, 'utf-8'))
                if len(list(iso_dict.keys()) ) > 0:
                    isocheck = False
                    for j in list(iso_dict.keys()):
                        g_j = nx.from_graph6_bytes( bytes(j, 'utf-8'))
                        if nx.is_isomorphic(g_i, g_j):
                            isocheck = True
                            break
                    if isocheck:
                        iso_list.append(i)
                        tmp_list = iso_dict[j]
                        tmp_list.append(i)
                        iso_dict[j] = tmp_list    
                    else:
                        iso_list.append(i)
                        iso_dict[i] = [i]
                else:
        	        iso_list.append(i)
        	        iso_dict[i] = [i]

        for i in list(iso_dict.keys()):
            current_iso = iso_dict[i]    
            counts = np.zeros(( 6 )) 		# total, rand, rand_sd, gi, gi_rand, gi_rand_sd

            for j in current_iso:
                counts_topo = np.zeros(( len(list_df) ))
                for d in range( len(list_df) ):
                    current_data = list_df[d]
                    if j in list(current_data['topo']):
                        counts_topo[d] = current_data[current_data['topo']==str(j)]['count'].item()
                    else:
                        counts_topo[d] = 0

                if j in list(result_topo['topo']):
                    counts_g6 = int( result_topo[result_topo['topo']==j]['count'].item() )
                else:
                    counts_g6 = 0

                # gi
                counts_topo_gi = np.zeros(( len(files_topogi) ))
                for dx, d in enumerate(files_topogi):
                    current_data = pd.read_csv(d, header=0, index_col=False, sep='\t')
                    if j in list(current_data['G6']):
                        counts_topo_gi[dx] = np.sum( np.array(list(current_data['G6']) ) == j )
                        counts_g6_gi = len( result_counts[result_counts['G6']==j])
                    else:
                        counts_topo_gi[dx] = 0
                        counts_g6_gi = 0
                counts += np.array([counts_g6, np.around(np.mean(counts_topo), 2), np.around(np.std(counts_topo), 2), counts_g6_gi, np.around(np.mean(counts_topo_gi),2), np.around(np.std(counts_topo_gi),2) ])

            resultDF.loc[len(resultDF)] = [i, k, counts[0], counts[1], counts[2], counts[3], counts[4], counts[5] ]
                 	
    resultDF.to_csv("../data/figures/figure1/motifDF.txt", header=True, index=False, sep='\t')

    # filter results for counts > random
    results2 = []
    for i in list(set(resultDF['k'])):
        current_df = resultDF[resultDF['k']==i]
        current_df = current_df[current_df['count'] > current_df['count_rand']]
        current_df.sort_values('count', ascending=False, inplace=True)
        current_df.sort_values('count_gi', ascending=False, inplace=True)
        if len(current_df) > 5:
            results2.append(current_df[0:7]) # only take the top of the sorted lists per k for the plot
        else:
            results2.append(current_df)
        resultDF2 = pd.concat(results2, ignore_index=True, sort=True)
    resultDF2.to_csv("../data/figures/figure1/motifDF2.txt", header=True, index=False, sep='\t')

    # isomorphic graphlets represented by first G6 name after sorting
    for i in glob.glob("../figures/Figure1/motifs/*.svg"):
        os.remove(i)

    for topo in resultDF2['topo']:
        print(topo)
        g = nx.from_graph6_bytes( bytes(topo, 'utf-8'))
        m = nx.to_numpy_array(g)
        np.savetxt('tmp.inp', m, fmt='%i')

        cmd_rand = "Rscript plot_motif.R"
        output = subprocess.run(cmd_rand, shell=True)
  
        dst = "../figures/Figure1/motifs/" + topo + ".svg"
        shutil.move("motif.svg", dst)
        os.remove("tmp.inp")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_gi(GI):
    """
    Function to compute summary statistics on genetic interactions
    in network motifs
    """

    #GI  = nx.read_gpickle("../data/processed/GI.nx")

    resultDF = pd.DataFrame(columns=['k', 'total', 'pos', 'neg', 'pct_pos', 'pct_neg'])

    for k in [3,4,5,6]:
        motif_data = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        current_columns = motif_data.columns
        current_orfs = list(current_columns[1:-1])
        motif_orfs = motif_data[current_orfs]
        
        for i in range(len(motif_data)): 
            current_orfs =  list(motif_orfs.loc[i]) 
            current_graph = nx.subgraph(GI, current_orfs)
            current_mat = nx.to_numpy_array(current_graph)

            current_mat[current_mat > 0] = 1        # discretize genetic interaction graph
            current_mat[current_mat < 0] = -1

            e_max = k*(k-1)/2.
            e_tot = np.sum(current_mat != 0)/2.
            e_pos = np.sum(current_mat > 0)/2.
            e_neg = np.sum(current_mat < 0)/2.

            pct_pos = e_pos / e_max
            pct_neg = e_neg / e_max

            resultDF.loc[len(resultDF)] = (int(k), int(e_tot), int(e_pos), int(e_neg), np.around(pct_pos,2), np.around(pct_neg,2))

            if i%10000==0:
                print(i, "lines processed")

    resultDF.to_csv("../data/figures/figure2/gi_summary.txt", header=True, index=False, sep='\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def histo_results():
    """
    Summary stats on genetic interactions in motifs
    Generates data for Figure S1A
    """

    results = []

    for i in [25, 50, 100]:
        for j in [3,4,5,6]:
            gi = np.loadtxt("../data/motif/result_d"+str(i)+"_k"+str(j)+"_gihist.txt")
            mat = np.zeros(( len(gi), 4 ))
            mat[:,0] = gi/np.sum(gi)
            mat[:,1] = np.arange(1,21)/20.
            mat[:,2] = np.repeat(j, 20)
            mat[:,3] = np.repeat(i, 20)

            df = pd.DataFrame(data=mat, columns=['value', 'bin', 'k', 'deg'])
            results.append(df)

    resultDF = pd.concat(results, ignore_index=True, sort=True)

    resultDF.to_csv("../data/figures/figure1/gihistogram.txt", header=True, index=False, sep='\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def rand_stats():
    """
    Computes some stats on randomizations
    """

    randmat = np.zeros(( 30, 4 )) *np.nan      # 30 rand trails, shows that sufficient

    for kx, k in enumerate([3,4,5,6]):
        for i in range(30):
            fIN = "../data/motif/rand/rand"+str(i)+"_result_rand_d50_k"+str(k)+"_counts.txt"
            if os.path.exists(fIN):
                data = np.loadtxt(fIN)
                randmat[i,kx] = data[0]
 
    randDF = pd.DataFrame(columns=['N', 'rand_mean', 'rand_std', 'k'])
    for kx, k in enumerate([3,4,5,6]):
        for i in range(0, 30):
            randmat_mean = np.around( np.mean( randmat[0:(i+1),kx]), 2) 
            randmat_std = np.around( np.std( randmat[0:(i+1),kx]), 2) 
            randDF.loc[len(randDF)] = (i+1, randmat_mean, randmat_std, k)

    randDF.to_csv("../data/figures/figure1/rand.txt", header=True, index=False, sep='\t')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cluster_motifs(G_ppi):
    """
    Ad hoc clustering of a large set of network motifs based on overlap.
    This function will identify both clusters of motifs and motifs that don't cluster. 
    - binary motif matrix [motifs, ORFs]
    - motif network from mapping motif ORFs onto PPI
    - graphlet enumerate within motif network 
    Returns: number of motifs per ORF
    """

    def rename_graph(G, labellist):
        label_dict = {}
        label_num = 0 
        for node in labellist:
            label_dict[node] = label_num
            label_num += 1
        Gnew = nx.relabel_nodes(G, label_dict, copy=True)

        return Gnew


    def getcluster(motprof, motifmat, enumotifs, list_orfs):
        """
        subroutine to identify and prune cluster
        """
        ka = np.shape(enumotifs)[1]
        motif_profile_current = motprof[np.array(enumotifs, dtype=int)]  
        motif_profile_min = np.min(motif_profile_current, axis=1)

        # get max minimum coverage
        current_max = np.max(motif_profile_min)
        current_max_idx = np.where(motif_profile_min == current_max)[0]
        if len(current_max_idx) > 1:            # if 
            tmp_idx = np.argmax( np.sum(motif_profile_current[current_max_idx,:], axis=1) )
            current_max_idx = current_max_idx[tmp_idx]
        sel_motif = enumotifs[int(current_max_idx),:]

        # select columns from motifmat of motif with highest 'coverage'
        s = np.ones(( np.shape(motifmat)[0] ))
        for i in range( len(sel_motif) ): 
            current_sel = motifmat[:,sel_motif[i]] == 1
            s *= current_sel
        tmp_motifmat = motifmat[:,sel_motif]
        sel_mm = np.sum(tmp_motifmat,1) == ka
        sel_orf = np.sum(motifmat[sel_mm,:], 0) > 0

        # extract motif and update motif profile
        clust_motif = []
        for i in np.where(sel_orf)[0]:
            clust_motif.append(list_orfs[i])
            motprof[i] = 0

        return motprof, clust_motif


    #resultDF = pd.DataFrame(columns=['k', 'total', 'pos', 'neg', 'pct_pos', 'pct_neg'])

    # loading yeast ORF - gene name mapping
    f_names = '../data/yeast/metab/SysBioChalmers-yeast-GEM-c11e8dd/ComplementaryData/databases/SGDgeneNames.tsv'
    yeast_names = {}
    for line in open(f_names, 'r'):
        current_line = line.split()
        orf = current_line[0]
        if len(current_line) > 1:
            name = current_line[1]
            yeast_names[orf] = name
        else:
            yeast_names[orf] = orf
    
    # LOAD DATA
    list_orfs = []
    N = 0
    for k in [3,4,5,6]:
        motif_data = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        current_columns = motif_data.columns
        current_orfs = list(current_columns[1:-1])
        N += len(motif_data)
        for i in current_orfs:
            list_orfs += list(set(list(motif_data[i])))
    list_orfs = sorted(list(set(list_orfs)))           # make ordered so that following results are deterministic


    # PARSE ALL MOTIFS INTO MATRIX   
    motifmat = np.zeros(( N, len(list_orfs) ), dtype=int)
    idx = 0
    for k in [3,4,5,6]:
        motif_data = pd.read_csv('../data/motif/result_d50_k'+str(k)+'.txt', header=0, index_col=False, sep='\t')
        current_columns = motif_data.columns
        current_columns = list(current_columns[1:-1])
        for i in range(len(motif_data)):
            current_motif = motif_data.loc[i]
            for col in current_columns:
                current_orf = str(current_motif[col])
                if current_orf in list_orfs:
                    current_orf_idx = list_orfs.index(current_orf)
                    motifmat[idx,current_orf_idx] = 1
                else:
                    print("there is a problem")
            idx += 1

    G = nx.subgraph(G_ppi, list_orfs)
    G2 = rename_graph(G, list_orfs)
    motif_profile = np.sum(motifmat, axis=0)


    # save network for plotting, some is redundant
    Gmat = nx.to_numpy_array(G2, nodelist=list(np.arange(len(list_orfs))))
    np.savetxt("../data/figures/figure2/network_FNM.txt", Gmat, fmt='%i')
    np.savetxt("../data/figures/figure2/profile_FNM.txt", motif_profile, fmt='%i')
    netw = pd.DataFrame({'ORF':list_orfs, 'N': list(motif_profile)})
    netw.to_csv("../data/figures/figure2/network_FNM_profile.txt", header=True, index=False, sep='\t')
   

    # CLUSTERING STARTS HERE 
    start_time = time.time() 
    cluster_results = []
    for k in [5,4,3,2]:         # start with largest graphlets that required biggest overlap
        Gmat = nx.to_numpy_array(G2, nodelist=list(np.arange(len(list_orfs))))
        motifsres = enumerate_kmotif(Gmat, k)
        L = 100
        while L > 0:
            motif_profile, cluster_motif = getcluster(motif_profile, motifmat, motifsres, list_orfs)
            L = len(cluster_motif)
            if L > 0:
                cluster_results.append(cluster_motif)
                print(k, L)
   
    allmotif = []
    allmotif_idx = []
    
    for ix, i in enumerate(cluster_results):
        gsub = nx.subgraph(G_ppi, i)
        gmat = nx.to_numpy_array(gsub, nodelist=list(i))
        np.savetxt("../data/figures/figure2/clusters/cluster."+str(ix)+".txt", gmat, fmt='%i')
        gsub_names = []
        for j in i:
            allmotif.append(j)
            allmotif_idx.append( list_orfs.index(j))
            gsub_names.append( yeast_names.get(j) )
        nodesDF = pd.DataFrame({'ORF':gsub_names})
        nodesDF.to_csv("../data/figures/figure2/clusters/cluster."+str(ix)+"_orfs.txt", header=True, index=False)

    rest = []
    rest_orf = []
    for i in list_orfs:
        if i not in allmotif:
            rest.append(list_orfs.index(i))
            rest_orf.append(i)
    rest = np.array(rest, dtype=int)

    restDF = pd.DataFrame({'ORF':list(rest_orf)})
    restDF.to_csv("../data/figures/figure2/list_orfs_noclust.txt", header=True, index=False)

    motifmat2 = motifmat[:,rest]
    motifmat2 = motifmat2[np.sum(motifmat2,1)>0,:]

    number_motifs_total = np.shape(motifmat)[0]
    number_motifs_noclust = np.shape(motifmat)[0] - np.sum(np.sum(motifmat[:,np.array(allmotif_idx, dtype=int)], 1) > 0)
    number_orfs_motif = len(list_orfs)
    number_orfs_motif_noclust = len(rest)

    summaryDF = pd.DataFrame({'cat': ['Total', 'Noclust'], 'orf': [number_orfs_motif, number_orfs_motif_noclust], 'motif': [number_motifs_total, number_motifs_noclust]})
    summaryDF.to_csv('../data/figures/figure2/motif_summary.txt', header=True, index=False, sep='\t')

    netw2 = netw.iloc[rest]
    netw2.to_csv("../data/figures/figure2/network_FNM_profile_noclust.txt", header=True, index=False, sep='\t')

    G_noclust = nx.subgraph(G_ppi, rest_orf)
    G_noclust = rename_graph(G_noclust, rest_orf)
    Gmat_noclust = nx.to_numpy_array(G_noclust, nodelist=list(np.arange(len(rest_orf))))
    np.savetxt("../data/figures/figure2/network_FNM_noclust.txt", Gmat_noclust, fmt='%i')
   
    motifmat_rest = motifmat[:,rest]
    motifmat_rest = motifmat_rest[np.sum( motifmat[:, np.array(allmotif_idx, dtype=int)], 1) == 0,:]
    motifmat_rest_profile = np.sum(motifmat_rest, 0)
    np.savetxt("../data/figures/figure2/profile_FNM_noclust.txt", motifmat_rest_profile, fmt='%i')

    print("--- %s seconds --- " % (time.time() - start_time))
    
  
    return netw


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def enumerate_kmotif(G_mat, k):
    """
    Subroutine of Kavosh motif enumeration
    """

    def partition(n):
        result = []
        result.append((n, ))
        for x in range(1, n):
            for y in partition(n - x):
                result.append( tuple(sorted((x, ) + y)))
        return result

    def composition(n):
        composition = []
        pat = partition(n)
        for p in pat:
            for comp in list(set(list(it.permutations(p)))):
                if comp not in composition:
                    composition.append(comp)
        return composition

    def parse_motif(motif_list):
        motif = []
        for i in motif_list:
            if np.size(i) == 1:
                motif.append(i)
            elif np.size(i) > 1:
                for j in i:
                    motif.append(j)
        return motif 
   
    result_mat = np.zeros(( 100000000, k ), dtype=int)
    count_motif = 0
    for S0 in range(np.shape(G_mat)[0]):
        G_mat[np.arange(S0), :] = 0
        G_mat[:, np.arange(S0)] = 0
        list_graphlet = composition(k-1)

        for g in range(len(list_graphlet)):
            current_graphlet = list_graphlet[g]
            num_levels = np.size(current_graphlet) 
            current_level = 0

            visited = np.zeros(( np.shape(G_mat)[0]  ), dtype=bool)
            visited[S0] = True

            queue = []
            queue.append([S0])
            queue2 = []                         # parallel copy of queue without nesting!
            queue2.append([S0])
            for i in range(np.size(current_graphlet)):
                queue.append([])
                queue2.append([])

            read_pointer = np.zeros(( num_levels + 1 ), dtype=int)
            len_queue = [len(queue[x]) for x in range(len(queue))]
            len_queue2 = [len(queue2[x]) for x in range(len(queue2))]

            current_motif = []
            s = 1
     
            while np.any(read_pointer < (len_queue) ): 
                current_funcmot = 0
                #print(read_pointer, len_queue)
                if np.size(current_graphlet) > 1: 
                    if current_level < np.size(current_graphlet):
                        s = current_graphlet[current_level]
                    else:
                        s = 1 #placeholder
                elif np.size(current_graphlet) == 1:
                    s = int(current_graphlet[0])

                if read_pointer[current_level] < len(queue[current_level]):  
                    current_node = queue[current_level][read_pointer[current_level]]
                    read_pointer[current_level] += 1
                    current_motif.append(current_node)

                    if current_level < num_levels:
                        if s == 1:
                            if np.size(current_node) == 1:
                                current_parent = current_node
                            elif np.size(current_node ) > 1:
                                current_parent = current_node[0] 
                            current_children = np.where(G_mat[current_parent,:] == 1)[0] 
                            current_children = current_children[ visited[current_children] == False ] 
                            visited[current_children] = True
                            queue[int(current_level+1)] += list(current_children)
                            queue2[int(current_level+1)] += list(current_children)
                            current_level += 1               

                        elif s > 1:
                            if np.size(current_node) == 1:
                                current_parent = current_node
                            elif np.size(current_node ) > 1:
                                current_parent = current_node[0] 
                            current_children = np.where( G_mat[current_parent,:] == 1)[0] 
                            current_children = current_children[ visited[current_children] == False ] 
                            if current_level < num_levels - 1:
                                current_combs = [list(x) for x in it.permutations(current_children, int(s) ) ]
                            elif current_level == num_levels -1:
                                current_combs = list(it.combinations(current_children,int(s) ))
                            queue[int(current_level+1)] += current_combs
                            visited[current_children] = True               
                            queue2[int(current_level+1)] += list(current_children)
                            current_level += 1               

                    elif current_level == num_levels:
                        result_motif = parse_motif(current_motif)
                        #print(result_motif)
                        result_mat[count_motif, :] = result_motif
                        count_motif += 1
                        current_motif.pop()

                elif read_pointer[current_level] == len(queue[current_level]):
                    list_current_level = np.array(queue2[current_level], dtype=int)
                    visited[list_current_level] = False 
                    queue[current_level] = []
                    queue2[current_level] = []
                    read_pointer[current_level] = 0
                    current_motif.pop()
                    current_level -= 1 
            
                len_queue = [len(queue[x]) for x in range(len(queue))]
                len_queue2 = [len(queue2[x]) for x in range(len(queue2))]

                if np.all(read_pointer == len_queue):   
                    break

    result_mat = result_mat[np.sum(result_mat,1)>0,:]

    return result_mat


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_stats(G, GI, countDF):
    """
    Adds some statistics of network properties for each ORF to input DF
    - degree in PPI and GI networks
    - betweenness centrality in PPI and GI networks
    """

    result_dict = {}

    deg_g = dict(nx.degree(G))
    deg_gi = dict(nx.degree(GI))

    betw_g = nx.betweenness_centrality(G)
    betw_gi = nx.betweenness_centrality(GI)

    list_g = []
    list_gi = []
    listb_g = []
    listb_gi = []
    for i in list(countDF['ORF']): 
        current_deg_g = deg_g.get(i, 0)
        current_deg_gi = deg_gi.get(i, 0)
        current_betw_g = betw_g.get(i,0)
        current_betw_gi = betw_gi.get(i,0)
 
        list_g.append(current_deg_g)
        list_gi.append(current_deg_gi)
        listb_g.append(current_betw_g)
        listb_gi.append(current_betw_gi)        

    countDF['degG'] = list_g
    countDF['degGI'] = list_gi
    countDF['betwG'] = listb_g
    countDF['betwGI'] = listb_gi

    countDF.to_csv("../data/figures/figure2/degree.txt", header=True, index=False, sep='\t')





    

if __name__ == '__main__':

    # load network data
    G_ppi = nx.read_gpickle("../data/processed/yeast_ppi.nx")
    G_gi  = nx.read_gpickle("../data/processed/GI.nx")


    # FNM analyses
    motif_summary()         # counts of motifs
    motif_gi(G_gi)          # stats on genetic interactions in motifs
    histo_results()         # stats on genetic interactions in motifs
    rand_stats()            # stats on randomizations


    # motif clustering
    orf_counts = cluster_motifs(G_ppi)      # clustering of motifs
    motif_stats(G_ppi, G_gi, orf_counts)


    


