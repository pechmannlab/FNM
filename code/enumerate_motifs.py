import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import itertools as it
import time
import multiprocessing as mp
import glob
import subprocess




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def partition(n):
    """
    partition of integer n
    """

    result = []
    result.append((n, ))
    for x in range(1, n):
        for y in partition(n - x):
            result.append( tuple(sorted((x, ) + y)))

    return result


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def composition(n):
    """
    computes combinations of partition of integer n
    """

    composition = []
    pat = partition(n)
    for p in pat:
        for comp in list(set(list(it.permutations(p)))):
            if comp not in composition:
                composition.append(comp)

    return composition


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_motif(motif_list):
    """
    parses motif, flattens potentially nested list
    """

    motif = []
    for i in motif_list:
        if np.size(i) == 1:
            motif.append(i)
        elif np.size(i) > 1:
            for j in i:
                motif.append(j)
    #motif = list(set(motif))

    return motif 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def to_g6(motif_array):
    """
    extract motif and convert to graph6 format
    """
    
    g_motif = nx.from_numpy_array(motif_array)
    g6 = nx.to_graph6_bytes(g_motif, nodes=None, header=False)

    g6 = g6.rstrip().decode('utf-8')

    return g6


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def rename_graph(G, labels=[]):
    """
    renames nodes in graph to integers
    returns renamed graph and dictionary with node mappings
    """

    if len(labels) == 0:
        label_dict = {}
        label_rev = {}
        label_num = 0 
        for node in G.nodes():
            label_dict[node] = label_num
            label_rev[label_num] = node
            label_num += 1
    else:
        label_dict = labels 
        label_rev = {}
        for i in list(label_dict.keys()):
            label_rev[label_dict[i]] = i

    Gnew = nx.relabel_nodes(G, label_dict, copy=True)

    return Gnew, label_dict, label_rev
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def check_complexes(cmplx, label_dict):
    """
    only one complex subunit per motif to focus on between complex motifs
    label_dict from renaming nodes to integers
    """
    
    complex_mat = np.zeros(( len(label_dict), len(label_dict) ))*np.nan	        # can be max all genes
    complex_list = np.zeros(( len(label_dict) ), dtype=bool)
    for i in list(label_dict.keys()) :
        if i in list(cmplx['ORF']):
            current_label = int( label_dict.get(i, -1) )  
            if current_label > -1:     # not needed but still
                complex_list[current_label] = True
                current_list_complex = []
                current_complexes = list(cmplx[cmplx['ORF']==str(i)]['Complex'])
                current_ids = [] 
                for c in cmplx[cmplx['ORF']==str(i)]['Complex']:
                    current_list_complex_orf = list( cmplx[cmplx['Complex']==str(c)]['ORF'] )
                    current_list_complex = np.array( [label_dict.get(i, -1) for i in current_list_complex_orf] )
                    current_list_complex = current_list_complex[current_list_complex > -1] 
                    current_ids += list(current_list_complex)
                complex_mat[current_label, 0:len(current_ids)] = np.array(current_ids)
               
    complex_mat = complex_mat[ :, np.sum(np.isnan(complex_mat), 0) != np.shape(complex_mat)[0] ]

    return complex_list, complex_mat


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def enumerate_kmotif(G_mat, GI_mat, S0, k, list_cmplx, mat_cmplx, combomat, labels_rev, theta_func, writeDF, TMPDIR):
    """
    Implementation of Kavosh motif enumeration (Kashani et al. BMC Bioinfo, 2009) with own modifications
    G_mat: input graph (PPI) in np.array form
    GI_mat: genetic interaction graph in np.array form
    S0: starting 'source' node as index for graph arrays
    k: size of motifs to be enumerated
    list_cmplx: list of protein complexes
    mat_cmplx: matrix of protein complex subunit pairings
    combomat: precomputed relationships of protein complex subunits for fast queries
    labels_rev: dictionary to link matrix indices to ORF names
    writeDF: logical to write extended output to disk
    TMPDIR: temporary working directory
    """

    if k > 6:
        print("do you really want this? check your hardware first")
        return None

    list_graphlet = composition(k-1)

    discovered_graphlets = {} 
    num_discovered = 0

    resultcols = ['ORF'+str(i+1) for i in range(k)] + ['pct_gi', 'G6']   
    resultDF = pd.DataFrame(columns=resultcols)

    giscore = np.zeros(( 20 ), dtype=int)
    num = 0
    fnum = 0

    for g in range(len(list_graphlet)):
        current_graphlet = list_graphlet[g]
        num_levels = np.size(current_graphlet) 
        current_level = 0

        visited = np.zeros(( np.shape(G_mat)[0]  ), dtype=bool)
        visited[S0] = True

        if list_cmplx[S0]:
            current_complex = mat_cmplx[S0]
            current_complex = np.array(current_complex[np.isnan(current_complex) == False], dtype=int)
            visited[current_complex] = True

        queue = []
        queue.append([S0])
        queue2 = []                    # parallel copy of queue without nesting!
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
            if np.size(current_graphlet) > 1: 
                if current_level < np.size(current_graphlet):
                    s = current_graphlet[current_level]
                else:
                    s = 1 #placeholder, though not needed
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

                        complex_child = current_children[ list_cmplx[current_children] ]                         
                        if len(complex_child) > 0:
                            current_subunits = mat_cmplx[complex_child,:]
                            current_subunits = current_subunits.flatten()
                            current_subunits = np.array(current_subunits[np.isnan(current_subunits)==False], dtype=int)
                            current_subunits = current_subunits[ visited[current_subunits] == False ] 
                            visited[current_subunits] = True
                            queue2[int(current_level+1)] += list(set(current_subunits))

                        current_level += 1               

                    elif s > 1:
                        if np.size(current_node) == 1:
                            current_parent = current_node
                        elif np.size(current_node ) > 1:
                            current_parent = current_node[0] 
                        current_children = np.where( G_mat[current_parent,:] == 1)[0] 
                        current_children = current_children[ visited[current_children] == False ] 
                        #visited[current_children] = True ## see below

                        if current_level < num_levels - 1:
                            current_combs = [list(x) for x in it.permutations(current_children, int(s) ) ]
                        elif current_level == num_levels -1:
                            current_combs = list(it.combinations(current_children,int(s) ))

                        ## omit combinations that have two or more orfs from same complex
                        current_combs_nocomplex = []
                        current_children_notsamecomplex = []
                        sel_current_combs_nocomplex = np.zeros(( len(current_combs) ), dtype=bool)
                        for ix_combo, combo in enumerate(current_combs):
                            current_combo = combomat[combo,:][:,combo]

                            if np.all(current_combo == False):
                                #current_combs_nocomplex.append(combo)		# retaining tuple/list structure
                                sel_current_combs_nocomplex[ix_combo] = True
                                current_children_notsamecomplex += combo        # flattened list for queue
                                current_children_notsamecomplex = list(set(current_children_notsamecomplex))

                        current_children_notsamecomplex = np.array(current_children_notsamecomplex, dtype=int)
                        current_combs = np.array(current_combs)
                        current_combs_nocomplex = current_combs[sel_current_combs_nocomplex]

                        current_combs = list(current_combs_nocomplex)
                        queue[int(current_level+1)] += current_combs
                        visited[current_children_notsamecomplex] = True 
                       
                        queue2[int(current_level+1)] += list(current_children)
                        queue2[int(current_level+1)] += list(current_children_notsamecomplex)

                        notsamecomplex_child = current_children_notsamecomplex[ list_cmplx[current_children_notsamecomplex] ] 
                        if len(notsamecomplex_child) > 0:
                            current_subunits = mat_cmplx[notsamecomplex_child,:]
                            current_subunits = current_subunits.flatten()
                            current_subunits = np.array(current_subunits[np.isnan(current_subunits)==False], dtype=int)

                            # only add new ones, not subunits that were already added at upper level (should not tho)
                            current_subunits = current_subunits[ visited[current_subunits] == False ] 
                            visited[current_subunits] = True
                            queue2[int(current_level+1)] += list(set(current_subunits))

                        current_level += 1               

                elif current_level == num_levels:
                    result_motif = parse_motif(current_motif)
                    current_motif.pop()
   
                    g_array = G_mat[result_motif,:][:,result_motif]
                    gi_array = GI_mat[result_motif,:][:,result_motif]

                    current_g6 = to_g6(g_array)
                    if current_g6 in list(discovered_graphlets.keys()):
                        discovered_graphlets[current_g6] += 1
                    else:
                        discovered_graphlets[current_g6] = 1


                    current_funcmot, gi_score = functional_motif(g_array, gi_array, current_graphlet, theta_func)
                    num += 1
                    fnum += int(current_funcmot)
                    gi_idx = np.min([int(gi_score * 20), 19])
                    giscore[gi_idx] += 1


                    if writeDF:
                        if current_funcmot == 1:        # write output (FIX!)
                            result_orf = []
                            for node in result_motif:
                                result_orf.append(labels_rev[node])
                            output_motif = [result_orf[x] for x in range(len(result_orf))] 
                            resultDF.loc[len(resultDF)] = output_motif + [np.around(gi_score,2), current_g6] 

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

    if writeDF:
        resultDF.to_csv(TMPDIR + str(S0) + '_motif.txt', header=True, index=False, sep='\t')

    # write to disk easier for multiprocessing, debugging, restarting, etc.
    np.savetxt(TMPDIR+str(S0)+"_counts.txt", np.array([num, fnum]), fmt='%i')
    np.savetxt(TMPDIR+str(S0)+"_gi.txt", giscore, fmt='%i')

    discovered_graphletsDF = pd.DataFrame({'topo': list(discovered_graphlets.keys()), 'count': list(discovered_graphlets.values() )}  )
    if len(discovered_graphletsDF) > 0:
        discovered_graphletsDF.to_csv(TMPDIR+str(S0)+"_topo.txt", header=True, index=False, sep='\t')
   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def functional_motif(gmat, gimat, graphlet, theta_func=0.5):
    """
    checks GI edges in PPI motif:
    functional if there is GI for all source to last layer
    and total fraction of GI edges > theta_func
    """

    me = 0
    gi = 0

    # could also try just gi edges, not equiv    
    gi = np.sum(gimat != 0) / (np.prod(np.shape(gimat))-np.shape(gimat)[0] )
    
    # interactions between source and last layer
    number_last_level = graphlet[-1]
    distant_edges = gimat[0, -number_last_level:]
    if np.all(distant_edges != 0) and gi >= theta_func:      
        me = 1

    return me, np.around(gi, 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter_degree(Gppi, Ggi, max_degree):
    """
    filter networks based on maximum degree in PPI
    rationale: omit 'nonspecifically' interacting nodes and reduce compute cost ...
    removes nodes and associated edges with degree > max_degree
    """

    list_keep = []
    for entry in list(Gppi.degree()):
        if entry[1] <= max_degree:
            list_keep.append(entry[0])

    Gppi = nx.subgraph(Gppi, list_keep)
    Ggi  = nx.subgraph(Ggi, list_keep)

    return Gppi, Ggi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def find_motifs(G, GI, cmplx, k, Fout_pre, randomized, randomized_gi, theta_func, writeDF=False):
    """
    wrapper function to enumerate motifs
    goes through each node as source (i.e. takes time to compute)
    G:  PPI network in nx format
    GI: GI network in nx format
    k:  size of motif
    Fout_pre: prefix for output files
    randomized: logical whether PPI network is randomized, using R package Birewire
    randomized_gi: logical whether GI network is randomized
    theta_func: threshold for fraction of GI edges in motif to consider 'functional'
    writeDF: write extended output to disk
    """

    TMPDIR = 'TMP_' + Fout_pre + '/'
    if not os.path.exists(TMPDIR):
        os.mkdir(TMPDIR)

    G, labels, labels_rev = rename_graph(G)
    GI, labels, labels_rev = rename_graph(GI, labels)

    G_mat = nx.to_numpy_array(G, nodelist=list(np.arange(len(labels))) )
    np.fill_diagonal(G_mat, 0)

    GI_mat = nx.to_numpy_array(GI, nodelist=list(np.arange(len(labels))) )
    GI_mat[GI_mat < 0] = -1       # discretize GI interactions, prefiltered for significance
    GI_mat[GI_mat > 0] = 1
    np.fill_diagonal(GI_mat, 0)
    GI_mat = np.array(GI_mat, dtype=int)

    list_cmplx, mat_cmplx = check_complexes(cmplx, labels)

    if randomized:
        G_mat_rand = randomized_network(G_mat)
        G_mat = np.copy(G_mat_rand)
    
    if randomized_gi:
        GI_mat[GI_mat < 0] = 1      # needs [0,1] matrix for randomization, i.e. no neg ints
        GI_mat = randomized_network(GI_mat)

    # precompute protein complex matrix for logical queries
    combomat = np.zeros(np.shape(G_mat), dtype=bool)
    for x in range(np.shape(G_mat)[0]):
        for y in range(x+1, np.shape(G_mat)[0] ):
            if (y in list(mat_cmplx[x,:])) or (x in list(mat_cmplx[y,:])):
                combomat[x,y] = combomat[y,x] = True 

    iter_max = np.shape(G_mat)[0]
    # for restart omit what already computed
    iter_list = []
    computed = [int(name.split('/')[1].split('_')[0]) for name in glob.glob(TMPDIR+'*_counts.txt')]
    for ii in np.arange(iter_max):
        if ii not in computed:
            iter_list.append(ii)

    pool = mp.Pool(processes=20)        # parallelize computation
    for s in iter_list:
        pool.apply_async(motif_parallel, [(s, G_mat, GI_mat, k, list_cmplx, mat_cmplx, combomat, labels_rev, theta_func, writeDF, TMPDIR)] ) 
    pool.close()
    pool.join()

    # COLLECT DATA FROM INDIVIDUAL FILES (SIMPLER FOR RESTART)  
    num = 0
    fnum = 0
    gihist = np.zeros(( 20 ))
    topodict = {}
    CLEANFLAG = True

    # combine counts
    for res in glob.glob(TMPDIR+'*_counts.txt'):
        current_counts = np.loadtxt(res)
        num += current_counts[0]
        fnum += current_counts[1]
        if CLEANFLAG:
            os.remove(res)
    np.savetxt('../data/motif/' + Fout_pre + "_counts.txt", np.array([num, fnum]), fmt='%i')

    # combine GI histograms
    for res in glob.glob(TMPDIR+'*_gi.txt'):
        current_gi = np.loadtxt(res)
        gihist += current_gi
        if CLEANFLAG:
            os.remove(res)
    np.savetxt('../data/motif/' + Fout_pre + "_gihist.txt", gihist, fmt='%i')

    # combine all functional motifs into result dataframe
    for res in glob.glob(TMPDIR+'*_topo.txt'):
        current_topodict = pd.read_csv(res, header=0, index_col=False, sep='\t')
        if len(current_topodict) > 0:
            for jx in range(len(current_topodict)):
                current_topo = current_topodict.loc[jx]['topo']
                current_count = current_topodict.loc[jx]['count']
                if current_topo in list(topodict.keys()):
                    topodict[current_topo] += int(current_count)
                else:
                    topodict[current_topo] = 1
        if CLEANFLAG:
            os.remove(res)
    topoDF = pd.DataFrame({'topo': list(topodict.keys()), 'count': list(topodict.values())}  )
    topoDF.to_csv('../data/motif/' + Fout_pre + '_topo.txt', header=True, index=False, sep='\t')
   
    if writeDF:
        results = []
        for f in glob.glob(TMPDIR + '*_motif.txt'):
            df = pd.read_csv(f, header=0, index_col=False, sep='\t')
            results.append(df)
            if CLEANFLAG:
                os.remove(f)
        resultDF = pd.concat(results, ignore_index=True, sort=True)
        resultDF.to_csv('../data/motif/' + Fout_pre +'.txt', header=True, index=False, sep='\t')

    if CLEANFLAG:
        os.rmdir(TMPDIR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_parallel(packaged_input):
    """
    import/export wrapper to run motifs on multi-cores
    """

    # pack/unpack vars for muliprocessing
    source = packaged_input[0]
    G_mat = packaged_input[1]
    GI_mat = packaged_input[2]
    k = packaged_input[3]
    list_cmplx = packaged_input[4]
    mat_cmplx = packaged_input[5]
    combomat = packaged_input[6]
    labels_rev = packaged_input[7]
    theta_func = packaged_input[8]
    writeDF = packaged_input[9]
    TMPDIR = packaged_input[10]

    # set all nodes preceding the source to 0 as already computed
    G_mat[np.arange(source), :] = 0
    G_mat[:, np.arange(source)] = 0
    GI_mat[np.arange(source), :] = 0
    GI_mat[:, np.arange(source)] = 0

    enumerate_kmotif(G_mat, GI_mat, source, k, list_cmplx, mat_cmplx, combomat, labels_rev, theta_func, writeDF, TMPDIR)

    print(source)
 
    return None


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def randomized_network(gmat):
    """
    generate randomized network with same degree distribution
    wrapper for BiRewire R package
    """

    np.savetxt("tmp.inp", gmat, fmt='%i')

    cmd_rand = "Rscript randomize_network.R"
    output = subprocess.run(cmd_rand, shell=True)

    if os.path.exists("tmp.rewired"):
        gmat_rand = np.loadtxt("tmp.rewired")
    else:
        print("it's fucked")

    os.remove("tmp.inp")
    os.remove("tmp.rewired")

    return gmat_rand 







if __name__ == '__main__':


    print("loading network data")
    G_ppi = nx.read_gpickle("../data/processed/yeast_ppi.nx")
    G_gi  = nx.read_gpickle("../data/processed/GI.nx")
    complexes = pd.read_csv("../data/yeast/MIPS.dat", header=0, index_col=False, sep=' ')

    print("filtering network data")
    G_ppi_25, G_gi_25 = filter_degree(G_ppi, G_gi, 25)
    G_ppi_50, G_gi_50 = filter_degree(G_ppi, G_gi, 50)
    G_ppi_100, G_gi_100 = filter_degree(G_ppi, G_gi, 100)
    

    start_time = time.time()
    print("enumerating motifs")

    # note that runnning all parameters takes some time (days/weeks)!
    for k in [3,4,5,6]: 

        print("-> enumerating functional network motifs of size k =", k)

        # discover motifs
        motif_stats = find_motifs(G_ppi_25, G_gi_25, complexes, k, "result_d25_k"+str(k), False, False, 0.5, writeDF=True)
        motif_stats = find_motifs(G_ppi_50, G_gi_50, complexes, k, "result_d50_k"+str(k), False, False, 0.5, writeDF=True)
        motif_stats = find_motifs(G_ppi_100, G_gi_100, complexes, k, "result_d100_k"+str(k), False, False, 0.5, writeDF=False)
  
        #randomization experiments
        for r in range(30):
            motif_stats = find_motifs(G_ppi_50, G_gi_50, complexes, k, "rand"+str(r)+"_result_rand_d50_k"+str(k), True, False, 0.5, writeDF=False)
            motif_stats = find_motifs(G_ppi_50, G_gi_50, complexes, k, "rand"+str(r)+"_result_randgi_d50_k"+str(k), False, True, 0.5, writeDF=True)

    print("--- %s seconds --- " % (time.time() - start_time))

    # only PPI motifs - the counts and randomized counts are computed above, but this function also writes out the enumerated motifs
    motif_stats = find_motifs(G_ppi_50, G_gi_50, complexes, k, "result_d50_all_k"+str(k), False, False, 0, writeDF=True)



    # generate random PPI motifs with some GI content as FNMs
    for k in [3,4,5,6]:
        motif_stats = find_motifs(G_ppi_50, G_gi_50, complexes, k, "revision_rand_result_d50_all_k"+str(k), True, False, 0.5, writeDF=True)

        for r in range(10):
            motif_stats = find_motifs(G_ppi_50, G_gi_50, complexes, k, "revision_rand"+str(r)+"_result_d50_all_k"+str(k), True, False, 0.5, writeDF=True)
