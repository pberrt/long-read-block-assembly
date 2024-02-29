# graph.py

from ..common.utils import bytes2numseq, numseq2bytes, check_args, num2seq, print_ref_unitigs, rev_comp

from time import time

import graph_tool.all as gt

def equal_numseq(ns1, ns2):
    if len(ns1)!=len(ns2):
        return False
    for c1,c2 in zip(ns1,ns2):
        if c1!=c2:
            return False
    return True

def numseq_in_list(ns, l):
    for ll in l:
        if equal_numseq(ns, ll):
            return True
    return False

def get_kmer_count_from_sequence(sequence, k=3, kmers = None, n_b=2, cyclic=False):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # dict to store kmers
    if kmers is None:
        kmers = []
    # count how many times each occurred in this sequence (treated as cyclic)
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]
        
        # for cyclic sequence get kmers that wrap from end to beginning
        length = len(kmer)
        if cyclic:
            if len(kmer) != k:
                kmer += sequence[:(k - length)]
        
        # if not cyclic then skip kmers at end of sequence
        else:
            if len(kmer) != k:
                continue
        
        # count occurrence of this kmer in sequence
        bkmer = numseq2bytes(kmer,n_b=n_b)
        bkmer_r = numseq2bytes(rev_comp(kmer),n_b=n_b)
        bkmer_min, bkmer_max = min(bkmer, bkmer_r), max(bkmer, bkmer_r)
        if (bkmer_min, bkmer_max) not in kmers:
            kmers.append((bkmer_min,bkmer_max))

        
    return kmers


def get_kmer_count_from_sequences(sequences, k=3, n_b=2, cyclic=True):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # list to store kmers
    kmers = []
    for sequence in sequences:
        kmers = get_kmer_count_from_sequence(sequence, k, kmers, n_b=n_b, cyclic = cyclic)
    return kmers


def add_to_edges(d, key, n):
    if key in d:
        d[key].append(n)
    else:
        d[key]=[n]
def get_debruijn_edges_from_kmers(kmers, n_b=2):
    """
    Every possible k-mer is assigned
    to a node, and we connect one node to another if the k-mers overlaps 
    another. Nodes are k-mers, edges are (k+1)-mers.
    """
    # edges = [(i1,i2,1) for i2, k2 in enumerate(kmers) for i1,k1 in enumerate(kmers) if k1[0][n_b:]==k2[0][:-n_b]]
    edges = {}
    for i1,k1 in enumerate(kmers):
        for i2, k2 in enumerate(kmers):
            if i2>=i1:
                # having both edges of type 1 and 2 is equivalent to k2=rev_comp(k2) so only keep one edge (only possible for even values of k)
                if k1[0][n_b:]==k2[0][:-n_b]:
                    add_to_edges(edges,(i1,i2),1) 
                if k1[0][n_b:]==k2[1][:-n_b]:
                    add_to_edges(edges,(i1,i2),2)
                # having both edges of type -1 and -2 is equivalent to k2=rev_comp(k2) so only keep one edge (only possible for even values of k)
                if k1[1][n_b:]==k2[1][:-n_b]:
                    add_to_edges(edges,(i1,i2),-1)
                if k1[1][n_b:]==k2[0][:-n_b]:
                    add_to_edges(edges,(i1,i2),-2)

    # for i, (i1,i2,r) in enumerate(edges):
    #     try:
    #         edges_r.remove((i1,i2,-1))
    #     except ValueError:
    #         pass
    #     else:
    #         edges[i]=(i1,i2,0)
    # for i, (i1,i2,r) in enumerate(edges_r):
    #     if i1!=i2:
    #         try:
    #             edges_r.remove((i2,i1,-1))
    #         except ValueError:
    #             pass
    # edges = edges + edges_r
    
    
    dbg={i:[[],[],0,0,0,False] for i in range(len(kmers))}
    # for i1,i2, r in edges:
    #     dbg[i1][1].append((i2,r))
    #     dbg[i2][0].append((i1,r))
    # for i in range(len(kmers)):
    #     dbg[i][2], dbg[i][3] = len(dbg[i][0]), len(dbg[i][1])
    #     dbg[i][5] = (dbg[i][2]!=1 or dbg[i][3]!=1)

    for i1,i2 in edges:
        for r in edges[(i1,i2)]:
            match r:
                case 1:
                    dbg[i1][1].append((i2,1))
                    dbg[i2][0].append((i1,1))
                case -1:
                    dbg[i1][0].append((i2,1))
                    dbg[i2][1].append((i1,1))
                case 2:
                    dbg[i1][1].append((i2,-1))
                    # avoid doubling the same edge
                    # if i1!=i2:
                    dbg[i2][1].append((i1,-1))
                case -2:
                    dbg[i1][0].append((i2,-1))
                    # avoid doubling the same edge
                    # if i1!=i2:
                    dbg[i2][0].append((i1,-1))
    for i in range(len(kmers)):
        dbg[i][2], dbg[i][3] = len(dbg[i][0]), len(dbg[i][1])

    d_edges = edges.copy()
    edges=[]
    for i1,i2 in d_edges:
        print(d_edges[(i1,i2)],(i1,i2),kmers[i1], kmers[i2])
        for r in d_edges[(i1,i2)]:
            edges.append((i1,i2,r))
    return edges, dbg

def switch_index(i, mode):
    # return i if mode==1 and the other element of odd/even pair if mode==-1
    # e.g if mode=-1, i=2 return 3 and i=3 return 2
    return 2*(i//2)+(mode==-1)+i%2*mode

def expand_unitigs(unitigs,dbg, kmers, node, mode, n_b=2, start=True):
    # print("expand",node,mode)
    if not start:
        if dbg[node][switch_index(2,mode)]>1 or dbg[node][4]!=0:
            return unitigs,dbg
        else:
            unitigs[-1] = unitigs[-1]+kmers[node][switch_index(0,mode)][-n_b:]
            dbg[node][4]=1
    if dbg[node][switch_index(3,mode)]==1:
        next_node, edge_mode = dbg[node][switch_index(1,mode)][0]
        next_mode = mode*edge_mode
        return expand_unitigs(unitigs,dbg, kmers, next_node, next_mode, n_b=n_b, start=False)
    else:
        return unitigs,dbg

def retrieve_unitigs(unitigs, dbg, kmers, node, mode, n_b=2, root=None):
    # print("retrieve",node,mode)
    if dbg[node][4]!=0:
        return unitigs, dbg
    if dbg[node][switch_index(2,mode)]==1:
        next_node, edge_mode = dbg[node][switch_index(0,mode)][0]
        next_mode = mode*edge_mode
        # print(next_node,next_mode)
        if dbg[next_node][switch_index(3,next_mode)]==1:
            if root is None:
                return retrieve_unitigs(unitigs, dbg, kmers, next_node, next_mode, n_b=n_b, root=(node,mode))
            elif root is not None and (next_node,next_mode)!=root:
                return retrieve_unitigs(unitigs, dbg, kmers, next_node, next_mode, n_b=n_b, root=root)
    unitigs.append(kmers[node][switch_index(0,mode)])
    dbg[node][4]=1
    return expand_unitigs(unitigs, dbg, kmers, node, mode, n_b=n_b)




def get_unitigs_from_dbg(dbg, kmers, n_b=2):
    unitigs = []
    for node in dbg:
        unitigs, dbg = retrieve_unitigs(unitigs,dbg,kmers, node, 1, n_b=n_b)
    for k, unitig in enumerate(unitigs):
        bunitig_r = numseq2bytes(rev_comp(bytes2numseq(unitig,n_b)),n_b=n_b)
        bunitig_min, bunitig_max = min(unitig, bunitig_r), max(unitig, bunitig_r)
        if (bunitig_min, bunitig_max) not in unitigs:
            unitigs[k]=(bunitig_min,bunitig_max)
    return unitigs

def get_compacted_dbg_edges_from_unitigs(unitigs, k,n_b=2):
    # c_edges = [(i1,i2) for i2, u2 in enumerate(unitigs) for i1,u1 in enumerate(unitigs) if len(u1)>=(k*n_b) and len(u2)>=(k*n_b) and u1[-n_b*(k-1):]==u2[:n_b*(k-1)]]
    c_edges = []
    for i1,u1 in enumerate(unitigs):
        for i2, u2 in enumerate(unitigs):
            if i2>=i1 and len(u1[0])>=(k*n_b) and len(u2[0])>=(k*n_b):
                print(i1,i2)
                # having both edges of type 1 and 2 is equivalent to k2=rev_comp(k2) so only keep one edge (only possible for even values of k)
                if u1[0][-n_b*(k-1):]==u2[0][:n_b*(k-1)]:
                    print(1)
                    c_edges.append((i1,i2,1)) 
                if u1[0][-n_b*(k-1):]==u2[1][:n_b*(k-1)]:
                    print(2)
                    print(u1,u2)
                    print(u1[0][-n_b*(k-1):],u2[1][:n_b*(k-1)])
                    c_edges.append((i1,i2,2))
                # having both edges of type -1 and -2 is equivalent to k2=rev_comp(k2) so only keep one edge (only possible for even values of k)
                if u1[1][-n_b*(k-1):]==u2[1][:n_b*(k-1)]:
                    print(-1)
                    c_edges.append((i1,i2,-1))
                if u1[1][-n_b*(k-1):]==u2[0][:n_b*(k-1)]:
                    print(-2)
                    print(u1,u2)
                    print(u1[1][-n_b*(k-1):],u2[0][:n_b*(k-1)])
                    c_edges.append((i1,i2,-2))
    return c_edges


def get_gt_graph(edges,sequences):
    g = gt.Graph()
    g.add_vertex(len(sequences))
    edge_source_type=[]
    edge_target_type=[]
    for i1,i2, r in edges:
        g.add_edge(i1,i2)
        match r:
            case 1:
                st, tt = "none","arrow"
            case -1:
                st, tt = "arrow","none"
            case 2:
                st, tt = "square","square"
            case -2:
                st, tt = "diamond","diamond"
        edge_source_type.append(st)
        edge_target_type.append(tt)
    est = g.new_edge_property("string", vals = edge_source_type)
    ett = g.new_edge_property("string", vals = edge_target_type)

    cytoscape_dict = {"arrow":"Arrow","none":"None","square":"Square","diamond":"Diamond"}

    est_cytoscape = g.new_edge_property("string", vals = [cytoscape_dict[e] for e in edge_source_type])
    ett_cytoscape = g.new_edge_property("string", vals = [cytoscape_dict[e] for e in edge_target_type])
    vlen=g.new_vp("int", vals=[len(s) for s in sequences])
    g.vp["len"] = vlen
    vname=g.new_vp("string", vals=[str(s) for s in sequences])
    g.vp["name"] = vname
    g.vp["n"] = g.vertex_index
    g.ep["est"] = est
    g.ep["ett"] = ett
    g.ep["est_cytoscape"] = est_cytoscape
    g.ep["ett_cytoscape"] = ett_cytoscape
    return g


def graph_multi_k(verbose=True, **kwargs):
    ref_seq, reads, kmin, kmax, bi_alphabet = check_args(verbose=verbose, **kwargs)
    unitigs = []
    n_b = 1
    for k in range(kmin, kmax+1):
        sequences = reads+unitigs
        kmers = get_kmer_count_from_sequences(sequences, k=k, n_b=n_b, cyclic=False)
        # print([(num2seq(bytes2numseq(seq[0],n_b), bi_alphabet),num2seq(bytes2numseq(seq[1],n_b), bi_alphabet)) for seq in kmers])
        edges, dbg  = get_debruijn_edges_from_kmers(kmers, n_b=n_b)
        # print(edges)
        # for k in dbg:
        #     print(k,dbg[k][0],dbg[k][1])
        unitigs = get_unitigs_from_dbg(dbg, kmers, n_b=n_b)
        c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k,n_b=n_b)
        unitigs = [bytes2numseq(u[0],n_b) for u in unitigs]
        if verbose:
            print("k={} unitigs: ".format(k), [num2seq(u, bi_alphabet) for u in unitigs])
        # c_g = get_gt_graph(c_edges, unitigs)
        # g = get_gt_graph(edges, [bytes2numseq(kmer[0],n_b) for kmer in kmers])
        # res.append((g,c_g))

    reads = [num2seq(seq, bi_alphabet) for seq in reads]
    unitigs = [num2seq(seq, bi_alphabet) for seq in unitigs]
    kmers = [(num2seq(bytes2numseq(seq[0],n_b), bi_alphabet),num2seq(bytes2numseq(seq[1],n_b), bi_alphabet)) for seq in kmers]    
    # print(kmers)
    g = get_gt_graph(edges,kmers)
    c_g = get_gt_graph(c_edges, unitigs)

    if verbose:
        print_ref_unitigs(ref_seq, unitigs)
        print("MULTI-K :")
        print("\t", ref_seq)
        print("\t", unitigs)
    return ref_seq, reads, kmers, g, unitigs, c_g

# def get_debruijn_edges_from_kmers(kmers, nbytes=1):
#     """
#     Every possible k-mer is assigned
#     to a node, and we connect one node to another if the k-mers overlaps 
#     another. Nodes are k-mers, edges are (k+1)-mers.
#     """
#     k = len(kmers[0])//nbytes
#     edges = [(i1,i2) for i2, k2 in enumerate(kmers) for i1,k1 in enumerate(kmers) if (i1!=i2 and k1[1:]==k2[:-1])]

#     dbg={i:[[],[],0,0,0,False] for i in range(len(kmers))}

#     for i1,i2 in edges:
#         dbg[i1][1].append(i2)
#         dbg[i2][0].append(i1)
#     for i in range(len(kmers)):
#         dbg[i][2], dbg[i][3] = len(dbg[i][0]), len(dbg[i][1])
#         dbg[i][5] = (dbg[i][2]!=1 or dbg[i][3]!=1)

#     unitigs = []
#     for node in dbg:
#         unitigs, dbg = retrieve_unitigs(unitigs,dbg,kmers, node)
    
#     c_edges = set()
#     for k1,u1 in enumerate(unitigs):
#         for k2,u2 in enumerate(unitigs):
#             if u1[-(k-1):] == u2[:(k-1)]:
#                 c_edges.add((k1,k2))
    
#     return edges, dbg, unitigs, c_edges