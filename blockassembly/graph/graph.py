# graph.py

from ..common.utils import bytes2numseq, numseq2bytes, check_args, num2seq, print_ref_unitigs, rev_comp

from time import time

import graph_tool.all as gt


def get_kmer_count_from_sequence(sequence, k=3, kmers = None, n_b=2, cyclic=False):
    """
    Returns list of all possible kmers in a sequence with its reverse complement
    """
    # list to store kmers
    if kmers is None:
        kmers = []
    # check all kmers in the sequence and compare to the list
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]
        
        # for cyclic sequence get kmers that wrap from end to beginning
        length = len(kmer)
        if len(kmer) != k:
            if cyclic:
                kmer += sequence[:(k - length)]
            else:
                continue
        
        # cast numseq kmer to bytes and sort it with it reverse complement
        bkmer = numseq2bytes(kmer,n_b=n_b)
        bkmer_r = numseq2bytes(rev_comp(kmer),n_b=n_b)
        bkmer_min, bkmer_max = min(bkmer, bkmer_r), max(bkmer, bkmer_r)
        # only keep the canonical pair
        if (bkmer_min, bkmer_max) not in kmers:
            kmers.append((bkmer_min,bkmer_max))
    return kmers


def get_kmer_count_from_sequences(sequences, k=3, n_b=2, cyclic=True):
    """
    Returns list of all possible kmers in a batch of sequences with its reverse complement
    """
    # list to store kmers
    # TODO change list to set
    kmers = []
    for sequence in sequences:
        kmers = get_kmer_count_from_sequence(sequence, k, kmers, n_b=n_b, cyclic = cyclic)
    return kmers


def add_to_edges(d, edge, edge_type):
    """
    helper function to add edge type edge_type to an edge between two nodes in a dictionnary d
    """
    if edge in d:
        d[edge].append(edge_type)
    else:
        d[edge]=[edge_type]

def get_debruijn_edges_from_kmers(kmers, n_b=2):
    """
    Every possible k-mer is assigned to a node, and we connect one node to another if the k-mers overlaps another.
    Nodes are k-mers, edges are (k-1)-mers.
    There is 4 types of overlaps for 2 given nodes:
    - the end of canonical kmer of node 1 overlaps with the beginning of canonical kmer of node 2 --> Value 1
    - the end of reverse complement of node 1 overlaps with the beginning of reverse complement of node 2 --> Value -1
    - the end of canonical kmer of node 1 overlaps with the beginning of reverse complement of node 2 --> Value 2
    - the end of reverse complement of node 1 overlaps with the beginning of canonical kmer of node 2 --> Value -2
    """
    edges = {}
    for i1,k1 in enumerate(kmers):
        for i2, k2 in enumerate(kmers):
            if i2>=i1:
                if k1[0][n_b:]==k2[0][:-n_b]:
                    add_to_edges(edges,(i1,i2),1) 
                # if self-edge 1 is equivalent to -1
                if k1[1][n_b:]==k2[1][:-n_b] and i1!=i2:
                    add_to_edges(edges,(i1,i2),-1)
                if k1[0][n_b:]==k2[1][:-n_b]:
                    add_to_edges(edges,(i1,i2),2)
                if k1[1][n_b:]==k2[0][:-n_b]:
                    add_to_edges(edges,(i1,i2),-2)
    # having both edges of type 1/-1 and 2/-2 is equivalent to k2=rev_comp(k2) or k1=rev_comp(k1)
    # so only keep the edge e with abs(e) = 1
    # Note : this is only possible for even values of k
    d_edges = edges.copy()
    edges=[]
    for i1,i2 in d_edges:
        es = d_edges[(i1,i2)]
        has1 = False
        for r in [-1,1]:
            if r in es:
                edges.append((i1,i2,r))
                has1 = True
        for r in [-2,2]:
            if r in es and not has1:
                edges.append((i1,i2,r))
    return edges

def create_dbg_from_edges(edges,kmers):
    dbg={i:[[],[],0,0,0,0] for i in range(len(kmers))}
    for i1,i2, r in edges:
        match r:
            case 1:
                dbg[i1][1].append((i2,1))
                dbg[i2][0].append((i1,1))
            case -1:
                dbg[i1][0].append((i2,1))
                dbg[i2][1].append((i1,1))
            case 2:
                dbg[i1][1].append((i2,-1))
                # avoid doubling the same self-mirroring edges
                if i1!=i2:
                    dbg[i2][1].append((i1,-1))
            case -2:
                dbg[i1][0].append((i2,-1))
                # avoid doubling the same self-mirroring edges
                if i1!=i2:
                    dbg[i2][0].append((i1,-1))
    for i in range(len(kmers)):
        dbg[i][2], dbg[i][3] = len(dbg[i][0]), len(dbg[i][1])

    
    
    return dbg

def switch_index(i, mode):
    # return i if mode==1 and the other element of odd/even pair if mode==-1
    # e.g if mode=-1, i=2 return 3 and i=3 return 2
    return 2*(i//2)+(mode==-1)+i%2*mode

def expand_unitigs(unitigs,dbg, kmers, node, mode, n_b=2, start=True):
    # print("expand",node,mode)
    if not start:
        if dbg[node][switch_index(2,mode)]>1 or dbg[node][switch_index(4,mode)]!=0:
            return unitigs,dbg
        else:
            unitigs[-1] = unitigs[-1]+kmers[node][switch_index(0,mode)][-n_b:]
            dbg[node][switch_index(4,mode)]=1
    if dbg[node][switch_index(3,mode)]==1:
        next_node, edge_mode = dbg[node][switch_index(1,mode)][0]
        next_mode = mode*edge_mode
        return expand_unitigs(unitigs,dbg, kmers, next_node, next_mode, n_b=n_b, start=False)
    else:
        return unitigs,dbg

def retrieve_unitigs(unitigs, dbg, kmers, node, mode, n_b=2, root=None):
    # print("retrieve",node,mode)
    if dbg[node][switch_index(4,mode)]!=0:
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
    dbg[node][switch_index(4,mode)]=1
    return expand_unitigs(unitigs, dbg, kmers, node, mode, n_b=n_b)




def get_unitigs_from_dbg(dbg, kmers, n_b=2):
    """
    Return unitigs present in the DBG
    """
    unitigs = []
    for node in dbg:
        unitigs, dbg = retrieve_unitigs(unitigs,dbg,kmers, node, 1, n_b=n_b)
    # Add reverse complement of each unitig and sort to have the canonical in first
        unitigs_set = set()
    for unitig in unitigs:
        bunitig_r = numseq2bytes(rev_comp(bytes2numseq(unitig,n_b)),n_b=n_b)
        bunitig_min, bunitig_max = min(unitig, bunitig_r), max(unitig, bunitig_r)
        unitigs_set.add((bunitig_min,bunitig_max))
    return list(unitigs_set)

def get_compacted_dbg_edges_from_unitigs(unitigs, k,n_b=2):
    # c_edges = [(i1,i2) for i2, u2 in enumerate(unitigs) for i1,u1 in enumerate(unitigs) if len(u1)>=(k*n_b) and len(u2)>=(k*n_b) and u1[-n_b*(k-1):]==u2[:n_b*(k-1)]]
    c_edges = {}
    for i1,u1 in enumerate(unitigs):
        for i2, u2 in enumerate(unitigs):
            if i2>=i1 and len(u1[0])>=(k*n_b) and len(u2[0])>=(k*n_b):
                if u1[0][-n_b*(k-1):]==u2[0][:n_b*(k-1)]:
                    # print(1)
                    add_to_edges(c_edges,(i1,i2),1) 
                if u1[0][-n_b*(k-1):]==u2[1][:n_b*(k-1)]:
                    # print(2)
                    add_to_edges(c_edges,(i1,i2),2)
                if u1[1][-n_b*(k-1):]==u2[1][:n_b*(k-1)] and i1!=i2: # for self edges, 1 is equivalent to -1
                    # print(-1)
                    add_to_edges(c_edges,(i1,i2),-1) 
                if u1[1][-n_b*(k-1):]==u2[0][:n_b*(k-1)]:
                    # print(-2)
                    add_to_edges(c_edges,(i1,i2),-2) 
    
    d_edges = c_edges.copy()
    c_edges=set()
    for e in d_edges:
        for r in d_edges[e]:
            c_edges.add((i1,i2,r))
    return list(c_edges)


def get_gt_graph(edges,sequences):
    g = gt.Graph()
    g.add_vertex(len(sequences))
    edge_source_type=[]
    edge_target_type=[]
    edge_type=[]
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
        edge_type.append(r)
    est = g.new_edge_property("string", vals = edge_source_type)
    ett = g.new_edge_property("string", vals = edge_target_type)
    edge_type = g.new_edge_property("int", vals = edge_type)

    cytoscape_dict = {"arrow":"Arrow","none":"None","square":"Square","diamond":"Diamond"}

    est_cytoscape = g.new_edge_property("string", vals = [cytoscape_dict[e] for e in edge_source_type])
    ett_cytoscape = g.new_edge_property("string", vals = [cytoscape_dict[e] for e in edge_target_type])
    vlen=g.new_vp("int", vals=[len(s[0]) for s in sequences])
    vname=g.new_vp("int",vals= [int(i) for i in g.vertices()])
    g.vp["len"] = vlen
    vseq=g.new_vp("string", vals=[str(s[0]) for s in sequences])
    # print(len(sequences[0][0][0]),sequences[0][0][0])
    sep = "~~~" if len(sequences[0][0][0])>1 else ""
    vseq_dot=g.new_vp("string", vals=[sep.join(s[0]) for s in sequences])
    vrevcomp=g.new_vp("string", vals=[str(s[1]) for s in sequences])
    g.vp["seq"] = vseq
    g.vp["seq_dot"] = vseq_dot
    g.vp["rev_comp_seq"] = vrevcomp
    g.vp["id"] = vname
    g.ep["est"] = est
    g.ep["ett"] = ett
    g.ep["edge_type"] = edge_type
    g.ep["est_cytoscape"] = est_cytoscape
    g.ep["ett_cytoscape"] = ett_cytoscape
    g.vp["out_degree"] = g.new_vp("int", vals=[v.out_degree() for v in g.vertices()])
    g.vp["in_degree"] = g.new_vp("int", vals=[v.in_degree() for v in g.vertices()])
    g.vp["all_degree"] = g.new_vp("int", vals=[v.in_degree()+v.out_degree() for v in g.vertices()])
    
    return g


def graph_multi_k(verbose=True, **kwargs):
    ref_seq, reads, kmin, kmax, bi_alphabet = check_args(verbose=verbose, **kwargs)
    unitigs = []
    n_b = 1
    for k in range(kmin, kmax+1):
        sequences = reads+[seq[0] for seq in unitigs]
        kmers = get_kmer_count_from_sequences(sequences, k=k, n_b=n_b, cyclic=False)
        print([(num2seq(bytes2numseq(seq[0],n_b), bi_alphabet),num2seq(bytes2numseq(seq[1],n_b), bi_alphabet)) for seq in kmers])
        edges = get_debruijn_edges_from_kmers(kmers, n_b=n_b)
        dbg = create_dbg_from_edges(edges, kmers)
        # print(edges)
        # for k in dbg:
        #     print(k,dbg[k][0],dbg[k][1])
        unitigs = get_unitigs_from_dbg(dbg, kmers, n_b=n_b)
        c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k,n_b=n_b)
        unitigs = [(bytes2numseq(u[0],n_b),bytes2numseq(u[1],n_b)) for u in unitigs]
        if verbose:
            print("k={} unitigs: ".format(k), [num2seq(u[0], bi_alphabet) for u in unitigs])
        # c_g = get_gt_graph(c_edges, unitigs)
        # g = get_gt_graph(edges, [bytes2numseq(kmer[0],n_b) for kmer in kmers])
        # res.append((g,c_g))

    reads = [num2seq(seq, bi_alphabet) for seq in reads]
    unitigs = [(num2seq(seq[0], bi_alphabet),num2seq(seq[1], bi_alphabet)) for seq in unitigs]
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


def dbg_tip_clipping(dbg, k, tip_length):
    for node in dbg:
        for mode in [-1,1]:
            if dbg[switch_index(2,mode)]==0:
                dbg = clip_node(dbg, node, mode, tip_length-k)
    return dbg

def clip_node(dbg, node, mode, l):
    if l<0 or dbg[node][switch_index(3,mode)]:
        return dbg
    if dbg[switch_index(2)]==0:
        pass
