# graph.py

from utils import bytes2numseq, numseq2bytes, check_args, num2seq, print_ref_unitigs
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

def get_kmer_count_from_sequence(sequence, k=3, kmers = None, cyclic=False):
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
        bkmer = numseq2bytes(kmer)
        if bkmer not in kmers:
            kmers.append(bkmer)
        #     kmers[bkmer] += 1
        # else:
        
    return kmers


def get_kmer_count_from_sequences(sequences, k=3, cyclic=True):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # list to store kmers
    kmers = []

    for sequence in sequences:
        kmers = get_kmer_count_from_sequence(sequence, k, kmers, cyclic = cyclic)
    return kmers
    


def expand_unitigs(unitigs,dbg, kmers, node, start=True):
    if not start:
        if dbg[node][2]>1 or dbg[node][4]!=0:
            return unitigs,dbg
        else:
            unitigs[-1] = unitigs[-1]+kmers[node][-1:]
            dbg[node][4]=1
    if dbg[node][3]==1:
        return expand_unitigs(unitigs,dbg, kmers, dbg[node][1][0], False)
    else:
        return unitigs,dbg

def retrieve_unitigs(unitigs, dbg, kmers, node, root=None):
    if dbg[node][4]!=0:
        return unitigs, dbg
    if root is None:
        if dbg[node][2]==1:
            return retrieve_unitigs(unitigs, dbg, kmers, dbg[node][0][0], root=node)
        else:
            unitigs.append(kmers[node])
            dbg[node][4]=1
            return expand_unitigs(unitigs, dbg, kmers, node)
    else:
        if dbg[node][2]==1 and dbg[dbg[node][0][0]][3]==1 and dbg[node][0][0]!=root:
            return retrieve_unitigs(unitigs, dbg, kmers, dbg[node][0][0], root=root)
        else:
            unitigs.append(kmers[node])
            dbg[node][4]=1
            return expand_unitigs(unitigs, dbg, kmers, node)


def get_debruijn_edges_from_kmers(kmers, nbytes=1):
    """
    Every possible k-mer is assigned
    to a node, and we connect one node to another if the k-mers overlaps 
    another. Nodes are k-mers, edges are (k+1)-mers.
    """
    k = len(kmers[0])//nbytes
    edges = [(i1,i2) for i2, k2 in enumerate(kmers) for i1,k1 in enumerate(kmers) if k1[1:]==k2[:-1]]
    
    
    t= time()
    dbg={i:[[],[],0,0,0,False] for i in range(len(kmers))}
    for i1,i2 in edges:
        dbg[i1][1].append(i2)
        dbg[i2][0].append(i1)
    for i in range(len(kmers)):
        dbg[i][2], dbg[i][3] = len(dbg[i][0]), len(dbg[i][1])
        dbg[i][5] = (dbg[i][2]!=1 or dbg[i][3]!=1)
    return edges, dbg

def get_unitigs_from_dbg(dbg, kmers):
    unitigs = []
    for node in dbg:
        unitigs, dbg = retrieve_unitigs(unitigs,dbg,kmers, node)
    return unitigs

def get_compacted_dbg_edges_from_unitigs(unitigs, k):
    c_edges = [(i1,i2) for i2, u2 in enumerate(unitigs) for i1,u1 in enumerate(unitigs) if u1[-(k-1):]==u2[:(k-1)]]
    return c_edges


def get_gt_graph(edges,sequences):
    g = gt.Graph()
    g.add_vertex(len(sequences))
    for e in edges:
        g.add_edge(*e)
    return g


def graph_multi_k(verbose=True, **kwargs):
    ref_seq, reads, kmin, kmax, bi_alphabet = check_args(verbose=verbose, **kwargs)
    unitigs = []
    for k in range(kmin, kmax+1):
        sequences = reads+unitigs
        kmers = get_kmer_count_from_sequences(sequences, k=k,cyclic=False)
        edges, dbg  = get_debruijn_edges_from_kmers(kmers, nbytes=1)
        unitigs = get_unitigs_from_dbg(dbg, kmers)
        c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k)
        if verbose:
            print("k={} unitigs: ".format(k), [num2seq(bytes2numseq(u), bi_alphabet) for u in unitigs])

    unitigs = [bytes2numseq(u) for u in unitigs]

    reads = [num2seq(seq, bi_alphabet) for seq in reads]
    unitigs = [num2seq(seq, bi_alphabet) for seq in unitigs]
    kmers = [num2seq(seq, bi_alphabet) for seq in kmers]    

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