# graph.py

from common.utils import bytes2numseq, numseq2bytes, check_args, num2seq, print_ref_unitigs, rev_comp

from time import time
import numpy as np

import graph_tool.all as gt


def get_kmer_count_from_sequence(sequence, k=3, kmers = None, n_b=2, cyclic=False):
    """
    Returns list of all possible kmers in a sequence with its reverse complement
    """
    # list to store kmers
    if kmers is None:
        kmers = {}
    # check all kmers in the sequence and compare to the list
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]
        
        # for cyclic sequence get kmers that wrap from end to beginning
        length = len(kmer)
        if length != k:
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
            kmers[(bkmer_min,bkmer_max)]=1
        else:
            kmers[(bkmer_min,bkmer_max)]+=1
    return kmers


def get_kmer_count_from_sequences(sequences, k=3, n_b=2, cyclic=True):
    """
    Returns list of all possible kmers in a batch of sequences with its reverse complement
    """
    # list to store kmers
    kmers = {}
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
                # having both edges of type 1/-1 and 2/-2 is equivalent to k2=rev_comp(k2) or k1=rev_comp(k1)
                # so only keep the edge e with abs(e) = 1
                # Note : this is only possible for even values of k
                has1 = False
                if k1[0][n_b:]==k2[0][:-n_b]:
                    has1 = True
                    add_to_edges(edges,(i1,i2),1) 
                # if self-edge 1 is equivalent to -1
                if k1[1][n_b:]==k2[1][:-n_b] and i1!=i2:
                    has1 = True
                    add_to_edges(edges,(i1,i2),-1)
                if k1[0][n_b:]==k2[1][:-n_b] and not has1:
                    add_to_edges(edges,(i1,i2),2)
                if k1[1][n_b:]==k2[0][:-n_b] and not has1:
                    add_to_edges(edges,(i1,i2),-2)
    d_edges = edges.copy()
    edges=[]
    for i1,i2 in d_edges:
        for r in d_edges[(i1,i2)]:
            edges.append((i1,i2,r))
    return edges

def create_dbg_from_edges(edges,kmers):
    dbg={i:[[],[],0,0,0,0, kmers[kmer]] for i,kmer in enumerate(kmers)}
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

def expand_unitigs(unitigs,dbg, kmers, node, mode, n_b=2, start=True, verbose = False):
    if verbose:
        print("expand",node,mode)
    if not start:
        if dbg[node][switch_index(2,mode)]>1 or dbg[node][switch_index(4,mode)]!=0:
            return unitigs,dbg
        else:
            kmer = list(kmers.keys())[node]
            unitigs[-1][0] = unitigs[-1][0]+kmer[switch_index(0,mode)][-n_b:]
            unitigs[-1][1].append(kmers[kmer])
            dbg[node][switch_index(4,mode)]=1
    if dbg[node][switch_index(3,mode)]==1:
        next_node, edge_mode = dbg[node][switch_index(1,mode)][0]
        next_mode = mode*edge_mode
        return expand_unitigs(unitigs,dbg, kmers, next_node, next_mode, n_b=n_b, start=False, verbose=verbose)
    else:
        return unitigs,dbg

def retrieve_unitigs(unitigs, dbg, kmers, node, mode, n_b=2, root=None, verbose = False):
    if verbose:
        print("retrieve",node,mode)
    if dbg[node][switch_index(4,mode)]!=0:
        return unitigs, dbg
    if dbg[node][switch_index(2,mode)]==1:
        next_node, edge_mode = dbg[node][switch_index(0,mode)][0]
        next_mode = mode*edge_mode
        if verbose:
            print(next_node,next_mode)
        if dbg[next_node][switch_index(3,next_mode)]==1:
            if root is None:
                return retrieve_unitigs(unitigs, dbg, kmers, next_node, next_mode, n_b=n_b, root=(node,mode), verbose=verbose)
            elif root is not None and (next_node,next_mode)!=root:
                return retrieve_unitigs(unitigs, dbg, kmers, next_node, next_mode, n_b=n_b, root=root, verbose=verbose)
    kmer = list(kmers.keys())[node]
    unitigs.append([kmer[switch_index(0,mode)],[kmers[kmer]]])
    dbg[node][switch_index(4,mode)]=1
    dbg[node][switch_index(4,mode*-1)]=1
    return expand_unitigs(unitigs, dbg, kmers, node, mode, n_b=n_b, verbose=verbose)

def get_unitigs_from_dbg(dbg, kmers, n_b=2, verbose = False):
    """
    Return unitigs present in the DBG
    """
    unitigs = []
    for node in dbg:
        # verbose = True
        unitigs, dbg = retrieve_unitigs(unitigs,dbg,kmers, node, 1, n_b=n_b, verbose = verbose)
        # unitigs, dbg = retrieve_unitigs(unitigs,dbg,kmers, node, -1, n_b=n_b, verbose = verbose)
    # Add reverse complement of each unitig and sort to have the canonical in first
    unitigs_dict = {}
    for unitig, unitig_counts in unitigs:
        bunitig_r = numseq2bytes(rev_comp(bytes2numseq(unitig,n_b)),n_b=n_b)
        bunitig_min, bunitig_max = min(unitig, bunitig_r), max(unitig, bunitig_r)
        if (bunitig_min,bunitig_max) not in unitigs_dict:
            unitig_counts.sort()
            # keep meadian values and lowest of the middle values if even number of values
            unitigs_dict[(bunitig_min,bunitig_max)]=unitig_counts[(len(unitig_counts)-1)//2]
    return unitigs_dict

class Unitig:
    def __init__(self, value):
        self.value=value
        self.kmer =value

class Bcalm_kmer:
    def __init__(self, kmer_b, a, n_b):
        self.kmer=kmer_b
        self.unitig = Unitig(kmer_b)
        self.abundance = a
        self.can_concatenate=[True,True]
        self.n_b = n_b
        self.start=True
    def __repr__(self) -> str:
        return str(self.kmer)+"  "+str(self.unitig.value)+"  "+str(self.start)
    def check_canonical_unitig(self):
        if self.start:
            u_rev = rev_comp(self.unitig.value,True,True, self.n_b)
            if u_rev<self.unitig.value:
                self.unitig.value=u_rev

def add_submer(l,e):
    if len(l)==0 or l[-1]!=e:
        l.append(e)

def expand_same_overlap(i, e,l, mode, kmers):
    l_i=[]
    n = len(l)
    while i<n:
        print("overlap",i,mode)
        e_bis, ei = l[i]
        if e_bis==e:
            if kmers[ei].can_concatenate[switch_index(0,mode)]:
                l_i.append(ei)
            else:
                print("cant")
            i+=1
        else:
            break
    return l_i, i

def remove_doubles(l_double):
    l = []
    if len(l_double)==0:
        return l 
    is_in_double = False
    (s0,i0) = l_double[0]
    for (s,i) in l_double[1:]:
        if s==s0:
            is_in_double=True
        else:
            if not is_in_double:
                l.append((s0,i0))
            s0,i0 = s,i
            is_in_double = False
    if not is_in_double:
        l.append((s0,i0))
    return l

def get_unitig_from_id_list(u_list, kmers, n_b):
    u = kmers[u_list[0]].kmer
    for i in range(1,len(u_list)):
        u += kmers[u_list[i]].kmer[-n_b:]
    return u

def get_unitigs_bcalm(kmers,n_b):
    list_lr=[[],[]]
    kmers_bcalm = [Bcalm_kmer(kmer_b[0], a, n_b) for (kmer_b, a) in kmers.items()]
    k = len(kmers_bcalm[0].kmer)//n_b
    print(kmers_bcalm)
    for i,kmer in enumerate(kmers_bcalm):
        if kmer.can_concatenate[0]:
            kp = kmer.kmer[:(k-1)*n_b]
            kp_r = rev_comp(kp,True,True,n_b)
            if kp <= kp_r:
                add_submer(list_lr[0],(kp,i))
            if kp_r <= kp:
                add_submer(list_lr[1],(kp_r,i))
        if kmer.can_concatenate[1]:
            ks = kmer.kmer[-((k-1)*n_b):]
            ks_r = rev_comp(ks,True, True,n_b)
            if ks <= ks_r:
                add_submer(list_lr[1],(ks,i))
            if ks_r <= ks:
                add_submer(list_lr[0],(ks_r,i))
    list_lr[0].sort()
    list_lr[1].sort()
    print(list_lr)
    list_lr = [remove_doubles(l) for l in list_lr]
    print(list_lr)
    il, ir = 0, 0
    ll ,lr = len(list_lr[0]),len(list_lr[1])

    while il<ll and ir<lr:
        le, ile = list_lr[0][il]
        re, ire = list_lr[1][ir]

        print(il, ir, le,re)
        if le<re:
            il+=1
        elif le>re:
            ir+=1
        else: # le == re
            # ils, il = expand_same_overlap(il, le, list_lr[0],1,kmers_bcalm)
            # irs, ir = expand_same_overlap(ir, re, list_lr[1],-1,kmers_bcalm)
            # print("equal", il, ir, ils,irs)
            # if len(ils)==1 and len(irs)==1:
                # kmers_bcalm[irs[0]].unitig = kmers_bcalm[irs[0]].unitig + kmers_bcalm[ils[0]].unitig
                # kmers_bcalm[ils[0]].unitig = [None,ils[0]]
            if ile!=ire:
                if kmers_bcalm[ire].unitig.value[-(n_b*(k-1)):]!=re:
                    kmers_bcalm[ire].unitig.value = rev_comp(kmers_bcalm[ire].unitig.value,True,True,n_b)
                if kmers_bcalm[ile].unitig.value[:(n_b*(k-1))]!=le:
                    kmers_bcalm[ile].unitig.value = rev_comp(kmers_bcalm[ile].unitig.value,True,True,n_b)               
                kmers_bcalm[ire].unitig.value += kmers_bcalm[ile].unitig.value[(n_b*(k-1)):]
                kmers_bcalm[ile].unitig = kmers_bcalm[ire].unitig
                if kmers_bcalm[ile].kmer != kmers_bcalm[ire].unitig.kmer:
                    kmers_bcalm[ile].start =False
                


            il+=1
            ir+=1
    list(map(lambda kmer: kmer.check_canonical_unitig(), kmers_bcalm))
    print(kmers_bcalm)
    # unitigs = [get_unitig_from_id_list(kmer.unitig, kmers_bcalm, n_b) for kmer in kmers_bcalm if kmer.unitig[0] is not None]   
    unitigs = [kmer.unitig.value for kmer in kmers_bcalm if kmer.start]
    # unitigs.sort()
    return unitigs
        

def get_compacted_dbg_edges_from_unitigs(unitigs, k,n_b=2):
    # c_edges = [(i1,i2) for i2, u2 in enumerate(unitigs) for i1,u1 in enumerate(unitigs) if len(u1)>=(k*n_b) and len(u2)>=(k*n_b) and u1[-n_b*(k-1):]==u2[:n_b*(k-1)]]
    c_edges = {}
    for i1,u1 in enumerate(unitigs):
        for i2, u2 in enumerate(unitigs):
            if i2>=i1 and len(u1[0])>=(k*n_b) and len(u2[0])>=(k*n_b):
                isequal = ( u1[0]==u1[1] or u2[0]==u2[1] )
                if u1[0][-n_b*(k-1):]==u2[0][:n_b*(k-1)]:
                    # print(1)
                    add_to_edges(c_edges,(i1,i2),1) 
                if u1[1][-n_b*(k-1):]==u2[1][:n_b*(k-1)] and i1!=i2: # for self edges, 1 is equivalent to -1
                    # print(-1)
                    add_to_edges(c_edges,(i1,i2),-1)
                if u1[0][-n_b*(k-1):]==u2[1][:n_b*(k-1)] and not isequal:
                    # print(2)
                    add_to_edges(c_edges,(i1,i2),2) 
                if u1[1][-n_b*(k-1):]==u2[0][:n_b*(k-1)] and not isequal:
                    # print(-2)
                    add_to_edges(c_edges,(i1,i2),-2) 
    
    d_edges = c_edges.copy()
    c_edges=set()
    for e in d_edges:
        for r in d_edges[e]:
            c_edges.add((*e,r))
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
    vlen=g.new_vp("int", vals=[len(s[0]) for s in sequences.values()])
    vname=g.new_vp("int",vals= [int(i) for i in g.vertices()])
    g.vp["len"] = vlen
    vseq=g.new_vp("string", vals=[str(s[0]) for s in sequences.values()])
    # print(len(sequences[0][0][0]),sequences[0][0][0])
    sep = "~~~" if len(sequences[list(sequences.keys())[0]][0][0])>1 else ""
    vseq_dot=g.new_vp("string", vals=[sep.join(s[0]) for s in sequences.values()])
    vrevcomp=g.new_vp("string", vals=[str(s[1]) for s in sequences.values()])
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
    g.vp["abundance"] = g.new_vp("int", vals = [s[2] for s in sequences.values()])
    
    return g


def graph_multi_k(verbose=True, **kwargs):
    ref_seq, reads, kmin, kmax, bi_alphabet = check_args(verbose=verbose, **kwargs)
    unitigs = {}
    n_b = 1
    for k in range(kmin, kmax+1):
        sequences = reads+[seq[0] for seq in unitigs.values()]
        kmers = get_kmer_count_from_sequences(sequences, k=k, n_b=n_b, cyclic=False)
        print([(num2seq(bytes2numseq(seq0,n_b), bi_alphabet),num2seq(bytes2numseq(seq1,n_b), bi_alphabet)) for (seq0, seq1), a  in kmers.items()])
        edges = get_debruijn_edges_from_kmers(kmers, n_b=n_b)
        dbg = create_dbg_from_edges(edges, kmers)
        print(dbg)
        # dbg = dbg_tip_clipping(dbg,k,5)
        print(dbg)
        # print(edges)
        # for k in dbg:
        #     print(k,dbg[k][0],dbg[k][1])
        unitigs = get_unitigs_bcalm(kmers, n_b)
        unitigs = [bytes2numseq(seq,n_b) for seq in unitigs]
        unitigs = [(num2seq(seq,bi_alphabet)) for seq in unitigs] # num2seq(rev_comp(seq),bi_alphabet)
        unitigs.sort()
        print(unitigs)
        unitigs = get_unitigs_from_dbg(dbg, kmers, n_b=n_b)
        c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k,n_b=n_b)
        print(c_edges)
        unitigs = {seq[0]: (bytes2numseq(seq[0],n_b),bytes2numseq(seq[1],n_b),a) for seq, a in unitigs.items()}
        if verbose:
            redable_unitigs = [num2seq(u[0], bi_alphabet) for u in unitigs.values()]
            redable_unitigs.sort()
            print("k={} unitigs: ".format(k), redable_unitigs)
        # c_g = get_gt_graph(c_edges, unitigs)
        # g = get_gt_graph(edges, [bytes2numseq(kmer[0],n_b) for kmer in kmers])
        # res.append((g,c_g))

    reads = [num2seq(seq, bi_alphabet) for seq in reads]
    unitigs = {ub: (num2seq(seq[0], bi_alphabet),num2seq(seq[1], bi_alphabet), seq[2]) for ub, seq in unitigs.items()}
    kmers = {seq[0]: (num2seq(bytes2numseq(seq[0],n_b), bi_alphabet),num2seq(bytes2numseq(seq[1],n_b), bi_alphabet),a) for seq, a in kmers.items()}    
    # print(kmers)
    g = get_gt_graph(edges,kmers)
    c_g = get_gt_graph(c_edges, unitigs)

    if verbose:
        print_ref_unitigs(ref_seq, unitigs)
        print("MULTI-K :")
        print("\t", ref_seq)
        print("\t", unitigs)
    return ref_seq, reads, kmers, g, unitigs, c_g


def dbg_tip_clipping(dbg, k, tip_length, n_rounds = 1):
    for _ in range(n_rounds):
        for node in dbg:
            for mode in [-1,1]:
                # print("clipping" ,node,mode,dbg[node][switch_index(2,mode)])
                if dbg[node][switch_index(2,mode)]==0:
                    verbose = False
                    dbg, _ = clip_node(dbg, node, mode, tip_length+k,verbose)
                    pass

    return dbg

def clip_node(dbg, node, mode, l,verbose):
    if l<0 or (dbg[node][switch_index(3,mode)]!=1 and dbg[node][switch_index(2,mode)]<=1):
        if verbose:
            print("clipping", node,mode,l, "type 1")
        return dbg, False
    if dbg[node][switch_index(2,mode)]>1:
        if verbose:
            print("clipping", node,mode,l, "type 2")
        return dbg, True
    else:
        if verbose:
            print("clipping", node,mode,l, "type 3")
        next_node, edge_mode = dbg[node][switch_index(1,mode)][0]
        next_mode = mode*edge_mode
        dbg, cut = clip_node(dbg, next_node, next_mode, l-1,verbose)
        if cut :
            dbg[node][switch_index(1,mode)].remove((next_node,edge_mode))
            dbg[node][switch_index(3,mode)]-=1
            # if self edges to reverse complement, there is only one edge that was already removed in previous step
            if node!=next_node or edge_mode!=-1:
                dbg[next_node][switch_index(0,next_mode)].remove((node, edge_mode))
                dbg[next_node][switch_index(2,next_mode)]-=1
        return dbg, False 
