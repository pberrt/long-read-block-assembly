# graph.py

from common.utils import bytes2numseq, numseq2bytes, check_args, num2seq, print_ref_unitigs, rev_comp

from data.data import  Unitig, Kmer

from time import time
import numpy as np

import graph_tool.all as gt

from collections import UserDict

class Graph(UserDict):
    def compute_order(self):
        for i,s in enumerate(self.data):
            s.order = i            
    # def edges(self, on_order):
    #     edges = []
    #     if on_order:
    #         self.compute_order()
    #     for s1 in self:
    #         for s2 in s1.link[0]:
    #             v1,v2 = (s1.order, s2.order) if on_order else (s1.id, s2.id)
    #             if v1<= v2:
    #                 edges.append((v1 , v2))
    def compute_edges(self, k):
        edges = {}
        # print("k=",k)
        for unitig in self.data:
            unitig.link=[{},{}]
        for unitig1 in self.data:
            for unitig2 in self.data:
                if unitig2.id>=unitig1.id and len(unitig1.seq)>=(k*unitig1.n_b) and len(unitig2.seq)>=(k*unitig2.n_b):
                    isequal = ( unitig1.seq==unitig1.rev_comp or unitig2.seq==unitig2.rev_comp )
                    show=False
                    # if (unitig1.id,unitig2.id) in [(485,682),(682,2310)]:
                        # print(unitig1.id,unitig2.id)
                        # print("\t",unitig1.num())
                        # print("\t",unitig1.num(canonical=False))
                        # print("\t",unitig2.num())
                        # print("\t",unitig2.num(canonical=False))
                        # print(unitig1.rev_comp[-unitig1.n_b*(k-1):]==unitig2.rev_comp[:unitig2.n_b*(k-1)] , unitig1.can_concatenate[0] , unitig2.can_concatenate[1] , unitig1.id!=unitig2.id)
                        # print(unitig1.rev_comp[-unitig1.n_b*(k-1):],unitig2.rev_comp[:unitig2.n_b*(k-1)] , unitig1.can_concatenate[0] , unitig2.can_concatenate[1] ,unitig1.id,unitig2.id)
                        # print(unitig1.seq[-unitig1.n_b*(k-1):]==unitig2.rev_comp[:unitig2.n_b*(k-1)] , unitig1.can_concatenate[1] , unitig2.can_concatenate[1] , not isequal)
                        # print(unitig1.seq[-unitig1.n_b*(k-1):],unitig2.rev_comp[:unitig2.n_b*(k-1)] , unitig1.can_concatenate[1] , unitig2.can_concatenate[1] , isequal)

                    if unitig1.seq[-unitig1.n_b*(k-1):]==unitig2.seq[:unitig2.n_b*(k-1)] and unitig1.can_concatenate[1] and unitig2.can_concatenate[0]:
                        # print(unitig1.num(),unitig2.num(),1)
                        # print(unitig1.id,unitig2.id,1)
                        if show:
                            print(1)
                        add_to_edges(edges,(unitig1.id,unitig2.id),1)
                        add_to_dict(unitig1.link[1],unitig2,1)
                        add_to_dict(unitig2.link[0],unitig1,1)
                    if unitig1.rev_comp[-unitig1.n_b*(k-1):]==unitig2.rev_comp[:unitig2.n_b*(k-1)] and unitig1.can_concatenate[0] and unitig2.can_concatenate[1] and unitig1.id!=unitig2.id: # for self edges, 1 is equivalent to -1
                        # print(unitig1.num(),unitig2.num(),-1)
                        # print(unitig1.id,unitig2.id,-1)
                        if show:
                            print(-1)
                        add_to_edges(edges,(unitig1.id,unitig2.id),-1)
                        add_to_dict(unitig1.link[0],unitig2,1)
                        add_to_dict(unitig2.link[1],unitig1,1)
                    if unitig1.seq[-unitig1.n_b*(k-1):]==unitig2.rev_comp[:unitig2.n_b*(k-1)] and unitig1.can_concatenate[1] and unitig2.can_concatenate[1] and not isequal:
                        # print(unitig1.num(),unitig2.num(),2)
                        # print(unitig1.id,unitig2.id,2)
                        if show:
                            print(2)
                        add_to_edges(edges,(unitig1.id,unitig2.id),2)
                        add_to_dict(unitig1.link[1],unitig2,-1)
                        # avoid doubling the same self-mirroring edges
                        if unitig1.id!=unitig2.id:
                            add_to_dict(unitig2.link[1],unitig1,-1)
                    if unitig1.rev_comp[-unitig1.n_b*(k-1):]==unitig2.seq[:unitig2.n_b*(k-1)] and unitig1.can_concatenate[0] and unitig2.can_concatenate[0] and not isequal:
                        # print(unitig1.num(),unitig2.num(),-2)
                        # print(unitig1.id,unitig2.id,-2)
                        if show:
                            print(-2)
                        add_to_edges(edges,(unitig1.id,unitig2.id),-2)
                        add_to_dict(unitig1.link[0],unitig2,-1)
                        # avoid doubling the same self-mirroring edges
                        if unitig1.id!=unitig2.id:
                            add_to_dict(unitig2.link[0],unitig1,-1)
        d_edges = edges.copy()
        edges=set()
        for e in d_edges:
            for r in d_edges[e]:
                edges.add((*e,r))
        return list(edges)

    def get_edges(self, on_order = False):
        edges = {}
        if on_order:
            self.compute_order()
        for s1 in self.data:
            v1 = s1.order if on_order else s1.id
            for s2, edge_modes in self.data[s1].link[1].items():
                v2 = s2.order if on_order else s2.id
                if v1<=v2:
                    for edge_mode in edge_modes:
                        add_to_edges(edges,(v1,v2),(edge_mode,1))
            for s2, edge_modes in self.data[s1].link[0].items():
                v2 = s2.order if on_order else s2.id
                if v1<=v2:
                    for edge_mode in edge_modes:
                        add_to_edges(edges,(v1,v2),(edge_mode,-1))
        d_edges = edges.copy()
        edges=[]
        for (i1,i2), l in d_edges.items():
            l = list(dict.fromkeys(l))
            for r in l:
                match r:
                    case (1,1):
                        edges.append((i1,i2,1))
                    case (1,-1):
                        if i1!=i2:
                            edges.append((i1,i2,-1))
                    case (-1,1):
                        edges.append((i1,i2,2))
                    case (-1,-1):
                        edges.append((i1,i2,-2))
        return edges
    def clip(self, k,  n=None, a=None, n_neighbors=1):
        while True:
            n_cut = 0
            # Tag tips
            for s in self.data:
                # print(len(s),n)
                if s.is_tip is None:
                    s.is_tip =False
                # if str(s)== "+mltB~~~+pncC~~~+recA~~~+recX~~~+alaS~~~+yqaB~~~+gshA~~~+luxS~~~-emrB3000074U0009612812615281415420.emrB~~~-emrA3000027AP009048128100822811255411.emrA~~~-emrR3000516U000963281076928113005286.emrR~~~-ygaH~~~-ygaZ~~~-ygaY~~~-proX~~~-proW~~~-proV~~~-nrdF~~~-nrdE~~~-nrdI~~~-ygaM~~~+ygaC~~~+stpA~~~-ygaP~~~-ygaV~~~+kbp~~~-csiR~~~-gabP~~~-gabT~~~-gabD~~~-lhgO~~~-glaH~~~-ygaQ~~~+ypjC~~~+aidA.ypjA":
                #     print(n,len(s), s.abundance, a)
                if (n is None or len(s)<=n) and (a is None or s.abundance<=a):
                    for mode in [1,-1]:
                        if len(s.link[switch_index(0,mode)])==0:
                            s.is_tip = True
                # Remove links if tip is linked to only 1 non-tip
            for s in self.data:
                # if str(s)== "+mltB~~~+pncC~~~+recA~~~+recX~~~+alaS~~~+yqaB~~~+gshA~~~+luxS~~~-emrB3000074U0009612812615281415420.emrB~~~-emrA3000027AP009048128100822811255411.emrA~~~-emrR3000516U000963281076928113005286.emrR~~~-ygaH~~~-ygaZ~~~-ygaY~~~-proX~~~-proW~~~-proV~~~-nrdF~~~-nrdE~~~-nrdI~~~-ygaM~~~+ygaC~~~+stpA~~~-ygaP~~~-ygaV~~~+kbp~~~-csiR~~~-gabP~~~-gabT~~~-gabD~~~-lhgO~~~-glaH~~~-ygaQ~~~+ypjC~~~+aidA.ypjA":
                #     pass
                if s.is_tip:
                    for mode in [1,-1]:
                        n_no_tip = 0
                        for n_s in s.link[switch_index(1,mode)]:
                            if not n_s.is_tip:
                                n_no_tip+=1
                        if n_no_tip<=n_neighbors and s not in s.link[switch_index(1,mode)] and len(s.link[switch_index(1,mode)])>0:
                            n_cut+=1
                            for n_s in s.link[switch_index(1,mode)]:
                                if s in  n_s.link[switch_index(0,mode)]:
                                    _ = n_s.link[switch_index(0,mode)].pop(s)
                            s.link[switch_index(1,mode)].clear()
                            s.can_concatenate[switch_index(1,mode)]=False
            self.data= get_unitigs_bcalm(self, k, on_unitig=True)
            self.compute_edges(k)
            if n_cut==0:
                break
        # print(n_cut)
    def detach_abundance(self, k, a=None):
        for s in self.data:
            if s.abundance<=a:
                for mode in [1,-1]:
                    for n_s in s.link[switch_index(1,mode)]:
                        if s in  n_s.link[switch_index(0,mode)]:
                            _ = n_s.link[switch_index(0,mode)].pop(s)
                    s.link[switch_index(1,mode)].clear()
                    s.can_concatenate[switch_index(1,mode)]=False
        self.compute_edges(k)

            
    

def get_kmer_count_from_sequence(sequence, k=3, kmers = None, cyclic=False,count_key=0):
    """
    Returns list of all possible kmers in a sequence with its reverse complement
    """
    # list to store kmers
    is_unitig = isinstance(sequence, Unitig)
        
    if kmers is None:
        kmers = Graph()
    # check all kmers in the sequence and compare to the list
    length = len(sequence.seq)
    kb = k*sequence.n_b
    for i in range(0, length, sequence.n_b):
        bkmer = sequence.seq[i:i + kb]
        n = len(bkmer)
        # for cyclic sequence get kmers that wrap from end to beginning
        if  n!= kb:
            if cyclic:
                bkmer += sequence.seq[:(kb - n)]
            else:
                continue
        
        # cast numseq kmer to bytes and sort it with it reverse complement
        bkmer_r = rev_comp(bkmer,True, True,n_b=sequence.n_b)
        bkmer_min = min(bkmer, bkmer_r)
        # only keep the canonical pair
        if bkmer_min not in kmers:
            kmer = Kmer(len(kmers),bkmer_min,0)
            kmers[kmer] = kmer
        if is_unitig:
            kmers[bkmer_min].a_unitigs.append(sequence.abundance)
        else:
            kmers[bkmer_min].a_reads += sequence.abundance
    return kmers


def get_kmer_count_from_sequences(sequences, k=3, cyclic=False, kmers=None):
    """
    Returns list of all possible kmers in a batch of sequences with its reverse complement
    """
    if kmers is None:
        kmers=Graph()
    for sequence in sequences:
        kmers = get_kmer_count_from_sequence(sequence, k, kmers, cyclic = cyclic)
    return kmers


def add_to_edges(d, edge, edge_type):
    """
    helper function to add edge type edge_type to an edge between two nodes in a dictionnary d
    """
    if edge in d:
        d[edge].append(edge_type)
    else:
        d[edge]=[edge_type]
def add_to_dict(d, key, e):
    """
    helper function to add element e to a list using its key in dictionnary d
    """
    if key in d:
        d[key].append(e)
    else:
        d[key]=[e]

def get_debruijn_edges_from_kmers(kmers):
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
    for k1 in kmers:
        for k2 in kmers:
            if k2.id>=k1.id:    
                # having both edges of type 1/-1 and 2/-2 is equivalent to k2=rev_comp(k2) or k1=rev_comp(k1)
                # so only keep the edge e with abs(e) = 1
                # Note : this is only possible for even values of k
                if k1.seq[k1.n_b:]==k2.seq[:-k2.n_b]:
                    add_to_edges(edges,(k1.id,k2.id),1)
                    add_to_dict(k1.link[1],k2,1)
                    add_to_dict(k2.link[0],k1,1)
                # if self-edge 1 is equivalent to -1
                if k1.rev_comp[k1.n_b:]==k2.rev_comp[:-k2.n_b] and k1.id!=k2.id:
                    add_to_edges(edges,(k1.id,k2.id),-1)
                    add_to_dict(k1.link[0],k2,1)
                    add_to_dict(k1.link[1],k1,1)
                if k1.seq[k1.n_b:]==k2.rev_comp[:-k2.n_b]:
                    add_to_edges(edges,(k1.id,k2.id),2)
                    add_to_dict(k1.link[1],k2,-1)
                    # avoid doubling the same self-mirroring edges
                    if k1.id!=k2.id:
                        add_to_dict(k2.link[1],k1,-1)
                if k1.rev_comp[k1.n_b:]==k2.seq[:-k2.n_b]:
                    add_to_edges(edges,(k1.id,k2.id),-2)
                    add_to_dict(k1.link[0],k2,-1)
                    # avoid doubling the same self-mirroring edges
                    if k1.id!=k2.id:
                        add_to_dict(k2.link[0],k1,-1)
                    
    d_edges = edges.copy()
    edges=[]
    for i1,i2 in d_edges:
        for r in d_edges[(i1,i2)]:
            edges.append((i1,i2,r))
    return edges

def create_dbg_from_edges(edges,kmers):
    dbg={i:[[],[],0,0,0,0, (kmer[0]==kmer[1]), kmer, kmers[kmer]] for i,kmer in enumerate(kmers)}
    print([dbg[node][7] for node in dbg if dbg[node][6]])
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
        # TODO: update neighbouring count of palindromes 
        dbg[i][2], dbg[i][3] = len(dbg[i][0]), len(dbg[i][1])
    return dbg

def create_edges_from_dbg(dbg):
    edges = {}
    for i1 in dbg:
        for i2, edge_mode in dbg[i1][1]:
            if i1<=i2:
                add_to_edges(edges,(i1,i2),(edge_mode,1))
        for i2, edge_mode in dbg[i1][0]:
            if i1<=i2:
                add_to_edges(edges,(i1,i2),(edge_mode,-1))
    d_edges = edges.copy()
    edges=[]
    for (i1,i2), l in d_edges.items():
        l = list(dict.fromkeys(l))
        for r in l:
            match r:
                case (1,1):
                    edges.append((i1,i2,1))
                case (1,-1):
                    if i1!=i2:
                        edges.append((i1,i2,-1))
                case (-1,1):
                    edges.append((i1,i2,2))
                case (-1,-1):
                    edges.append((i1,i2,-2))
    return edges





def switch_index(i, mode):
    # return i if mode==1 and the other element of odd/even pair if mode==-1
    # e.g if mode=-1, i=2 return 3 and i=3 return 2
    return 2*(i//2)+(mode==-1)+i%2*mode

def expand_unitigs(unitigs,dbg, kmers, node, mode, n_b=2, start=True, verbose = False):
    # TODO : prevent navigating in reverse side of palindromes
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
            dbg[node][switch_index(4,mode*-1)]=1
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
    unitigs = {}
    for i, (u,a) in enumerate(unitigs_dict.items()):
        unitig = Unitig(i,u[0],a)
        unitigs[u]=unitig
    return unitigs


def add_submer(l, e, on_unitig):
    if on_unitig or len(l)==0 or l[-1]!=e:
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
    u = kmers[u_list[0]].seq
    for i in range(1,len(u_list)):
        u += kmers[u_list[i]].seq[-n_b:]
    return u

def get_unitigs_bcalm(kmers, k, on_unitig = False):
    list_lr=[[],[]]
    list(map(lambda kmer: kmer.create_unitig(k),kmers))
    # for kmer in kmers_bcalm:
    #     print(kmer)
    # print([kmer.unitig.abundance for kmer in kmers])
    for kmer in kmers:
        # print(kmer.unitig.a_list)
        # print(kmer)
        # print(kmer.can_concatenate)
        if kmer.can_concatenate[0]:
            # print(0,kmer)
            kp = kmer.seq[:(k-1)*kmer.n_b]
            kp_r = rev_comp(kp,True,True,kmer.n_b)
            # print(kp,kp_r)
            if kp <= kp_r:
                add_submer(list_lr[0],(kp,kmer), on_unitig)
                # print(0)
            if kp_r <= kp:
                add_submer(list_lr[1],(kp_r,kmer), on_unitig)
                # print(1)
        else:
            pass
            # print("not 0")
        if kmer.can_concatenate[1]:
            # print(1,kmer)
            ks = kmer.seq[-((k-1)*kmer.n_b):]
            ks_r = rev_comp(ks,True, True,kmer.n_b)
            # print(ks,ks_r)
            if ks <= ks_r:
                add_submer(list_lr[1],(ks,kmer), on_unitig)
                # print(0)
            if ks_r <= ks:
                add_submer(list_lr[0],(ks_r,kmer), on_unitig)
                # print(1)
        else:
            pass
        # print(list_lr)
            # print("not 1")
    # for kmer in kmers:
    #     print(kmer.num())
        
    # print(list_lr)
    list_lr[0].sort(key=lambda x: (x[0],x[1].id))
    list_lr[1].sort(key=lambda x: (x[0],x[1].id))
    list_lr = [remove_doubles(l) for l in list_lr]
    il, ir = 0, 0
    ll ,lr = len(list_lr[0]),len(list_lr[1])
    while il<ll and ir<lr:
        le, kmer_l = list_lr[0][il]
        re, kmer_r = list_lr[1][ir]

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
            if kmer_l.id!=kmer_r.id and kmer_l.unitig.id!=kmer_r.unitig.id:
                # print(ire,ile)
                if kmer_r.unitig.seq[-(kmer_r.n_b*(k-1)):]!=re:
                    kmer_r.unitig.switch()
                if kmer_l.unitig.seq[:(kmer_l.n_b*(k-1))]!=le:
                    kmer_l.unitig.switch()             
                kmer_r.unitig.seq += kmer_l.unitig.seq[(kmer.n_b*(k-1)):]
                kmer_r.unitig.a_list += kmer_l.unitig.a_list
                kmer_r.unitig.kmers+= kmer_l.unitig.kmers
                kmer_r.unitig.compute_revcomp()
                for kmer in kmer_l.unitig.kmers:
                    kmer.unitig = kmer_r.unitig
                # for kmer in kmers_bcalm:
                #     print(kmer)
            il+=1
            ir+=1

    # for kmer in kmers_bcalm:
    #     print(kmer)
    list(map(lambda kmer: kmer.check_start(), kmers))

    # for kmer in kmers_bcalm:
    #     print(kmer)
    # unitigs = [get_unitig_from_id_list(kmer.unitig, kmers_bcalm, n_b) for kmer in kmers_bcalm if kmer.unitig[0] is not None]
    # print([kmer for kmer in kmers])
    # unitigs = {(kmer.unitig.seq,rev_comp(kmer.unitig.seq,True,True,kmer.n_b)): kmer.unitig for kmer in kmers if kmer.start}
    
    unitigs = Graph({kmer.unitig: kmer.unitig for kmer in kmers if kmer.start})
    list(map(lambda unitig: unitig.compute_revcomp(), unitigs))
    list(map(lambda unitig: unitig.compute_can_concatenate(), unitigs))
    list(map(lambda unitig: unitig.compute_abundance(), unitigs))
    # print(unitigs)
    # unitigs.sort()
    for i,unitig in enumerate(unitigs):
        unitig.id = i
    return unitigs
        

def get_compacted_dbg_edges_from_unitigs(unitigs, k,n_b=2):
    # c_edges = [(i1,i2) for i2, u2 in enumerate(unitigs) for i1,u1 in enumerate(unitigs) if len(u1)>=(k*n_b) and len(u2)>=(k*n_b) and u1[-n_b*(k-1):]==u2[:n_b*(k-1)]]
    c_edges = {}
    for unitig1 in unitigs:
        for unitig2 in unitigs:
            if unitig2.id>=unitig1.id and len(unitig1.seq)>=(k*n_b) and len(unitig2.seq)>=(k*n_b):
                isequal = ( unitig1.seq==unitig1.rev_comp or unitig2.seq==unitig2.rev_comp )
                if unitig1.seq[-n_b*(k-1):]==unitig2.seq[:n_b*(k-1)] and unitig1.can_concatenate[1] and unitig2.can_concatenate[0]:
                    # print(1)
                    add_to_edges(c_edges,(unitig1.id,unitig2.id),1)
                    add_to_dict(unitig1.link[1],unitig2,1)
                    add_to_dict(unitig2.link[0],unitig1,1)
                if unitig1.rev_comp[-n_b*(k-1):]==unitig2.rev_comp[:n_b*(k-1)] and unitig1.can_concatenate[0] and unitig2.can_concatenate[1] and unitig1.id!=unitig2.id: # for self edges, 1 is equivalent to -1
                    # print(-1)
                    add_to_edges(c_edges,(unitig1.id,unitig2.id),-1)
                    add_to_dict(unitig1.link[0],unitig2,1)
                    add_to_dict(unitig1.link[1],unitig1,1)
                if unitig1.seq[-n_b*(k-1):]==unitig2.rev_comp[:n_b*(k-1)] and unitig1.can_concatenate[1] and unitig2.can_concatenate[1] and not isequal:
                    # print(2)
                    add_to_edges(c_edges,(unitig1.id,unitig2.id),2)
                    add_to_dict(unitig1.link[1],unitig2,-1)
                    # avoid doubling the same self-mirroring edges
                    if unitig1.id!=unitig2.id:
                        add_to_dict(unitig2.link[1],unitig1,-1)
                if unitig1.rev_comp[-n_b*(k-1):]==unitig2.seq[:n_b*(k-1)] and unitig1.can_concatenate[0] and unitig2.can_concatenate[0] and not isequal:
                    # print(-2)
                    add_to_edges(c_edges,(unitig1.id,unitig2.id),-2)
                    add_to_dict(unitig1.link[0],unitig2,-1)
                    # avoid doubling the same self-mirroring edges
                    if unitig1.id!=unitig2.id:
                        add_to_dict(unitig2.link[0],unitig1,-1)
    d_edges = c_edges.copy()
    c_edges=set()
    for e in d_edges:
        for r in d_edges[e]:
            c_edges.add((*e,r))
    return list(c_edges)


def get_gt_graph(sequences):
    edges = sequences.get_edges(on_order=True)
    g = gt.Graph()
    g.add_vertex(len(sequences))
    # print(len(sequences))

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
                st, tt = "none","none"
            case -2:
                st, tt = "arrow","arrow"
        edge_source_type.append(st)
        edge_target_type.append(tt)
        edge_type.append(r)
    est = g.new_edge_property("string", vals = edge_source_type)
    ett = g.new_edge_property("string", vals = edge_target_type)
    edge_type = g.new_edge_property("int", vals = edge_type)

    cytoscape_dict = {"arrow":"Arrow","none":"None","square":"Square","diamond":"Diamond"}

    est_cytoscape = g.new_edge_property("string", vals = [cytoscape_dict[e] for e in edge_source_type])
    ett_cytoscape = g.new_edge_property("string", vals = [cytoscape_dict[e] for e in edge_target_type])
    vlen=g.new_vp("int", vals=[len(s) for s in sequences])
    vname=g.new_vp("int",vals= [int(s.id) for s in sequences])
    vstability=g.new_vp("int",vals= [int(s.stability) for s in sequences])
    vorder=g.new_vp("int",vals= [s.order for s in sequences])
    vtopoorder=g.new_vp("int",vals= [s.topo_order for s in sequences])
    g.vp["len"] = vlen
    vseq=g.new_vp("string", vals=[str(s.num()) for s in sequences])
    # print(len(sequences[0][0][0]),sequences[0][0][0])
    vseq_dot=g.new_vp("string", vals=[str(s) for s in sequences])
    vrevcomp=g.new_vp("string", vals=[str(s.num(canonical=False)) for s in sequences])
    g.vp["debug"] = g.new_vp("string", vals=[str(s.debug) for s in sequences])
    g.vp["z_seq"] = vseq
    g.vp["z_seq_dot"] = vseq_dot
    g.vp["z_rev_comp_seq"] = vrevcomp
    g.vp["id"] = vname
    g.vp["stability"] = vstability
    g.vp["order"] = vorder
    g.vp["topo_order"] = vtopoorder
    g.ep["est"] = est
    g.ep["ett"] = ett
    g.ep["edge_type"] = edge_type
    g.ep["est_cytoscape"] = est_cytoscape
    g.ep["ett_cytoscape"] = ett_cytoscape
    g.vp["out_degree"] = g.new_vp("int", vals=[v.out_degree() for v in g.vertices()])
    g.vp["in_degree"] = g.new_vp("int", vals=[v.in_degree() for v in g.vertices()])
    g.vp["all_degree"] = g.new_vp("int", vals=[v.in_degree()+v.out_degree() for v in g.vertices()])
    g.vp["abundance"] = g.new_vp("int", vals = [s.abundance for s in sequences])
    
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
        unitigs = get_unitigs_bcalm(kmers, n_b, bi_alphabet)
        # unitigs = [bytes2numseq(seq,n_b) for seq in unitigs]
        # unitigs = [(num2seq(seq,bi_alphabet)) for seq in unitigs] # num2seq(rev_comp(seq),bi_alphabet)
        # unitigs.sort()
        print(unitigs)
        # unitigs = get_unitigs_from_dbg(dbg, kmers, n_b=n_b)
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


# TODO Clip tip if all but one neighbourgs are tips : first tag the small enough tips ; then discard the one that match these restrictions (only 1 non-tip neighbourg)
def dbg_tip_clipping(dbg, k, tip_length, n_rounds = 1, kmers_bcalm=None):
    for _ in range(n_rounds):
        for node in dbg:
            for mode in [-1,1]:
                # print("clipping" ,node,mode,dbg[node][switch_index(2,mode)])
                if dbg[node][switch_index(2,mode)]==0:
                    verbose = False
                    dbg, _ , kmers_bcalm= old_clip_node(dbg, node, mode, tip_length, kmers_bcalm,verbose)

    return dbg, kmers_bcalm

def old_clip_node(dbg, node, mode, l, kmers_bcalm,verbose):
    if l<0 or (dbg[node][switch_index(3,mode)]!=1 and dbg[node][switch_index(2,mode)]<=1):
        if verbose:
            print("clipping", node,mode,l, "type 1")
        return dbg, False, kmers_bcalm
    if dbg[node][switch_index(2,mode)]>1:
        if verbose:
            print("clipping", node,mode,l, "type 2")
        return dbg, True, kmers_bcalm
    else:
        if verbose:
            print("clipping", node,mode,l, "type 3")
        next_node, edge_mode = dbg[node][switch_index(1,mode)][0]
        next_mode = mode*edge_mode
        dbg, cut, kmers_bcalm = old_clip_node(dbg, next_node, next_mode, l-1,kmers_bcalm,verbose)
        if cut :
            if kmers_bcalm is not None:
                kmers_bcalm[node].can_concatenate[switch_index(1,mode)]=False
            dbg[node][switch_index(1,mode)].remove((next_node,edge_mode))
            dbg[node][switch_index(3,mode)]-=1
            # if self edges to reverse complement, there is only one edge that was already removed in previous step
            if node!=next_node or edge_mode!=-1:
                dbg[next_node][switch_index(0,next_mode)].remove((node, edge_mode))
                dbg[next_node][switch_index(2,next_mode)]-=1
        return dbg, False, kmers_bcalm
    



