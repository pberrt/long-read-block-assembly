# script_graph_multi_k.py

import os
import sys
from time import time

import json

import matplotlib.pyplot as plt
import numpy as np

import gfapy

sys.path.append('..')

from blockassembly.common.utils import bytes2numseq, numseq2bytes, seq2num, num2seq
from blockassembly.graph.graph import graph_multi_k
from blockassembly.graph.graph import get_kmer_count_from_sequences, get_debruijn_edges_from_kmers, get_unitigs_from_dbg, get_compacted_dbg_edges_from_unitigs, get_gt_graph, create_dbg_from_edges
from blockassembly.data.visu import plot_debruijn_graph_gt, plot_compacted_debruijn_graph_gt
from blockassembly.data.inout import create_gfa

if __name__ == '__main__':
    
    PROJECT_FOLDER = os.getcwd()
    RES_OUTPUT_FOLDER = os.path.join(PROJECT_FOLDER,"..","res","ecoli")
    
    
    sys. setrecursionlimit(10000)
    
    ref_file = "../input/truth_data/GCA_027944875.1_ASM2794487v1_genomic.truth_genes.json"
    with open(ref_file, 'r') as f:
        ref_data = json.load(f)
    
    read_file = "../input/truth_data/SRR23044204_1.subset.truth_genes.json"
    # read_file = "input/real_data/SRR23044204_1.subset.pandora_block_calls.json"
    with open(read_file, 'r') as f:
        read_data = json.load(f)

    blocks = set()
    for k,g in ref_data.items():
        for block in g:
            blocks.add(block[1:])
    for k,g in read_data.items():
        for block in g:
            blocks.add(block[1:])

    # blocks = dict([(p1,p2+1) for p1,p2 in zip(list(blocks), list(range(len(blocks))))])
    alphabet = [("+"+p1,"-"+p1) for p1 in list(blocks)]
    bi_alphabet = (alphabet,{k:((i+1)*((-1)**j)) for i, ks in enumerate(alphabet) for j, k in enumerate(ks)})

    ref_seqs = [seq2num(seq,bi_alphabet) for seq in ref_data.values()]
    # for k,g in ref_data.items():
    #     a = np.zeros(len(g),dtype = np.int16)
    #     for i,block in enumerate(g):
    #         a[i]=np.array(int(str(block[0])+str(blocks[block[1:]]))).astype(a.dtype)
    #     ref_seqs.append(a)
    
    read_data_trimmed = [["_".join(block.split("_")[:-1]) for block in seq] for seq in read_data.values()]
    read_seqs = [seq2num(seq,bi_alphabet) for seq in read_data_trimmed]
    # seqs = []
    # for k,g in read_data.items():
    #     a = np.zeros(len(g),dtype = np.int16)
    #     for i,block in enumerate(g):
    #         a[i]=np.array(int(str(block[0])+str(blocks["_".join(block.split("_")[:-1])[1:]]))).astype(a.dtype)
    #     seqs.append(a)

    read_length=[len(a) for a in read_seqs]
    plt.hist(read_length,bins=np.arange(0.5,max(read_length)+0.5))


    ref_read_sets = [set(),set()]
    seq_sets = []
    for rs in ref_seqs:
        s = set()
        for c in rs:
            s.add(abs(c))
            ref_read_sets[0].add(abs(c))
        seq_sets.append(s)
    for i, s1 in enumerate(seq_sets):
        print("Sequence {} contains {} different element(s)".format(i,len(s1)))
        for j, s2 in enumerate(seq_sets):
            if i!=j:
                inter, diff1, diff2 = s1.intersection(s2), s1.difference(s2), s2.difference(s1)
                linter = len(inter)
                if linter > 0:
                    print("\tShares {} elements with sequence {}: {}".format(linter,j,inter))
                if len(diff1) == 0:
                    print("{} is included in {}".format(i,j)) 
                if len(diff2) == 0:
                    print("{} is included in {}".format(j,i))
        if len(s1)<10:
            print("\t",s1)
    for rs in read_seqs:
        for c in rs:
            ref_read_sets[1].add(abs(c))

    print("The {} reference sequences contains {} unique elements.".format(len(ref_seqs),len(ref_read_sets[0])))
    print("The {} read sequences contains {} unique elements.".format(len(read_seqs),len(ref_read_sets[1])))
    inter = ref_read_sets[0].intersection(ref_read_sets[1])
    diff1 = ref_read_sets[0].difference(ref_read_sets[1])
    diff2 = ref_read_sets[1].difference(ref_read_sets[0])
    print("{} elements are common to both reference sequences and reads.".format(len(inter)))
    print("{} elements are only present in reference sequences:".format(len(diff1)))
    print("\t",*diff1)
    print("{} elements are only present in reads:".format(len(diff2)))
    print("\t",*diff2)

    match read_seqs[0].dtype:
        case np.int8:
            n_b = 1
        case np.int16:
            n_b = 2

    subseq = read_seqs[:]
    kmin, kmax = 11,11
    unitigs = []
    prev_unitigs = []

    res = []

    for k in range(kmin, kmax+1):
        t1 = time()
        sequences = subseq+unitigs
        kmers = get_kmer_count_from_sequences(sequences, k=k, n_b=n_b, cyclic=False)
        edges  = get_debruijn_edges_from_kmers(kmers, n_b=n_b)
        dbg = create_dbg_from_edges(edges, kmers)
        unitigs = get_unitigs_from_dbg(dbg, kmers, n_b=n_b)
        c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k, n_b=n_b)

        # unitigs_total = unitigs+prev_unitigs
        for i,u in enumerate(unitigs):
            if len(u[0])==(n_b*k) or k ==kmax:
                # print(u, k)
                foundInEdges = False
                if k!=kmax:
                    for i1, i2, _ in c_edges:
                        if i1==i or i2==i:
                            foundInEdges =True
                            # print(i,i1,i2)
                            break
                if not foundInEdges:
                    # print("not found",u)
                    prev_unitigs.append(u)
        
        
        unitigs_readable = [(num2seq(bytes2numseq(u[0],n_b),bi_alphabet), num2seq(bytes2numseq(u[1],n_b), bi_alphabet)) for u in unitigs]
        unitigs = [u[0] for u in unitigs]

        c_g = get_gt_graph(c_edges, unitigs_readable)
        g = get_gt_graph(edges, [(num2seq(bytes2numseq(kmer[0],n_b),bi_alphabet), num2seq(bytes2numseq(kmer[1],n_b),bi_alphabet)) for kmer in kmers])

                
        
        res.append((g,c_g))

        t2 = time()
        t=t2-t1
        h,m,s = t//3600, t%3600//60,t%60
        print("Step k = {:{}} processed in {:0=2.0f}:{:0=2.0f}:{:0=2.2f}, there are {} kmers and {} unitigs with {} unitigs selected" .format(k,len(str(kmax+1)),h,m,s,len(kmers),len(unitigs),len(prev_unitigs)))
    prev_unitigs = [bytes2numseq(u[0],n_b) for u in prev_unitigs]


    for k, (g,c_g) in zip(range(kmin, kmax+1),res):
        g.save(os.path.join(RES_OUTPUT_FOLDER,"graphtest_revcomp_k_{}_{}.graphml".format(kmin,k)))
        c_g.save(os.path.join(RES_OUTPUT_FOLDER,"c_graphtest_revcomp_k_{}_{}.graphml".format(kmin,k)))
        gfa_g = create_gfa(g,k)    
        gfa_g.to_file(os.path.join(RES_OUTPUT_FOLDER,"graphtest_k_{}_{}.gfa".format(kmin,k)))
        gfa_c_g = create_gfa(c_g,k)    
        gfa_c_g.to_file(os.path.join(RES_OUTPUT_FOLDER,"c_graphtest_k_{}_{}.gfa".format(kmin,k)))
