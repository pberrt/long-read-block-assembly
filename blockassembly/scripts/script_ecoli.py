# script_graph_multi_k.py

import os
import sys
from time import time

import json
import argparse

import matplotlib.pyplot as plt
import numpy as np
import csv

import gfapy
import paired

sys.path.append('..')

from blockassembly.common.utils import bytes2numseq, numseq2bytes, seq2num, num2seq
from blockassembly.graph.graph import graph_multi_k, dbg_tip_clipping
from blockassembly.graph.graph import get_kmer_count_from_sequences, get_debruijn_edges_from_kmers, get_unitigs_from_dbg, get_compacted_dbg_edges_from_unitigs, get_gt_graph, create_dbg_from_edges
from blockassembly.data.visu import plot_debruijn_graph_gt, plot_compacted_debruijn_graph_gt
from blockassembly.data.inout import create_gfa_csv

def get_args():
    parser = argparse.ArgumentParser(
                    prog='DBG multi-k',
                    description='')
    parser.add_argument('--exp', default=None)
    parser.add_argument('--kmin',type=int,default=3)
    parser.add_argument('--kmax', type=int, default=None)
    parser.add_argument('--clipping', action="store_true")
    args = parser.parse_args()
    if args.kmax == None:
        args.kmax = args.kmin
    return args

if __name__ == '__main__':
    
    args = get_args()

    PROJECT_FOLDER = os.getcwd()
    RES_OUTPUT_FOLDER = os.path.join(PROJECT_FOLDER,"..","res","ecoli","{}_data".format(args.exp))
    if not os.path.isdir(RES_OUTPUT_FOLDER):
        os.mkdir(RES_OUTPUT_FOLDER)

    sys. setrecursionlimit(10000)

    ref_file = "../input/truth_data/GCA_027944875.1_ASM2794487v1_genomic.truth_genes.json"
    with open(ref_file, 'r') as f:
        ref_data = json.load(f)
    args.exp = "truth"
    match args.exp:
        case "truth":
            read_file = "../input/truth_data/SRR23044204_1.subset.truth_genes.json"
        case "simulated":
            read_file = "../input/simulated_data/SRR23044204_1.subset.simulated_gene_dropout.json"
        case "real":
            read_file = "../input/real_data/SRR23044204_1.subset.pandora_gene_calls.json"
    
    with open(read_file, 'r') as f:
        read_data = json.load(f)

    blocks = {}
    for k,g in read_data.items():
        for block in g:
            blocks[block[1:]]=None
    for k,g in ref_data.items():
        for block in g:
            blocks[block[1:]]=None
    blocks = list(blocks.keys())
    # blocks = dict([(p1,p2+1) for p1,p2 in zip(list(blocks), list(range(len(blocks))))])
    alphabet = [("+"+p1,"-"+p1) for p1 in list(blocks)]
    alphabet.sort()
    print(alphabet[:10])
    bi_alphabet = (alphabet,{k:((i+1)*((-1)**j)) for i, ks in enumerate(alphabet) for j, k in enumerate(ks)})

    ref_seqs = [seq2num(seq,bi_alphabet) for seq in ref_data.values()]
    # for k,g in ref_data.items():
    #     a = np.zeros(len(g),dtype = np.int16)
    #     for i,block in enumerate(g):
    #         a[i]=np.array(int(str(block[0])+str(blocks[block[1:]]))).astype(a.dtype)
    #     ref_seqs.append(a)
    
    if args.exp != "real":
        read_data_trimmed = [["_".join(block.split("_")[:-1]) for block in seq] for seq in read_data.values()]  
    else:
        read_data_trimmed = read_data.values()
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
    kmin, kmax = args.kmin, args.kmax
    unitigs = []
    prev_unitigs = []

    res = []


    for k in range(kmin, kmax+1):
        t1 = time()
        sequences = subseq+[u[0] for u in unitigs]
        kmers = get_kmer_count_from_sequences(sequences, k=k, n_b=n_b, cyclic=False)
        edges  = get_debruijn_edges_from_kmers(kmers, n_b=n_b)
        dbg = create_dbg_from_edges(edges, kmers)

        if args.clipping:
            dbg = dbg_tip_clipping(dbg,k,50,3)

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
        

        unitigs = [(bytes2numseq(u[0],n_b), bytes2numseq(u[1],n_b)) for u in unitigs]

        gap_score = -5
        match_score = 1
        mismatch_score = -1

        u_ref = []
        for both_u in unitigs:
            u_ref_scores = []
            for seq in ref_seqs:
                u_ref_scores.append([])
                for u in both_u:
                    # r = paired.align(seq,u, match_score=match_score, mismatch_score=mismatch_score, gap_score=gap_score)
                    # s = 0
                    # i_start = None
                    # i_end = None
                    # for i,(r1,r2) in enumerate(r):
                    #     if (r1!=None and r2!= None) or i_start is not None:
                    #         # print(i)
                    #         if i_start is None:
                    #             # print("step{}".format(1))
                    #             i_start=i
                    #         if r1 is None or r2 is None:
                    #             # print("step{}".format(2))
                    #             s+=gap_score
                    #         else:
                    #             i_end=i
                    #             if seq[r1]==u[r2]:
                    #                 # print("step{}".format(3))
                    #                 s+=match_score
                    #             else:
                    #                 # print("step{}".format(4))
                    #                 s+= mismatch_score
                    #         # print(i_end, matching_string)
                    # s -= (len(r)-i_end-1)*gap_score
                    s = len(set(seq).intersection(set(u)))
                    u_ref_scores[-1].append(s)
                
            u_ref_scores = [(max(urs),i) for i, urs in enumerate(u_ref_scores)]
            u_ref_scores.sort()
            # print(u_ref_scores)
            ref_string_tuples = []
            for s,i in u_ref_scores[::-1]:
                if s>0:
                    # and s>u_ref_scores[-1][0]/2:
                    ref_string_tuples.append((i,s))
            # ref_string_tuples.sort()
            if len(ref_string_tuples)==0:
                ref_string = "?"
            elif len(ref_string_tuples)==1:
                ref_string = str(ref_string_tuples[0][0])
            else:
                ref_string = ""
                for i,s in  ref_string_tuples:
                    ref_string+="{}({}) | ".format(i,s)
                ref_string=ref_string[0:-3]
            # print(ref_string)
            u_ref.append(ref_string)


        unitigs_readable = [(num2seq(u[0],bi_alphabet), num2seq(u[1], bi_alphabet)) for u in unitigs]
        c_g = get_gt_graph(c_edges, unitigs_readable)
        c_g.vp["ref"] = c_g.new_vp("string", vals=u_ref)

        g = get_gt_graph(edges, [(num2seq(bytes2numseq(kmer[0],n_b),bi_alphabet), num2seq(bytes2numseq(kmer[1],n_b),bi_alphabet)) for kmer in kmers])

                
        
        res.append((g,c_g))

        t2 = time()
        t=t2-t1
        h,m,s = t//3600, t%3600//60,t%60
        print("Step k = {:{}} processed in {:0=2.0f}:{:0=2.0f}:{:0=2.2f}, there are {} kmers and {} unitigs with {} unitigs selected" .format(k,len(str(kmax+1)),h,m,s,len(kmers),len(unitigs),len(prev_unitigs)))
    prev_unitigs = [bytes2numseq(u[0],n_b) for u in prev_unitigs]

    if args.clipping:
        exp = "clipped_"
    else:
        exp = ""
    for k, (g,c_g) in zip(range(kmin, kmax+1),res):
        g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}.graphml".format(exp,kmin,k)))
        c_g.save(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}.graphml".format(exp,kmin,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}{{}}".format(exp,kmin,k)),g,k)
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}{{}}".format(exp,kmin,k)),c_g,k, vp = ["id","ref"])    
