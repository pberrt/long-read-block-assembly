# script_ecoli.py

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

from common.utils import bytes2numseq, numseq2bytes, seq2num, num2seq
from graph.graph import graph_multi_k, dbg_tip_clipping, get_unitigs_bcalm, Bcalm_kmer, create_edges_from_dbg, add_to_dict
from graph.graph import get_kmer_count_from_sequences, get_debruijn_edges_from_kmers, get_unitigs_from_dbg, get_compacted_dbg_edges_from_unitigs, get_gt_graph, create_dbg_from_edges
from data.visu import plot_debruijn_graph_gt, plot_compacted_debruijn_graph_gt
from data.inout import create_gfa_csv

def get_args():
    parser = argparse.ArgumentParser(
                    prog='DBG multi-k',
                    description='')
    parser.add_argument('--exp', default=None)
    parser.add_argument('--kmin',type=int,default=3)
    parser.add_argument('--kmax', type=int, default=None)
    parser.add_argument('--clipping', action="store_true")
    parser.add_argument('--no-multi-k', dest="multi_k", action="store_false")
    parser.add_argument('--mode', choices = ["normal","bcalm"], default="bcalm")
    parser.add_argument('--kmer-abundance', choices = ["reads", "median_unitigs", "max_unitigs_reads"], default="max_unitigs_reads")
    args = parser.parse_args()
    if args.kmax == None:
        args.kmax = args.kmin
    return args

if __name__ == '__main__':
    
    args = get_args()

    # args.exp = "test"
    # args.kmin=3
    # args.kmax=4
    # args.clipping = True
    # args.mode = "normal"
    PROJECT_FOLDER = os.getcwd()
    RES_OUTPUT_FOLDER = os.path.join(PROJECT_FOLDER,"..","res","ecoli{}".format(("_"+args.mode) if args.mode!="normal" else ""),"{}_data".format(args.exp))
    if not os.path.isdir(RES_OUTPUT_FOLDER):
        os.mkdir(RES_OUTPUT_FOLDER)

    exp = args.kmer_abundance+"_"
    if args.clipping:
        exp+="clippedTEST_"
    if not args.multi_k:
        exp+="nmk_"
    sys. setrecursionlimit(10000)

    ref_file = "../input/truth_data/GCA_027944875.1_ASM2794487v1_genomic.truth_genes.json"
    with open(ref_file, 'r') as f:
        ref_data = json.load(f)
    match args.exp:
        case "truth":
            read_file = "../input/truth_data/SRR23044204_1.subset.truth_genes.json"
        case "simulated":
            read_file = "../input/simulated_data/SRR23044204_1.subset.simulated_gene_dropout.json"
        case "real":
            read_file = "../input/real_data/SRR23044204_1.subset.pandora_gene_calls.json"
        case "test":
            read_file = "../input/test.json"
    
    with open(read_file, 'r') as f:
        read_data = json.load(f)

    blocks2reads = {}
    for k,g in read_data.items():
        for block in g:
            if args.exp in ["truth","simulated"]:
                add_to_dict(blocks2reads,"_".join(block[1:].split("_")[:-1]),k)
            else:
                add_to_dict(blocks2reads,block[1:],k)
    for k,g in ref_data.items():
        for block in g:
            add_to_dict(blocks2reads,block[1:],k)
 
    # len_blocks = [(len(l),block) for block, l in blocks2reads.items()]
    # len_blocks.sort(key = lambda x: (x[0], x[1]))
    # print(len(len_blocks),len_blocks[:5],len_blocks[-5:])
    # c,b,_ = plt.hist([l[0] for l in len_blocks], bins = np.arange(len_blocks[0][0],len_blocks[-1][0]+1)-0.5)
    # plt.show()
    # plt.plot(np.cumsum(c))
    # plt.show()
    
    blocks = list(blocks2reads.keys())
    blocks.sort()

    # blocks = dict([(p1,p2+1) for p1,p2 in zip(list(blocks), list(range(len(blocks))))])
    alphabet = [("+"+p1,"-"+p1) for p1 in blocks]
    print(alphabet[:10])
    bi_alphabet = (alphabet,{k:((i+1)*((-1)**j)) for i, ks in enumerate(alphabet) for j, k in enumerate(ks)})

    ref_seqs = [seq2num(seq,bi_alphabet) for seq in ref_data.values()]
    # for k,g in ref_data.items():
    #     a = np.zeros(len(g),dtype = np.int16)
    #     for i,block in enumerate(g):
    #         a[i]=np.array(int(str(block[0])+str(blocks[block[1:]]))).astype(a.dtype)
    #     ref_seqs.append(a)
    
    if args.exp in ["truth","simulated"]:
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
    with open('reads.npy', 'wb') as f:
        np.save(f,np.array(read_seqs,dtype=object))
    b= {}
    for i, (p1,p2) in enumerate(alphabet):
        b[i+1]= p1
        b[-(i+1)]=p2
    j = json.dumps(b)
    
    with open('blocks.json', 'w') as f:
        f.write(j)
    
    read_length=[len(a) for a in read_seqs]
    print(sum(read_length),len(read_length))
    plt.hist(read_length,bins=np.arange(0.5,max(read_length)+0.5))
    plt.savefig(os.path.join(RES_OUTPUT_FOLDER,"read_hist.png"))


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
    unitigs = {}
    prev_unitigs = []

    res = []
    kmers_prev = None
    kmer_sets = []
    kmer_count_check =[ ]

    for k in range(kmin, kmax+1):
        t1 = time()
        ### Count kmers
        
        if args.multi_k:
            klow = kmin
            kmers = get_kmer_count_from_sequences([(u[0],u[2]) for u in unitigs.values()], k=k, n_b=n_b, cyclic=False,mode=args.kmer_abundance,count_key=0)
        else:
            klow = k
            kmers = {}

        kmers = get_kmer_count_from_sequences(subseq, k=k, n_b=n_b, cyclic=False, kmers=kmers,mode=args.kmer_abundance,count_key=1)

        for kmer in kmers:
            match args.kmer_abundance:
                case "reads":
                    kmers[kmer]=sum(kmers[kmer][:])
                case "median_unitigs":
                    kmers[kmer]=kmers[kmer][0] if kmers[kmer][0]!=0 else kmers[kmer][1]
                case  "max_unitigs_reads":
                    kmers[kmer]=max(kmers[kmer])
        
        kmers_num = {kmer[0]:(bytes2numseq(kmer[0],n_b), bytes2numseq(kmer[1],n_b), a) for kmer, a  in kmers.items()}
        
        if kmers_prev is not None:
            kmers_km1 = get_kmer_count_from_sequences([(u[0],u[2]) for u in kmers_num.values()], k=k-1, n_b=n_b, cyclic=False)
            kmer_sets.append((set(kmers_prev.keys()),set(kmers_km1.keys()),set(kmers_kp1.keys()),set(kmers.keys())))
            s_prev, s_next, s_prev_kp1, s_next_kp1 = kmer_sets[-1]
            kmer_count_check.append([len(s_prev&s_next), len(s_prev-s_next), len(s_next-s_prev), len(s_prev_kp1&s_next_kp1), len(s_prev_kp1-s_next_kp1), len(s_next_kp1-s_prev_kp1)])
            print(kmer_count_check[-1])
            if kmer_count_check[-1][4]>0:
                print([("~~~".join(num2seq(bytes2numseq(k1,n_b),bi_alphabet)),"~~~".join(num2seq(bytes2numseq(k2,n_b),bi_alphabet))) for k1,k2 in (kmer_sets[-1][2]-kmer_sets[-1][3])])
        kmers_prev = kmers.copy()

        if args.multi_k:
            kmers_kp1 = get_kmer_count_from_sequences([(u[0],u[2]) for u in unitigs.values()], k=k+1, n_b=n_b, cyclic=False,mode=args.kmer_abundance,count_key=0)
        else:
            kmers_kp1 = {}
        kmers_kp1 = get_kmer_count_from_sequences(subseq, k=k+1, n_b=n_b, cyclic=False, kmers=kmers_kp1,mode=args.kmer_abundance,count_key=1)
        for kmer in kmers_kp1:
            kmers_kp1[kmer]=kmers_kp1[kmer][0]

        ### Create DBG on kmers
        edges  = get_debruijn_edges_from_kmers(kmers, n_b=n_b)
        dbg = create_dbg_from_edges(edges, kmers)
        
        ### Save DBG on kmers
        kmers_readables = {kmer: (num2seq(k1,bi_alphabet), num2seq(k2,bi_alphabet),a) for kmer,(k1,k2,a) in kmers_num.items()}
        g = get_gt_graph(edges, kmers_readables)
        g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}{{}}".format(exp,klow,k)),g,k, ["id","abundance"])
        
        ### Graph cleaning
        kmers_bcalm = [Bcalm_kmer(i,kmer_b[0], a, n_b, bi_alphabet) for i,(kmer_b, a) in enumerate(kmers.items())] if args.mode=="bcalm" else None
        if args.clipping:
            dbg , kmers_bcalm = dbg_tip_clipping(dbg,k,1,3,kmers_bcalm)

        ### Save cleaned graph
        edges_cleaned = create_edges_from_dbg(dbg)
        g_cleaned = get_gt_graph(edges_cleaned, kmers_readables)
        g_cleaned.save(os.path.join(RES_OUTPUT_FOLDER,"graph_clean_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"graph_clean_{}k_{}_{}{{}}".format(exp,klow,k)),g_cleaned,k, ["id","abundance"])

        ### Retrieve unitigs
        match args.mode:
            case "normal":
                unitigs = get_unitigs_from_dbg(dbg, kmers, n_b=n_b)
            case "bcalm":
                unitigs = get_unitigs_bcalm(kmers_bcalm, n_b=n_b, on_kmer=True)

        ### Create compacted DBG on unitigs
        c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k, n_b=n_b)
        
        ### Compute unitig ref

        unitigs = {u[0]:(bytes2numseq(u[0],n_b), bytes2numseq(u[1],n_b),unit.abundance) for u,unit  in unitigs.items()}
        gap_score = -5
        match_score = 1
        mismatch_score = -1

        u_ref = []
        for u_bytes in unitigs:
            both_u = unitigs[u_bytes][:2]
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

        ### Save compacted DBG
        unitigs_readable = {u: (num2seq(u0,bi_alphabet), num2seq(u1, bi_alphabet),a) for u, (u0,u1,a) in unitigs.items()}
        c_g = get_gt_graph(c_edges, unitigs_readable)
        c_g.vp["ref"] = c_g.new_vp("string", vals=u_ref)
        c_g.save(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}{{}}".format(exp,klow,k)),c_g,k, vp = ["id","ref","abundance"])    

        # # unitigs_total = unitigs+prev_unitigs
        # for i,u in enumerate(unitigs):
        #     if len(u[0])==(n_b*k) or k ==kmax:
        #         # print(u, k)
        #         foundInEdges = False
        #         if k!=kmax:
        #             for i1, i2, _ in c_edges:
        #                 if i1==i or i2==i:
        #                     foundInEdges =True
        #                     # print(i,i1,i2)
        #                     break
        #         if not foundInEdges:
        #             # print("not found",u)
        #             prev_unitigs.append(u)
        t2 = time()
        t=t2-t1
        h,m,s = t//3600, t%3600//60,t%60
        print("Step k = {:{}} processed in {:0=2.0f}:{:0=2.0f}:{:0=2.2f}, there are {} kmers and {} unitigs with {} unitigs selected" .format(k,len(str(kmax+1)),h,m,s,len(kmers),len(unitigs),len(prev_unitigs)))
    prev_unitigs = [bytes2numseq(u[0],n_b) for u in prev_unitigs]


if len(kmer_count_check)>0:
    width = 0.5 
    kmer_count_check = np.array(kmer_count_check).T
    types = ["common k-1","lost k-1","added k-1","common k+1","lost k+1","added k+1"]
    weight_counts = { t:k for t,k in zip(types,kmer_count_check)}
    print(weight_counts)
    weight_counts = {key:weight_counts[key] for key in weight_counts if key in ["lost k-1","added k-1","lost k+1","added k+1"]}
    fig, ax = plt.subplots()
    bottom = np.zeros(kmax-kmin)
    print(weight_counts)
    for boolean, weight_count in weight_counts.items():
        p = ax.bar(np.arange(kmin,kmax)+0.5, weight_count, width, label=boolean, bottom=bottom)

        bottom += weight_count

    ax.set_title("")
    ax.legend(loc="upper right")

    plt.savefig(os.path.join(RES_OUTPUT_FOLDER,"kmer_count_{}{}_{}.png".format(exp,kmin,kmax)))