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

from rapidfuzz.distance import Levenshtein

from common.utils import bytes2numseq, numseq2bytes, seq2num, num2seq, compute_unitig_ref, print_ops
from data.data import Sequence
from graph.graph import graph_multi_k, dbg_tip_clipping, get_unitigs_bcalm, create_edges_from_dbg, add_to_dict
from graph.graph import get_kmer_count_from_sequences, get_debruijn_edges_from_kmers, get_unitigs_from_dbg, get_compacted_dbg_edges_from_unitigs, get_gt_graph, create_dbg_from_edges
from graph.graph import Graph
from data.visu import plot_debruijn_graph_gt, plot_compacted_debruijn_graph_gt
from data.inout import create_gfa_csv, save_sequences, load_sequences

def get_args():
    parser = argparse.ArgumentParser(
                    prog='DBG multi-k',
                    description='')
    parser.add_argument('--exp', default=None)
    parser.add_argument('--kmin',type=int,default=3)
    parser.add_argument('--kmax', type=int, default=None)
    parser.add_argument('--clipping', action="store_true")
    parser.add_argument('--no-multi-k', dest="multi_k", action="store_false")
    parser.add_argument('--kmer-abundance', choices = ["reads", "median_unitigs", "max_unitigs_reads"], default="max_unitigs_reads")
    args = parser.parse_args()
    if args.kmax == None:
        args.kmax = args.kmin
    return args

if __name__ == '__main__':
    
    args = get_args()

    # args.exp = "real"
    # args.kmin=23
    # args.kmax=23
    # args.clipping = True
    # args.kmer_abundance = "max_unitigs_reads"

    PROJECT_FOLDER = os.getcwd()
    RES_OUTPUT_FOLDER = os.path.join(PROJECT_FOLDER,"..","res","ecoli_bcalm","{}_data".format(args.exp))
    
    # filename = os.path.join(RES_OUTPUT_FOLDER,"res_unitigs_max_unitigs_reads_clippedTEST_k_2_35.pkl")
    # unitigs = load_sequences(filename)
    unitigs = Graph()

    
    
    
    if not os.path.isdir(RES_OUTPUT_FOLDER):
        os.mkdir(RES_OUTPUT_FOLDER)
    exp = args.kmer_abundance+"_"
    if args.clipping:
        exp+="clippedTEST_"
    if not args.multi_k:
        exp+="nmk_"
    sys. setrecursionlimit(10000)

    ref_file = "../input/truth_data/GCA_027944875.1_ASM2794487v1_genomic.truth_genes.json"
    # if args.exp == "real":
    #     ref_file = "../input/new_data/truth_data/GCA_027944875.1_ASM2794487v1_genomic.truth_genes.json"
    with open(ref_file, 'r') as f:
        ref_data = json.load(f)
    match args.exp:
        case "truth":
            read_file = "../input/truth_data/SRR23044204_1.subset.truth_genes.json"
        case "simulated":
            read_file = "../input/simulated_data/SRR23044204_1.subset.simulated_gene_dropout.json"
        case "new":
            read_file = "../input/new_data/SRR23044204_1.pandora_gene_calls.json"
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
    l_alphabet = len(alphabet)
    if l_alphabet < 2**7:
        Sequence.n_b = 1
    elif l_alphabet < 2**15:
        Sequence.n_b = 2
    elif l_alphabet < 2**31:
        Sequence.n_b = 4
    else:
        Sequence.n_b = 8
    Sequence.bi_alphabet = bi_alphabet

    ref_seqs = [Sequence(i,numseq2bytes(seq2num(seq,Sequence.bi_alphabet),Sequence.n_b),1) for i,seq in enumerate(ref_data.values())]
    # for k,g in ref_data.items():
    #     a = np.zeros(len(g),dtype = np.int16)
    #     for i,block in enumerate(g):
    #         a[i]=np.array(int(str(block[0])+str(blocks[block[1:]]))).astype(a.dtype)
    #     ref_seqs.append(a)
    
    if args.exp in ["truth","simulated"]:
        read_data_trimmed = [["_".join(block.split("_")[:-1]) for block in seq] for seq in read_data.values()]  
    else:
        read_data_trimmed = read_data.values()
    read_seqs = [Sequence(i, numseq2bytes(seq2num(seq,Sequence.bi_alphabet), Sequence.n_b), 1) for i,seq in enumerate(read_data_trimmed)]
    # seqs = []
    # for k,g in read_data.items():
    #     a = np.zeros(len(g),dtype = np.int16)
    #     for i,block in enumerate(g):
    #         a[i]=np.array(int(str(block[0])+str(blocks["_".join(block.split("_")[:-1])[1:]]))).astype(a.dtype)
    #     seqs.append(a)
    with open('reads.npy', 'wb') as f:
        np.save(f,np.array([seq.seq for seq in read_seqs],dtype=object))
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

    sims = []
    ss = []
    n=1
    for r in read_seqs[:]:
        sim=[]
        found = False
        for i in range(len(ref_seqs)):
            a = ref_seqs[i].num()
            scores = []
            for b in (r.num(),r.num(canonical=False)):
                ops = Levenshtein.opcodes(a,b)
                s = [op.dest_end-op.dest_start for op in ops if op.tag=="equal"]
                scores.append(sum(s)/len(b)* (len(s)==1))
                if sum(s)/len(b)==1 and len(s)!=1 and not found:
                    n+=1
                    # print(i,len(s))
                    # print(len(b))
                    found=True
                    ss.append(ops)
                    # print_ops(new_func(a),b,ops)
            #         break
            # if found:
            #     break
            sim.append(max(scores))
        sims.append(sim)
        # plt.hist(sim,bins=np.arange(0,1.1,0.1)-0.05)
    # for s in ss:
        # print([(o.src_end-o.src_start,o.src_start) for o in s if o.tag=="delete"])
    sim_a = np.array(sims)
    for i in range(len(ref_seqs)):
        plt.clf()
        plt.hist(sim_a[:,i])
        plt.savefig(os.path.join(RES_OUTPUT_FOLDER,"read_ref_{}.png".format(i+1)))
    sim_a_m = sim_a.max(axis=1)
    # print(sum(sim_a_m==1), sim_a_m.shape[0],sum(sim_a_m==1)/sim_a_m.shape[0], n)

    sim_t = [tuple(sim_a[k,:]) for k in range(sim_a.shape[0])]
    sim_t.sort()
    print(sim_t[-10:])
    print(sim_t[:10])
    plt.clf()
    fig, ax = plt.subplots()
    p = plt.imshow(sim_t, cmap='Greens', vmin=0, vmax=0.1, interpolation='none')
    ax.set_aspect(sim_a.shape[1]/sim_a.shape[0])
    fig.savefig(os.path.join(RES_OUTPUT_FOLDER,"read_ref_hm.png"))

    ref_read_sets = [set(),set()]
    seq_sets = []
    for rs in ref_seqs:
        s = set()
        for c in rs.num():
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
        for c in rs.num():
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

 
    subseq = read_seqs[:]
    kmin, kmax = args.kmin, args.kmax
    prev_unitigs = []

    res = []
    kmers_prev = None
    kmer_sets = []
    kmer_count_check =[]

    n_clip = 5

    for k in range(kmin, kmax+1):
        t1 = time()
        ### Count kmers
        
        if args.multi_k:
            klow = kmin
            kmers = get_kmer_count_from_sequences(unitigs, k=k, cyclic=False)
        else:
            klow = k
            kmers = Graph()

        kmers = get_kmer_count_from_sequences(subseq, k=k, cyclic=False, kmers=kmers)

        [kmer.compute_abundance(args.kmer_abundance) for kmer in kmers]
        
        
        if kmers_prev is not None:
            kmers_km1 = get_kmer_count_from_sequences(kmers, k=k-1, cyclic=False)
            kmer_sets.append((set(kmers_prev.keys()),set(kmers_km1.keys()),set(kmers_kp1.keys()),set(kmers.keys())))
            s_prev, s_next, s_prev_kp1, s_next_kp1 = kmer_sets[-1]
            kmer_count_check.append([len(s_prev&s_next), len(s_prev-s_next), len(s_next-s_prev), len(s_prev_kp1&s_next_kp1), len(s_prev_kp1-s_next_kp1), len(s_next_kp1-s_prev_kp1)])
            # print(kmer_count_check[-1])
            # if kmer_count_check[-1][4]>0:
            #     print([str(k1) for k1 in (kmer_sets[-1][2]-kmer_sets[-1][3])])
        kmers_prev = kmers.copy()

        if args.multi_k:
            # print(unitigs)
            kmers_kp1 = get_kmer_count_from_sequences(unitigs, k=k+1, cyclic=False)
        else:
            kmers_kp1 = Graph()
        kmers_kp1 = get_kmer_count_from_sequences(subseq, k=k+1, cyclic=False)
        [kmer.compute_abundance(args.kmer_abundance) for kmer in kmers_kp1]

        ### Create DBG on kmers
        kmers.compute_edges(k)
        # dbg = create_dbg_from_edges(edges, kmers)

        ### Save DBG on kmers
        g = get_gt_graph(kmers)
        g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}{{}}".format(exp,klow,k)),g,k, ["id","abundance"])

        ### Retrieve unitigs
        unitigs = get_unitigs_bcalm(kmers, k, on_unitig=False)
        unitigs.compute_edges(k)
        # l_u = [len(u) for u in unitigs]
        # l_u.sort()
        # print(l_u[:10])
        # print(l_u[-10:])
        ### Create compacted DBG on unitigs
        # c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k)
        
        
        
        u_ref_string, u_ref_id = compute_unitig_ref(unitigs, ref_seqs)
        # for kmer in kmers:
        #     print(kmer.__repr__(), kmer.num())
        # for unitig in unitigs:
        #     print(unitig.seq, [kmer.id for kmer in unitig.kmers], unitig.num())
        ### Save compacted DBG
        c_g = get_gt_graph(unitigs)
        c_g.vp["ref"] = c_g.new_vp("string", vals=u_ref_string)
        print(len(u_ref_id))
        for i, u_ref in enumerate(u_ref_id):
            c_g.vp["ref_ {}".format(i)] = c_g.new_vp("float", vals=u_ref)
        c_g.save(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}{{}}".format(exp,klow,k)),c_g,k, vp = ["id","ref","abundance"])    

        ### Graph cleaning
        if args.clipping:
            # dbg , kmers_bcalm = dbg_tip_clipping(dbg,k,1,3,kmers)
            ### Save cleaned graph
            # edges_cleaned = create_edges_from_dbg(dbg)
            unitigs.clip(n=k-1+n_clip)
            unitigs= get_unitigs_bcalm(unitigs, k, on_unitig=True)
            unitigs.compute_edges(k)
            unitigs.clip(a=10)
            unitigs= get_unitigs_bcalm(unitigs, k, on_unitig=True)
            unitigs.compute_edges(k)
            if k==kmax:
                unitigs.clip(a=10, n_neighbors=5)
                unitigs= get_unitigs_bcalm(unitigs, k, on_unitig=True)
                unitigs.compute_edges(k)
            
            g_cleaned = get_gt_graph(unitigs)
            g_cleaned.save(os.path.join(RES_OUTPUT_FOLDER,"c_graph_clean_{}k_{}_{}.graphml".format(exp,klow,k)))
            create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"c_graph_clean_{}k_{}_{}{{}}".format(exp,klow,k)),g_cleaned,k, ["id","abundance"])

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
        filename = os.path.join(RES_OUTPUT_FOLDER,"res_unitigs_{}k_{}_{}.pkl".format(exp,klow,k))
        save_sequences(unitigs, filename)

    prev_unitigs = [u for u in prev_unitigs]

    filename = os.path.join(RES_OUTPUT_FOLDER,"res_unitigs_{}k_{}_{}.pkl".format(exp,klow,k))
    # TEST LOAD SAV
    t1 = time()
    save_sequences(unitigs, filename)

    unitigs_neww = load_sequences(filename)
    t2 = time()
    print("Saving and loading final unitigs took {:0=2.2f} seconds".format(t2-t1))
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