# script_ecoli.py

import os
import sys
from time import time
import glob

import json
import argparse

import numpy as np

from common.utils import numseq2bytes, seq2num, hms
from data.data import Sequence, get_num
from graph.graph import get_unitigs_bcalm, add_to_dict, switch_index
from graph.graph import get_kmer_count_from_sequences, get_gt_graph
from graph.graph import Graph
from graph.poa import POAGraph, generate_poa_graph
from graph.components import add_unitigs_sets, check_cycles
from graph.dag import get_gt_graph_poa, get_max_block, get_all_path, add_coordinate, needleman_wunsch_inclusion
from data.inout import create_gfa_csv, save_sequences

import numpy as np

from Bio import SeqIO, SeqRecord
import pyfastx


def get_args():
    parser = argparse.ArgumentParser(
                    prog='DBG multi-k',
                    description='')
    parser.add_argument('--sample', default="GCA_027944575.1_ASM2794457v1_genomic")
    parser.add_argument('--exp', default=None)
    parser.add_argument('--kmin',type=int,default=2)
    parser.add_argument('--kmax', type=int, default=35)
    parser.add_argument('--clipping', action="store_true")
    parser.add_argument('--kmer-abundance', choices = ["reads", "median_unitigs", "max_unitigs_reads"], default="reads")
    args = parser.parse_args()
    return args
       
if __name__ == '__main__':
    t_start = time()
    args = get_args()
    
    exp=args.exp
    exp = "TEST"
    sample_name = args.sample
    # args.exp = "real"
    # args.kmin=23
    # args.kmax=23
    # args.clipping = True
    # args.kmer_abundance = "max_unitigs_reads"

    PROJECT_FOLDER = os.getcwd()

    # TODO Change result folder according to exp
    RES_OUTPUT_FOLDER = os.path.join(PROJECT_FOLDER,"..","res","ecoli_bcalm","{}_data".format(exp),sample_name)
    if not os.path.isdir(RES_OUTPUT_FOLDER):
        os.mkdir(RES_OUTPUT_FOLDER)
    RES_FASTA_FOLDER = os.path.join(RES_OUTPUT_FOLDER,"coordinates_fasta")
    if not os.path.exists(RES_FASTA_FOLDER):
        os.mkdir(RES_FASTA_FOLDER)

    sys. setrecursionlimit(10000)

    t= time()
    print("\n"+"-"*50+"\n"+"\t\tLOAD INPUT\n"+"-"*50)
    INPUT_FOLDER = os.path.join("..","input","all_data",sample_name)
    fastq_file = os.path.join(INPUT_FOLDER,(glob.glob(INPUT_FOLDER+"/*.fastq")[0]).split("/")[-1])
    read_file = os.path.join(INPUT_FOLDER,"gene_calls_with_gene_filtering.json")
    read_pos_file = os.path.join(INPUT_FOLDER,"gene_positions_with_gene_filtering.json")
    # read_file = os.path.join(INPUT_FOLDER,"gene_calls_without_gene_filtering.json")
    # read_pos_file = os.path.join(INPUT_FOLDER,"gene_positions_without_gene_filtering.json")


    # TODO add ref parsing and ref data
    # ref_file = "../input/truth_data/GCA_027944875.1_ASM2794487v1_genomic.truth_genes.json"
    # read_file = "../input/new_data/SRR23044204_1.pandora_gene_calls.json"
    # read_pos_file = "../input/new_data/SRR23044204_1.pandora_gene_positions.json"
    # fastq_file ="../input/new_data/SRR23044204_1.fastq"
    print("Loading FASTQ file...: ", fastq_file)
    # fastq = SeqIO.to_dict(SeqIO.parse(fastq_file,"fastq"))
    fastq = pyfastx.Fastq(fastq_file)
    print("Loading gene calls file...:... ", read_file)
    with open(read_file, 'r') as f:
        read_data = json.load(f)
    print("Loading gene positions file...: ", read_pos_file)
    with open(read_pos_file, 'r') as f:
        gene_positions_reads = json.load(f)
    # print("\tloading... ", ref_file)
    # with open(ref_file, 'r') as f:
    #     ref_data = json.load(f)
    print("Input loading time: {:0=2.0f}:{:0=2.0f}:{:0=2.2f}".format(*hms(time()-t)))
    
    t= time()
    print("\n"+"-"*50+"\n"+"\t\tINPUT PROCESSING\n"+"-"*50)
    # filename = os.path.join(RES_OUTPUT_FOLDER,"res_unitigs_max_unitigs_reads_clippedTEST_k_2_35.pkl")
    # unitigs = load_sequences(filename)
    # unitigs = Graph()

    blocks2reads = {}
    for k,g in read_data.items():
        for block in g:
            add_to_dict(blocks2reads,block[1:],k)
    # for k,g in ref_data.items():
    #     for block in g:
    #         add_to_dict(blocks2reads,block[1:],k)
     
    blocks = list(blocks2reads.keys())
    blocks.sort()
    
    alphabet = [("+"+p1,"-"+p1) for p1 in blocks]
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

    # ref_seqs = [Sequence(i,numseq2bytes(seq2num(seq,Sequence.bi_alphabet),Sequence.n_b),1) for i,seq in enumerate(ref_data.values())]
    
    read_data_trimmed = read_data.values()
    read_seqs = [Sequence(i, numseq2bytes(seq2num(seq,Sequence.bi_alphabet), Sequence.n_b), 1) for i,seq in enumerate(read_data_trimmed)]
    print("The sample contains {} reads and {} different blocks.".format(len(read_seqs), len(blocks)))
    # with open('reads.npy', 'wb') as f:
    #     np.save(f,np.array([seq.seq for seq in read_seqs],dtype=object))
    # b= {}
    # for i, (p1,p2) in enumerate(alphabet):
    #     b[i+1]= p1
    #     b[-(i+1)]=p2
    # j = json.dumps(b)
    
    # with open('blocks.json', 'w') as f:
    #     f.write(j)
    
    # read_length=[len(a) for a in read_seqs]
    # print(sum(read_length),len(read_length))
    # plt.hist(read_length,bins=np.arange(0.5,max(read_length)+0.5))
    # plt.savefig(os.path.join(RES_OUTPUT_FOLDER,"read_hist.png"))

    # sims = []
    # ss = []
    # n=1
    # for r in read_seqs[:]:
    #     sim=[]
    #     found = False
    #     for i in range(len(ref_seqs)):
    #         a = ref_seqs[i].num()
    #         scores = []
    #         for b in (r.num(),r.num(canonical=False)):
    #             ops = Levenshtein.opcodes(a,b)
    #             s = [op.dest_end-op.dest_start for op in ops if op.tag=="equal"]
    #             scores.append(sum(s)/len(b)* (len(s)==1))
    #             if sum(s)/len(b)==1 and len(s)!=1 and not found:
    #                 n+=1
    #                 # print(i,len(s))
    #                 # print(len(b))
    #                 found=True
    #                 ss.append(ops)
    #                 # print_ops(new_func(a),b,ops)
    #         #         break
    #         # if found:
    #         #     break
    #         sim.append(max(scores))
    #     sims.append(sim)
    #     # plt.hist(sim,bins=np.arange(0,1.1,0.1)-0.05)
    # # for s in ss:
    #     # print([(o.src_end-o.src_start,o.src_start) for o in s if o.tag=="delete"])
    # sim_a = np.array(sims)
    # for i in range(len(ref_seqs)):
    #     plt.clf()
    #     plt.hist(sim_a[:,i])
    #     plt.savefig(os.path.join(RES_OUTPUT_FOLDER,"read_ref_{}.png".format(i+1)))
    # sim_a_m = sim_a.max(axis=1)
    # # print(sum(sim_a_m==1), sim_a_m.shape[0],sum(sim_a_m==1)/sim_a_m.shape[0], n)

    # sim_t = [tuple(sim_a[k,:]) for k in range(sim_a.shape[0])]
    # sim_t.sort()
    # print(sim_t[-10:])
    # print(sim_t[:10])
    # plt.clf()
    # fig, ax = plt.subplots()
    # p = plt.imshow(sim_t, cmap='Greens', vmin=0, vmax=0.1, interpolation='none')
    # ax.set_aspect(sim_a.shape[1]/sim_a.shape[0])
    # fig.savefig(os.path.join(RES_OUTPUT_FOLDER,"read_ref_hm.png"))

    # ref_read_sets = [set(),set()]
    # seq_sets = []
    # for rs in ref_seqs:
    #     s = set()
    #     for c in rs.num():
    #         s.add(abs(c))
    #         ref_read_sets[0].add(abs(c))
    #     seq_sets.append(s)
    # for i, s1 in enumerate(seq_sets):
    #     print("Sequence {} contains {} different element(s)".format(i,len(s1)))
    #     for j, s2 in enumerate(seq_sets):
    #         if i!=j:
    #             inter, diff1, diff2 = s1.intersection(s2), s1.difference(s2), s2.difference(s1)
    #             linter = len(inter)
    #             if linter > 0:
    #                 print("\tShares {} elements with sequence {}: {}".format(linter,j,inter))
    #             if len(diff1) == 0:
    #                 print("{} is included in {}".format(i,j)) 
    #             if len(diff2) == 0:
    #                 print("{} is included in {}".format(j,i))
    #     if len(s1)<10:
    #         print("\t",s1)
    # for rs in read_seqs:
    #     for c in rs.num():
    #         ref_read_sets[1].add(abs(c))

    # print("The {} reference sequences contains {} unique elements.".format(len(ref_seqs),len(ref_read_sets[0])))
    # print("The {} read sequences contains {} unique elements.".format(len(read_seqs),len(ref_read_sets[1])))
    # inter = ref_read_sets[0].intersection(ref_read_sets[1])
    # diff1 = ref_read_sets[0].difference(ref_read_sets[1])
    # diff2 = ref_read_sets[1].difference(ref_read_sets[0])
    # print("{} elements are common to both reference sequences and reads.".format(len(inter)))
    # print("{} elements are only present in reference sequences:".format(len(diff1)))
    # print("\t",*diff1)
    # print("{} elements are only present in reads:".format(len(diff2)))
    # print("\t",*diff2)

    print("Input processing time: {:0=2.0f}:{:0=2.0f}:{:0=2.2f}".format(*hms(time()-t)))
 
    t= time()
    print("\n"+"-"*50+"\n"+"\tDe Bruijn graph backbone construction\n"+"-"*50)
    subseq = read_seqs[:]
    kmin, kmax = args.kmin, args.kmax
    # prev_unitigs = []
    # unitigs = []

    res = []
    kmers_prev = None
    kmer_sets = []
    kmer_count_check =[]
    need_break=False
    n_clip = 5

    for k in range(kmin, kmax+1):
        t1 = time()
        print("Start k={}".format(k))
        ### Count kmers
        
        klow = kmin
        # kmers = get_kmer_count_from_sequences(unitigs, k=k, cyclic=False)
        kmers = Graph()
        kmers = get_kmer_count_from_sequences(subseq, k=k, cyclic=False, kmers=kmers)

        [kmer.compute_abundance(args.kmer_abundance) for kmer in kmers]
        print("\tThere are {} kmers".format(len(kmers)))
        
        # if kmers_prev is not None:
        #     kmers_km1 = get_kmer_count_from_sequences(kmers, k=k-1, cyclic=False)
        #     kmer_sets.append((set(kmers_prev.keys()),set(kmers_km1.keys()),set(kmers_kp1.keys()),set(kmers.keys())))
        #     s_prev, s_next, s_prev_kp1, s_next_kp1 = kmer_sets[-1]
        #     kmer_count_check.append([len(s_prev&s_next), len(s_prev-s_next), len(s_next-s_prev), len(s_prev_kp1&s_next_kp1), len(s_prev_kp1-s_next_kp1), len(s_next_kp1-s_prev_kp1)])
        #     # print(kmer_count_check[-1])
        #     # if kmer_count_check[-1][4]>0:
        #     #     print([str(k1) for k1 in (kmer_sets[-1][2]-kmer_sets[-1][3])])
        # kmers_prev = kmers.copy()

        # kmers_kp1 = get_kmer_count_from_sequences(unitigs, k=k+1, cyclic=False)
        # kmers_kp1 = get_kmer_count_from_sequences(subseq, k=k+1, cyclic=False)
        # [kmer.compute_abundance(args.kmer_abundance) for kmer in kmers_kp1]

        ### Create DBG on kmers
        kmers.compute_edges(k)
        # dbg = create_dbg_from_edges(edges, kmers)

        ### Save DBG on kmers
        g = get_gt_graph(kmers)
        g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"graph_{}k_{}_{}{{}}".format(exp,klow,k)),g,k, ["id","abundance"])
        
        # kmers.detach_abundance(k,1)
        
        # g = get_gt_graph(kmers)
        # g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_clean_{}k_{}_{}.graphml".format(exp,klow,k)))
        # create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"graph_clean_{}k_{}_{}{{}}".format(exp,klow,k)),g,k, ["id","abundance"])
        
        ### Retrieve unitigs
        unitigs = get_unitigs_bcalm(kmers, k, on_unitig=False)
        unitigs.compute_edges(k)
        # l_u = [len(u) for u in unitigs]
        # l_u.sort()
        # print(l_u[:10])
        # print(l_u[-10:])
        ### Create compacted DBG on unitigs
        # c_edges = get_compacted_dbg_edges_from_unitigs(unitigs,k)
        
        
        
        # u_ref_string, u_ref_id = compute_unitig_ref(unitigs, ref_seqs)
        # for kmer in kmers:
        #     print(kmer.__repr__(), kmer.num())
        # for unitig in unitigs:
        #     print(unitig.seq, [kmer.id for kmer in unitig.kmers], unitig.num())
        ### Save compacted DBG
        c_g = get_gt_graph(unitigs)
        # c_g.vp["ref"] = c_g.new_vp("string", vals=u_ref_string)
        # for i, u_ref in enumerate(u_ref_id):
        #     c_g.vp["ref_ {}".format(i)] = c_g.new_vp("float", vals=u_ref)
        c_g.save(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}k_{}_{}{{}}".format(exp,klow,k)),c_g,k, vp = ["id","abundance"])    

        ### Graph cleaning
        # if args.clipping:
            # dbg , kmers_bcalm = dbg_tip_clipping(dbg,k,1,3,kmers)
            ### Save cleaned graph
            # edges_cleaned = create_edges_from_dbg(dbg)
        unitigs.clip(k,n=k-1+n_clip)
        unitigs.clip(k,a=5, n_neighbors=5)
        # print('huge clipping')
        # print(len(unitigs))
        # if k==kmax:
        #     unitigs.clip(k, a=10, n_neighbors=5)
        # print(len(unitigs))
        
        g_cleaned = get_gt_graph(unitigs)
        g_cleaned.save(os.path.join(RES_OUTPUT_FOLDER,"c_graph_clean_{}k_{}_{}.graphml".format(exp,klow,k)))
        create_gfa_csv(os.path.join(RES_OUTPUT_FOLDER,"c_graph_clean_{}k_{}_{}{{}}".format(exp,klow,k)),g_cleaned,k, ["id","abundance"])
        components = []
        for u in unitigs:
            found = False
            for key in range(len(components)):
                if u in components[key]:
                    found=True
                    break
            if not found:
                components.append(Graph())
                key = len(components)-1
                components[key][u]=u
                components = add_unitigs_sets({uu:uu for uu in u.link[0]}, components, key=key)
                components = add_unitigs_sets({uu:uu for uu in u.link[1]}, components, key=key)
        print("\tThere are {} connected components".format(len(components)))
        # components = [s for s in components if len(s)>1]
        component_length = [sum([len(u) for u in component])-(len(component)-1)*(k-1) for component in components]
        components = [(component, cl) for (cl,component) in sorted(zip(component_length,components), key=lambda pair: pair[0])][::-1]
        for i,(component,cl) in enumerate(components):
            sep = "\t  -->\t" if i==0 else "\t\t"
            print("{}#{}: {} unitigs (~{} blocks)".format(sep,i+1,len(component),cl))
        chromosome = components[0][0]
        # g = get_gt_graph(unitigs)
        # g.save(os.path.join("graph_on_unitigs.graphml"))
        # create_gfa_csv("graph_on_unitigs{}",g,k)

        cyclic, seen = check_cycles(chromosome, verbose=False)
        if cyclic:
            print("\tAt least one cycle found, trying to decycle...")
            best = max(chromosome, key = len)
            best = chromosome.pop(best)
            n = len(best)
            
            # [0  (split_s2  [split_s1    (n
            # split_s1-split_s2=k-2
            split_s2 = (n-k+2)//2
            split_s1 = split_s2+k-2
            
            # s1 = s[0:split_s1(exclu)]
            # s2 = s[split_s2:n(exclu)]
            s1 = Sequence(best.id,best.seq[:Sequence.n_b*split_s1], best.abundance)
            s2 = Sequence(len(unitigs),best.seq[Sequence.n_b*split_s2:], best.abundance)
            s1.stability, s2.stability = -2,-2
            chromosome.data[s1]=s1
            chromosome.data[s2]=s2
            _ = chromosome.compute_edges(k)
            # print([[u.id for u in u_l] for u_l in s1.link])
            # print([[u.id for u in u_l] for u_l in s2.link])
            # chromosome.compute_order()
            # for s in chromosome:
            #     if s.id in [682, 2435, 2310,1802,316,485,]:
            #         print(s.order,[[u.order for u in u_l] for u_l in s.link], s.num(),s.num(canonical=False))
            starts= [s1,s2]
        else:
            starts=[]
            for u in chromosome:
                l1,l2 = len(u.link[0]) , len(u.link[1])
                if l1*l2==0 and l1+l2!=0:
                    starts.append(u)
        if not check_cycles(chromosome)[0]:
            start, end = starts[1],starts[0]
            print("\tAcyclic graph found")
            need_break=True
        else:
            print("\tStill cycles, larger k.")
        t2 = time()
        t_tot=t2-t1
        print("\tStep k = {:{}} processed in {:0=2.0f}:{:0=2.0f}:{:0=2.2f}, there are {} kmers and {} unitigs" .format(k,len(str(kmax+1)),*hms(t_tot),len(kmers),len(unitigs)))
        filename = os.path.join(RES_OUTPUT_FOLDER,"res_unitigs_{}k_{}_{}.pkl".format(exp,klow,k))
        save_sequences(unitigs, filename)
            
        if need_break:
            break

    print("De Bruijn graph backbone construction time: {:0=2.0f}:{:0=2.0f}:{:0=2.2f}".format(*hms(time()-t)))
 
    t= time()
    print("\n"+"-"*50+"\n"+"\tDAG processing\n"+"-"*50)
    
    coordinates_1 = dict()
    coordinates_2 = dict()
    s = start
    c = 0

    coordinates_1 = add_coordinate(start, 1, coordinates_1,0,k)
    coordinates_2 = add_coordinate(end, -1, coordinates_2,0,k)

    max_1 = max([x[0][1] for x in coordinates_1.values()])
    max_2 = max([x[0][1] for x in coordinates_2.values()])

    coordinates_2 = {s:([max_2-x[0][1], max_2-x[0][0]],-x[1]) for s,x in coordinates_2.items()}

    coordinates = dict()
    
    g = get_gt_graph(chromosome)
    g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_stability.graphml"))
    
    for s in chromosome:
        if s in coordinates_1 and s in coordinates_2:
            coordinates[s] = (np.array(coordinates_1[s][0])+np.array(coordinates_2[s][0]),coordinates_1[s][1])
            s.debug=str(coordinates_1[s]),str(coordinates_2[s]),str(coordinates[s])
    max_3 = max([coordinates[s][0][1] for s in coordinates])
    
    coordinates_diff = {s:coordinates_2[s][0][0]-coordinates_1[s][0][0] for s in coordinates}
    for u in coordinates:
        u.stability = coordinates_diff[u]
    
    g = get_gt_graph(chromosome)
    g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_stability_1.graphml"))
    for s in coordinates:
        if (coordinates_2[s][0][0]-coordinates_1[s][0][0])>6:
            print(str(s))
    
    consensus_set_coordinates = [set() for _ in range(max_3)]
    for s in coordinates:
        c, mode = coordinates[s]
        for i,b in enumerate(s.num(canonical=(mode==1))):
            consensus_set_coordinates[c[0]+2*i].add(b)
            consensus_set_coordinates[c[0]+2*i+1].add(b)
    ## Compute consensus from local POA
    # sequences = ["abcdef","abcddef"]
    # sequences = [[12301,12302,12303,12304,12305,12306,12307,81230],[12301,12302,12302,12303,12304,12305,12306,12307,81230]]
    # sequences = [read.num() for read in read_seqs[:100]]
    # # test_order_of_alignment_case1()
    # g = generate_poa_graph(sequences)
    # # g.generateAlignmentStrings()
    # Create topological order from coordinates
    topological_order = [d+ (s,) for s,d  in coordinates.items()]
    topological_order.sort(key = lambda x: x[0][0])
    for i,s in enumerate(topological_order):
        s[2].topo_order = i

    # Compute topological order edges
    topo_dict = {s[2]:i for i,s in enumerate(topological_order)}
    e_topo = []
    for i, (_, mode , s) in enumerate(topological_order):
        for n in s.link[switch_index(1,mode)]:
            e_topo.append([topo_dict[ss] for ss in [s,n] if n in topo_dict])

    # Get stables node in the DAG
    stables = []
    n_nodes = set()
    for i, (_, mode , s) in enumerate(topological_order):
        if i in n_nodes:
            n_nodes.remove(i)
        if len(n_nodes) == 0:
            stables.append(i)
        for n in s.link[switch_index(1,mode)]:
            if n in topo_dict:
                n_nodes.add(topo_dict[n])
    for s in chromosome:
        s.stability=0
    for s in stables:
        topological_order[s][2].stability=1
    
    bubbles = [[s[2] for s in topological_order[stables[i]+1:stables[i+1]]] for i in range(len(stables)-1)]
    # while True:
    #     mode = coordinates[u][1]
    #     bubble=[]
    #     next = []
    #     while True:
    #         for uu in u.link[switch_index(1,mode)]:
                
    #             print(str(uu), coordinates[uu][1])
    #         break
    #     break

    all_path_bubbles = []
    for i in range(len(stables)-1):
        s1,s2 = topological_order[stables[i]][2], topological_order[stables[i+1]][2]
        all_path_bubbles.append(get_all_path(s1,s2,coordinates[s1][1],set(),[],coordinates,[]))
    for all_path in all_path_bubbles:
        # print("------new--------")
        for p,path in enumerate(all_path):
            fusion = list(get_num(*path[0])[-(k-1):])
            # print(len(list(get_num(*path[0])[-(k-1):])))
            for i in range(1,len(path)-1):
                fusion = fusion+list(get_num(*path[i])[(k-1):])
                # print(len(list(get_num(*path[i])[(k-1):])))
            all_path[p]=[all_path[p],fusion]
            # print(len(fusion))
    
    g = get_gt_graph(chromosome)
    g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_stability_2.graphml"))

    # Compute POA and set list consensus in each bubbles and merge theme with stables unitigs
    poa_list = []
    consensus_set_poa = [set([c]) for c in get_num(*(topological_order[stables[0]][1:][::-1]))]
    for i,bubble in enumerate(all_path_bubbles):
        sequences = [s[1] for s in bubble]
        print([s.id for s in sequences])
        g = generate_poa_graph(sequences)
        poa_list.append(g)
        # g.generateAlignmentStrings()
        file = "res_all_path"
        if i in [0,3,7,12,17,19,21,23,27,29,30,31,33,35,37,38]:
            file = os.path.join(file,"best_poa_all_path_{}.html".format(i))
        else:
            file = os.path.join(file,"simple_poa_all_path_{}.html".format(i))
        with open(file,"w") as f:
            g.htmlOutput(f)
        l = g.generateAlignment()
        for j,ll in enumerate(l):
            if ll == set([-1570, -542, None]):
                print(i)
        consensus_set_poa = consensus_set_poa + l[(k-1):]
        consensus_set_poa = consensus_set_poa + [set([c]) for c in get_num(*topological_order[stables[i+1]][1:][::-1])][(k-1):]
    consensus_set_poa = consensus_set_poa[:-(k-2)]
    # print(len(consensus_set_poa))

    # g = poa_list[38]
    # for a in g.generateAlignmentStrings():
    #     print(a[1])
    # print(g.generateAlignment())
    # poa_list = []
    # for k,bubble in enumerate(bubbles):
    #     print(k)
    #     sequences = []
    #     for unitig in bubble:
    #         s = unitig.num() if coordinates[unitig][1]==1 else unitig.num(canonical=False)
    #         sequences.append(s) 
    #     # test_order_of_alignment_case1()
    #     g = generate_poa_graph(sequences)
    #     poa_list.append(g)
    #     # g.generateAlignmentStrings()
    #     with open("poa_{}.html".format(k),"w") as f:
    #         g.htmlOutput(f)
    get_num(*topological_order[stables[0]][1:][::-1])

    g = POAGraph(seq=get_num(*topological_order[stables[0]][1:][::-1]))
    pn = g._nextnodeID-1
    for i,poa_g in enumerate(poa_list):
        start=True
        edges = []
        alignement_class = []
        old_new_id = {}
        # l = poa_g.generateAlignment()
        for cn in poa_g.nodeidlist[(k-2):]:
            node = poa_g.nodedict[cn]
            if start:
                nid = pn
                start=False
            else:
                nid = g.addNode(node.base)
            for e in node.outEdges.values():
                # if i==37:
                #     print(e)
                edges.append((e.outNodeID,e.inNodeID))
            if len(node.alignedTo)>0:
                found=False
                for c in alignement_class:
                    if cn in c:
                        found=True
                        break
                    for alignedID in node.alignedTo:
                        if alignedID in c:
                            found=True
                            break
                if not found:
                    alignement_class.append(set(node.alignedTo+[cn]))
                else:
                    c.update(node.alignedTo+[cn])
            old_new_id[cn] = nid
        # print(i)
        # if i==37:
        #     print(edges)
        for e in edges:
            g.addEdge(*[old_new_id[i] for i in e],None)
        for c in alignement_class:
            for cn in c:
                g.nodedict[old_new_id[cn]].alignedTo = [old_new_id[alignedID] for alignedID in c if alignedID!=cn]
                # print("ALIGN",old_new_id[cn], g.nodedict[old_new_id[cn]].alignedTo)
        pn=nid
        if i==37 or i==36:
            print(get_num(*topological_order[stables[i+1]][1:][::-1]))
        for c in get_num(*topological_order[stables[i+1]][1:][::-1])[(k-1):]:
            nid = g.addNode(c)
            g.addEdge(pn, nid,None)
            pn = nid
    g.toposort()

    g_gt = get_gt_graph_poa(g)
    g_gt.save(os.path.join(RES_OUTPUT_FOLDER,"full_poa.graphml"))


    file = os.path.join(RES_OUTPUT_FOLDER,"full_poa.html")
    with open(file,"w") as f:
        g.htmlOutput(f)
    # topo_order = np.zeros(g._nextnodeID)
    consensus_set_poa_2=[]
    for pnode in g._simplified_graph_rep():
        consensus_set_poa_2.append(set())
        for nid in pnode.node_ids:
            consensus_set_poa_2[-1].add(g.nodedict[nid].base)
    for i,(s1,s2) in enumerate(zip(consensus_set_poa,consensus_set_poa_2)):
        if s1!=s2:
            if s1.symmetric_difference(s2) != set([None]):
                print(s1,s2,"NON")
                break
    ## Align read on consensus
    ### Rough localisation
    l_thresh, u_thresh = 0.5,2

    read_sets_cano = [set(read.num()) for read in read_seqs]
    read_sets_ncano = [set(read.num(canonical=False)) for read in read_seqs]
    # consensus_reads_coordinates = np.zeros((len(read_seqs),len(consensus_set_coordinates)))
    # consensus_reads_poa = np.zeros((len(read_seqs),len(consensus_set_poa)))
    consensus_reads_coordinates = []
    consensus_reads_poa = []
    main_cano = []
    for consensus_reads, consensus_set, norm in zip([consensus_reads_coordinates,consensus_reads_poa],[consensus_set_coordinates,consensus_set_poa],[2,1]):
        main_cano_list = []
        for i, (rs, rsn) in enumerate(zip(read_sets_cano,read_sets_ncano)):
            a = [len(cs.intersection(rs))>0 for cs in consensus_set]
            an = [len(cs.intersection(rsn))>0 for cs in consensus_set]
            iaf, af = max(enumerate([a,an]),key= lambda x: sum(x[1]))
            # saf= get_max_block(af)[0]/norm
            lread = len(read_seqs[i])
            stored=False
            if True: #saf>(l_thresh*lread) and saf<(u_thresh*lread):
                consensus_reads.append(af)
                stored = True
            main_cano_list.append([iaf,sum(a),sum(an),stored])
        main_cano.append(main_cano_list)
    consensus_reads_coordinates = np.array(consensus_reads_coordinates)
    consensus_reads_poa = np.array(consensus_reads_poa)

    ratio_total_poa = [sum(consensus_reads_poa[i])/len(s) for i, s in enumerate(read_seqs)]
    max_block2_poa = [get_max_block(c,2,cyclic=True) for c in consensus_reads_poa]

    consensus_reads_poa_block2 = np.zeros(consensus_reads_poa.shape)
    id_max = consensus_reads_poa_block2.shape[1]
    median_maxblock = []
    len_maxblock = []
    for i,(_,s,e,l) in enumerate(max_block2_poa):
        median_maxblock.append((s+l/2,l,s,e,i))
        len_maxblock.append(l)
        consensus_reads_poa_block2[i,s:(e+1)]=1
        if s+l>id_max:
            consensus_reads_poa_block2[i,0:(id_max-(s-l))]=1

    consensus_reads_poa_block2_l = [x for _, x in sorted(zip(median_maxblock, consensus_reads_poa_block2))]

    ### Classification (in_chr vector)
    ratio_max_block2_poa = [mb[0]/len(s) for mb, s in zip(max_block2_poa, read_seqs)]
    y_block2_poa = np.array(ratio_max_block2_poa)>0.4

    ### Align on sub-sequence of the consensus POA
    in_chr_list = (y_block2_poa==1)
    start_end_block_list = [(s,e,l) for _,s,e,l in max_block2_poa]
    is_canonical_list = np.array([mc[0] for mc in main_cano[1]])==0

    
    id_list_of_parts = [[] for _ in range(len(consensus_set_poa_2))]
    id_max = len(id_list_of_parts)
    for read_id, (in_chr, (start,end,length), is_canonical,s) in enumerate(zip(in_chr_list,start_end_block_list,is_canonical_list, read_seqs)):
        if in_chr:
            a = get_num(s,1*is_canonical)
            b = consensus_set_poa_2[start:(start+length)]
            if start+length>id_max:
                b+=consensus_set_poa_2[0:(id_max-(start+length))]
            alignment, alignment_score, _ = needleman_wunsch_inclusion(a,b)
            i, j = 0,0
            show = len(a)!= len(b)
            show = False
            if show:
                print(a)
                print(b)
                print(alignment)
            for a in alignment:
                ma = True
                if a[0]=="*":
                    j+=1
                    ma = False
                if a[1]=="*":
                    i+=1
                    ma = False
                if ma:
                    if show:
                        print(read_id,is_canonical,j, i)
                    id_list_of_parts[(start+j)%id_max].append((read_id,is_canonical,i))
                    i+=1
                    j+=1            
            # if read_id>1000:
            #     break
    coverage = [len(l) for l in id_list_of_parts]
    
    gene_margin_funs = {"short":lambda x,i: x[1-i], "long": lambda x,i: x[i], 'middle': lambda x,i: (x[0]+x[1])//2}

    t = time()
    margin_mode = ("next_gene","middle")
    if margin_mode is not None and margin_mode[0]=="next_gene":
        gene_margin_fun = gene_margin_funs[margin_mode[1]]
    # read_ids = {read_name:i for i, read_name in enumerate(read_data)}
    read_names = [read_name for read_name in read_data]
    marginmode = "_".join([str(x) for x in margin_mode])
    for i, block_list in enumerate(id_list_of_parts):
        block=[]
        for read_id, is_canonical, i in block_list:
            read_name = read_names[read_id]
            s_read = fastq[read_name].seq
            read, read_pos, read_seq = read_data[read_name], gene_positions_reads[read_name], read_seqs[read_id]
            # print(len(read),len(read_seq))
            if is_canonical^read_seq.original_order:
                rev = True
                pos = -i-1
            else:
                pos = i
                rev = False
            block_b, block_e = read_pos[pos]
            if margin_mode is not None:
                match margin_mode[0]:
                    case "nucleotides":
                        p_b, p_e= max(0,block_b - margin_mode[1]), min(block_e + margin_mode[1],len(s_read)-1)
                    case "next_gene":
                        lrs = len(read_seq)
                        if pos==0 or pos==-lrs:
                            p_b = block_b
                        else:
                            p_b = gene_margin_fun(read_pos[pos-1],0)
                        if pos==(lrs-1) or pos==[-1]:
                            p_e = block_e
                        else:
                            p_e = gene_margin_fun(read_pos[pos+1],1)
                        
                # print(block_b,b,block_e,e)
                block_b, block_e = min(block_b,p_b), max(block_e,p_e)
            # print(block_b, block_e)            
            # print(read[pos],read[-pos-1],read_seq.num()[pos], bi_alphabet[0][abs(read_seq.num()[pos])-1],rev)
            # print(block_b,block_e, len(s_read))
            s = s_read[block_b:(block_e+1)]
            if is_canonical^read_seq.original_order:
                s = pyfastx.reverse_complement(s)
            block.append(SeqRecord.SeqRecord(SeqRecord.Seq(s),id=read_name,name=read_name, description=read_name))
        with open(os.path.join(RES_FASTA_FOLDER,"{}_{}.fasta".format(marginmode,i)), "w") as output_handle:
            SeqIO.write(block, output_handle, "fasta")
    print("Time for MSA fasta: {:0=2.0f}:{:0=2.0f}:{:0=2.2f}".format(*hms(time()-t)))
    
    print("\nTOTAL TIME: {:0=2.0f}:{:0=2.0f}:{:0=2.2f}".format(*hms(time()-t_start)))