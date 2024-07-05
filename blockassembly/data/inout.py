import gfapy
import csv
import pickle

from data.data import Sequence
from graph.graph import Graph

## EXPORT

def create_gfa_csv(output_file, g,k,vp=["id"],ep=[]):
    gfa = gfapy.Gfa(vlevel=0)
    with open(output_file.format(".csv"), 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Name"]+vp)
        for v in g.vertices():
            gfa.add_line("S\t{}\t{}\tLN:i:{}\tdp:i:{}\tid:i:{}".format(g.vp["z_seq_dot"][v],"*",g.vp["len"][v],g.vp["abundance"][v],g.vp["id"][v]))
            # gfa.add_line("S\t{}\t{}\tLN:i:{}".format(g.vp["z_seq_dot"][v],"a"*g.vp["len"][v],g.vp["len"][v]))
            # gfa.add_line("S\t{}\t{}\tLN:i:{}".format(int(v),"a"*g.vp["len"][v],g.vp["len"][v]))
            csvwriter.writerow([g.vp["z_seq_dot"][v],*[g.vp[vpp][v] for vpp in vp]])
    for e in g.edges():
        r = g.ep["edge_type"][e]
        vs,vt = e.source(), e.target()
        # s,t = g.vp["z_seq"][vs], g.vp["z_seq"][vt]
        s,t = int(vs), int(vt)
        s,t = g.vp["z_seq_dot"][vs], g.vp["z_seq_dot"][vt]
        match r:
            case 1:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"+",t,"+",k))
            case -1:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"-",t,"-",k))
            case 2:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"+",t,"-",k))
            case -2:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"-",t,"+",k))
    gfa.to_file(output_file.format(".gfa"))

def save_sequences(sequences,filename):
    with open(filename, 'wb') as f:
        pickle.dump([(s.seq, s.abundance,[[(ss.id, edge) for ss,edge in d.items()] for d in s.link], s.can_concatenate) for s in sequences], f)

def load_sequences(filename):
    with open(filename, 'rb') as f:
        sequence_tuples = pickle.load(f)
    sequences = [Sequence(i,seq,a) for i, (seq,a,_,_) in enumerate(sequence_tuples)]
    for k,s in enumerate(sequences):
        for lr in [0,1]:
            s.link[lr] = {sequences[i]:edge for i, edge in sequence_tuples[k][2][lr]}
        s.can_concatenate = sequence_tuples[k][2]
    return Graph({s:s for s in sequences})
    