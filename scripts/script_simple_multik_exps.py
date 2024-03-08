# script_graph_multi_k.py

import sys
import argparse

sys.path.append('..')

from  blockassembly.graph.graph import graph_multi_k
from blockassembly.data.visu import plot_debruijn_graph_gt, plot_compacted_debruijn_graph_gt
from blockassembly.data.inout import create_gfa

import os

import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(
                    prog='DBG multi-k',
                    description='')
    parser.add_argument('exp')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    PROJECT_FOLDER = os.getcwd()
    RES_OUTPUT_FOLDER = os.path.join(PROJECT_FOLDER,"..","res","simple_multik")
    
    args = get_args()

    spades = {"ref_seq":'ggacatcagata',
              "reads":["acat", "catc", "atca", "tcag", "caga", "agat", "gata", "tagg", "ggac"]}
    multik = {"ref_seq":"aaaaatcgatctcatcgaatt",
              "reads":["aaaaatcgatctc","tctcatcgaatt"]}
    artifacts = {
        "ref_seq":"cactatagtt",
        "reads":["cac","aac", "act", "cta", "ata"]
    }
    
    kwargs = {}
    match args.exp:
        case "exp_artifacts":
            kwargs = {"ref_seq":artifacts["ref_seq"],"alphabet": [("a","t"),("c","g")], "reads":artifacts["reads"], "shuffle":True}
            kmins, kmaxs = [3], [3]
        case "exp_spades":
            kwargs = {"ref_seq":spades["ref_seq"], "reads":spades["reads"], "shuffle":True}
            kmins, kmaxs = [3,4,3], [3,4,4]
        case "exp_spades_comp":
            kwargs = {"ref_seq":spades["ref_seq"], "alphabet": [("a","t"),("c","g")],"reads":spades["reads"], "shuffle":True}
            kmins, kmaxs = [3,4,3], [3,4,4]
        case "exp_multik":
            kwargs = {"ref_seq":multik["ref_seq"], "reads":multik["reads"], "shuffle":True}
            kmins, kmaxs = [5,7,5], [5,7,7]
        case "exp_multik_comp":
            kwargs = {"ref_seq":multik["ref_seq"], "alphabet": [("a","t"),("c","g")],"reads":multik["reads"], "shuffle":True}
            kmins, kmaxs = [5,7,5], [5,7,7]
        case "exp1":
            kwargs = {"ref_seq":'cddcdbdcbabbbdaadcccdccddaccdbcdbcdbdabcdadd',
                      "read_length":15, "mean_coverage":5}
            kmins, kmaxs = [4,5,6,7,8,4], [4,5,6,7,8,8]
        case "exp2":
            kwargs = {"ref_seq":"jeabjabijaigedaegibiebdgfdjgjbjecghiijcaghibhbeaifehiicgciiggfgaagjbijbcijjfachd",
                          "read_length":10, "mean_coverage":20}
            kmins, kmaxs = [3,4,5,3], [3,4,5,5]
        case "exp3":
            kwargs = {"ref_seq":"abcdefghituvwxyzjkldefmnop",
                          "read_length":10, "mean_coverage":20}
            kmins, kmaxs = [3,4,5], [3,4,5]
        case "exp4":
            kwargs = {"ref_seq":"abcdefghituvwxyzjkldefmntuvwxyzop",
                          "read_length":10, "mean_coverage":20}
            kmins, kmaxs = [3,4,5], [3,4,5]
        case "exp5":
            kwargs = {"ref_seq":"jeabjabijaigedaegibibonsaiigeebrepeatdgfbjbonsaiabrepeatidjgjbigbonsaiejebjabicg",
                          "read_length":10, "mean_coverage":20}
            kmins, kmaxs = [3,4,5,6,7,8,9,3], [3,4,5,6,7,8,9,9]
        case "exprandom":
            kwargs = {"ref_length":80, "alphabet":[(chr(97+i), chr(65+i)) for i in range(10)],
                          "read_length":10, "mean_coverage":20, "seed":10}
            kmins, kmaxs = [3,4,5,3], [3,4,5,5]
    # kwargs = {"ref_seq":"abcdefghijkl", "reads":["abcd","cdef","efgh", "ijkl"], "shuffle":True}
    # kmins, kmaxs = [3,4,3], [3,4,4] 
    if kwargs:
        nrow, ncol = len(kmins), 2
        plt.switch_backend("cairo")
        fig, axs = plt.subplots(nrow, ncol, squeeze=False, figsize=(8*ncol,5*nrow))
        for i, (kmin, kmax) in enumerate(zip(kmins, kmaxs)):
            ref_seq, reads, kmers, g, unitigs, c_g = graph_multi_k(**kwargs, kmin = kmin, kmax = kmax)
            print(kmers, unitigs)
            # plot_debruijn_graph_gt(g, [seq[0] for seq in kmers], os.path.join(RES_OUTPUT_FOLDER, args.exp+"graph.png"), ref_seq=ref_seq, axext=axs[i,0])
            axs[i,0].set_title("DB graph (kmin={}, kmax={})".format(kmin,kmax), fontweight="bold")
            # plot_compacted_debruijn_graph_gt(c_g, [seq[0] for seq in unitigs], os.path.join(RES_OUTPUT_FOLDER, args.exp+"compacted_graph.png"), ref_seq=ref_seq, axext=axs[i,1])
            axs[i,1].set_title("Compacted DB graph (kmin={}, kmax={})".format(kmin,kmax), fontweight="bold")
            g.save(os.path.join(RES_OUTPUT_FOLDER,"graph_{}_k_{}_{}.graphml".format(args.exp,kmin,kmax)))
            c_g.save(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}_k_{}_{}.graphml".format(args.exp,kmin,kmax)))
            gfa_g = create_gfa(g,kmax)    
            gfa_g.to_file(os.path.join(RES_OUTPUT_FOLDER,"graph_{}_k_{}_{}.gfa".format(args.exp,kmin,kmax)))
            gfa_c_g = create_gfa(c_g,kmax)    
            gfa_c_g.to_file(os.path.join(RES_OUTPUT_FOLDER,"c_graph_{}_k_{}_{}.gfa".format(args.exp,kmin,kmax)))
        fig.tight_layout()
        plt.axis("off")
        fig.savefig(os.path.join(RES_OUTPUT_FOLDER, args.exp+"_full.png"))
        