# visu.py

import numpy as np

import matplotlib.pyplot as plt
import cairo

import networkx as nx
import graph_tool.all as gt



# def plot_debruijn_graph_toyplot(edges, kmers, file = None, width=1000, height=1000):
#     "returns a toyplot graph from an input of edges"
#     graph = toyplot.graph(
#         np.array(list(edges)),
#         [k for k in kmers],
#         width=width,
#         height=height,
#         tmarker=">", 
#         vsize=25,
#         vstyle={"stroke": "none", "stroke-width": 2, "fill": "none"},
#         vlstyle={"font-size": "11px","fill":"white"},
#         estyle={"stroke": "blue", "stroke-width": 2},
#         layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
#     graph[0].style= {"background-color":"dimgrey"}
#     if file is not None:
#         toyplot.html.render(graph[0], file)
#     return graph



# def plot_compacted_debruijn_graph_toyplot(edges, unitigs, file = None, width=500, height=500):
#     "returns a toyplot graph from an input of edges"
#     if len(edges)==0:
#         n=len(unitigs)
#         edges = [[i,i] for i in range(n)]
#     graph = toyplot.graph(
#         np.array(list(edges)),
#         [k for k in range(len(unitigs))],
#         width=width,
#         height=height,
#         tmarker=">", 
#         vsize=25,
#         vstyle={"stroke": "none", "stroke-width": 1, "fill": "black"},
#         vlstyle={"font-size": "11px","fill":"white"},
#         estyle={"stroke": "red", "stroke-width": 2},
#         layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
#     graph[0].style= {"background-color":"dimgrey"}
#     if file is not None:
#         toyplot.html.render(graph[0], file)
#     return graph

def plot_dbg_axis(G,pos,labels, ref_seq, ax):
    font = {"color": "k", "fontsize": 12}
    # nx.draw_networkx_nodes(G, pos, ax=ax, node_size=800)
    nx.draw_networkx_labels(G, pos, labels=labels, font_color="#1f78b4", font_weight="bold", ax=ax)
    nx.draw_networkx_edges(
        G, pos,
        # connectionstyle="arc3,rad=0.1",  # <-- THIS IS IT
        ax=ax
    )
    if ref_seq is not None:
        max_char=63
        if len(ref_seq)>63:
            ref_lab=ref_seq[:30]+"..."+ref_seq[-30:]
        else:
            ref_lab = ref_seq
        ax.text(
                    0,
                    0.98,
                    "ref: "+ref_lab,
                    horizontalalignment="left",
                    transform=ax.transAxes,
                    fontdict=font
                )
    ax.margins(0, 0)

def plot_debruijn_graph(edges, kmers, file = None, ref_seq = None, save=True, axext=None):
    "returns a toyplot graph from an input of edges"
   # Compute position of nodes
    print("visu: ",kmers, edges)
    G=nx.DiGraph()
    for k, kmer_name in enumerate(kmers):
        G.add_nodes_from([(k,{"name":kmer_name,"nb":k+1})])
    G.add_edges_from(edges)
    pos = nx.kamada_kawai_layout(G)
    labels = nx.get_node_attributes(G,"name")
    print(labels)
    # labels = {k:("..."+ll[-1]) for k,ll in nx.get_node_attributes(G,"name").items()}

    if save:
        # Draw nodes and edges
        fig, ax = plt.subplots(figsize=(20, 15))
        plot_dbg_axis(G,pos,labels, ref_seq, ax)
        fig.tight_layout()
        plt.axis("off")
        fig.savefig(file)
    if axext is not None:
        plot_dbg_axis(G,pos,labels, ref_seq, axext)
        axext.axis("off")

    return ax

def plot_cdbg_axis(G, pos, labels, unitigs, ref_seq, ax):
    
    font = {"color": "k", "fontsize": 12}

    # Draw nodes and edges
    nx.draw_networkx_nodes(G, pos, ax=ax)
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax)
    nx.draw_networkx_edges(
        G, pos,
        connectionstyle="arc3,rad=0.1",  # <-- THIS IS IT
        ax=ax
    )

    max_unitigs = 10
    n = min(len(unitigs),max_unitigs)

    for k, u in enumerate(unitigs[0:max_unitigs]):
        if k==(max_unitigs-1) and n!=len(unitigs):
            lab = "..."
        elif len(u)>23:
            lab = u[:10]+"..."+u[-10:]
        else:
            lab = u
        ax.text(
            1,
            0.04*(n-k-0.5),
            str(k+1)+": "+lab,
            ha="right",
            transform=ax.transAxes,
            fontdict=font
        )
    if ref_seq is not None:
        max_char=63
        if len(ref_seq)>63:
            ref_lab=ref_seq[:30]+"..."+ref_seq[-30:]
        else:
            ref_lab = ref_seq
        ax.text(
                    0,
                    1,
                    "ref: "+ref_lab,
                    ha="left", va="top",
                    transform=ax.transAxes,
                    fontdict=font
                )
    # ax.margins(0, 0)


def plot_compacted_debruijn_graph(edges, unitigs, file = None, ref_seq = None, save=True, axext=None):
    # Compute position of nodes
    G=nx.DiGraph()
    for k, unitig in enumerate(unitigs):
        G.add_nodes_from([(k,{"name":unitig,"nb":k+1})])
    G.add_edges_from(edges)
    pos = nx.kamada_kawai_layout(G)
    labels = nx.get_node_attributes(G,"nb")

    if save:
        # Draw nodes and edges
        fig, ax = plt.subplots(figsize=(20, 15))
        plot_cdbg_axis(G,pos,labels, unitigs, ref_seq, ax)
        fig.tight_layout()
        plt.axis("off")
        fig.savefig(file)
    if axext is not None:
        plot_cdbg_axis(G,pos,labels, unitigs, ref_seq, axext)
        axext.axis("off")

    return ax
    
def plot_both_graphs(kmers, edges, unitigs, c_edges, name, ref_seq = None):
    graph = plot_debruijn_graph(edges, kmers, name+"graph.png", ref_seq=ref_seq)
    plot_compacted_debruijn_graph(c_edges, unitigs, name+"compacted_graph.png", ref_seq=ref_seq)

    return graph

def plot_dbg_axis_gt(g,pos,kmers, ref_seq, ax):
    font = {"color": "k", "fontsize": 12}
    # nx.draw_networkx_nodes(G, pos, ax=ax, node_size=800)
    vlabs = g.new_vp("string", vals=kmers) 
    gt.graph_draw(g, pos, vertex_text=vlabs, mplfig=ax,
                  vertex_size=0.3, vertex_anchor=0, vertex_aspect=max(1,1+(len(kmers[0])-3)*0.1), vertex_shape="none", edge_mid_marker="arrow", edge_end_marker="none", edge_color=[0, 0, 0, 0.3],vertex_color=[0,0,0,0], vertex_font_weight= cairo.FONT_WEIGHT_BOLD,vertex_text_color="red", vertex_fill_color=[0,0,0,0])
    if ref_seq is not None:
        max_char=63
        if len(ref_seq)>63:
            ref_lab=ref_seq[:30]+"..."+ref_seq[-30:]
        else:
            ref_lab = ref_seq
        ax.text(
                    0,
                    0.98,
                    "ref: "+ref_lab,
                    horizontalalignment="left",
                    transform=ax.transAxes,
                    fontdict=font
                )
    ax.margins(0, 0)

def plot_debruijn_graph_gt(g, kmers, file = None, ref_seq = None, save=True, axext=None):
    "returns a toyplot graph from an input of edges"
    # Compute position of nodes
    pos = gt.fruchterman_reingold_layout(g)

    if save:
        # Draw nodes and edges
        fig, ax = plt.subplots(figsize=(20, 15))
        plot_dbg_axis_gt(g,pos,kmers, ref_seq, ax)
        fig.tight_layout()
        plt.axis("off")
        fig.savefig(file)
    if axext is not None:
        plot_dbg_axis_gt(g,pos,kmers, ref_seq, axext)
        axext.axis("off")

    return ax

def plot_cdbg_axis_gt(g, pos, unitigs, ref_seq, ax, plot_graph=True):
    
    font = {"color": "k", "fontsize": 12}
    
    plot_graph = plot_graph and g.num_vertices()>1

    if plot_graph:
        vsize = g.new_vp("double") 
        for k, u in enumerate(unitigs):
            vsize[k] = len(u)
        vmin,vmax = min(vsize), max(vsize)
        smin, smax = 0.02, 0.32
        for i,vs in enumerate(vsize):
            # print(vs, vmin)
            vsize[i] = smin+np.sqrt(vs-vmin)/np.sqrt(vmax-vmin+1e-6)*(smax-smin)
        vlabs=g.new_vp("int", vals=[l+1 for l in g.vertex_index])
        gt.graph_draw(g, pos, vertex_size=vsize, vertex_font_size=0.15, vertex_text=vlabs, mplfig=ax)
        

    max_unitigs = 10
    n = min(len(unitigs),max_unitigs)

    for k, u in enumerate(unitigs[0:n]):
        max_char = 23 if plot_graph else 63
        if k==(max_unitigs-1) and n!=len(unitigs):
            lab = "..."
        elif len(u)>max_char:
            lab = u[:((max_char-3)//2)]+"..."+u[-((max_char-3)//2):]
        else:
            lab = u
        if plot_graph:
            ax.text(
                1,
                0.04*(n-k-0.5),
                str(k+1)+": "+lab,
                ha="right",
                transform=ax.transAxes,
                fontdict=font
            ) 
        else:
            ax.text(
                0.1,
                0.5+0.04*(n-k-0.5),
                str(k+1)+": "+lab,
                ha="left",
                transform=ax.transAxes,
                fontdict=font
            )
    if ref_seq is not None:
        max_char=63
        if len(ref_seq)>max_char:
            ref_lab=ref_seq[:((max_char-3)//2)]+"..."+ref_seq[-((max_char-3)//2):]
        else:
            ref_lab = ref_seq
        ax.text(
                    0,
                    1,
                    "ref: "+ref_lab,
                    ha="left", va="top",
                    transform=ax.transAxes,
                    fontdict=font
                )
    # ax.margins(0, 0)


def plot_compacted_debruijn_graph_gt(c_g, unitigs, file = None, ref_seq = None, save=True, axext=None):
    # Compute position of nodes
    
    pos = gt.fruchterman_reingold_layout(c_g, grid=False)
    if save:
        # Draw nodes and edges
        fig, ax = plt.subplots(figsize=(20, 15))
        plot_cdbg_axis_gt(c_g,pos, unitigs, ref_seq, ax)
        fig.tight_layout()
        plt.axis("off")
        fig.savefig(file)
    if axext is not None:
        plot_cdbg_axis_gt(c_g,pos, unitigs, ref_seq, axext)
        axext.axis("off")

    return ax
