import graph_tool as gt
import warnings
import numpy as np
from graph.graph import switch_index
from itertools import product

def sort_by_median(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        m = np.median(np.argwhere(x).flatten())
        if np.isnan(m):
            m=-1
    return m
def get_max_block(x, wait=0, cyclic=False):
    c = [0,-1,-1]
    n_c = [0,-1,-1]
    new= True
    stop = False
    w = wait
    lx = len(x)
    i=0
    while True:
        xx= x[i]
        if xx==1 or (w>0 and not new):
            if xx==1:
                if n_c[1]==i:
                    break
                else:
                    n_c[2]=i
                    n_c[0]+=1
                w=wait
            else:
                w-=1
            if new:
                n_c[1]=i
                new=False
        else:
            if stop:
                break
            if not new:
                c = max(c,n_c,key=lambda x: x[0])
                n_c = [0,-1,-1]
            new = True
        i+=1
        if i>=lx:
            if cyclic:
                i=0
                stop = True
            else:
                break

    c = max(c,n_c,key=lambda x: x[0])
    c.append(c[2]+1-c[1])
    if c[3]<=0:
        c[3]+=lx
    return c
def get_max_block_mean(x, wait=0, cyclic=False):
    _,s,e, l = get_max_block(x, wait, cyclic)
    return s+(l-1)/2

def get_gt_graph_poa(g_poa):
        g = gt.Graph()
        g.add_vertex(g_poa._nextnodeID)
        islink = g.new_edge_property("int")
        est = g.new_edge_property("string")
        ett = g.new_edge_property("string")
        elt = g.new_edge_property("string")
        ect = g.new_edge_property("string")
        names =[g_poa.nodedict[nID].base for nID in range(g_poa._nextnodeID)]
        for nID in g_poa.nodeidlist:
            node = g_poa.nodedict[nID]
            edges = node.outEdges
            for nextID in edges:
                g.add_edge(nID, nextID)
            for alignedID in node.alignedTo:
                if nID<alignedID:
                    g.add_edge(nID, alignedID)
                    islink[(nID, alignedID)]=1
                    est[(nID,alignedID)] = "None"
                    ett[(nID,alignedID)] = "None"
                    elt[(nID,alignedID)] = "Equal-Dash"
                    ect[(nID,alignedID)] = "#FF0000"
        for e in g.edges():
            if not islink[e]:
                islink[e], est[e], ett[e], elt[e], ect[e] = 0, "None", "Arrow", "Solid", "#000000"
        topo_order = g.new_vertex_property("int")
        for pnode in g_poa._simplified_graph_rep():
            for nid in pnode.node_ids:
                topo_order[nid]=pnode.pnode_id
        vid=g.new_vp("int",vals = [nID for nID in range(g_poa._nextnodeID)])
        vname=g.new_vp("int",vals= names)
        g.vp["id"] = vid
        g.vp["base"] = vname
        g.vp["out_degree"] = g.new_vp("int", vals=[v.out_degree() for v in g.vertices()])
        g.vp["in_degree"] = g.new_vp("int", vals=[v.in_degree() for v in g.vertices()])
        g.vp["topo_order"] = g.new_vp("int", vals=topo_order)
        g.ep["islink"] = islink
        g.ep["est"] = est
        g.ep["ett"] = ett
        g.ep["elt"] = elt
        g.ep["ect"] = ect
        return g

def get_all_path(s1,s2, mode, visited, path, coordinates, res):
    path.append((s1,mode))
    visited.add(s1)
    if s1 == s2:
        res.append(path[::])
    else:
        for n_s1 in s1.link[switch_index(1,mode)]:
            if n_s1 not in visited and n_s1 in coordinates:
                n_mode = coordinates[n_s1][1]
                res = get_all_path(n_s1, s2, n_mode, visited,path, coordinates, res)
    
    path.pop()
    visited.remove(s1)
    return res

def add_coordinate(s, mode, coordinates,c, k):
    # print(s.id)
    if s not in coordinates or c> coordinates[s][0][0]:
        # if all([prev_s in coordinates for prev_s in s.link[switch_index(0,mode)]]) or c==0:
            # if len(s.link[switch_index(0,mode)])==0 or c==0:
            # if c==0:
            #     new_c = 0
            # else:
            #     new_c = max([coordinates[prev_s][0][1] - (k-1) for prev_s in s.link[switch_index(0,mode)]])
            # coordinates[s]=([new_c,new_c+len(s)],mode)
        coordinates[s]=([c,c+len(s)],mode)
        # print(s.id)
        # c = c+1
        # for next_s, edge_modes in s.link[switch_index(0,mode)].items():
        #     for edge_mode in edge_modes:
        #         c = max(c, coordinates[next_s][0]+1)
        # if s.id in [112,132]:
        #     print([coordinates[prev_s][0] for prev_s in s.link[switch_index(0,mode)]])
        #     print(s.id, c, new_c)
        for next_s, edge_modes in s.link[switch_index(1,mode)].items():
            for edge_mode in edge_modes:
                next_mode = mode*edge_mode
                coordinates = add_coordinate(next_s, next_mode, coordinates, c + len(s) - (k-1)+1,k)
            
    return coordinates

def needleman_wunsch_inclusion(x, y):
    # get the lengths of x and y
    N, M = len(x), len(y)
    gap_score = 0.5
    def score(a, b):
        # Scoring function: returns 1 if elements are equal, 0 otherwise
        if a in b:
            return 1, a, a
        else:
            return -1, a, "M"

    # Direction constants for traceback
    DIAG, LEFT, UP = (-1, -1), (-1, 0), (0, -1)
    # Initialize score (F) and pointer (Ptr) matrices
    F, Ptr, Elem = {}, {}, {}
    F[-1, -1] = 0
    # Initial scoring for gaps along x
    for i in range(N):
        # F[i, -1] = -i
        F[i, -1] = 0
    # Initial scoring for gaps along y
    for j in range(M):
        # F[-1, j] = -j
        F[-1, j] = 0
    # Option for Ptr to trace back alignment
    option_Ptr = DIAG, LEFT, UP
    # Fill F and Ptr tables
    for i, j in product(range(N), range(M)):
        # Score options: match/mismatch, gap in x, gap in y
        gap_x = 0 if j == M-1 else gap_score
        gap_y = 0 if i == N-1 else gap_score
        s, e1, e2 = score(x[i], y[j])
        option_F = (
            F[i - 1, j - 1] + s,  # Match/mismatch
            F[i - 1, j] - gap_x,  # Gap in x
            F[i, j - 1] - gap_y,  # Gap in y
        )
        option_Elem = [(e1,e2),(e1,"*"),("*",e2)]
        # Choose best option for F and Ptr
        F[i, j], Ptr[i, j], Elem[i,j] = max(zip(option_F, option_Ptr, option_Elem))
    # Trace back to get the alignment
    alignment_score = F[N-1,M-1]
    alignment = []
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        # Add aligned elements or gaps based on direction
        if direction == DIAG:
            element = x[i]
        elif direction == LEFT:
            element = x[i] # Insert gap in y
        elif direction == UP:
            element = "*"  # Insert gap in x
        alignment.append(Elem[i,j])
        di, dj = direction
        i, j = i + di, j + dj
    # Add remaining gaps if any
    while i >= 0:
        alignment.append((x[i],"*"))  # Gap in y
        i -= 1
    while j >= 0:
        alignment.append(("*","M"))  # Gap in x
        j -= 1
    return alignment[::-1], alignment_score,F