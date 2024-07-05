from graph.graph import switch_index
def add_unitigs_sets(unitigs,components,key):
    for u in unitigs:
        if u not in components[key]:
            components[key][u]=u
            components = add_unitigs_sets({uu:uu for uu in u.link[0]}, components, key=key)
            components = add_unitigs_sets({uu:uu for uu in u.link[1]}, components, key=key)
    return components

def go_through(seen, rp, n, mode, verbose=False):
    if verbose:
        print(n.id,mode)
    if (n,mode) in rp:
        return True
    if (n,mode) in seen:
        return False
    rp[(n,mode)]=None
    seen[(n,mode)]=None
    for next_n, edge_list in n.link[switch_index(1,mode)].items():
        for edge_mode in edge_list:
            next_mode = mode*edge_mode
            if go_through(seen,rp,next_n, next_mode, verbose=verbose):
                return True
    _ = rp.pop((n,mode))
    return False
def check_cycles(graph, verbose=False):
    seen = {}
    rec_pile = {}
    for u in graph:
        if go_through(seen,rec_pile,u,1, verbose=verbose):
            return True, seen
        if go_through(seen,rec_pile,u,-1, verbose=verbose):
            return True, seen
    return False, seen