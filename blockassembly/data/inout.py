import gfapy







## EXPORT

def create_gfa(g,k):
    gfa = gfapy.Gfa(vlevel=0)
    for v in g.vertices():
        gfa.add_line("S\t{}\t{}\tLN:i:{}".format(g.vp["seq_dot"][v],"*",g.vp["len"][v]))
        # gfa.add_line("S\t{}\t{}\tLN:i:{}".format(g.vp["seq_dot"][v],"a"*g.vp["len"][v],g.vp["len"][v]))
        # gfa.add_line("S\t{}\t{}\tLN:i:{}".format(int(v),"a"*g.vp["len"][v],g.vp["len"][v]))
    for e in g.edges():
        r = g.ep["edge_type"][e]
        vs,vt = e.source(), e.target()
        # s,t = g.vp["seq"][vs], g.vp["seq"][vt]
        s,t = int(vs), int(vt)
        s,t = g.vp["seq_dot"][vs], g.vp["seq_dot"][vt]
        match r:
            case 1:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"+",t,"+",k))
            case -1:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"-",t,"-",k))
            case 2:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"+",t,"-",k))
            case -2:
                gfa.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s,"-",t,"+",k))
    return gfa
        