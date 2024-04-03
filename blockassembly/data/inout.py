import gfapy
import csv






## EXPORT

def create_gfa_csv(output_file, g,k,vp=["id"],ep=[]):
    gfa = gfapy.Gfa(vlevel=0)
    with open(output_file.format(".csv"), 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Name"]+vp)
        for v in g.vertices():
            gfa.add_line("S\t{}\t{}\tLN:i:{}\tdp:i:{}\tid:i:{}".format(g.vp["seq_dot"][v],"*",g.vp["len"][v],g.vp["abundance"][v],g.vp["id"][v]))
            # gfa.add_line("S\t{}\t{}\tLN:i:{}".format(g.vp["seq_dot"][v],"a"*g.vp["len"][v],g.vp["len"][v]))
            # gfa.add_line("S\t{}\t{}\tLN:i:{}".format(int(v),"a"*g.vp["len"][v],g.vp["len"][v]))
            csvwriter.writerow([g.vp["seq_dot"][v],*[g.vp[vpp][v] for vpp in vp]])
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
    gfa.to_file(output_file.format(".gfa"))
        