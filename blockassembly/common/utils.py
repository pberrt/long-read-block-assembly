# utils.py

import argparse

import numpy as np
import numpy.random as rd

from rapidfuzz.distance import Levenshtein

def get_args():
    parser = argparse.ArgumentParser(
                    prog='DBG multi-k',
                    description='')
    parser.add_argument('exp')
    args = parser.parse_args()
    return args


def rev_comp(seq, in_bytes = False, out_bytes =False, n_b=None):
    numseq = bytes2numseq(seq, n_b) if in_bytes else seq
    rc = -numseq[::-1]
    seq_out = numseq2bytes(rc, n_b) if out_bytes else rc
    return seq_out


def check_args(verbose=False,ref_seq = None, alphabet = None, ref_length = None, cyclic=False,
                  reads = None, read_length = None, n_reads = None, mean_coverage = None,
                  kmin = 4, kmax = None, seed = 4, shuffle = False):
    np.random.seed(seed)
    # Check ref_seq and its related variables
    if ref_seq is None:
        if (alphabet is None or ref_length is None):
            raise ValueError("Empty reference sequence requires alphabet and ref_length")
        alphabet = [("a","t"),("c","g")]
        ref_seq = create_ref(ref_length,alphabet)
        if verbose:
            print("Generate ref_seq with {} characters: ".format(ref_length), end="")
    else:
        ref_length = len(ref_seq)
        if alphabet is None:
            letters = set([l for l in ref_seq])
            if reads is not None:
                for r in reads:
                    letters |= set([l for l in r])
            letters = list(letters)
            letters.sort()
            print(letters)
            alphabet = [(l,chr(65+ord(l)-97)) for l in letters]
        if verbose:
            print("Reference sequence as input:",end="")
    if verbose:
        if ref_length>20:
            print(ref_seq[:10],"...",ref_seq[-10:],sep="")
        else:
            print(ref_seq)
        print("Used alphabet: ",alphabet)

    # Create bi_alphabet
    bi_alphabet = (alphabet,{k:((i+1)*((-1)**j)) for i, ks in enumerate(alphabet) for j, k in enumerate(ks)})

    # Check reads and its related variables
    if reads is None:
        if (read_length is None and n_reads is not None and mean_coverage is not None):
            read_length = round(ref_length*mean_coverage/n_reads)
        elif (read_length is not None and n_reads is None and mean_coverage is not None):
            n_reads = round(ref_length*mean_coverage/read_length)
        elif (read_length is not None and n_reads is not None and mean_coverage is None):
            mean_coverage = round(read_length*n_reads/ref_length) 
        else:
            raise ValueError("Empty reads requires at least two of read_length, n_reads, mean_coverage")
        reads = create_reads(seq2num(ref_seq,bi_alphabet), n_reads, read_length, cyclic=cyclic)
        if verbose:
            print([len(read) for read in reads])
            print("Generate {} reads with {} characters resulting in a mean coverage of {}.".format(n_reads, read_length, mean_coverage))
    else:
        reads = [seq2num(read, bi_alphabet) for read in reads]
        n_reads = len(reads)
        read_length = sum([len(r) for r in reads])/len(reads)
        mean_coverage = round(read_length*n_reads/ref_length)
        read_length=round(read_length)
        if verbose:
            print("{} reads are given as input with a mean coverage of {} and containing on average {} characters.".format(n_reads, mean_coverage, read_length))
    if shuffle:
        np.random.shuffle(reads)
    
    # Check kmin/kmax
    if kmax is None:
        kmax = kmin
    if kmax<kmin:
        raise ValueError("kmin should be lower than kmax")
    if verbose:
        print("Performing multi-k from k={} to k={}.".format(kmin,kmax))


    return ref_seq, reads, kmin, kmax, bi_alphabet

def create_ref(ref_length, alphabet):
    return "".join(list(rd.choice(alphabet,ref_length,replace=True)))
    
def create_read(ref_seq, read_length, cyclic=False):
    n = len(ref_seq)
    # k = rd.randint(len(ref_seq)-read_length+1)
    if cyclic:
        k = rd.randint(n)
        read = ref_seq[k:k+read_length]
        if k+read_length>n:
            read = read+ref_seq[0:min(k+read_length-n,k)]
    else:
        k = rd.randint(n)-read_length//2
        read = ref_seq[max(0,k):k+read_length]
    return read
def create_reads(ref_seq, n_reads, read_length, cyclic=False):
    return [create_read(ref_seq, read_length, cyclic) for _ in range(n_reads)]
    

def print_ref_unitigs(ref_seq,unitigs):
    n = len(ref_seq)
    dn = int(np.log10(len(unitigs)))+1
    print("R"+" "*(dn)+ref_seq)
    for k,u in enumerate(unitigs):
        l = len(u)
        found = False
        for i in range(n-l+1):
            if ref_seq[i:i+l]==u:
                found=True
                print(str(k+1)+" "*(dn-(int(np.log10(k+1))))+"|"*i+u+"|"*(n-l-i))

def unitig_classification(u1, u2):
    compact_unitigs={}
    for u in u1:
        if u in compact_unitigs:
            compact_unitigs[u].add(1)
        else:
            compact_unitigs[u]=set([1])
    for u in u2:
        if u in compact_unitigs:
            compact_unitigs[u].add(2)
        else:
            compact_unitigs[u]=set([2])

    sorted_unitigs = {"common":[], "only1":[], "only2": []}
    for u in compact_unitigs:
        if compact_unitigs[u]==set([1,2]):
            sorted_unitigs["common"].append(u)
        if compact_unitigs[u]==set([1]):
            sorted_unitigs["only1"].append(u)
        if compact_unitigs[u]==set([2]):
            sorted_unitigs["only2"].append(u)
    return sorted_unitigs

def seq2num(seq, bi_alphabet):
    dtype = np.int8 if len(bi_alphabet[0])<128 else np.int16
    num = np.array([bi_alphabet[1][l] for l in seq], dtype=dtype)
    return num

def num2seq(num, bi_alphabet, in_str = False):
    # print(bi_alphabet)
    # print(num)
    # print([(abs(n)-1,(-np.sign(n)+1)//2) for n in num])
        
    # seq = sep.join(np.array([bi_alphabet[0][abs(n)-1][(-np.sign(n)+1)//2] for n in num]))
    seq = [bi_alphabet[0][abs(n)-1][(-np.sign(n)+1)//2] for n in num]
    if in_str:
        sep = "~~~" if len(bi_alphabet[0][0])>1 else ""
        seq = sep.join(seq)
    return seq

def numseq2bytes(ns , n_b):
    # print("nb numseq2b ", ns, n_b)
    bs = b''
    for n in ns:
        bs = bs + int(n).to_bytes(n_b, "big", signed=True)
    return bs

def bytes2numseq(bs, n_b):
    # print("nb b2numseq ", n_b)
    match n_b:
        case 1:
            dtype = np.int8
        case 2:
            dtype = np.int16
        case 4:
            dtype = np.int32
        case 8:
            dtype = np.int64
    ns = np.array([int.from_bytes(bs[k*n_b:(k+1)*n_b],"big",signed=True) for k in range(len(bs)//n_b)], dtype=dtype)
    # print(ns)
    return ns
    
def compute_unitig_ref(unitigs, ref_seqs):
    ### Compute unitig ref

    gap_score = -5
    match_score = 1
    mismatch_score = -1

    u_ref = []
    u_ref_by_id = []
    for unitig in unitigs:
        both_u = (unitig.num(), unitig.num(canonical=False))
        u_ref_scores = []
        for seq in ref_seqs:
            u_ref_scores.append([])
            for u in both_u:
                if False:
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
                u_ref_scores[-1].append(Levenshtein.similarity(seq.num(),u)/len(u))
        u_ref_scores = [max(urs) for urs in u_ref_scores]
        u_ref_by_id.append(u_ref_scores)
        u_ref_scores_id = [(urs,i) for i, urs in enumerate(u_ref_scores)]
        u_ref_scores_id.sort()
        # print(u_ref_scores)
        ref_string_tuples = []
        for s,i in u_ref_scores_id[::-1]:
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
    return u_ref, list(np.array(u_ref_by_id)[:,:3].T)


def print_ops(source,dest,ops):
    s1,s2 = "",""
    for op in ops:
        match op.tag:
            case "equal":
                n = op.src_end-op.src_start
                s1+=",".join(n*["{:5d}"]).format(*list(source[op.src_start:op.src_end]))
                s2+=",".join(n*["{:5d}"]).format(*list(dest[op.dest_start:op.dest_end]))
            case "delete":
                n = op.src_end-op.src_start
                s1+=",".join(["{:5d},..({:5d})..,{:5d}"]).format(source[op.src_start],n-2,source[op.src_end-1])
                s2+="."*23
    print(s1)
    print(s2)

def hms(t):
    return t//3600, t%3600//60,t%60