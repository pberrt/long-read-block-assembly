# utils.py

import argparse

import numpy as np
import numpy.random as rd

def get_args():
    parser = argparse.ArgumentParser(
                    prog='DBG multi-k',
                    description='')
    parser.add_argument('exp')
    args = parser.parse_args()
    return args


def rev_comp(numseq):
    return -numseq[::-1]


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

def num2seq(num, bi_alphabet):
    # print(bi_alphabet)
    # print(num)
    # print([(abs(n)-1,(-np.sign(n)+1)//2) for n in num])
    sep = "~~~" if len(bi_alphabet[0][0])>1 else ""
        
    # seq = sep.join(np.array([bi_alphabet[0][abs(n)-1][(-np.sign(n)+1)//2] for n in num]))
    seq = [bi_alphabet[0][abs(n)-1][(-np.sign(n)+1)//2] for n in num]
    return seq

def numseq2bytes(ns , n_b=2):
    # print("nb numseq2b ", ns, n_b)
    bs = b''
    for n in ns:
        bs = bs + int(n).to_bytes(n_b, "big", signed=True)
    return bs

def bytes2numseq(bs, n_b=2):
    # print("nb b2numseq ", n_b)
    match n_b:
        case 1:
            dtype = np.int8
        case 2:
            dtype = np.int16
    ns = np.array([int.from_bytes(bs[k*n_b:(k+1)*n_b],"big",signed=True) for k in range(len(bs)//n_b)], dtype=dtype)
    # print(ns)
    return ns
    