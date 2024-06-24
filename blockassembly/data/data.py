from common.utils import bytes2numseq, num2seq, rev_comp

import numpy as np

import sys


class Sequence:
    bi_alphabet = None
    n_b = None

    def __init__(self, i, seq, abundance, name = None):
        self.seq=seq
        self.len = len(self.seq)//self.n_b
        self.id=i
        self.order = None
        self.name = name
        self.topo_order =-1
        self.abundance = abundance
        self.can_concatenate = [True,True]
        self.link = [{},{}]
        self.is_tip=None
        self.stability=-1
        self.original_order = True
        self.compute_revcomp()
        
    def __hash__(self) -> int:
        return hash(self.seq)
    def __len__(self) -> int:
        return len(self.seq)//self.n_b
    def __str__(self):
        return num2seq(bytes2numseq(self.seq, self.n_b), self.bi_alphabet, in_str=True)
    def __repr__(self):
        return str(self.id)+" : "+str(self.seq)
    def str(self, canonical=True):
        if canonical:
            return num2seq(bytes2numseq(self.seq, self.n_b), self.bi_alphabet, in_str=True)
        else:
            return num2seq(bytes2numseq(self.rev_comp, self.n_b), self.bi_alphabet, in_str=True)

    def num(self, canonical=True):
        if canonical:
            return bytes2numseq(self.seq, self.n_b)
        else:
            return bytes2numseq(self.rev_comp, self.n_b)

    def blocks(self, canonical=True):
        if canonical:
            return num2seq(bytes2numseq(self.seq, self.n_b), self.bi_alphabet, in_str=False)
        else:
            return num2seq(bytes2numseq(self.rev_comp, self.n_b), self.bi_alphabet, in_str=False)

    def switch(self):
        self.seq, self.rev_comp = self.rev_comp, self.seq
        self.can_concatenate = self.can_concatenate[::-1]
        self.can_concatenate = self.link[::-1]
        self.original_order = not self.original_order
        
    def compute_revcomp(self):
        self.rev_comp = rev_comp(self.seq,True,True, self.n_b)
        if self.rev_comp<self.seq:
            self.switch()
        
    def __eq__(self,s):
        return self.seq ==s or (isinstance(s,Sequence) and self.seq ==s.seq)
    def check_start(self):
        self.start = (self.id == self.unitig.kmers[0].id)
    def create_unitig(self,k):
        self.unitig = Unitig(self,k)


class Kmer(Sequence):
    def __init__(self, i, kmer_b, a):
        self.start=None
        self.a_reads = 0
        self.a_unitigs = []
        super().__init__(i ,kmer_b, a)
        self.unitig = None
    # def __repr__(self) -> str:
    #     return "({},{})   ".format(self.id,self.unitig.kmers)+str(bytes2numseq(self.seq,self.n_b))+"  "+str(bytes2numseq(self.unitig.seq,self.n_b))+"  "+str(self.start)
    def compute_abundance(self, mode):
        match mode:
            case "reads":
                self.abundance=self.a_reads + len(self.a_unitigs)
            case "median_unitigs":
                self.abundance= sum(self.a_unitigs)
                if self.abundance == 0:
                    self.abundance = self.a_reads
            case  "max_unitigs_reads":
                self.abundance=max(self.a_reads + len(self.a_unitigs), sum(self.a_unitigs))


class Unitig(Kmer):
    def __init__(self, kmer,k=None):
        k = len(kmer) if k is None else k
        self.kmers=[kmer]
        self.a_list = [kmer.abundance]*(len(kmer)-k+1)
        super().__init__(kmer.id,kmer.seq,0)
    def switch(self,val=None):
        super().switch()
        self.kmers = self.kmers[::-1]
    def compute_abundance(self):
        self.abundance = np.median(self.a_list).astype(int)
    def compute_can_concatenate(self):
        self.can_concatenate = [self.kmers[0].can_concatenate[0],self.kmers[-1].can_concatenate[1]]
    