# /bin/bash

# Align consecutive MSA consensus

# for exp in Raw nucleotides_250 next_gene_short next_gene_long next_gene_middle ; do
for exp in next_gene_middle; do
    ls consensus_fasta/"$exp"_*.fasta |wc
    END=4071
    for i in $(seq 0 $(($END-1))); do
        fasta1="consensus_fasta/$exp"_"$i"_consensus.fasta
        fasta2="consensus_fasta/$exp"_"$((i+1))"_consensus.fasta
        echo $fasta1 $fasta2;
        minimap2 -x ava-ont  "$fasta1" "$fasta2" > minimap_paf/"$exp"_"$i"_"$((i+1))"_minimap.paf
    done
done
    # for filename in consensus_fasta/"$exp"_*.fasta"consensus_fasta/$(basename "$filename" .fasta)_consensus.fasta"; do