# /bin/bash

# Compute MSA on each fasta file (each block)

N=6

# for exp in Raw nucleotides_250 next_gene_short next_gene_long next_gene_middle ; do
for exp in next_gene_middle; do
    ls coordinates_fasta/"$exp"_*.fasta |wc
    # for filename in coordinates_fasta/"$exp"_*.fasta; do
    for filename in coordinates_fasta/next_gene_middle_3067.fasta coordinates_fasta/next_gene_middle_1122.fasta; do
        echo "starting task $exp $filename.."
        (
            # .. do your stuff here
            # echo "starting task $filename.."
            spoa "$filename" > "consensus_fasta/$(basename "$filename" .fasta)_consensus.fasta"
        ) &

        # allow to execute up to $N jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
            # now there are $N jobs already running, so wait here for any job
            # to be finished so there is a place to start next one.
            wait -n
        fi

    done
done
# for filename in coordinates_fasta/*.fasta; do
#     echo $filename
#     spoa "$filename" > "consensus_fasta/$(basename "$filename" .fasta)_consensus.fasta" &
# done


# for filename in coordinates_fasta/*.fasta; do
#     echo $filename Do more stuff here
#     echo "spoa "$filename" > "consensus_fasta/$(basename "$filename" .fasta)_consensus.fasta""
#     sem -j+0 spoa "$filename" > "consensus_fasta/$(basename "$filename" .fasta)_consensus.fasta"
# done
# sem --wait
