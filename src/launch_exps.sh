# /bin/bash

for EXP in exp1 exp2 exp3 exp4 exp5 exp_spades exp_multik
do
    # EXP=exp5
    echo $EXP
    python3 script_graph_multi_k.py $EXP > res/$EXP.txt
    # exit
done
