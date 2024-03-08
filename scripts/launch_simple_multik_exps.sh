# /bin/bash

# for EXP in exp1 exp2 exp3 exp4 exp5 exp_spades exp_multik exp_spades_comp exp_multik_comp exp_artifacts
# # for EXP in exp_spades exp_multik
# do
#     # EXP=exp5
#     echo $EXP
#     python3 script_simple_multik_exps.py $EXP > ../res/simple_multik/$EXP.txt
#     # exit
# done


for EXP in truth simulated
do
    python3 script_ecoli.py $EXP --kmin 3 --kmax 11 > ../res/ecoli/$EXP.txt
    # exit
done