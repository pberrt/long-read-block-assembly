# /bin/bash

# for EXP in exp1 exp2 exp3 exp4 exp5 exp_spades exp_multik exp_spades_comp exp_multik_comp exp_artifacts
# # for EXP in exp_spades exp_multik
# do
#     # EXP=exp5
#     echo $EXP
#     python3 script_simple_multik_exps.py $EXP > ../res/simple_multik/$EXP.txt
#     # exit
# done


for EXP in truth real
do
    python3 script_ecoli.py $EXP --kmin 3 --kmax 23 --clipping > ../res/ecoli/$EXP_3_23_clipping.txt
    python3 script_ecoli.py $EXP --kmin 7 --kmax 7 --clipping > ../res/ecoli/$EXP_7_7_clipping.txt
    python3 script_ecoli.py $EXP --kmin 11 --kmax 11 --clipping > ../res/ecoli/$EXP_11_11_clipping.txt
    python3 script_ecoli.py $EXP --kmin 23 --kmax 23 --clipping > ../res/ecoli/$EXP_23_23_clipping.txt

    python3 script_ecoli.py $EXP --kmin 3 --kmax 23 > ../res/ecoli/$EXP_3_23.txt
    python3 script_ecoli.py $EXP --kmin 7 --kmax 7 > ../res/ecoli/$EXP_7_7.txt
    python3 script_ecoli.py $EXP --kmin 11 --kmax 11 > ../res/ecoli/$EXP_11_11.txt
    python3 script_ecoli.py $EXP --kmin 23 --kmax 23 > ../res/ecoli/$EXP_23_23.txt
    # exit
done