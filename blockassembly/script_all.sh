# /bin/bash


# for file in GCA_027944575.1_ASM2794457v1_genomic GCA_027944595.1_ASM2794459v1_genomic GCA_027944615.1_ASM2794461v1_genomic GCA_027944635.1_ASM2794463v1_genomic GCA_027944655.1_ASM2794465v1_genomic GCA_027944675.1_ASM2794467v1_genomic GCA_027944695.1_ASM2794469v1_genomic GCA_027944715.1_ASM2794471v1_genomic GCA_027944735.1_ASM2794473v1_genomic GCA_027944775.1_ASM2794477v1_genomic GCA_027944795.1_ASM2794479v1_genomic GCA_027944815.1_ASM2794481v1_genomic GCA_027944835.1_ASM2794483v1_genomic GCA_027944855.1_ASM2794485v1_genomic GCA_027944875.1_ASM2794487v1_genomic GCA_027944895.1_ASM2794489v1_genomic GCA_027944915.1_ASM2794491v1_genomic GCA_027944935.1_ASM2794493v1_genomic GCA_027944955.1_ASM2794495v1_genomic GCA_027944995.1_ASM2794499v1_genomic GCA_027945015.1_ASM2794501v1_genomic GCA_027945035.1_ASM2794503v1_genomic GCA_027945055.1_ASM2794505v1_genomic GCA_028551585.1_ASM2855158v1_genomic
# LIGHT FILES
# for file in GCA_027944575.1_ASM2794457v1_genomic GCA_027944615.1_ASM2794461v1_genomic GCA_027944695.1_ASM2794469v1_genomic GCA_027944715.1_ASM2794471v1_genomic GCA_027944735.1_ASM2794473v1_genomic GCA_027944775.1_ASM2794477v1_genomic GCA_027944795.1_ASM2794479v1_genomic GCA_027944815.1_ASM2794481v1_genomic GCA_027944855.1_ASM2794485v1_genomic GCA_027944875.1_ASM2794487v1_genomic GCA_027944895.1_ASM2794489v1_genomic GCA_027944995.1_ASM2794499v1_genomic GCA_027945015.1_ASM2794501v1_genomic GCA_027945035.1_ASM2794503v1_genomic GCA_028551585.1_ASM2855158v1_genomic
# HUGE FILES
# for file in GCA_027944595.1_ASM2794459v1_genomic GCA_027944635.1_ASM2794463v1_genomic GCA_027944655.1_ASM2794465v1_genomic GCA_027944675.1_ASM2794467v1_genomic GCA_027944835.1_ASM2794483v1_genomic GCA_027944915.1_ASM2794491v1_genomic GCA_027944935.1_ASM2794493v1_genomic GCA_027944955.1_ASM2794495v1_genomic GCA_027945055.1_ASM2794505v1_genomic 
# for file in GCA_027944575.1_ASM2794457v1_genomic
for file in GCA_027944595.1_ASM2794459v1_genomic GCA_027944615.1_ASM2794461v1_genomic GCA_027944635.1_ASM2794463v1_genomic GCA_027944655.1_ASM2794465v1_genomic GCA_027944675.1_ASM2794467v1_genomic GCA_027944695.1_ASM2794469v1_genomic GCA_027944715.1_ASM2794471v1_genomic GCA_027944735.1_ASM2794473v1_genomic GCA_027944775.1_ASM2794477v1_genomic GCA_027944795.1_ASM2794479v1_genomic GCA_027944815.1_ASM2794481v1_genomic GCA_027944835.1_ASM2794483v1_genomic GCA_027944855.1_ASM2794485v1_genomic GCA_027944875.1_ASM2794487v1_genomic GCA_027944895.1_ASM2794489v1_genomic GCA_027944915.1_ASM2794491v1_genomic GCA_027944935.1_ASM2794493v1_genomic GCA_027944955.1_ASM2794495v1_genomic GCA_027944995.1_ASM2794499v1_genomic GCA_027945015.1_ASM2794501v1_genomic GCA_027945035.1_ASM2794503v1_genomic GCA_027945055.1_ASM2794505v1_genomic GCA_028551585.1_ASM2855158v1_genomic
# for file in GCA_027944815.1_ASM2794481v1_genomic #Error topo_order GCA_027945035.1_ASM2794503v1_genomic
# for file in GCA_027945035.1_ASM2794503v1_genomic
# for file in GCA_027944875.1_ASM2794487v1_genomic
do
    # python -u script_local_msa.py --sample "$file" --kmin 2
    python -u script_local_msa.py --sample "$file" > ../res/ecoli_bcalm/TEST_nac_data/"$file"_log.out 2>&1 
    # python -u script_local_msa.py --abundance_cleaning --sample "$file" > ../res/ecoli_bcalm/TEST_ac_data/"$file"_log.out 2>&1 
done