##
for i in $(cat orthohmm_gene_count.txt.out1.deal1.group)
do
echo "singularity run omm_macse.v10.02.sif --in_seq_file $PWD/out/${i}/${i}.best.fa --out_dir $PWD/out/${i}/omm_macse --out_file_prefix ${i} --save_details --java_mem 10000m --min_percent_NT_at_ends 0.7" >> all_OMM_MACSE.sh
done

##
parallel --jobs 48 bash -c {} < all_OMM_MACSE.sh

# Or split all_OMM_MACSE.sh and run
