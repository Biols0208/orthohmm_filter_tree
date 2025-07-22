##
mkdir out
for i in $(cat orthohmm_gene_count.txt.out1.deal1.group)
do 
echo "sh get_best_from_family.sh $PWD/cds/${i}.cds -o ${i} -d $PWD/out/${i} -t 8 --no-timestamp --force --log-level DEBUG " >> get_best.sh
done

##
parallel --jobs 48 bash -c {} < get_best.sh

# Or split get_best.sh and run
