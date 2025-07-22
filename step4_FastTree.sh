##
for i in $(cat orthohmm_gene_count.txt.out1.deal1.group)
do
cp constructing_gene_tree.sh $PWD/out/${i}/omm_macse
done

##
for i in $(cat orthohmm_gene_count.txt.out1.deal1.group)
do
cd $PWD/out/${i}/omm_macse && sh constructing_gene_tree.sh && cd ../../../
done
