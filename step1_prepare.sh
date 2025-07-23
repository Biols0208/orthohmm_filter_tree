## calculate
python deal_orthohmm_gene_count.py orthohmm_gene_count.txt orthohmm_gene_count.txt.out1

## StdDev<1 && Species_Present_Percent>90
awk '$69<1 && $(NF-2)>90' orthohmm_gene_count.txt.out1 > orthohmm_gene_count.txt.out1.deal1
awk '{print $1}' orthohmm_gene_count.txt.out1.deal1 | sed 's/://g' > orthohmm_gene_count.txt.out1.deal1.group

## get cds sequence from all cds file
[ -d cds ] || mkdir cds

## all_species_in_one.cds.fa format
## >gene_name species
## AATTGGCCC

input=all_species_in_one.cds.fa
thread=48
cat orthohmm_gene_count.txt.out1.deal1.group | parallel --jobs ${thread} '
    grep -f <(echo {}) orthohmm_orthogroups.txt \
    | sed "s/ /\n/g" | sed "1d" \
    | seqtk subseq ${input} - > cds/{}.cds
'

#for i in $(cat orthohmm_gene_count.txt.out1.deal1.group); do echo ${i} | grep -f - orthohmm_orthogroups.txt | sed 's/ /\n/g' | sed '1d' | seqtk subseq all.deal.cds - > cds/${i}.cds; done
