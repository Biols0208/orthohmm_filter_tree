for i in $(ls *_final_align*.aln); do python fa2phy.py ${i} ${i}.phy; done

# https://github.com/inab/trimal
for i in $(ls *_AA.aln.phy); do trimal -in ${i} -out ${i}.trimal1.phy -phylip_paml -gt 0.1 -st 0.0001 -cons 80 ; done
for i in $(ls *_AA.aln.phy); do trimal -in ${i} -out ${i}.trimal2.phy -phylip_paml -gt 0.1 -st 0.0001 ; done
for i in $(ls *_AA.aln.phy); do trimal -in ${i} -out ${i}.trimal3.phy -phylip_paml -automated1 ; done

for i in $(ls *_NT.aln.phy); do trimal -in ${i} -out ${i}.trimal4.phy -phylip_paml -gt 0.1 -st 0.0001 -cons 80 ; done
for i in $(ls *_NT.aln.phy); do trimal -in ${i} -out ${i}.trimal5.phy -phylip_paml -gt 0.1 -st 0.0001 ; done
for i in $(ls *_NT.aln.phy); do trimal -in ${i} -out ${i}.trimal6.phy -phylip_paml -automated1 ; done

for i in $(ls *_final_unmask_align_NT.aln); do sed -i "s/\!/-/g;s/\*/-/g" ${i}; done
for i in $(ls *_final_unmask_align_AA.aln); do sed -i "s/\!/-/g;s/\*/-/g" ${i}; done

for i in $(ls *_final_unmask_align*.aln); do python fa2phy.py ${i} ${i}.phy; done

for i in $(ls *_final_unmask_align_AA.aln.phy); do trimal -in ${i} -out ${i}.trimal7.phy -phylip_paml -gt 0.1 -st 0.0001 -cons 80 ; done
for i in $(ls *_final_unmask_align_AA.aln.phy); do trimal -in ${i} -out ${i}.trimal8.phy -phylip_paml -gt 0.1 -st 0.0001 ; done
for i in $(ls *_final_unmask_align_AA.aln.phy); do trimal -in ${i} -out ${i}.trimal9.phy -phylip_paml -automated1 ; done

for i in $(ls *_final_unmask_align_NT.aln.phy); do trimal -in ${i} -out ${i}.trimal10.phy -phylip_paml -gt 0.1 -st 0.0001 -cons 80 ; done
for i in $(ls *_final_unmask_align_NT.aln.phy); do trimal -in ${i} -out ${i}.trimal11.phy -phylip_paml -gt 0.1 -st 0.0001 ; done
for i in $(ls *_final_unmask_align_NT.aln.phy); do trimal -in ${i} -out ${i}.trimal12.phy -phylip_paml -automated1 ; done

# https://morgannprice.github.io/fasttree/
for i in $(ls *NT.aln.phy*phy); do FastTree -gamma -pseudo -nt ${i} > ${i}.gtr.tree ; done
for i in $(ls *NT.aln.phy*phy); do FastTree -gamma -pseudo -nt ${i} > ${i}.JCC.tree ; done
for i in $(ls *AA.aln.phy*phy); do FastTree -gamma -pseudo ${i} > ${i}.JTT.tree ; done
for i in $(ls *AA.aln.phy*phy); do FastTree -gamma -pseudo -wag ${i} > ${i}.wag.tree ; done
for i in $(ls *AA.aln.phy*phy); do FastTree -gamma -pseudo -lg ${i} > ${i}.lg.tree ; done

for i in $(ls *_NT.aln.phy); do python get_4fold_sites.py --output-format fasta --include-gaps ${i} ${i}.gaps.4d; done
for i in $(ls *_NT.aln.phy); do python get_4fold_sites.py --output-format fasta ${i} ${i}.nogaps.4d; done
for i in $(ls *.4d); do python fa2phy.py ${i} ${i}.phy; done

cat *tree > all.tree.nwk
tar -zcf all.tree.gz *.tree --remove-files
tar -zcf all.phy.gz *.phy *.4d *_stats.txt --remove-files --exclude='*_final_align_AA.aln.phy' --exclude='*_final_align_NT.aln.phy' --exclude='*_final_align_NT.aln.phy.gaps.4d.phy' --exclude='*_final_align_NT.aln.phy.nogaps.4d.phy' --exclude='*_final_align_NT.aln.phy.nogaps.4d' --exclude='*_final_unmask_align_NT.aln.phy.nogaps.4d'
