grep -i "fail\|error" *.e >> ~/fox_repo/data/sorterrors.txt
sed -n 's/^\(.*\)bwa\.sorted\.e.*/\1bwa.bam/p' sorterrors.txt | sort | uniq >> rerun_sort.txt

#anders' suggestions for taking a look at errors
#tail *.e | less -S
#wc -l *e 
#grep -i "error" / grep -i "fail" *.e
