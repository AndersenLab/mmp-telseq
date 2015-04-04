cat `ls *100bp.telseq_elegans.TTAGGC.noreadgroup.txt` | grep -v 'Read' | cut -f 1,7 

for i in `ls *100bp.telseq_elegans.TTAGGC.noreadgroup.txt`; do
    awk -v i=${i} '{ print i "\t" $0}' $i | tr '.' '\t' | grep -v 'ReadGroup' >> telseq_stats.txt
done;

sbatch mmp_telseq.py 0
for i in `seq 1 300`; do sbatch mmp_telseq.py $i; done;
for i in `seq 301 601`; do sbatch mmp_telseq.py $i; done;


# START HERE
for i in `seq 602 2529`; do sbatch mmp_telseq.py $i --nodelist=$(rand_element 2 3 4 5 6); done;


rand() {
	printf $((  $1 *  RANDOM  / 32767   ))
}
rand_element () {
    local -a th=("$@")
    unset th[0]
    printf $'%s\n' "${th[$(($(rand "${#th[*]}")+1))]}"
}
