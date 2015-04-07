cat `ls *100bp.telseq_elegans.TTAGGC.noreadgroup.txt` | grep -v 'Read' | cut -f 1,7 

rm telseq_mmp.txt
paste <(echo "strain_length" | tr '_' '\t') <(head -n 1 VC30134.76bp.telseq_elegans.TTAGGC.noreadgroup.txt) >> telseq_mmp.txt
for i in `ls *.telseq_elegans.TTAGGC.noreadgroup.txt`; do
    in=`echo ${i} | cut -f 1,2 -d '.'`
    awk -v i=${in} '{ print i "\t" $0}' $i | tr '.' '\t' | grep -v 'ReadGroup' >> telseq_mmp.txt
done;



# START HERE
for i in `seq 1 2600`; do 
    line=$(sed -n "${i}p" "strain_info.txt")
    strain=`echo $line | cut -f 3 -d ' '`
    size=`echo $line | cut -f 5 -d ' '`
    if [ ! -f telseq/${strain}.${size}.telseq_elegans.TTAGGC.noreadgroup.txt ]; then
        echo $i
        sbatch --nodelist=node$(rand_element 2 3 4 5) mmp_telseq.py $i;
    fi
done;


rand() {
    printf $((  $1 *  RANDOM  / 32767   ))
}
rand_element () {
    local -a th=("$@")
    unset th[0]
    printf $'%s\n' "${th[$(($(rand "${#th[*]}")+1))]}"
}
