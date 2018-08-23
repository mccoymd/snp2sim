#!/bin/bash

printf "1\n2\n3\n4\n5\n6\n7\n8\n9" > rank.txt

echo -e 'rank\tprotein\tligand\tvariant\tscaffold\taffinity\trmsd_lb\trmsd_ub' > vinaSummary.txt
for f in `ls */*pdbqt`
do
    grep RESULT $f > tempResults.txt
    name1=${f#*/}
    name2=${name1%.*}
    protein=${name2%%.*}
    name3=${name2#*.}
    variant=${name3%%.*}
    name4=${name3#*.}
    input=${name4%%.*}
    lig=${name4##*.}
    id="$protein $lig $variant $input"
#    echo $id
    sed -i -e "s/REMARK VINA RESULT:/\"$id\"/g" tempResults.txt
    sed -i -e "s/\"//g" tempResults.txt
    sed -i -e 's/ \+/\t/g' tempResults.txt
    paste rank.txt tempResults.txt >> vinaSummary.txt
done

rm rank.txt tempResults.txt
