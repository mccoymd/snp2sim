#!/bin/bash

#Argument: path to directory with vina pdbqt results
#Argument: path to output directory

cd $2

printf "1\n2\n3\n4\n5\n6\n7\n8\n9" > rank.txt

echo -e 'rank\tprotein\tligand\tlibrary\tvariant\tscaffold\taffinity\trmsd_lb\trmsd_ub' > vinaSummary.txt
for f in $1/*.pdbqt
do
    grep RESULT $f > tempResults.txt

    f=$(basename $f)
    protein=$(echo $f | cut -d "." -f1)
    variant=$(echo $f | cut -d "." -f2)
    input=$(echo $f | cut -d "." -f3)
    library=$(echo $f | cut -d "." -f4)
    lig=$(echo $f | cut -d "." -f5)

    id="$protein $lig $library $variant $input"
    sed -i -e "s/REMARK VINA RESULT:/\"$id\"/g" tempResults.txt
    sed -i -e "s/\"//g" tempResults.txt
    sed -i -e 's/ \+/\t/g' tempResults.txt
    paste rank.txt tempResults.txt >> vinaSummary.txt
done

mv vinaSummary.txt "$2/vinaSummary_$protein.txt"
rm rank.txt tempResults.txt*
