cd $(dirname $1)
echo "$(tail -n +2 $1)" > "$1"
awk '{gsub(/} {/,"\n")}1' $1 > "$1_temp"
awk '{gsub(/}/,"")}1' "$1_temp" > $1
awk '{gsub(/{/,"")}1' $1 > "$1_temp"
awk '{gsub(/ /,",")}1' "$1_temp" > $1
rm "$1_temp"