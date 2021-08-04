files=$1
for file in `cat $files | cut -f1 -d '.'`
do
	echo $file	
	zgrep 6872 $file.gz | grep 6872 -w > taf1.$file
done