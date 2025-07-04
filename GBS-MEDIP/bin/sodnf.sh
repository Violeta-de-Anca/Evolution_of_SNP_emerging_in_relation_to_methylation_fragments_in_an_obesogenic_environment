cd /proj/naiss2023-23-55/GBS_violeta/trimmed
for i in $(ls *1.trimmed.fq.gz | cut -d "." -f 1); do
	input2=$i.2.trimmed.fq.gz
	echo "$input2"
done

