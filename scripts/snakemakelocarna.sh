#!/usr/bin/env bash
# take calls from mlocarna script and make a snakemake file, then deploy to the
# cluster for parallel computation at the last minute

# fake writing the output files so locarna will continue
echo $* | egrep -o '[^=]+intermediate[0-9]+\.[a-z]+' | xargs touch

# find the name and number of the pp file we've been asked to write
out=$(echo $* | sed -r 's|.*--pp=(intermediates/intermediate[0-9]+\.pp).*|\1|')
i=$(echo $out | egrep -o "[0-9]+")

# for the first file, write prereqs table header
(( i == 1 )) && echo "i,in1,in2" >smlocarna.csv

# find the two prerequisite files
in1=$(echo $1 | sed 's|intermediates/|snakeintermediates/|')
shift
in2=$(echo $1 | sed 's|intermediates/|snakeintermediates/|')

# write prereqs to the table
echo "$i,$in1,$in2" >>smlocarna.csv

# find the total number of sequences in the input file
n=$(grep -c "^[^#].* 60$" input/input.fa)
# the total number of nodes is one less than the number of sequences
# on the last one, do all the work
# assume that all parameters except input and output files are the same in each
# run
if (( i == n - 1 )); then
	loc_command="locarna {input.in1} {input.in2}"
	sm_command="snakemake"
	output="output:"
	# process the command line
	while shift; do
		case $1 in
			--conda) conda="conda: \"$2\"
    "; sm_command="$sm_command --use-conda "; shift ;;
			--profile) sm_command="$sm_command --profile $2 "; shift;;
			--clustal=*)
				loc_command="$loc_command --clustal={output.clustal}"
				output="$output
        clustal=\"snakeintermediates/intermediate{i,\\d+}.aln\"," ;;
			--stockholm=*)
				loc_command="$loc_command --stockholm={output.stockholm}"
				output="$output
        stockholm=\"snakeintermediates/intermediate{i,\\d+}.stk\",";;
			--pp=*)
				loc_command="$loc_command --pp={output.pp}"
				output="$output
        pp=\"snakeintermediates/intermediate{i,\\d+}.pp\",";;
			*) loc_command="$loc_command $1"
		esac
	done
	# remove trailing comma
	output=$(echo $output | sed 's/,$//')
	echo "writing Snakefile"
	cat <<-XXX >Snakefile
		import pandas as pd
		queue = pd.read_csv("smlocarna.csv").set_index("i")
		def getinputs(wildcards):
		    return { 'in1': queue.in1[int(wildcards.i)], 'in2': queue.in2[int(wildcards.i)] }
		rule locarna:
		    $output
		    input: unpack(getinputs)
		    ${conda}shell: 
		        "$loc_command"
		XXX
	echo "executing..."
	rm intermediates/intermediate*.pp
	rm intermediates/intermediate*.stk
	rm intermediates/intermediate*.aln
        echo "#!/bin/env sh" >snakemake.sh
	echo "$sm_command -r snakeintermediates/intermediate$i.pp" >>snakemake.sh
        chmod +x snakemake.sh
        snakemake -t snakeintermediate$i.pp
	for f in $(ls snakeintermediates)
	  do
	  ln -s "../snakeintermediates/$f" "intermediates/$f"
	done
fi
