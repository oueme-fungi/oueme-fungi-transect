#!/usr/bin/env bash
# take calls from mlocarna script and make a snakemake file, then deploy to the
# cluster for parallel computation at the last minute

# fake writing the output files so locarna will finish
echo $* | egrep -o '[^=]+intermediate[0-9]+\.[a-z]+' | xargs touch

out=$(echo $* | sed -r 's|.*--pp=(intermediates/intermediate[0-9]+\.pp).*|\1|')
echo "out=$out"
i=$(echo $out | egrep -o "[0-9]+")
echo "i=$i"
(( i == 1 )) && [ -e makelocarna.make ] && rm makelocarna.make

echo "$out: $1 $2" >>makelocarna.make
echo "	locarna $1 $2 ${@:4}" >>makelocarna.make

n=$(grep -c "^[^#].* 60$" input/input.fa)
echo "n=$n"
(( i == n - 1 )) && echo "executing..."
(( i == n - 1 )) && rm intermediates/intermediate*.pp && make -f makelocarna.make -j $3 $out

exit 0
