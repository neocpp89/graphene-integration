#!/bin/bash

QLIST=$(seq 25 25 400)
QPLIST=$(seq 25 25 400)

prefix=B
gridfile=${prefix}.scan
echo -ne > $gridfile
for q in $QLIST; do
    for qp in $QPLIST; do
        qfile=${prefix}q${q}_qp${qp} 
        qpfile=${prefix}qp${qp}
        outputdir=${qfile}.output
        if [ ! -d $outputdir ]; then
            echo $qfile $qpfile
            \time -o timing ./run_integration.sh $q $qp $qfile $qpfile
            mv timing ${outputdir}/
        fi
        k=$(tail -n1 ${outputdir}/info | cut -d' ' -f 5 | rev | cut -c 2- | rev)
        echo -ne "$k " >> $gridfile
    done
    echo >> $gridfile
done
