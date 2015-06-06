#!/bin/bash

#
# Uncomment this block to have grids with differing numbers of points in x and y.
#
# if  [[ $# -ne 6 ]]; then
#     echo "need arguments: N_Q M_Q N_QPRIME M_QPRIME QGRID_FILENAME QPRIMEGRID_FILENAME";
#     exit;
# fi
# 
# nq="$1"
# mq="$2"
# nqprime="$3"
# mqprime="$4"
# qgrid_filename="$5"
# qprimegrid_filename="$6"

if  [[ $# -ne 4 ]]; then
    echo "need arguments: N_Q N_QPRIME QGRID_FILENAME QPRIMEGRID_FILENAME";
    exit;
fi

nq="$1"
mq=$nq
nqprime="$2"
mqprime=$nqprime
qgrid_filename="$3"
qprimegrid_filename="$4"

calculate_scatterings() {
    # $1 is the file containing the q points to process
    # $2 is the file containing the q-prime points
    q_grid="$1"
    qprime_grid="$2"
    if hash parallel 2>/dev/null; then
        time parallel --progress --colsep=' ' ./scattering "$qprime_grid" {1} {2}  :::: "$q_grid"
    else
        time while IFS='' read -r line
        do
            IFS=' ' read -a qp <<< "$line";
            qx="${qp[0]}"
            qy="${qp[1]}"
            ./scattering "$qprime_grid" "$qx" "$qy";
        done < "$q_grid";
    fi
}

# generate both grids across entire BZ
./gengrid $qprimegrid_filename $nqprime $mqprime;
./gengrid $qgrid_filename $nq $mq;

# filter q points to IBZ
python FilterToIBZ.py $qgrid_filename > ibz_$qgrid_filename;

numqppts=`wc -l < $qprimegrid_filename`;
echo $numqppts;
echo "qprime grid has $numqppts points.";

numqpts=`wc -l < ibz_$qgrid_filename`;
echo $numqpts;
echo "filtered q grid has $numqpts points.";

mkdir output;
calculate_scatterings ibz_$qgrid_filename $qprimegrid_filename;
outdir=$qgrid_filename.output;
while [ -d "$outdir" ]
do
    outdir=$outdir.new;
done
mv output $outdir;

conductivity=`./integrate $outdir/*`;
echo "Calculated thermal conductivity is $conductivity.";

echo "qprime grid has $numqppts points." >> $outdir/info;
echo "filtered q grid has $numqpts points." >> $outdir/info;
echo "Calculated thermal conductivity is $conductivity." >> $outdir/info;
