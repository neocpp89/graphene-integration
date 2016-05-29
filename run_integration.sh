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
    prefix="$3"
    if hash parallel 2>/dev/null; then
        time parallel --slf machines.slf --eta --progress --colsep=' ' "cd $(pwd); ./scattering \"$qprime_grid\" {1} {2} \"$prefix\" "  :::: "$q_grid"
    else
        time while IFS='' read -r line
        do
            IFS=' ' read -a qp <<< "$line";
            qx="${qp[0]}"
            qy="${qp[1]}"
            ./scattering "$qprime_grid" "$qx" "$qy" "$prefix";
        done < "$q_grid";
    fi
}

outdir=$qgrid_filename.output;
while [ -d "$outdir" ]
do
    outdir=$outdir.new;
done
prefix=$outdir/tauinv_q

mkdir -p $outdir

# generate both grids across entire BZ
./gengrid $outdir/$qprimegrid_filename $nqprime $mqprime;
./gengrid $outdir/$qgrid_filename $nq $mq;

# filter q points to IBZ
python FilterToFiniteIBZ.py $outdir/$qgrid_filename > $outdir/ibz_$qgrid_filename;

numqppts=`wc -l < $outdir/$qprimegrid_filename`;
echo $numqppts;
echo "qprime grid has $numqppts points.";

numqpts=`wc -l < $outdir/ibz_$qgrid_filename`;
echo $numqpts;
echo "filtered q grid has $numqpts points.";

calculate_scatterings $outdir/ibz_$qgrid_filename $outdir/$qprimegrid_filename $prefix;

conductivity=`./integrate ${prefix}*`;
echo "Calculated thermal conductivity is $conductivity.";

echo "qprime grid has $numqppts points." >> $outdir/info;
echo "filtered q grid has $numqpts points." >> $outdir/info;
echo "Calculated thermal conductivity is $conductivity." >> $outdir/info;
