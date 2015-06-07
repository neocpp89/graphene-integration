# graphene-integration
Conversion of Sandeep's integration code for thermal conductivity of graphene.

The basic steps to compile are the usual ones:

    make clean
    make

This will create a set of executables. For ease of use, I have included the
script `./run_integration.sh` which tries to take care of everything. You
just need to provide it with Q point density, Q prime point density, file to
save Q points, and file to save Q prime points (in order).

After running, the results will be saved in the folder
`"Q Point filename".output`, and the summary will be given in the
info file in this directory.
