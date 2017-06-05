#! /bin/sh -x
# a script to automate running FormCalc-generated code
# handy for testing/debugging
# this file is part of FormCalc
# last modified 21 Feb 17 th

# THIS IS A TEMPLATE AND NOT A FINAL PRODUCT.
# YOU MUST ADAPT THIS TO THE PROCESS YOU'RE LOOKING AT.

make || exit 1

arg=
polarizations=${1:-uuuu}
energy=${2:-00500}

rm -fr run$arg.$polarizations.$energy

CUBACORES=0 LTVERSION=43690 ./run $arg $polarizations $energy

less run$arg.$polarizations.$energy/0000001

