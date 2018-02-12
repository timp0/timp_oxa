#!/bin/bash

smaugdir=/atium/Data/Nanopore/Analysis/170104_reassemble

if [ $1 == gfa2dp ]; then
    for i in 1 2 4 6 7 8 9 10 12 ;
    do
	tiggfa=$smaugdir/170104_reassemble/canu/$i.nanopore/$i.nanopore.contigs.gfa
	cp $tiggfa ~/Dropbox/Lab/oxa/canu_gfa
    done
fi


if [ $1 == circ2smaug ]; then
    for i in 1 2 4 6 7 8 9 10 12 ;
    do
	piloncirc=~/work/oxa/$i.circlate/$i.pilon_circ
	canucirc=~/work/oxa/$i.circlate/$i.canu_circ
	scp -r $piloncirc smaug:$smaugdir/circlate/
	scp -r $canucirc smaug:$smaugdir/circlate/
    done
fi
