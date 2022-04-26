#!/bin/bash
MASTER='thunder-master'

for i in include Makefile MACHINES a.GLOBAL b.FUNCTIONS c.SYSTEM f.MPI g.XC_FUNCTIONALS p.THEORY
do
    if [ -e $i ]
    then
	rm $i
    fi
    ln -s ../../${MASTER}/src/$i
done
