#!/bin/bash

for i in 4 5 6
do
    python snp2sim.py --mode varMDsim --protein dev1 --simLength 1 > wtRun.$i.log
    python snp2sim.py --mode varMDsim --protein dev1 --simLength 1 --varResID 68 --varAA L > varRun.$i.log 
done
