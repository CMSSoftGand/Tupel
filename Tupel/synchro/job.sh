#!/bin/sh
#ntuplesynchro.root
#Lyon/Andrey_v1.root
eval `scram runtime -sh` 
root -l -b <<EOF 
.x simpleReader.cc+ ("ntuplesynchro_v4.root","outputntuple") 
.q; 
EOF
