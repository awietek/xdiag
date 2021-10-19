#!/bin/bash

liladir="/Users/awietek/Research/Software/lila"
limedir="/Users/awietek/Research/Software/lime"
claradir="/Users/awietek/Research/Software/Clara"
hydradir="/Users/awietek/Research/Software/hydra"
includes="-I$hydradir -I$liladir -I$limedir -I$claradir/include -I/usr/local/Cellar/open-mpi/4.0.5/include/"

# sources=(`find $hydradir/hydra -name *.h -not -path "*old*"` `find $hydradir/hydra -name *.cpp -not -path "*old*"`)

sources=(`find $hydradir/hydra -name *.h -not -path "*old*"`)


target_dir=hyde
logfile=creade_hyde_markdown.log

rm -f $logfile

for src in ${sources[@]}; do
    echo Compiling $src
    hyde -hyde-yaml-dir=$target_dir -hyde-update --access-filter-public $src -- -x c++ -std=c++17 -c -DLILA_USE_ACCELERATE -Wno-pragma-once-outside-header $includes >> $logfile 2>&1

#    2>&1 | grep -v 'extraneous file' | grep -v '.md@' 
done
