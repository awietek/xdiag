#!/bin/bash

normal_headers=`grep -rl XDIAG_API ../xdiag  --exclude=all.hpp  --exclude-dir=old --exclude=*.cpp --exclude *distributed*`
distributed_headers=`grep -rl XDIAG_API ../xdiag  --exclude=all.hpp  --exclude-dir=old --exclude=*.cpp --include *distributed*`

# sort lexicographically
normal_headers=$(echo $normal_headers | xargs -n1 | sort | xargs)
distributed_headers=$(echo $distributed_headers | xargs -n1 | sort | xargs)

echo "#pragma once"

echo
echo "#include <xdiag/armadillo.hpp>"


for header in ${normal_headers[@]}; do
    header=`echo $header | sed "s/..\///"`
    echo "#include <$header>"
done

echo

# echo "#ifdef XDIAG_USE_MPI"
# for header in ${distributed_headers[@]}; do
#     header=`echo $header | sed "s/..\///"`
#     echo "#include <$header>"
# done
# echo "#endif"

# undef the XDIAG macros
macronames=`grep -rhoP '#define\s+\KXDIAG_[A-Za-z0-9_]+' ../xdiag | sort -u`
for macro in ${macronames[@]}; do
    if [[ ! "$macro" == "XDIAG_SHOW" ]]; then
	echo "#undef $macro"
    fi
done
