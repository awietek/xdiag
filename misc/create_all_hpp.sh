#!/bin/bash

normal_headers=`grep -rl XDIAG_API ../xdiag  --exclude=../xdiag/common.hpp --exclude=../xdiag/all.hpp  --exclude=*.cpp --exclude *distributed*`
distributed_headers=`grep -rl XDIAG_API ../xdiag  --exclude=../xdiag/common.hpp --exclude=../xdiag/all.hpp --exclude=*.cpp --include *distributed*`

# sort lexicographically
normal_headers=$(echo $normal_headers | xargs -n1 | sort | xargs)
distributed_headers=$(echo $distributed_headers | xargs -n1 | sort | xargs)

echo "#pragma once"

echo
echo "#include <xdiag/common.hpp>"
echo

for header in ${normal_headers[@]}; do
    header=`echo $header | sed "s/..\///"`
    echo "#include <$header>"
done

echo



echo "#ifdef XDIAG_USE_MPI"
for header in ${distributed_headers[@]}; do
    header=`echo $header | sed "s/..\///"`
    echo "#include <$header>"
done
echo "#endif"

echo

echo "#undef XDIAG_THROW"
echo "#undef XDIAG_RETHROW"
echo "#undef XDIAG_API"
echo "#undef XDIAG_OFFSET"
