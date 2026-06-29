#!/bin/bash

normal_headers=`grep -rl XDIAG_API ../xdiag  --exclude=all.hpp  --exclude-dir=old --exclude=*.cpp --exclude *distributed*`
distributed_headers=`grep -rl XDIAG_API ../xdiag  --exclude=all.hpp  --exclude-dir=old --exclude=*.cpp --include *distributed*`

# sort lexicographically
normal_headers=$(echo $normal_headers | xargs -n1 | sort | xargs)
distributed_headers=$(echo $distributed_headers | xargs -n1 | sort | xargs)

# add license on top
year=`date +%Y`
author="Alexander Wietek"

echo "/*
Copyright $year $author

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
"

echo "#pragma once"

echo
echo "#include <xdiag/armadillo.hpp>"


for header in ${normal_headers[@]}; do
    header=`echo $header | sed "s/..\///"`
    echo "#include <$header>"
done

echo

echo "#ifdef XDIAG_DISTRIBUTED"
for header in ${distributed_headers[@]}; do
    header=`echo $header | sed "s/..\///"`
    echo "#include <$header>"
done
echo "#endif"
echo

# undef the XDIAG macros (portable: BSD grep has no -P/\K)
macronames=`grep -rhE '^[[:space:]]*#[[:space:]]*define[[:space:]]+XDIAG_' ../xdiag | sed -E 's/^[[:space:]]*#[[:space:]]*define[[:space:]]+(XDIAG_[A-Za-z0-9_]+).*/\1/' | sort -u`
for macro in ${macronames[@]}; do
    if [[ ! "$macro" == "XDIAG_SHOW" ]]; then
	echo "#undef $macro"
    fi
done
