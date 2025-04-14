#!/bin/bash

license="Apache-2.0"
author="Alexander Wietek <awietek@pks.mpg.de>"

cppfiles=`find ../xdiag -type f -name \*.cpp -not -path "../xdiag/extern/*"`
hppfiles=`find ../xdiag -type f -name \*.hpp -not -path "../xdiag/extern/*"` 
tomlfiles=`find .. -type f -name \*.toml -not -path "../xdiag/extern/*"` 
mdfiles=`find .. -type f -name \*.md -not -path "../xdiag/extern/*"` 
pyfiles=`find .. -type f -name \*.py -not -path "../xdiag/extern/*"` 
jlfiles=`find .. -type f -name \*.py -not -path "../xdiag/extern/*"` 

for f in ${hppfiles[@]}; do
    reuse annotate --copyright="Alexander Wietek <awietek@pks.mpg.de>" --license="Apache-2.0" $f
done

for f in ${cppfiles[@]}; do
    reuse annotate --copyright="Alexander Wietek <awietek@pks.mpg.de>" --license="Apache-2.0" $f
done

