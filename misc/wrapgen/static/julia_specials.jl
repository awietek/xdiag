# SPDX-License-Identifier: Apache-2.0
# STATIC hand-written special -- copied verbatim by generate.sh.

# make_complex modifies State -> put a ! at the end of the name in Julia
make_complex!(obj::State) = mth_make_complex(obj.cxx_object)

# Utility to have the ^-operator in Julia
Base.:^(p::Permutation, power::Int64)::Permutation = Permutation(fun_pow(p.cxx_object, power))
