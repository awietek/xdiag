// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/representation.hpp>

// Apply H*v via matrix multiply and via apply(); check H^2*v and e0 match
void test_apply(int64_t N, xdiag::OpSum ops);

// Same as test_apply but with a multi-column matrix instead of a single vector
void test_apply_mat(int64_t N, xdiag::OpSum ops);

// Verify that onsite op12 equals op1*op1 for all sites and nup sectors
void test_onsite(std::string op1, std::string op12);

// Check the lowest eigenvalue per irrep against a reference for the
// Kitaev-Gamma model on the N=8 honeycomb lattice
void test_kitaev_gamma(double K, double G,
                       std::vector<std::pair<std::string, double>> irrep_names_e0);

// Check apply agrees with matrix multiply for each (nup, irrep) sector
void test_spinhalf_symmetric_apply(xdiag::OpSum ops, int64_t nsites,
                                   std::vector<xdiag::Representation> const &irreps);

// Same check without fixing nup (full Hilbert space with symmetry)
void test_spinhalf_symmetric_apply_no_sz(
    xdiag::OpSum ops, int64_t nsites,
    std::vector<xdiag::Representation> const &irreps);

// Verify that symmetry-resolved spectra reconstruct the full spectrum
void test_spinhalf_symmetric_spectra(xdiag::OpSum ops, int64_t nsites,
                                     std::vector<xdiag::Representation> irreps,
                                     std::vector<int64_t> multiplicities);

// Same verification without fixing nup
void test_spinhalf_symmetric_spectra_no_sz(
    xdiag::OpSum ops, int64_t nsites,
    std::vector<xdiag::Representation> irreps,
    std::vector<int64_t> multiplicities);
