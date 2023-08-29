#pragma once

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

namespace hydra {

double norm(State const &v);

// dot
double dot(State const &v, State const &w);
complex dotC(State const &v, State const &w);

double inner(BondList const &bonds, State const &v);
complex innerC(BondList const &bonds, State const &v);
double inner(Bond const &bonds, State const &v);
complex innerC(Bond const &bonds, State const &v);

State &operator*=(State &X, complex alpha);
State &operator*=(State &X, double alpha);
State &operator/=(State &X, complex alpha);
State &operator/=(State &X, double alpha);

// // inner
// double inner(Bond const &bond, State const &v);
// double inner(State const &w, BondList const &bonds, State const &v);
// double inner(State const &w, Bond const &bond, State const &v);
// complex innerC(BondList const &bonds, State const &v);
// complex innerC(Bond const &bond, State const &v);
// complex innerC(State const &w, BondList const &bonds, State const &v);
// complex innerC(State const &w, Bond const &bond, State const &v);

// State &operator*=(State &X, complex alpha);
// State &operator*=(State &X, double alpha);

} // namespace hydra
