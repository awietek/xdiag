#pragma once
namespace hydra {

template <typename coeff_t, typename block_t, class fill_f>
void generic_fill(block_t const &block, coeff_t *ptr, fill_f fill) try {

} catch (...) {
  rethrow(__func__, "cannot perform generic fill");
}

} // namespace hydra
