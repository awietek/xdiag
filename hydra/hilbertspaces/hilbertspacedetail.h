#ifndef HYDRA_HILBERTSPACES_DETAIL_
#define HYDRA_HILBERTSPACES_DETAIL_

namespace hydra { namespace hilbertspaces { namespace detail {

      template <int localdim> 
      inline constexpr int n_bits_for_localdim(int localdim);
      template <> inline constexpr int n_bits_for_localdim<2>() { return 1; }
      template <> inline constexpr int n_bits_for_localdim<3>() { return 2; }
      template <> inline constexpr int n_bits_for_localdim<4>() { return 2; }
      template <> inline constexpr int n_bits_for_localdim<5>() { return 3; }
      template <> inline constexpr int n_bits_for_localdim<6>() { return 3; }
      template <> inline constexpr int n_bits_for_localdim<7>() { return 3; }
      template <> inline constexpr int n_bits_for_localdim<8>() { return 3; }

      template<int localdim, class state_t>
      state_t compress_state(const state_t& state)
      {
	state_t state_compressed = 0;
	state_t workstate = state;
	state_t basis = 1;
	const state_t site_zero_mask = ((state_t)1 
					<< n_bits_for_localdim<localdim>()) - 1;
	while (workstate != 0)
	  {
	    state_compressed += (workstate & site_zero_mask)*basis; 
	    basis *= localdim;
	    workstate >>= n_bits_for_localdim<localdim>();
	  }
	return state_compressed;
      }

      template<int localdim, class state_t>
      state_t decompress_state(const state_t& state)
      {
	state_t state_decompressed = 0;
	state_t workstate = state;
	int shift = 0;
	while (workstate != 0)
	  {
	    state_decompressed |= ((state_t)(workstate%localdim)) 
	      << (n_bits_for_localdim<localdim>()*shift);
	    ++shift;
	    workstate /= localdim;
	  }
	return state_decompressed;
      }
      

    }
  }
}


#endif
