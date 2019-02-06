#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/utils/typedefs.h>

#include "siteoperations.h"

namespace hydra { namespace hilbertspaces {
    
    template <class state_t>
    int get_site_val_single(const state_t& state, const int& site)
    { return (state >> site) & 1; }
    
    template <class state_t>
    void set_site_val_single(state_t* state, const int& site, const int& val)
    { *state = (*state & ~((state_t)1 << site)) | ((state_t)val << site);}


    // Template instantiations
    template 
    int get_site_val_single<uint16>(const uint16& state, const int& site);
    template 
    int get_site_val_single<uint32>(const uint32& state, const int& site);
    template 
    int get_site_val_single<uint64>(const uint64& state, const int& site);

    template 
    void set_site_val_single<uint16>(uint16* state, const int& site, const int& val);
    template 
    void set_site_val_single<uint32>(uint32* state, const int& site, const int& val);
    template 
    void set_site_val_single<uint64>(uint64* state, const int& site, const int& val);


    // Template specializations
    template <>
    int get_site_val<Spinhalf<uint16>>(const uint16& state, const int& site)
    { return get_site_val_single<uint16>(state, site); }
    template <>
    int get_site_val<Spinhalf<uint32>>(const uint32& state, const int& site)
    { return get_site_val_single<uint32>(state, site); }
    template <>
    int get_site_val<Spinhalf<uint64>>(const uint64& state, const int& site)
    { return get_site_val_single<uint64>(state, site); }
    template <>
    void set_site_val<Spinhalf<uint16>>(uint16* state, const int& site, const int& val)
    { return set_site_val_single<uint16>(state, site, val); }
    template <>
    void set_site_val<Spinhalf<uint32>>(uint32* state, const int& site, const int& val)
    { return set_site_val_single<uint32>(state, site, val); }
    template <>
    void set_site_val<Spinhalf<uint64>>(uint64* state, const int& site, const int& val)
    { return set_site_val_single<uint64>(state, site, val); }

    
  }  // namespace hilbertspaces
}  // namespace hydra
