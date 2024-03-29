
#ifndef HYDRA_EXPORT_H
#define HYDRA_EXPORT_H

#ifdef HYDRA_STATIC_DEFINE
#  define HYDRA_EXPORT
#  define HYDRA_NO_EXPORT
#else
#  ifndef HYDRA_EXPORT
#    ifdef hydra_EXPORTS
        /* We are building this library */
#      define HYDRA_EXPORT 
#    else
        /* We are using this library */
#      define HYDRA_EXPORT 
#    endif
#  endif

#  ifndef HYDRA_NO_EXPORT
#    define HYDRA_NO_EXPORT 
#  endif
#endif

#ifndef HYDRA_DEPRECATED
#  define HYDRA_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef HYDRA_DEPRECATED_EXPORT
#  define HYDRA_DEPRECATED_EXPORT HYDRA_EXPORT HYDRA_DEPRECATED
#endif

#ifndef HYDRA_DEPRECATED_NO_EXPORT
#  define HYDRA_DEPRECATED_NO_EXPORT HYDRA_NO_EXPORT HYDRA_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef HYDRA_NO_DEPRECATED
#    define HYDRA_NO_DEPRECATED
#  endif
#endif

#endif /* HYDRA_EXPORT_H */
