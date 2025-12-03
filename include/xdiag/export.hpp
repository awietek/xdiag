
#ifndef XDIAG_EXPORT_H
#define XDIAG_EXPORT_H

#ifdef XDIAG_STATIC_DEFINE
#  define XDIAG_EXPORT
#  define XDIAG_NO_EXPORT
#else
#  ifndef XDIAG_EXPORT
#    ifdef xdiag_EXPORTS
        /* We are building this library */
#      define XDIAG_EXPORT 
#    else
        /* We are using this library */
#      define XDIAG_EXPORT 
#    endif
#  endif

#  ifndef XDIAG_NO_EXPORT
#    define XDIAG_NO_EXPORT 
#  endif
#endif

#ifndef XDIAG_DEPRECATED
#  define XDIAG_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef XDIAG_DEPRECATED_EXPORT
#  define XDIAG_DEPRECATED_EXPORT XDIAG_EXPORT XDIAG_DEPRECATED
#endif

#ifndef XDIAG_DEPRECATED_NO_EXPORT
#  define XDIAG_DEPRECATED_NO_EXPORT XDIAG_NO_EXPORT XDIAG_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef XDIAG_NO_DEPRECATED
#    define XDIAG_NO_DEPRECATED
#  endif
#endif

#endif /* XDIAG_EXPORT_H */
