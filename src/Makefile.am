lib_LTLIBRARIES = libwlc.la
pkginclude_HEADERS = wlc.h
libwlc_la_SOURCES = wlc.c utils.c\
		    f_root.c f_root.h\
		    f_min.c f_min.h\
		    f_deriv.c f_deriv.h\
		    cavity_macros.h\
		    cavity.c cavity.h\
		    cavity_alloc.c cavity_alloc.h\
		    cavity_init.c cavity_init.h\
		    cavity_integrate.c cavity_integrate.h\
		    cavity_gradient_alloc.c cavity_gradient_alloc.h\
		    cavity_gradient.c cavity_gradient.h\
		    cavity_gradient_init.c cavity_gradient_init.h\
		    cavity_gradient_scalar.c cavity_gradient_scalar.h\
		    cavity_scalar.c cavity_scalar.h

libwlc_la_LIBADD = @GSL_LIBS@
