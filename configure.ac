AC_INIT([Kepler], [1.0], [jasonfrowe@gmail.com])
AC_PROG_FC([gfortran ifort ifx])
AC_PROG_F77([gfortran ifort ifx])
AC_PROG_CC

LIBDIRS="$LIBDIRS \
         $HOME/Software/pgplot/build \
         /Users/jasonrowe/Software/pgplot/build \
         /opt/homebrew/lib \
         /opt/homebrew/opt/cfitsio/lib \
         /opt/X11/lib \
         /usr/lib \
         /usr/local/lib \
         /usr/local/pgplot \
         /usr/local/cfitsio \
         /opt/local/lib"

AC_MSG_NOTICE([Checking library directories...])
for LIBDIR in $LIBDIRS ; do
  if test -d "$LIBDIR"; then
    LDFLAGS="$LDFLAGS -L$LIBDIR"
  fi
done

INCDIRS="/opt/homebrew/include /opt/homebrew/opt/cfitsio/include /opt/X11/include /usr/local/include"
for INCDIR in $INCDIRS ; do
  if test -d "$INCDIR"; then
    FCFLAGS="$FCFLAGS -I$INCDIR"
    FFLAGS="$FFLAGS -I$INCDIR"
  fi
done

AC_MSG_NOTICE([  -- PGPLOT -- ])
AC_CHECK_LIB([pgplot],[pgopen_],[PGPLOTLIB="-lpgplot -lpng -lz"],[
  AC_MSG_WARN([PGPLOT library not found in standard search paths. Will attempt build with stub/custom path if available.])
  PGPLOTLIB=""
],[-lpng -lz])

AC_CHECK_LIB([png],[png_error],[PGPLOTLIB="$PGPLOTLIB -lpng -lz"],[],[$PGPLOTLIB])

AC_MSG_NOTICE([  -- CFITSIO -- ])
AC_CHECK_LIB([cfitsio],[ffopen],[CFITSIOLIB="-lcfitsio"],[
  AC_MSG_WARN([CFITSIO library not found. kfitsread target may require CFITSIO.])
  CFITSIOLIB=""
])

AC_PATH_X
if test -n "$ac_x_libraries"; then
  LFLAGS="$LDFLAGS -L$ac_x_libraries"
else
  LFLAGS="$LDFLAGS"
fi
XFLAGS="-lX11"

AC_SUBST(LFLAGS)
AC_SUBST(XFLAGS)
AC_SUBST(PGPLOTLIB)
AC_SUBST(CFITSIOLIB)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
