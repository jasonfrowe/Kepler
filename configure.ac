AC_INIT([jwst], [0.1], [jasonfrowe@gmail.com])
AC_PROG_F77([ifort gfortran])
AC_PROG_FC([ifort gfortran])
LIBDIRS="$LIBDIRS           \
         /usr/lib           \
         /usr/local/lib     \
         /usr/local/pgplot  \
         /usr/local/cfitsio \
         /local/lib         \
         /local/pgplot      \
         /local/cfitsio     \
         /opt/local/lib"
AC_MSG_NOTICE(LIBDIRS)
for LIBDIR in $LIBDIRS ; do
  AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
done
AC_MSG_NOTICE([  -- PGPLOT -- ])
AC_CHECK_LIB([pgplot],[pgopen_],[PGPLOTLIB="-lpgplot"],[AC_MSG_ERROR([PGPLOT NOT FOUND])],)
AC_CHECK_LIB([png],     [png_error], [PGPLOTLIB="-lpng $PGPLOTLIB"], 
  	     [], [$PGPLOTLIB])
#AC_MSG_NOTICE([$PGPLOTLIB])
#AC_MSG_NOTICE([$LDFLAGS])
AC_MSG_NOTICE([  -- CFITSIO -- ])
AC_CHECK_LIB([cfitsio], [ffopen],[CFITSIOLIB="-lcfitsio"],[AC_MSG_ERROR([CFITSIO NOT FOUND])],)
AC_PATH_X
LFLAGS="$LDFLAGS -L$ac_x_libraries"
XFLAGS="-lX11"
AC_MSG_NOTICE([$LFLAGS])
AC_SUBST(LFLAGS)
AC_SUBST(XFLAGS)
AC_SUBST(PGPLOTLIB)
AC_SUBST(CFITSIOLIB)
AC_CONFIG_FILES([Makefile]) ## Create Makefile
AC_OUTPUT
