if PYTHONINTERFACE

AM_CPPFLAGS += $(PYTHON_INCLUDE)
BH_PATH_OBJ = ../src

pkglib_LTLIBRARIES = 

if ! PUBLIC
PythonInterface.d:	PythonInterface.i 
	swig $(AM_CPPFLAGS) -M -Wall -w321 -python -c++ -o PythonInterface.cpp $(srcdir)/PythonInterface.i 

@am__include@ ./PythonInterface.d

PythonInterface.cpp:	PythonInterface.d 
	swig $(AM_CPPFLAGS) -Wall -w321 -python -c++ -o PythonInterface.cpp.old $(srcdir)/PythonInterface.i 
	python $(srcdir)/fixSwig.py
	cp PythonInterface.cpp.new PythonInterface.cpp

PythonInterface.py:	PythonInterface.cpp
	@


pkglib_LTLIBRARIES += libPythonInterface.la

libPythonInterface_la_SOURCES =  PythonInterface.cpp
libPythonInterface_la_CXXFLAGS =  -fpermissive $(AM_CXXFLAGS)

endif

if BHLHPY

if HAS_SWIG
BH_LH_Interface.d:	BH_LH_Interface.i 
	swig $(AM_CPPFLAGS) -M -Wall -python -c++ -o BH_LH_Interface.cpp $(srcdir)/BH_LH_Interface.i 

@am__include@ ./BH_LH_Interface.d

BH_LH_Interface.cpp:	BH_LH_Interface.d 
	swig $(AM_CPPFLAGS) -Wall -w321 -python -c++ -o BH_LH_Interface.cpp $(srcdir)/BH_LH_Interface.i

BHLH.py:	BH_LH_Interface.cpp
	@
endif

EXTRA_DIST = BH_LH_Interface.cpp BHLH.py BH_LH_Interface.i


pkglib_LTLIBRARIES +=  libBH_LH_Interface.la
libBH_LH_Interface_la_SOURCES =  BH_LH_Interface.cpp
libBH_LH_Interface_la_CXXFLAGS =  -fpermissive $(AM_CXXFLAGS)

endif



LIBSTOLINK = \
	$(BH_PATH_OBJ)/Interface/.libs/libInterface.la \
	$(BH_PATH_OBJ)/.libs/libassembly.la \
	$(BH_PATH_OBJ)/.libs/libOLA.la \
	$(BH_PATH_OBJ)/.libs/libRatPart.la \
	$(BH_PATH_OBJ)/.libs/libCutPart.la \
	$(BH_PATH_OBJ)/ratext/.libs/libratext.la \
	$(BH_PATH_OBJ)/.libs/libBH.la \
	$(BH_PATH_OBJ)/trees_eval/.libs/libAmpl_eval.la \
	$(BH_PATH_OBJ)/rational_eval/.libs/libRateval.la \
	$(BH_PATH_OBJ)/cut/.libs/libCut_wCI.la \
	$(BH_PATH_OBJ)/cut_eval/.libs/libCuteval.la \
	$(BH_PATH_OBJ)/.libs/libBG.la       \
	$(BH_PATH_OBJ)/.libs/libBHcore.la    \
	$(BH_PATH_OBJ)/.libs/libIntegrals.la    \
	$(BH_PATH_OBJ)/.libs/libSpinors.la           
if BH_USE_GMP
LIBSTOLINK += $(BH_PATH_OBJ)/.libs/libmpreal.la    
endif
if ! PUBLIC
LIBSTOLINK += $(BH_PATH_OBJ)/../FeynDiagram/.libs/libfd.la

endif

LINKER_OPTS = \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/Interface/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/cut/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/cut_eval/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/ratext/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/rational_eval/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/../qd/src/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/../FeynDiagram/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/trees/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/trees_eval/.libs

LDFLAGS += $(LINKER_OPTS)




if ! PUBLIC
libPythonInterface_la_LIBADD = $(LIBSTOLINK) $(GMP_MPFR_LIBS) $(LIBS)
endif

if BHLHPY
libBH_LH_Interface_la_LIBADD = $(LIBSTOLINK) $(GMP_MPFR_LIBS) $(LIBS)
endif

install-exec-hook:
if ! PUBLIC
	if test -f ${prefix}/lib/blackhat/libPythonInterface.so; then ln -fs ${prefix}/lib/blackhat/libPythonInterface.so ${prefix}/lib/blackhat/_BH.so; else ln -fs ${prefix}/lib/blackhat/libPythonInterface.dylib ${prefix}/lib/blackhat/_BH.so; fi 
	cp $(builddir)/BH.py ${prefix}/lib/blackhat/
endif

if BHLHPY
	if test -f ${prefix}/lib/blackhat/libBH_LH_Interface.so; then ln -fs ${prefix}/lib/blackhat/libBH_LH_Interface.so ${prefix}/lib/blackhat/_BHLH.so; else ln -fs ${prefix}/lib/blackhat/libBH_LH_Interface.dylib ${prefix}/lib/blackhat/_BHLH.so; fi 
	cp $(builddir)/BHLH.py ${prefix}/lib/blackhat
endif

uninstall-hook:
if ! PUBLIC
	rm -rf ${prefix}/lib/blackhat/_BH.so 
	rm -rf ${prefix}/lib/blackhat/BH.py
endif

if BHLHPY
	rm -rf ${prefix}/lib/blackhat/_BHLH.so
	rm -rf ${prefix}/lib/blackhat/BHLH.py
endif

endif
