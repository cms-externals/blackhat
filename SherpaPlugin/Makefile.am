if SHERPAPLUGIN

AM_CPPFLAGS += @CONDITIONAL_SHERPAINCS@

BH_PATH_OBJ = ../src

pkglib_LTLIBRARIES = libBH_Sherpa_Interface.la

libBH_Sherpa_Interface_la_SOURCES = BH_Sherpa_Interface.cpp

LIBSTOLINK = \
	$(BH_PATH_OBJ)/Interface/.libs/libInterface.la \
	$(BH_PATH_OBJ)/.libs/libassembly.la \
	$(BH_PATH_OBJ)/.libs/libOLA.la \
	$(BH_PATH_OBJ)/.libs/libCutPart.la \
	$(BH_PATH_OBJ)/.libs/libRatPart.la \
	$(BH_PATH_OBJ)/ratext/.libs/libratext.la \
	$(BH_PATH_OBJ)/.libs/libBH.la \
	$(BH_PATH_OBJ)/trees_eval/.libs/libAmpl_eval.la \
	$(BH_PATH_OBJ)/rational_eval/.libs/libRateval.la \
	$(BH_PATH_OBJ)/cut_eval/.libs/libCuteval.la \
	$(BH_PATH_OBJ)/cut/.libs/libCut_wCI.la \
	$(BH_PATH_OBJ)/.libs/libBG.la       \
	$(BH_PATH_OBJ)/.libs/libBHcore.la    \
	$(BH_PATH_OBJ)/.libs/libIntegrals.la    \
	$(BH_PATH_OBJ)/.libs/libSpinors.la  

if ! PUBLIC
LIBSTOLINK += $(BH_PATH_OBJ)/../FeynDiagram/.libs/libfd.la    
endif

libBH_Sherpa_Interface_la_LIBADD = $(LIBSTOLINK) -lqd	

endif