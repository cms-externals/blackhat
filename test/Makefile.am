
datadir= data


BH_PATH_OBJ = $(abs_top_builddir)/src/

PRIMITIVE_TESTS=parents_2q2l.dat parents_1g2q2l.dat parents_2g2q2l.dat parents_4q2l.dat parents_3g2q2l.dat parents_1g4q2l.dat parents_4g2q2l.dat parents_2g4q2l.dat parents_6q2l.dat parents_1g2q1y.dat parents_2g2q1y.dat parents_4q1y.dat parents_3g2q1y.dat parents_1g4q1y.dat parents_4g2q1y.dat parents_2g4q1y.dat parents_6q1y.dat parents_4q.dat parents_2g2q.dat parents_4g.dat parents_1g4q.dat parents_3g2q.dat parents_5g.dat parents_6q.dat parents_2g4q.dat parents_4g2q.dat parents_6g.dat 

LONG_PRIMITIVE_TESTS=parents_1g6q2l.dat parents_3g4q2l.dat parents_5g2q2l.dat parents_1g6q1y.dat parents_3g4q1y.dat parents_5g2q1y.dat parents_1g6q.dat parents_3g4q.dat parents_5g2q.dat parents_7g.dat


if PUBLIC
J_TESTS= 2j_ME2_check 3j_ME2_check 4j_ME2_check 4j_loop_ME2_check
Z_TESTS= Z1j_ME2_check Z2j_ME2_check Z3j_ME2_check 
Y_TESTS= Y1j_ME2_check Y2j_ME2_check Y3j_ME2_check 
YY_TESTS=  
W_TESTS= W_ME2_check W1j_ME2_check  W2j_ME2_check W3j_partial_and_ME2_check 
PS_TESTS= Wm_PS Wm1_PS Wm2_PS Wm3_PS Z_PS Z1_PS Z2_PS Z3_PS Y1_PS Y2_PS Y3_PS 2j_PS 3j_PS
else
J_TESTS= 2j_ME2_check G3_3j_ME2_check 3j_ME2_check 4j_ME2_check 4j_loop_ME2_check 5j_ME2_check
Z_TESTS= Z1j_ME2_check Z2j_ME2_check Z3j_ME2_check Z4j_ME2_check Z5j_ME2_check
YY_TESTS= YY_ME2_check YY1j_ME2_check YY2j_ME2_check YY3j_ME2_check YY4j_ME2_check 
Y_TESTS= Y1j_ME2_check Y2j_ME2_check Y3j_ME2_check Y4j_ME2_check Y5j_ME2_check
W_TESTS= W_ME2_check W1j_ME2_check W2j_ME2_check W3j_partial_and_ME2_check W4j_ME2_check W4j_ME2_check_auto W5j_ME2_check
PS_TESTS= Wm_PS Wm1_PS Wm2_PS Wm3_PS Z_PS Z1_PS Z2_PS Z3_PS Y1_PS Y2_PS Y3_PS 2j_PS 3j_PS
endif

SHORT_TESTS= $(W_TESTS) gkm_check sub_leading_test Z4q_test Z2q2G_test Z2qNg_test 2q2GNg_test photon_test ME_5parton_2lept_check ME_6parton_2lepton_lc_check  $(Y_TESTS) $(Z_TESTS) $(J_TESTS) Vnj_bornME2_check higgs_check

LONG_TESTS= $(W_TESTS) gkm_check sub_leading_test Z4q_test Z2q2G_test Z2qNg_test 2q2GNg_test photon_test ME_5parton_2lept_check ME_6parton_2lepton_lc_check  $(Y_TESTS) $(Z_TESTS) $(J_TESTS) Vnj_bornME2_check higgs_check

if ! PUBLIC
COLLINEAR_TESTS = SplittingTestTree SplittingTestnf SplittingTestPartial SplittingTestPartial4q
endif




#generate_data compare hist.sh

if ! PUBLIC
check_PROGRAMS = \
		gkm_check \
		sub_leading_test \
		Z4q_test Z2q2G_test Z2qNg_test 2q2GNg_test \
		photon_test \
		BH2Sherpa_test \
		ME_5parton_2lept_check ME_6parton_2lepton_lc_check  \
		$(W_TESTS) \
		$(Y_TESTS) \
		$(Z_TESTS) \
		$(YY_TESTS) \
		$(J_TESTS) \
		Vnj_bornME2_check higgs_check \
		$(PS_TESTS) \
		check_parent_diagrams

TESTS= \
	short_test primitive_test

else
check_PROGRAMS = \
		$(W_TESTS) \
		$(Y_TESTS) \
		$(Z_TESTS) \
		$(YY_TESTS) \
		$(J_TESTS) \
		$(PS_TESTS) 

TESTS= \
	$(W_TESTS) \
		$(Y_TESTS) \
		$(Z_TESTS) \
		$(YY_TESTS) \
		$(J_TESTS) \
		$(PS_TESTS) 

endif

if !PUBLIC
check_PROGRAMS += check_parent_diagrams

if BH_USE_GMP
check_PROGRAMS += CollinearLimit
endif
endif

W_ME2_check_SOURCES = W_ME2_check.cpp
W1j_ME2_check_SOURCES = W1j_ME2_check.cpp
W2j_ME2_check_SOURCES = W2j_ME2_check.cpp
W3j_partial_and_ME2_check_SOURCES = W3j_partial_and_ME2_check.cpp 

Y1j_ME2_check_SOURCES = Y1j_ME2_check.cpp
Y2j_ME2_check_SOURCES = Y2j_ME2_check.cpp
Y3j_ME2_check_SOURCES = Y3j_ME2_check.cpp

Z1j_ME2_check_SOURCES = Z1j_ME2_check.cpp
Z2j_ME2_check_SOURCES = Z2j_ME2_check.cpp
Z3j_ME2_check_SOURCES = Z3j_ME2_check.cpp


Wm_PS_SOURCES = Wm_PS.cpp
Wm1_PS_SOURCES = Wm1_PS.cpp
Wm2_PS_SOURCES = Wm2_PS.cpp
Wm3_PS_SOURCES = Wm3_PS.cpp
Z_PS_SOURCES = Z_PS.cpp PStest.hpp theconstants.h
Z1_PS_SOURCES = Z1_PS.cpp
Z2_PS_SOURCES = Z2_PS.cpp
Z3_PS_SOURCES = Z3_PS.cpp
Y1_PS_SOURCES = Y1_PS.cpp
Y2_PS_SOURCES = Y2_PS.cpp
Y3_PS_SOURCES = Y3_PS.cpp
2j_PS_SOURCES = 2j_PS.cpp
3j_PS_SOURCES = 3j_PS.cpp

2j_ME2_check_SOURCES = 2j_ME2_check.cpp
3j_ME2_check_SOURCES = 3j_ME2_check.cpp
4j_ME2_check_SOURCES = 4j_ME2_check.cpp
4j_loop_ME2_check_SOURCES = 4j_loop_ME2_check.cpp


if ! PUBLIC
nodist_G3_3j_ME2_check_SOURCES = G3_3j_ME2_check.cpp
nodist_gkm_check_SOURCES = gkm_check.cpp 
nodist_sub_leading_test_SOURCES = sub_leading_test.cpp  
nodist_Z4q_test_SOURCES = Z4q_test.cpp 
nodist_Z2q2G_test_SOURCES = Z2q2G_test.cpp 
nodist_Z2qNg_test_SOURCES = Z2qNg_test.cpp 
nodist_2q2GNg_test_SOURCES = 2q2GNg_test.cpp
nodist_photon_test_SOURCES = photon_test.cpp 
nodist_BH2Sherpa_test_SOURCES = BH2Sherpa_test.cpp 
nodist_ME_5parton_2lept_check_SOURCES = ME_5parton_2lept_check.cpp 
nodist_ME_6parton_2lepton_lc_check_SOURCES = ME_6parton_2lepton_lc_check.cpp
nodist_W4j_ME2_check_SOURCES = W4j_ME2_check.cpp
nodist_W4j_ME2_check_auto_SOURCES = W4j_ME2_check_auto.cpp
nodist_W5j_ME2_check_SOURCES = W5j_ME2_check.cpp
nodist_Y4j_ME2_check_SOURCES = Y4j_ME2_check.cpp
nodist_Y5j_ME2_check_SOURCES = Y5j_ME2_check.cpp
nodist_Z4j_ME2_check_SOURCES = Z4j_ME2_check.cpp
nodist_Z5j_ME2_check_SOURCES = Z5j_ME2_check.cpp
nodist_YY_ME2_check_SOURCES = YY_ME2_check.cpp
nodist_YY1j_ME2_check_SOURCES = YY1j_ME2_check.cpp
nodist_YY2j_ME2_check_SOURCES = YY2j_ME2_check.cpp
nodist_YY3j_ME2_check_SOURCES = YY3j_ME2_check.cpp
nodist_YY4j_ME2_check_SOURCES = YY4j_ME2_check.cpp
nodist_5j_ME2_check_SOURCES = 5j_ME2_check.cpp
nodist_Vnj_bornME2_check_SOURCES = Vnj_bornME2_check.cpp
nodist_higgs_check_SOURCES = higgs_check.cpp
nodist_check_parent_diagrams_SOURCES =	check_parent_diagrams.cpp
if BH_USE_GMP
CollinearLimit_SOURCES = CollinearLimit.cpp
endif
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
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/trees/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/trees_eval/.libs 

if ! PUBLIC
LINKER_OPTS += \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/../FeynDiagram/.libs 
endif

LDFLAGS += $(LINKER_OPTS)


if BH_USE_GMP
LIBSTOLINK += -lmpfr -lgmpxx -lgmp 
endif
if BH_USE_GMP
LIBSTOLINK += $(BH_PATH_OBJ)/.libs//libmpreal.la    
endif

#LIBS = $(LIBSTOLINK) $(LIBS) 


W3j_partial_and_ME2_check_LDADD =  $(LIBSTOLINK)
W_ME2_check_LDADD =  $(LIBSTOLINK)
W1j_ME2_check_LDADD =  $(LIBSTOLINK)
W2j_ME2_check_LDADD =  $(LIBSTOLINK)

Z1j_ME2_check_LDADD =  $(LIBSTOLINK)
Z2j_ME2_check_LDADD =  $(LIBSTOLINK)
Z3j_ME2_check_LDADD =  $(LIBSTOLINK)

Y1j_ME2_check_LDADD =  $(LIBSTOLINK)
Y2j_ME2_check_LDADD =  $(LIBSTOLINK)
Y3j_ME2_check_LDADD =  $(LIBSTOLINK)

Wm_PS_LDADD =  $(LIBSTOLINK)
Wm1_PS_LDADD =  $(LIBSTOLINK)
Wm2_PS_LDADD =  $(LIBSTOLINK)
Wm3_PS_LDADD =  $(LIBSTOLINK)
Z_PS_LDADD =  $(LIBSTOLINK)
Z1_PS_LDADD =  $(LIBSTOLINK)
Z2_PS_LDADD =  $(LIBSTOLINK)
Z3_PS_LDADD =  $(LIBSTOLINK)
Y1_PS_LDADD =  $(LIBSTOLINK)
Y2_PS_LDADD =  $(LIBSTOLINK)
Y3_PS_LDADD =  $(LIBSTOLINK)
2j_PS_LDADD =  $(LIBSTOLINK)
3j_PS_LDADD =  $(LIBSTOLINK)

2j_ME2_check_LDADD =  $(LIBSTOLINK)
3j_ME2_check_LDADD =  $(LIBSTOLINK)
4j_ME2_check_LDADD =  $(LIBSTOLINK)
G3_3j_ME2_check_LDADD =  $(LIBSTOLINK)
4j_loop_ME2_check_LDADD =  $(LIBSTOLINK)

if ! PUBLIC
gkm_check_LDADD = $(LIBSTOLINK)
sub_leading_test_LDADD = $(LIBSTOLINK)
Z4q_test_LDADD =  $(LIBSTOLINK)
Z2q2G_test_LDADD =  $(LIBSTOLINK)
Z2qNg_test_LDADD =  $(LIBSTOLINK)
2q2GNg_test_LDADD =  $(LIBSTOLINK)
photon_test_LDADD =  $(LIBSTOLINK)
BH2Sherpa_test_LDADD =  $(LIBSTOLINK)
ME_5parton_2lept_check_LDADD =  $(LIBSTOLINK)
ME_6parton_2lepton_lc_check_LDADD =  $(LIBSTOLINK)
W4j_ME2_check_LDADD =  $(LIBSTOLINK)
W4j_ME2_check_auto_LDADD =  $(LIBSTOLINK)
W5j_ME2_check_LDADD =  $(LIBSTOLINK)
Y4j_ME2_check_LDADD =  $(LIBSTOLINK)
Y5j_ME2_check_LDADD =  $(LIBSTOLINK)
Z4j_ME2_check_LDADD =  $(LIBSTOLINK)
Z5j_ME2_check_LDADD =  $(LIBSTOLINK)
YY_ME2_check_LDADD =  $(LIBSTOLINK)
YY1j_ME2_check_LDADD =  $(LIBSTOLINK)
YY2j_ME2_check_LDADD =  $(LIBSTOLINK)
YY3j_ME2_check_LDADD =  $(LIBSTOLINK)
YY4j_ME2_check_LDADD =  $(LIBSTOLINK)
5j_ME2_check_LDADD =  $(LIBSTOLINK)
Vnj_bornME2_check_LDADD = $(LIBSTOLINK)
higgs_check_LDADD = $(LIBSTOLINK)
check_parent_diagrams_LDADD = $(LIBSTOLINK)	
if BH_USE_GMP
CollinearLimit_LDADD = $(LIBSTOLINK)
endif
endif

short_test:	$(SHORT_TESTS)
	perl $(srcdir)/do_tests.pl $(SHORT_TESTS)

j_test:	$(J_TESTS)
	perl $(srcdir)/do_tests.pl $(J_TESTS)

z_test:	$(Z_TESTS)
	perl $(srcdir)/do_tests.pl $(Z_TESTS)

yy_test:	$(YY_TESTS)
	perl $(srcdir)/do_tests.pl $(YY_TESTS)

y_test:	$(Y_TESTS)
	perl $(srcdir)/do_tests.pl $(Y_TESTS)

w_test:	$(W_TESTS)
	perl $(srcdir)/do_tests.pl $(W_TESTS)

ps_test:	$(PS_TESTS)
	perl $(srcdir)/do_tests.pl $(PS_TESTS)

primitive_test:	check_parent_diagrams
	perl -s $(srcdir)/do_tests.pl -primitive $(PRIMITIVE_TESTS)

long_primitive_test: check_parent_diagrams
	perl -s $(srcdir)/do_tests.pl -longprimitive $(LONG_PRIMITIVE_TESTS)

collinear_test:	$(COLLINEAR_TESTS)
	perl $(srcdir)/do_tests.pl $(COLLINEAR_TESTS)

bsub_test:      $(SHORT_TESTS)
	for t in $(SHORT_TESTS) ; do \
                echo "$$PWD/.libs/$$t > $$PWD/$${t}.log.$$(date +%d%m%y_%H%M%S) " > script_$${t}.sh ;\
                chmod +x  script_$${t}.sh;   bsub -q 1nd -J $$t script_$${t}.sh; \
        done
if  PUBLIC
check_DATA = dataInstall.log

dataInstall.log:	$(top_srcdir)/BHdata.tar.gz 
	sh $(abs_top_builddir)/dataInstall $(top_srcdir)/BHdata.tar.gz > dataInstall.log

endif

EXTRA_DIST = do_tests.pl testAll.py singleTest.sh CollTest.py PStest

#dist_data_DATA = hist.sh compare_mpmpmp.grace.in compare_mmpmpp.grace.in compare_mmpmmp.grace.in
