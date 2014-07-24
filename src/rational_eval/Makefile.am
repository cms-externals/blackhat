

CXXCOMPILESLOW = $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) -g


%.slow:	%.cpp
	$(CXXCOMPILESLOW) -MT $*.o -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $*.o $<
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
	touch $@ 
#	source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
#	$(CXXCOMPILE) -c -o $@ $<

slow:	R_2q2G1g2l_eval.slow R_2q3g2l_eval.slow R_2q2g2l_eval.slow	



pkglib_LTLIBRARIES = libRateval.la 

		     		
libRateval_la_SOURCES = R_4g_eval.cpp R_4g_eval.h \
					R_5g_eval.cpp  R_5g_eval.h\
					R_6g_eval.cpp R_6g_eval.h \
					R_7g_eval.cpp  R_7g_eval.h\
					R_2q2g_eval.cpp R_2q2g_eval.h \
					R_2g2q_eval.cpp R_2g2q_eval.h \
					R_2q3g_eval.cpp R_2q3g_eval.h \
					R_2q3g_r_eval.cpp R_2q3g_r_eval.h \
					R_2q4g_eval.cpp R_2q4g_eval.h \
					R_4q_eval.cpp R_4q_eval.h \
					R_2q2l_eval.cpp R_2q2l_eval.h \
					R_2q1g1y_eval.cpp R_2q1g1y_eval.h \
					R_2q1g2l_eval.cpp R_2q1g2l_eval.h \
					R_2q2g2l_eval.cpp R_2q2g2l_eval.h \
					R_2q3g2l_eval.cpp R_2q3g2l_eval.h \
					R_2q3g2l_place_holder_eval.cpp \
					R_2q2g1y_eval.cpp R_2q2g1y_eval.h \
					R_2q2Q1y_eval.cpp R_2q2Q1y_eval.h \
					R_8g_eval.cpp R_8g_eval.h \
					R_9g_eval.cpp R_9g_eval.h \
					R_galln_eval.cpp R_alln_eval.h \
					R_2q2Q2l_eval.cpp R_2q2Q2l_eval.h \
					R_2q2Gl2l_eval.cpp R_2q2Gl2l_eval.h \
					R_2q2Q_eval.cpp R_2q2Q_eval.h \
					R_2q2Q1g_eval.cpp R_2q2Q1g_eval.h \
					R_2q2G1g2l_eval.cpp R_2q2G1g2l_eval.h \
					R_3g1ph_eval.cpp R_3g1ph_eval.h \
					R_4g1ph_eval.cpp R_4g1ph_eval.h \
					R_2q2g1ph_eval.cpp R_2q2g1ph_eval.h \
					R_2q2Q1ph_eval.cpp R_2q2Q1ph_eval.h \
					amplitudes_rat_eval.h 


