
pkglib_LTLIBRARIES = libAmpl_eval.la		     		
		     		

include A0_2q3g2l_eval.inc

SPLITFILES = $(ALL_A0_2q3g2l_eval)

libAmpl_eval_la_SOURCES = A0_3g_eval.cpp A0_4g_eval.cpp A0_5g_eval.cpp A0_6g_eval.cpp \
		     		A0_2q1g_eval.cpp A0_2q2g_eval.cpp A0_2q3g_eval.cpp \
		     		A0_2q2l_eval.cpp A0_2q1g2l_eval.cpp A0_2q2g2l_eval.cpp A0_2q4g2l_eval.cpp A0_2q5g2l_eval.cpp \
		     		A0_2q2Q_eval.cpp A0_4q_eval.cpp A0_4q1g_eval.cpp A0_4q2g_eval.cpp A0_2s_eval.cpp\
		     		A0_2q1g2Q_eval.cpp A0_2q2g2Q_eval.cpp A0_2q1_2q2_2q3_eval.cpp\
		     		A0_4q1_2q2_eval.cpp A0_2q_massive_eval.cpp A0_2q_massive_lepton_eval.cpp\
		     		A0_2q2Q2l_eval.cpp A0_2q2Q1g2l_eval.cpp A0_2q2Q2g2l_eval.cpp A0_2q2Q3g2l_eval.cpp\
				A0_2q2l2G_eval.cpp A0_2q1g1y_eval.cpp A0_2q2g1y_eval.cpp A0_2q2Q1y_eval.cpp A0_2q2Gl1y_eval.cpp\
		     		amplitudes_tree_eval.h amplitudes_2q_eval.cpp\
		     		A0_1phi_eval.cpp A0_1phi_scalar.cpp A0_1phi_scalar_2.cpp A0_1phi_Q_massive.cpp\
				A0_sc_massless.cpp\
				A0_2Gsc_eval.cpp A0_1phi_scalar_massive.cpp\
				A0_G3_eval.cpp\
				$(SPLITFILES)

if ! PUBLIC
nodist_libAmpl_eval_la_SOURCES = \
				amplitudes_tree_grass.h A_Nglu_Tree_Ptr_Grass.cpp A_2qNg_Tree_Ptr_Grass.cpp A_2qNg2l_Tree_Ptr_Grass.cpp\
				A_4qNg_Tree_Ptr_Grass.cpp A_6qNg_Tree_Ptr_Grass.cpp 

endif

#include_HEADERS = amplitudes_tree_eval.h
