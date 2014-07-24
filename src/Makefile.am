

pkglib_LTLIBRARIES = libBH.la libBG.la libBHcore.la libassembly.la libOLA.la libCutPart.la libRatPart.la libSpinors.la libIntegrals.la 
if ! PUBLIC
pkglib_LTLIBRARIES += libtools.la 
endif

if BH_USE_GMP
pkglib_LTLIBRARIES += libmpreal.la
endif


libBH_la_SOURCES =  \
					BH_typedefs.h mode_dependent_typedefs.h BH_error.h BH_utilities.h BH_particle_names.h BH_forward.h BH_A0.h\
					partitions.h partitions_inline.hpp partitions_fwd.h amplitudes.h tree_amplitudes.h factorization.h \
					Series.h polylog.h polylog.hpp str_constants.h\
					tree.cpp \
		     		histogram.h histogram.cpp \
		     		rec_BB.h rec_BB.cpp \
		     		rat_amp.h tree_amp.h \
		     		known_rec_tree.cpp  \
				cut_part.h OneLoopHelAmpl.h \
		     		computable.h computable.cpp\
		     		extend.cpp \
		     		cut_part_base.h cut_part_base.cpp \
		     		scheme.h scheme.cpp\
		     		process_utils_basic.cpp process_utils_basic.h\
		     		process_utils_massive.cpp process_utils_massive.h\
		     		iterators.h \
		     		settings_reader.cpp settings_reader.h \
					rational_eval.h rational_eval.cpp rec_tree_eval.h rec_tree_eval.cpp \
					worker_utils.h worker_utils.cpp \
					worker_tree.h worker_tree.cpp \
					rational.cpp rational.h  \
					Get3ptType.cpp \
					rec_tree.h process_utils.h

if ! PUBLIC
nodist_libBH_la_SOURCES =  \
		     		partitions.cpp \
		     		FHZ_D.h FHZ_D.cpp 		     		\
		     		intexp.cpp intexp.h IntExpGD.cpp IntExpGD.h  \
		     		options.h options.cpp\
		     		FHZbox_D.cpp FHZbox_D.h FHZtri_D.cpp FHZtri_D.h FHZ_utils.cpp FHZ_utils.h\
		     		darrenD_new.h darrenD_spur.h bubbles_new.cpp bubbles_spur.cpp darren_num_cnst.h \
		     		pole_strategy.cpp\
				grass_tree.h grass_tree.cpp \
		     		inf_rat.h inf_rat.cpp inf_rat_amp.h \
		     		spurious.cpp spurious_ep.cpp spurious_triangle.cpp spurious_box.cpp spurious.h\
					rec_tree.cpp \
					shifts.h shifts.cpp \
		     		partial_order.cpp partial_order.h ordering.cpp ordering.h \
		     		bubbles_new_ep.cpp bubbles_spur_ep.cpp \
					coupling_summary.h coupling_summary.cpp \
					cut_part_spurious.h cut_part_spurious.cpp \
					is_zero.h is_zero.cpp \
					BG_tree.h BG_vertices.h BG_tree.cpp BG_tree_V.cpp BG_tree_cache.cpp BG_vertex.cpp BG_gluon_vertex.cpp BG_gluon_vertex_V.cpp \
					rec_tree_new.h rec_tree_new.cpp rec_tree_cache.cpp grass_tree.h grass_tree.cpp \
					BH_grass.h BH_grass.cpp grass_contour6.cpp grass_contour7.cpp grass_contour8.cpp grass_contour9.cpp \
					process_utils.cpp process_utils.h \
					cut_part.cpp cut_part_ep.cpp 

endif


libIntegrals_la_SOURCES = \
	     		Series.cpp polylog.cpp polylog_HP.hpp polylog_VHP.hpp \
	     		integrals.cpp  integrals.h \
	     		integrals_ep.cpp integrals_ep.h \
	     		cached_integral.cpp cached_integral.h  

libCutPart_la_SOURCES = \
		     		cut_Darren.h cut_Darren.cpp\
		     		cut_Darren.hpp cut_Darren_worker.cpp \
		     		cut_Darren_ep.hpp \
		     		cut_worker.h cut_worker.cpp \
					cut_part_factory.cpp cut_part_factory.h \
		     		cut_part_worker.h cut_part_worker.cpp cut_part_worker_ep.cpp \
		     		standard_cut_part.h standard_cut_part.hpp standard_cut_part.cpp standard_cut_part_ep.hpp standard_cut_part_ep.cpp \
		     		IR_checked_cut_part.h IR_checked_cut_part.cpp \
		     		momenta_evaluator.h momenta_evaluator.hpp momenta_evaluator.cpp \
		     		bubble_specification.h triangle_specification.h triangle_subtraction.h

if ! PUBLIC
nodist_libCutPart_la_SOURCES = \
			cut_Darren_normal.cpp \
			cut_part_normal.h cut_part_normal.cpp cut_part_normal_ep.cpp 
endif

libRatPart_la_SOURCES = 	\
				rational_factory.h rational_factory.cpp \
				rec_rational_factory.h rec_rational_factory.cpp \
				known_rational.h known_rational.cpp known_rational_ep.cpp

if ! PUBLIC				
nodist_libRatPart_la_SOURCES = rec_rational.cpp rec_rational_ep.cpp rec_rational.h \
				spurious_test.h spurious_test.cpp \
				inf_num.cpp inf_num.h  \
				info.h info.cpp
endif

libOLA_la_SOURCES = \
				cached_OLHA.h cached_OLHA.cpp \
				cached_THA.h cached_THA.cpp \
				IR_checked_amplitudes.cpp IR_checked_amplitudes_ep.cpp \
				IR_checked.h  \
				OneLoopHelAmpl.cpp \
				splitamp_tree.cpp splitamp_loop.cpp fact_momenta.cpp \
				cached_EP.h cached_EP.cpp 	

if ! PUBLIC
nodist_libOLA_la_SOURCES = \
		splittingtest.h splittingtest.cpp \
		IR_checked_pole.h IR_checked_pole.cpp
endif     		 

libBHcore_la_SOURCES = \
				spinor.h spinor_inline.h splib.cpp  \
				particles.cpp particles.h \
	     		process.cpp process.h \
				mom_conf.h mom_conf.hpp mom_conf_inline.h  mom_conf.cpp \
				std_momenta.hpp \
		     	amplitudes.cpp  \
		     	genkey.cc genkey.h \
		     	settings.h settings.cpp \
	     		external_momenta.h external_momenta.cpp \
	     		qd_suppl.h qd_suppl.cpp \
	     		index_vector.cpp index_vector.h \
	     		color_structure.cpp \
				color_algebra.h color_algebra.cpp \
	     		arrange_flavors.cpp arrange_flavors.h \
	     		counted.h counted.cpp \
	     		eval_param.h eval_param.cpp\
	     		eval_param_reader.h eval_param_reader.cpp \
	     		constants.h constants.cpp \
	     		BH_debug.h BH_debug.cpp \
				timing.h timing.cpp \
                helcode.h helcode.cpp NoAmpl.h coeff_param.h \
                path.h path.cpp \
                key.cpp key.h \
                from_file.h from_file.cpp \
                multi_precision_constant.h \
				mc_storage.h mc_storage.hpp mc_storage.cpp environment.hpp

if ! PUBLIC
nodist_libBHcore_la_SOURCES = RawAmplitude.cpp \
	     		parent_diagrams_from_file.h  parent_diagrams_from_file.cpp \
	     		parent_diagrams.h  parent_diagrams.cpp 

endif

if BH_USE_GMP
libBHcore_la_SOURCES+= splib_gmp.cpp gmp_r.h gmp_r.cpp	GMP_random.h GMP_random.cpp
libmpreal_la_SOURCES =  mpreal.cpp mpreal.h
endif

libSpinors_la_SOURCES = \
	spinor.h spinor_inline.h splib.cpp  

libassembly_la_SOURCES = \
		     		assembly.h assembly.cpp \
					primitive_ampl_map.h primitive_ampl_map.cpp \
					ME2_from_file.h  ME2_from_file.cpp \
		     		partial_amplitude.cpp 

if ! PUBLIC
libtools_la_SOURCES = \
			print_cutD.cpp print_cutD.h
nodist_libassembly_la_SOURCES = \
			ME2.h ME2.cpp \
			print_cutD.cpp print_cutD.h
endif




libBG_la_SOURCES = tree1.cc tree2.cc tree3.cc util.h util.cc zero.h Tree.h Tree_impl.h treebase.cc BerendsGiele.cpp BerendsGiele.h BerendsGiele_impl.h particleid.h 


#dist_data_DATA = std_momenta.cpp partitions_inline.cpp polylog_HP.cpp polylog_VHP.cpp IntExpGD.cpp  

pkginclude_HEADERS = setable.h BH_error.h

#pkginclude_HEADERS = spinor.h mom_conf.h mom_conf_inline.h spinor_inline.h qd_suppl.h \
					 BH_typedefs.h particles.h BH_error.h BH_utilities.h tree_amplitudes.h counted.h process.h cut_part.h setable.h
