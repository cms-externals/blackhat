/*
 * momenta_evaluator.h
 *
 *  Created on: 10 Aug 2009
 *      Author: daniel
 */

#ifndef MOMENTA_EVALUATOR_H_
#define MOMENTA_EVALUATOR_H_

#include "BH_typedefs.h"
#include <vector>
#if BH_USE_GMP
#include "gmp_r.h"
#endif
#include "BH_debug.h"

namespace BH {

template <class T> class momentum_configuration;
template <class T, int N> class coeffparam;
template <class T> class eval_param;
template <class T> class Cmom;

namespace cut {

namespace Darren {

class MomentumEvaluator {
public:
	MomentumEvaluator(){};
	~MomentumEvaluator(){};
};

template <int TPOINTSBUB, int YPOINTS> class bubble_Darren_eval_points : public MomentumEvaluator {
	static C _circpos[TPOINTSBUB];
	static CHP _circpos_HP[TPOINTSBUB];
	static CVHP _circpos_VHP[TPOINTSBUB];

	static C _circpos_y[YPOINTS];
	static CHP _circpos_y_HP[YPOINTS];
	static CVHP _circpos_y_VHP[YPOINTS];

	static C _circpos_y_matrix[YPOINTS*YPOINTS];
	static CHP _circpos_y_matrix_HP[YPOINTS*YPOINTS];
	static CVHP _circpos_y_matrix_VHP[YPOINTS*YPOINTS];

#if BH_USE_GMP
	static std::complex<RGMP> _circpos_GMP[TPOINTSBUB];
	static std::complex<RGMP> _circpos_y_GMP[YPOINTS];
	static std::complex<RGMP> _circpos_y_matrix_GMP[YPOINTS*YPOINTS];
	static int s_GMP_precision;
#endif

	static long int _nbr_instances;

public:

	typedef int MomentaGridType[YPOINTS][TPOINTSBUB];

	enum { Ypoints = YPOINTS};
	enum { Tpointsbub = TPOINTSBUB};
	bubble_Darren_eval_points(){init();};
    virtual ~bubble_Darren_eval_points(){}

    // Get the points on the circle for evaluating the t extraction in the bubble
    static void get_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos;}
    static void get_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_HP;}
    static void get_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_VHP;}

    static void get_yeval_pts(const C*& ptr_circpos){ptr_circpos=_circpos_y;}
    static void get_yeval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_y_HP;}
    static void get_yeval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_y_VHP;}

    static void get_y_matrix_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos_y_matrix;}
    static void get_y_matrix_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_y_matrix_HP;}
    static void get_y_matrix_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_y_matrix_VHP;}

#if BH_USE_GMP
    static void get_eval_pts(const std::complex<RGMP>*& ptr_circpos){if (RGMP::get_current_precision()> s_GMP_precision){ BH_DEBUG_MESSAGE2("recomputing y points for precision ",RGMP::get_current_precision());get_y_points(); }; ptr_circpos=_circpos_GMP;}
    static void get_yeval_pts(const std::complex<RGMP>*& ptr_circpos){if (RGMP::get_current_precision()> s_GMP_precision); ptr_circpos=_circpos_y_GMP;}
    static void get_y_matrix_eval_pts(const std::complex<RGMP>*& ptr_circpos){if (RGMP::get_current_precision()> s_GMP_precision); ptr_circpos=_circpos_y_matrix_GMP;}
#endif

    template <class T,int CPOINTS,class cutDbase> void compute_momenta(momentum_configuration<T>& mc,const std::vector<int>& ind ,MomentaGridType (&result)[2],MomentaGridType (&resultm)[2] ,coeffparam<T,CPOINTS>& tp,const cutDbase& cb);
	template <class T,int CPOINTS,class cutDbase> void compute_momenta(const eval_param<T>& ep,Cmom<T> (&momenta)[2][YPOINTS][TPOINTSBUB],Cmom<T> (&momenta_m)[2][YPOINTS][TPOINTSBUB] ,coeffparam<T,CPOINTS>& tp,const cutDbase& cb);
private:
    void init();

    static const bool get_y_points();
};


template <int YPOINTS,int TPOINTSBUB,class MomentumEvaluator> class General_Bubble_Combiner {
public:
	template <class T> std::complex<T> combine(
			std::complex<T> (&ampl)[YPOINTS],
			std::complex<T> (&ampl_err)[YPOINTS],
			std::complex<T>& trisuby,
			std::complex<T>& amp_error,
			const RVHP& tri_sub_acc,
			RVHP& accuracy);
};

template <int YPOINTS,int TPOINTSBUB,class MomentumEvaluator> class Normal_Bubble_Combiner {
public:
	template <class T> std::complex<T> combine(
			std::complex<T> (&ampl)[YPOINTS],
			std::complex<T> (&ampl_err)[YPOINTS],
			std::complex<T>& trisuby,
			std::complex<T>& amp_error,
			const RVHP& tri_sub_acc,
			RVHP& accuracy);
};


template <class MomentumEvaluator,class CutType,int CPOINTS> class Normal_Corner_Tree_Strategy : public MomentumEvaluator {
public:
	template <class T> void fill_trees(
			BH::momentum_configuration<T>&,
			std::vector<int> const&,
		    std::vector<std::complex<T> >& trees_result_1,
		    std::vector<std::complex<T> >& trees_result_2,
		    coeffparam<T,CPOINTS>& tp,
		    CutType& self);
	template <class T> void fill_trees(
			BH::momentum_configuration<T>&,
			std::vector<int> const&,
		    std::vector<std::complex<T> >& trees_result_1,
		    std::vector<std::complex<T> >& trees_result_2,
		    std::vector<std::complex<T> >& trees_result_3,
		    coeffparam<T,CPOINTS>& tp,
		    CutType& self,int k1leg,int k2leg,int k3leg,int masslessleg_type,int reverse);
	template <class T> void fill_trees(
				const eval_param<T>&,
                eval_param<T>**& epc,
			    std::vector<std::complex<T> >& trees_result_1,
			    std::vector<std::complex<T> >& trees_result_2,
			    coeffparam<T,CPOINTS>& tp,
			    CutType& self);
};

template <class MomentumEvaluator,class CutType,int CPOINTS> class Skeleton_Corner_Tree_Strategy : public MomentumEvaluator {
public:
	template <class T> void fill_trees(
			BH::momentum_configuration<T>&,
			std::vector<int> const&,
			std::vector<std::complex<T> >& trees_result_1,
			std::vector<std::complex<T> >& trees_result_2,
			coeffparam<T,CPOINTS>& tp,
			CutType& self);
	template <class T> void fill_trees(
			BH::momentum_configuration<T>&,
			std::vector<int> const&,
		    std::vector<std::complex<T> >& trees_result_1,
		    std::vector<std::complex<T> >& trees_result_2,
		    std::vector<std::complex<T> >& trees_result_3,
		    coeffparam<T,CPOINTS>& tp,
		    CutType& self,int k1leg,int k2leg,int k3leg,int masslessleg_type,int reverse);
	template <class T> void fill_trees(
			const eval_param<T>&,
            eval_param<T>**&,
			std::vector<std::complex<T> >& trees_result_1,
			std::vector<std::complex<T> >& trees_result_2,
			coeffparam<T,CPOINTS>& tp,
			CutType& self);
};



template <int CPOINTS, int TPOINTSTRI> class triangle_Darren_eval_points : public MomentumEvaluator {
	static C coeffext[CPOINTS*TPOINTSTRI];
	static CHP coeffext_HP[CPOINTS*TPOINTSTRI];
	static CVHP coeffext_VHP[CPOINTS*TPOINTSTRI];

#if BH_USE_GMP
	static CGMP coeffext_GMP[CPOINTS*TPOINTSTRI];
	static int s_GMP_precision;
#endif


	static long int _nbr_instances;

public:
	enum { Cpoints = CPOINTS};
	enum { Tpointstri = TPOINTSTRI};
	typedef int MomentaGridType[3*TPOINTSTRI];

	triangle_Darren_eval_points(){init();};
    virtual ~triangle_Darren_eval_points(){};

    // Get the points on the circle for extracting the 7 coeffs of C
    static void get_coeff_pts(const C*& ptr_circpos){ptr_circpos=coeffext;}
    static void get_coeff_pts(const CHP*& ptr_circpos){ptr_circpos=coeffext_HP;}
    static void get_coeff_pts(const CVHP*& ptr_circpos){ptr_circpos=coeffext_VHP;}
#if BH_USE_GMP
    static void get_coeff_pts(const CGMP*& ptr_circpos){if (RGMP::get_current_precision()> s_GMP_precision){get_matrix_points();};ptr_circpos=coeffext_GMP;}
#endif

    template <class T,class cutDbase> void compute_momenta(momentum_configuration<T>& mc,const std::vector<int>& ind ,MomentaGridType& result,MomentaGridType& resultm,coeffparam<T,CPOINTS>& tp,const cutDbase& cb,int k1leg,int k2leg,int k3leg,int masslessleg_type);

private:
    void init();

    static const bool get_matrix_points();
};



template <int TPOINTSTRI> class box_Darren_eval_points : public MomentumEvaluator {
	//The points we evaluate the circle in t around
	static C circpos[TPOINTSTRI];
	static CHP circpos_HP[TPOINTSTRI];
	static CVHP circpos_VHP[TPOINTSTRI];
#if BH_USE_GMP
	static CGMP circpos_GMP[TPOINTSTRI];
	static int s_GMP_precision;
#endif

	static long int _nbr_instances;

public:
	box_Darren_eval_points(){init();};
    virtual ~box_Darren_eval_points(){}

    // Get the points on the circle for evaluating the t extraction
    static void get_eval_pts(const C*& ptr_circpos){ptr_circpos=circpos;}
    static void get_eval_pts(const CHP*& ptr_circpos){ptr_circpos=circpos_HP;}
    static void get_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=circpos_VHP;}
#if BH_USE_GMP
    static void get_eval_pts(const CGMP*& ptr_circpos){if (RGMP::get_current_precision()> s_GMP_precision){get_t_points();}; ptr_circpos=circpos_GMP;}
#endif

private:
    void init();

    static const bool get_t_points();
};



} /* Darren */

} /* cut */

} /* BH */

#endif /* MOMENTA_EVALUATOR_H_ */
