/*
 * bubble_ratext.h
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef BUBBLE_RATEXT_H_
#define BUBBLE_RATEXT_H_

#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "ratext/rat_ext.h"
#include "eval_param.h"
#include "momenta_evaluator_ratext.h"
#include "BH_debug.h"


namespace BH {
template <class T> class momentum_configuration;
}



namespace BH {

namespace ratext {

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class box_rat_eval_points; // Declared in box_ratext.h

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class bub_rat_eval_points : public MomentumEvaluator {
	static C _circpos[CBUBPOINTS_T];
	static CHP _circpos_HP[CBUBPOINTS_T];
	static CVHP _circpos_VHP[CBUBPOINTS_T];

	static C _circpos_matrix[CBUBPOINTS_T*CBUBPOINTS_T];
	static CHP _circpos_matrix_HP[CBUBPOINTS_T*CBUBPOINTS_T];
	static CVHP _circpos_matrix_VHP[CBUBPOINTS_T*CBUBPOINTS_T];

	static C _y_circpos[YPOINTS_T];
	static CHP _y_circpos_HP[YPOINTS_T];
	static CVHP _y_circpos_VHP[YPOINTS_T];

	static C _y_circpos_matrix[YPOINTS_T*YPOINTS_T];
	static CHP _y_circpos_matrix_HP[YPOINTS_T*YPOINTS_T];
	static CVHP _y_circpos_matrix_VHP[YPOINTS_T*YPOINTS_T];

	// For now we only store the results for at most a mu^4 bubble, this is good for a max power of l^5
	static const C _rat_int[4];
	static const CHP _rat_int_HP[4];
	static const CVHP _rat_int_VHP[4];

	// For now we only store the results for at most a mu^4 bubble, this is good for a max power of l^5
	static const C _extra_fac[4];
	static const CHP _extra_fac_HP[4];
	static const CVHP _extra_fac_VHP[4];

	static long int _nbr_instances;

#if BH_USE_GMP
	static CGMP _circpos_GMP[CBUBPOINTS_T];
	static CGMP _circpos_matrix_GMP[CBUBPOINTS_T*CBUBPOINTS_T];
	static CGMP _y_circpos_GMP[YPOINTS_T];
	static CGMP _y_circpos_matrix_GMP[YPOINTS_T*YPOINTS_T];
	static  CGMP _rat_int_GMP[4];
	static  CGMP _extra_fac_GMP[4];  // precision can change
	static int s_GMP_precision;
#endif

public:
	enum { MUPOINTS = MUPOINTS_T};
	enum { MUTRIPOINTS = MUTRIPOINTS_T};
	enum { MUBUBPOINTS = MUBUBPOINTS_T};
	enum { CPOINTS = CPOINTS_T};
	enum { CBUBPOINTS = CBUBPOINTS_T};
	enum { YPOINTS = YPOINTS_T};
	bub_rat_eval_points(){init();};
    virtual ~bub_rat_eval_points(){}

    // Get the points on the circle for evaluating the t extraction in the bubble
    static void get_t_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos;}
    static void get_t_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_HP;}
    static void get_t_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_VHP;}

    // Get the points on the circle for evaluating the t extraction in the bubble
    static void get_t_matrix_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos_matrix;}
    static void get_t_matrix_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_matrix_HP;}
    static void get_t_matrix_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_matrix_VHP;}

    // Get the points on the circle for evaluating the y extraction in the bubble
    static void get_y_eval_pts(const C*& ptr_circpos){ptr_circpos=_y_circpos;}
    static void get_y_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_y_circpos_HP;}
    static void get_y_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_y_circpos_VHP;}

    // Get the points on the circle for evaluating the y extraction in the bubble
    static void get_y_matrix_eval_pts(const C*& ptr_circpos){ptr_circpos=_y_circpos_matrix;}
    static void get_y_matrix_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_y_circpos_matrix_HP;}
    static void get_y_matrix_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_y_circpos_matrix_VHP;}

    // Get the rational integral coefficient
    static const C get_rat_integral(int element, const C& s1);
    static const CHP get_rat_integral(int element, const CHP& s1);
    static const CVHP get_rat_integral(int element, const CVHP& s1);

    // Get the y^(2+n) for n>=0 for the extra contribution from the y^(2+n) terms to the bubble
    static const C get_extra_fac(int element, const C& s1);
    static const CHP get_extra_fac(int element, const CHP& s1);
    static const CVHP get_extra_fac(int element, const CVHP& s1);

#if BH_USE_GMP
    static void get_t_eval_pts(const CGMP*& ptr_circpos){ if (RGMP::get_current_precision()> s_GMP_precision){ recomuteGMP(); }; ptr_circpos=_circpos_GMP;}
    static void get_t_matrix_eval_pts(const CGMP*& ptr_circpos){ if (RGMP::get_current_precision()> s_GMP_precision){ recomuteGMP(); };ptr_circpos=_circpos_matrix_GMP;}
    static void get_y_eval_pts(const CGMP*& ptr_circpos){ if (RGMP::get_current_precision()> s_GMP_precision){ recomuteGMP(); };ptr_circpos=_y_circpos_GMP;}
    static void get_y_matrix_eval_pts(const CGMP*& ptr_circpos){ if (RGMP::get_current_precision()> s_GMP_precision){ recomuteGMP(); };ptr_circpos=_y_circpos_matrix_GMP;}
    static const CGMP get_rat_integral(int element, const CGMP& s1);
    static const CGMP get_extra_fac(int element, const CGMP& s1);
    static void recomuteGMP();
#endif

    static box_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> BoxPoints;

	template <class T,class cutDbase> void compute_momenta(const eval_param<T>& ep,Cmom<T> (&momenta)[2],Cmom<T> (&momenta_m)[2],const cutDbase& cb);

	static void init();
private:
    static const bool get_points();
};


//template <class MomentumEvaluator> class General_RatBub_Combiner {
//public:
//	template <class T> std::complex<T> combine(
//			std::complex<T>* ampl,
//			std::complex<T>* ampl_err);
//};
//
//template <class MomentumEvaluator> class Normal_RatBub_Combiner {
//public:
//	template <class T> std::complex<T> combine(
//			std::complex<T>* ampl,
//			std::complex<T>* ampl_err);
//};


template <class MomentumEvaluator,class CutType> class Normal_RatBub_Corner_Tree_Strategy : public MomentumEvaluator {
public:
	template <class T> void fill_trees(
			const eval_param<T>&,
			std::vector<std::complex<T> > (&trees_result)[2],
			std::vector<size_t>& cut_mass,
			CutType& self);
};


template <class CutType> struct Normal_RatBub_Specification {
public:
	enum { MUPOINTS=3 }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=2 }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=2 }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=7 }; // Number of C coefficients
	enum { CBUBPOINTS=3}; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=3}; // Number of points in y used to evaluate for the bubble

	typedef bub_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatBub_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
//	typedef Normal_RatBub_Combiner<MomentaEvaluator> Combiner;
};

template <class CutType, int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> struct General_RatBub_Specification {
public:
	enum { MUPOINTS=MUPOINTS_T }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=MUTRIPOINTS_T }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=MUBUBPOINTS_T }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=CPOINTS_T }; // Number of C coefficients
	enum { CBUBPOINTS=CBUBPOINTS_T }; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=YPOINTS_T}; // Number of points in y used to evaluate for the bubble

	typedef bub_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatBub_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
//	typedef General_RatBub_Combiner<MomentaEvaluator> Combiner;
};

template <class CutType> struct Higgs_RatBub_Specification : public  General_RatBub_Specification<CutType,MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS> {};



//! Class for Darren's rational bubble diagrams
template <class BASE,  class RatBubSpecs/* = Normal_RatBub_Specification<BASE> */> class bubble_Rat : public BASE, public RatBubSpecs::CornerTreeStrategy /*,  public RatBubSpecs::Combiner*/ {
	long int mcID;
	unsigned long int epID;
	std::vector<int> indID;

	C B0coeff;
	CHP B0coeff_HP;
	CVHP B0coeff_VHP;

#if BH_USE_GMP
	CGMP B0coeff_GMP;
#endif
	double d_accuracy;

protected:
	//This vector is used for storing the indexes that get passed to the cut trees
	std::vector<int> indlst[2];

	eval_param<R>* _ep[2];
	eval_param<RHP>* _ep_HP[2];
	eval_param<RVHP>* _ep_VHP[2];

#if BH_USE_GMP
	eval_param<RGMP>* _ep_GMP[2];
#endif

	std::vector<size_t> _cut_mass; // Store the unique mass_params of the cut legs


	void set_accuracy(double acc){d_accuracy=acc;};

public:
	mass_param_coll* _leg_masses; // Store the mass_params for each cut leg
	template <class T> bubble_Rat(T& t): BASE(t), RatBubSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(2)), mcID(-1){init();};
	template <class T> bubble_Rat(T* t): BASE(t), RatBubSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(2)), mcID(-1){ init();};
	virtual ~bubble_Rat();

	virtual C eval(momentum_configuration<R>& mc,const std::vector<int>& ind);
	virtual CHP eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind);
	virtual CVHP eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind);

	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);

#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
	virtual CGMP eval(const eval_param<RGMP>& ep);
#endif

	inline void set_eval(long int mcset, const std::vector<int>& indset){mcID=mcset;indID=indset;};
	inline void set_eval(unsigned long int mcset){epID=mcset;};
	inline bool is_eval(long int mcset, const std::vector<int>& indset){if((mcset==mcID)&&(indset==indID)){return true;} return false;};
	inline bool is_eval(unsigned long int mcset){if(mcset==epID){return true;} return false;};

	double get_accuracy(){return d_accuracy;}; // Returns the estimated accuracy of the computation

	// Add a cutD to the list
	void add_mass(const cutD& pcd);
	void add_mass(int m1,int m2);

   	typedef typename RatBubSpecs::MomentaEvaluator MomentaEvaluatorType;
   	typedef typename daughter_info<bubble_Rat<BASE,RatBubSpecs> >::daughter_type* dau_type_p;
protected:
	void set_B0coeff(C value){B0coeff=value;};
	void set_B0coeff(CHP value){B0coeff_HP=value;B0coeff=to_double(value);};
	void set_B0coeff(CVHP value){B0coeff_VHP=value;B0coeff_HP=to_HP(value);B0coeff=to_double(value);};

	void get_B0coeff(C& value){value=B0coeff;};
	void get_B0coeff(CHP& value){value=B0coeff_HP;};
	void get_B0coeff(CVHP& value){value=B0coeff_VHP;};

	void get_ep(eval_param<R>**& ep){ep=_ep;}
	void get_ep(eval_param<RHP>**& ep){ep= _ep_HP;}
	void get_ep(eval_param<RVHP>**& ep){ep= _ep_VHP;}

#if BH_USE_GMP
	void set_B0coeff(const CGMP& value){B0coeff_GMP=value;B0coeff_HP=to_HP(value);B0coeff_VHP=to_VHP(value);B0coeff=to_double(value);};
	void get_B0coeff(CGMP& value){value=B0coeff_GMP;};
	void get_ep(eval_param<RGMP>**& ep){if (_ep_GMP[0]->ref()->E().real().get_precision() < RGMP::get_current_precision() ){ init(); }; ep= _ep_GMP;}
#endif

private:
	//Computes the actual bubble coefficients
	template <class T> std::complex<T> get_coeffs(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T> get_coeffs(const eval_param<T>& ep);
	void init();
};

}

}
#endif /* BUBBLE_RATEXT_H_ */
