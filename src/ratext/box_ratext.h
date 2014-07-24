/*
 * box_ratext.h
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef BOX_RATEXT_H_
#define BOX_RATEXT_H_

#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "ratext/rat_ext.h"
#include "eval_param.h"
#include "momenta_evaluator_ratext.h"
#include "triangle_ratext.h"
#include "BH_debug.h"

namespace BH {

namespace ratext {

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class box_rat_eval_points : public MomentumEvaluator {
	static C _circpos[MUPOINTS_T];
	static CHP _circpos_HP[MUPOINTS_T];
	static CVHP _circpos_VHP[MUPOINTS_T];

	static C _circpos_matrix[MUPOINTS_T*MUPOINTS_T];
	static CHP _circpos_matrix_HP[MUPOINTS_T*MUPOINTS_T];
	static CVHP _circpos_matrix_VHP[MUPOINTS_T*MUPOINTS_T];

	// For now we only store the results for at most a mu^6 box, this is good for a max power of l^7
	static const C _rat_int[4];
	static const CHP _rat_int_HP[4];
	static const CVHP _rat_int_VHP[4];

#if BH_USE_GMP
	static CGMP _circpos_GMP[MUPOINTS_T];
	static CGMP _circpos_matrix_GMP[MUPOINTS_T*MUPOINTS_T];
	static CGMP _rat_int_GMP[4];
	static int s_GMP_precision;
#endif

	static long int _nbr_instances;

public:
	enum { MUPOINTS = MUPOINTS_T};
	enum { MUTRIPOINTS = MUTRIPOINTS_T};
	enum { MUBUBPOINTS = MUBUBPOINTS_T};
	enum { CPOINTS = CPOINTS_T};
	enum { CBUBPOINTS = CBUBPOINTS_T};
	enum { YPOINTS = YPOINTS_T};
	box_rat_eval_points(){init();};
    virtual ~box_rat_eval_points(){}

    // Get the points on the circle for evaluating the t extraction in the bubble
    static void get_mu_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos;}
    static void get_mu_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_HP;}
    static void get_mu_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_VHP;}

    // Get the points on the circle for evaluating the t extraction in the bubble
    static void get_mu_matrix_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos_matrix;}
    static void get_mu_matrix_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_matrix_HP;}
    static void get_mu_matrix_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_matrix_VHP;}

    // Get the rational integral coefficient
    static const C get_rat_integral(int element, const C& s1);
    static const CHP get_rat_integral(int element, const CHP& s1);
    static const CVHP get_rat_integral(int element, const CVHP& s1);

#if BH_USE_GMP
    static void get_mu_eval_pts(const CGMP*& ptr_circpos){  if (RGMP::get_current_precision()> s_GMP_precision){ recomputeGMP(); };  ptr_circpos=_circpos_GMP;}
    static void get_mu_matrix_eval_pts(const CGMP*& ptr_circpos){ptr_circpos=_circpos_matrix_GMP;}
    static const CGMP get_rat_integral(int element, const CGMP& s1);
	static void recomputeGMP();
#endif

    static tri_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> TriPoints;

	template <class T,class cutDbase> void compute_momenta(const eval_param<T>& ep,Cmom<T> (&momenta)[4],Cmom<T> (&momenta_m)[4], const cutDbase& cb);

    static void init();

private:
    static const bool get_points();
};


//template <class MomentumEvaluator> class General_RatBox_Combiner {
//public:
//	template <class T> std::complex<T> combine(
//			std::complex<T>* ampl,
//			std::complex<T>* ampl_err);
//};
//
//template <class MomentumEvaluator> class Normal_RatBox_Combiner {
//public:
//	template <class T> std::complex<T> combine(
//			std::complex<T>* ampl,
//			std::complex<T>* ampl_err);
//};


template <class MomentumEvaluator,class CutType> class Normal_RatBox_Corner_Tree_Strategy : public MomentumEvaluator {
public:
	template <class T> void fill_trees(
			const eval_param<T>&,
			std::vector<std::complex<T> > (&trees_result)[4],
			std::vector<size_t>& cut_mass,
			CutType& self);
};


template <class CutType> struct Normal_RatBox_Specification {
public:
	enum { MUPOINTS=3 }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=2 }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=2 }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=7 }; // Number of C coefficients
	enum { CBUBPOINTS=3}; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=3}; // Number of points in y used to evaluate for the bubble

	typedef box_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatBox_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
//	typedef Normal_RatBox_Combiner<MomentaEvaluator> Combiner;
};

template <class CutType, int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> struct General_RatBox_Specification {
public:
	enum { MUPOINTS=MUPOINTS_T }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=MUTRIPOINTS_T }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=MUBUBPOINTS_T }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=CPOINTS_T }; // Number of C coefficients
	enum { CBUBPOINTS=CBUBPOINTS_T }; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=YPOINTS_T}; // Number of points in y used to evaluate for the bubble

	typedef box_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatBox_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
//	typedef General_RatBox_Combiner<MomentaEvaluator> Combiner;
};

template <class CutType> struct Higgs_RatBox_Specification : public  General_RatBox_Specification<CutType,MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS> {};


// A class to store all the parameters to be passed between the triangle and the box
template <class T, int MUPOINTS, int CPOINTS> class box_param
{
public:
	momentum<std::complex<T> > K1, K2, K4; // The momenta of the corners
	Cmom<T> K1flatc, K2flatc; // The momenta of the basis vectors
	momentum<std::complex<T> > vec1c, vec2c;
	std::complex<T> alp1, alp2; // Basis momenta factors

	std::complex<T> gamma; // Basis factor

	std::complex<T> tri_sub_ret[(CPOINTS+1)*MUPOINTS]; // The box subtraction term
	std::complex<T> boxcoeffs[2*MUPOINTS]; // Store the box results for later use
	std::complex<T> boxpoles[2*MUPOINTS];
	std::complex<T> denfac1;

	size_t box_corner; // Stores the corner of the box that was opened to get the pentagon

	double accuracy; // Store the accuracy of our computation

	long int mcID; // Used in the mom_conf code
	size_t vertex_ref; // Used in the mom_conf code
};



//template <class BASE> class Normal_RatTri_Specification; // Declared in triangle_ratext.h

//! Class for box cut diagrams
template <class BASE,  class RatBoxSpecs/* = Normal_RatBox_Specification<BASE> */> class box_Rat : public BASE,  public RatBoxSpecs::CornerTreeStrategy /*,  public RatBoxSpecs::Combiner*/  {
	unsigned long int mcID, epID; // Used for tracking whether we have evaluated anything or not
	std::vector<int> indID;

	C coeffs[2*RatBoxSpecs::MUPOINTS];
	CHP coeffs_HP[2*RatBoxSpecs::MUPOINTS];
	CVHP coeffs_VHP[2*RatBoxSpecs::MUPOINTS];

	C D0coeff;
    CHP D0coeff_HP;
    CVHP D0coeff_VHP;

    // Stores the information we need to relate the d_+ and d_- solutions to the correct pole
    Cmom<R> K1fsmom[RatBoxSpecs::MUPOINTS], K2fsmom[RatBoxSpecs::MUPOINTS];
    momentum<std::complex<R> > K4fsmom[RatBoxSpecs::MUPOINTS];
    Cmom<RHP> K1fsmom_HP[RatBoxSpecs::MUPOINTS], K2fsmom_HP[RatBoxSpecs::MUPOINTS];
    momentum<std::complex<RHP> > K4fsmom_HP[RatBoxSpecs::MUPOINTS];
    Cmom<RVHP> K1fsmom_VHP[RatBoxSpecs::MUPOINTS], K2fsmom_VHP[RatBoxSpecs::MUPOINTS];
    momentum<std::complex<RVHP> > K4fsmom_VHP[RatBoxSpecs::MUPOINTS];
    C Tnorm[RatBoxSpecs::MUPOINTS];
    CHP Tnorm_HP[RatBoxSpecs::MUPOINTS];
    CVHP Tnorm_VHP[RatBoxSpecs::MUPOINTS];
#if BH_USE_GMP
	CGMP coeffs_GMP[2*RatBoxSpecs::MUPOINTS];
    CGMP D0coeff_GMP;
    Cmom<RGMP> K1fsmom_GMP[RatBoxSpecs::MUPOINTS], K2fsmom_GMP[RatBoxSpecs::MUPOINTS];
    momentum<std::complex<RGMP> > K4fsmom_GMP[RatBoxSpecs::MUPOINTS];
    CGMP Tnorm_GMP[RatBoxSpecs::MUPOINTS];
#endif
protected:
    //This vector is used for storing the indexes that get passed to the cut trees
	std::vector<int> indlst[4];

	eval_param<R>* _ep[4];
	eval_param<RHP>* _ep_HP[4];
	eval_param<RVHP>* _ep_VHP[4];

	std::vector<size_t> _cut_mass; // Store the unique mass_params of the cut legs


	bool _k1massive, _k2massive; // Keeps track of wether the K1 and K2 legs are massive or not

#if BH_USE_GMP
	eval_param<RGMP>* _ep_GMP[4];
#endif

public:
	mass_param_coll* _leg_masses; // Store the mass_params for each cut leg
	template <class T> box_Rat(T& t);
	template <class T> box_Rat(T* t);
	virtual ~box_Rat();

	//! Evaluation of box coefficients
	/** \param mc momentum configuration to be used \param ind indices of the momenta to be used \returns box coefficient*/
	virtual C eval(mom_conf& mc,const std::vector<int>& ind);
	virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind);
	virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind);

	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);

	// Computes the box subtraction terms for the triangle
	template <class T> void get_sub_terms(momentum_configuration<T>& mc,const std::vector<int>& ind, box_param<T,RatBoxSpecs::MUTRIPOINTS,RatBoxSpecs::CPOINTS>& bp);
	template <class T> void get_sub_terms(const eval_param<T>& ep, box_param<T,RatBoxSpecs::MUTRIPOINTS,RatBoxSpecs::CPOINTS>& bp);

	inline void set_eval(long int mcset, const std::vector<int>& indset){mcID=mcset;indID=indset;};
	inline void set_eval(unsigned long int mcset){epID=mcset;};
	inline bool is_eval(long int mcset, const std::vector<int>& indset){if((mcset==mcID)&&(indset==indID)){return true;} return false;};
	inline bool is_eval(unsigned long int mcset){if(mcset==epID){return true;} return false;};

	// Add a cutD to the list
	void add_mass(const cutD& pcd);
	void add_mass(int m1,int m2,int m3,int m4);
	bool is_k1massive() const { return _k1massive;}
	bool is_k2massive() const { return _k2massive;}

#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
	virtual CGMP eval(const eval_param<RGMP>& ep);
#endif

   	typedef typename RatBoxSpecs::MomentaEvaluator MomentaEvaluatorType;
   	typedef typename daughter_info<box_Rat<BASE,RatBoxSpecs> >::daughter_type* dau_type_p;
protected:
	// Computes the box coefficients
	virtual void get_coeffs(momentum_configuration<R>& mc, const std::vector<int>& ind, double& accuracy , const long int original_mcID, const size_t vertex_ref){get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};
	virtual void get_coeffs(momentum_configuration<RHP>& mc, const std::vector<int>& ind, double& accuracy , const long int original_mcID, const size_t vertex_ref){get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};
	virtual void get_coeffs(momentum_configuration<RVHP>& mc, const std::vector<int>& ind, double& accuracy , const long int original_mcID, const size_t vertex_ref){ get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};
	template <class T> void get_coeffs_fn(momentum_configuration<T>& mc, const std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref);

	virtual void get_coeffs(const eval_param<R>& ep,double& accuracy){get_coeffs_fn(ep,accuracy);};
	virtual void get_coeffs(const eval_param<RHP>& ep,double& accuracy){get_coeffs_fn(ep,accuracy);};
	virtual void get_coeffs(const eval_param<RVHP>& ep,double& accuracy){ get_coeffs_fn(ep,accuracy);};
	template <class T> void get_coeffs_fn(const eval_param<T>& ep,double& accuracy);
#if BH_USE_GMP
	virtual void get_coeffs(momentum_configuration<RGMP>& mc, const std::vector<int>& ind, double& accuracy , const long int original_mcID, const size_t vertex_ref){ get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};
	virtual void get_coeffs(const eval_param<RGMP>& ep,double& accuracy){ get_coeffs_fn(ep,accuracy);};
#endif

protected:
	void set_coeff(size_t i, size_t m0pole, const C& value){coeffs[i+2*m0pole]=value;};
	void set_coeff(size_t i, size_t m0pole, const CHP& value){coeffs_HP[i+2*m0pole]=value;};
	void set_coeff(size_t i, size_t m0pole, const CVHP& value){coeffs_VHP[i+2*m0pole]=value;};

	void get_coeff(size_t m0pole, C& value1, C& value2){value1=coeffs[2*m0pole]; value2=coeffs[1+2*m0pole];};
	void get_coeff(size_t m0pole, CHP& value1, CHP& value2){value1=coeffs_HP[2*m0pole]; value2=coeffs_HP[1+2*m0pole];};
	void get_coeff(size_t m0pole, CVHP& value1, CVHP& value2){value1=coeffs_VHP[2*m0pole]; value2=coeffs_VHP[1+2*m0pole];};

	void set_Kif_coeffs(size_t ratpoints, const Cmom<R>& K1fmom, const Cmom<R>& K2fmom, const momentum<std::complex<R> >& K4fmom){K1fsmom[ratpoints]=K1fmom;K2fsmom[ratpoints]=K2fmom;K4fsmom[ratpoints]=K4fmom;};
	void set_Kif_coeffs(size_t ratpoints, const Cmom<RHP>& K1fmom, const Cmom<RHP>& K2fmom, const momentum<std::complex<RHP> >& K4fmom){K1fsmom_HP[ratpoints]=K1fmom;K2fsmom_HP[ratpoints]=K2fmom;K4fsmom_HP[ratpoints]=K4fmom;};
	void set_Kif_coeffs(size_t ratpoints, const Cmom<RVHP>& K1fmom, const Cmom<RVHP>& K2fmom, const momentum<std::complex<RVHP> >& K4fmom){K1fsmom_VHP[ratpoints]=K1fmom;K2fsmom_VHP[ratpoints]=K2fmom;K4fsmom_VHP[ratpoints]=K4fmom;};

	void get_Kif_coeffs(size_t ratpoints, Cmom<R>& K1f, Cmom<R>& K2f, momentum<std::complex<R> >& K4f){K1f=K1fsmom[ratpoints];K2f=K2fsmom[ratpoints];K4f=K4fsmom[ratpoints];};
	void get_Kif_coeffs(size_t ratpoints, Cmom<RHP>& K1f, Cmom<RHP>& K2f, momentum<std::complex<RHP> >& K4f){K1f=K1fsmom_HP[ratpoints];K2f=K2fsmom_HP[ratpoints];K4f=K4fsmom_HP[ratpoints];};
	void get_Kif_coeffs(size_t ratpoints, Cmom<RVHP>& K1f, Cmom<RVHP>& K2f, momentum<std::complex<RVHP> >& K4f){K1f=K1fsmom_VHP[ratpoints];K2f=K2fsmom_VHP[ratpoints];K4f=K4fsmom_VHP[ratpoints];};

	void set_denfac(int m0pole, const C& value1){Tnorm[m0pole]=value1;};
	void set_denfac(int m0pole, const CHP& value1){Tnorm_HP[m0pole]=value1;};
	void set_denfac(int m0pole, const CVHP& value1){Tnorm_VHP[m0pole]=value1;};

	void get_denfac(int m0pole, C& value1){value1=Tnorm[m0pole];};
	void get_denfac(int m0pole, CHP& value1){value1=Tnorm_HP[m0pole];};
	void get_denfac(int m0pole, CVHP& value1){value1=Tnorm_VHP[m0pole];};

	void set_D0coeff(const C& value){D0coeff=value;};
	void set_D0coeff(const CHP& value){D0coeff_HP=value;D0coeff=to_double(value);};
	void set_D0coeff(const CVHP& value){D0coeff_VHP=value;D0coeff_HP=to_HP(value);D0coeff=to_double(value);};

	void get_D0coeff(C& value){value=D0coeff;};
	void get_D0coeff(CHP& value){value=D0coeff_HP;};
	void get_D0coeff(CVHP& value){value=D0coeff_VHP;};

	void get_ep(eval_param<R>**& ep){ep=_ep;}
	void get_ep(eval_param<RHP>**& ep){ep= _ep_HP;}
	void get_ep(eval_param<RVHP>**& ep){ep= _ep_VHP;}

#if BH_USE_GMP
	void set_coeff(size_t i, size_t m0pole, const CGMP& value){coeffs_GMP[i+2*m0pole]=value;};
	void get_coeff(size_t m0pole, CGMP& value1, CGMP& value2){value1=coeffs_GMP[2*m0pole]; value2=coeffs_GMP[1+2*m0pole];};
	void set_Kif_coeffs(size_t ratpoints, const Cmom<RGMP>& K1fmom, const Cmom<RGMP>& K2fmom, const momentum<std::complex<RGMP> >& K4fmom){K1fsmom_GMP[ratpoints]=K1fmom;K2fsmom_GMP[ratpoints]=K2fmom;K4fsmom_GMP[ratpoints]=K4fmom;};
	void get_Kif_coeffs(size_t ratpoints, Cmom<RGMP>& K1f, Cmom<RGMP>& K2f, momentum<std::complex<RGMP> >& K4f){K1f=K1fsmom_GMP[ratpoints];K2f=K2fsmom_GMP[ratpoints];K4f=K4fsmom_GMP[ratpoints];};
	void set_denfac(int m0pole, const CGMP& value1){Tnorm_GMP[m0pole]=value1;};
	void get_denfac(int m0pole, CGMP& value1){value1=Tnorm_GMP[m0pole];};
	void set_D0coeff(const CGMP& value){D0coeff_GMP=value;D0coeff_HP=to_HP(value);D0coeff=to_double(value);};
	void get_D0coeff(CGMP& value){value=D0coeff_GMP;};
	void get_ep(eval_param<RGMP>**& ep){if (_ep_GMP[0]->ref()->E().real().get_precision() < RGMP::get_current_precision() ){ init(); }; ep= _ep_GMP;}
#endif


	void init();
};



}
}
#endif /* BOX_RATEXT_H_ */
