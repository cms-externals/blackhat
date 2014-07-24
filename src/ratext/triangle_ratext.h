/*
 * triangle_ratext.h
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef TRIANGLE_RATEXT_H_
#define TRIANGLE_RATEXT_H_

#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "ratext/rat_ext.h"
#include "eval_param.h"
#include "momenta_evaluator_ratext.h"
#include "BH_debug.h"
namespace BH {

namespace ratext {

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class box_rat_eval_points; // Declared in box_ratext.h
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class bub_rat_eval_points; // Declared in bubble_ratext.h

// The rational integral function pointer typedefs
typedef C (*Rat_Int_Fac_Fn_Ptr)(const C& s1, const C& s2, const C& s3);
typedef CHP (*Rat_Int_Fac_Fn_Ptr_HP)(const CHP& s1, const CHP& s2, const CHP& s3);
typedef CVHP (*Rat_Int_Fac_Fn_Ptr_VHP)(const CVHP& s1, const CVHP& s2, const CVHP& s3);	

#if BH_USE_GMP
typedef CGMP (*Rat_Int_Fac_Fn_Ptr_GMP)(const CGMP& s1, const CGMP& s2, const CGMP& s3);
#endif

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class tri_rat_eval_points : public MomentumEvaluator {
	static C _circpos[CPOINTS_T];
	static CHP _circpos_HP[CPOINTS_T];
	static CVHP _circpos_VHP[CPOINTS_T];

	static C _circpos_matrix[CPOINTS_T*CPOINTS_T];
	static CHP _circpos_matrix_HP[CPOINTS_T*CPOINTS_T];
	static CVHP _circpos_matrix_VHP[CPOINTS_T*CPOINTS_T];

	// For now we only store the results for at most a mu^4 triangle, this is good for a max power of l^5
	static Rat_Int_Fac_Fn_Ptr _rat_int[4];
	static Rat_Int_Fac_Fn_Ptr_HP _rat_int_HP[4];
	static Rat_Int_Fac_Fn_Ptr_VHP _rat_int_VHP[4];
#if BH_USE_GMP
	static CGMP _circpos_GMP[CPOINTS_T];
	static CGMP _circpos_matrix_GMP[CPOINTS_T*CPOINTS_T];
	static Rat_Int_Fac_Fn_Ptr_GMP _rat_int_GMP[4];
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
	tri_rat_eval_points(){init();};
    virtual ~tri_rat_eval_points(){}

    // Get the points on the circle for evaluating the t extraction in the bubble
    static void get_t_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos;}
    static void get_t_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_HP;}
    static void get_t_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_VHP;}

    // Get the points on the circle for evaluating the t extraction in the bubble
    static void get_t_matrix_eval_pts(const C*& ptr_circpos){ptr_circpos=_circpos_matrix;}
    static void get_t_matrix_eval_pts(const CHP*& ptr_circpos){ptr_circpos=_circpos_matrix_HP;}
    static void get_t_matrix_eval_pts(const CVHP*& ptr_circpos){ptr_circpos=_circpos_matrix_VHP;}

    // Get the rational integral coefficient
    static const C get_rat_integral(int element, const C& s1, const C& s2, const C& s3);
    static const CHP get_rat_integral(int element, const CHP& s1, const CHP& s2, const CHP& s3);
    static const CVHP get_rat_integral(int element, const CVHP& s1, const CVHP& s2, const CVHP& s3);

#if BH_USE_GMP
    static void get_t_eval_pts(const CGMP*& ptr_circpos){if (RGMP::get_current_precision()> s_GMP_precision){ BH_DEBUG_MESSAGE2("recomputing ratext triangle points for precision ",RGMP::get_current_precision());get_points(); }; ptr_circpos=_circpos_GMP;}
    static void get_t_matrix_eval_pts(const CGMP*& ptr_circpos){if (RGMP::get_current_precision()> s_GMP_precision){ BH_DEBUG_MESSAGE2("recomputing ratext triangle points for precision ",RGMP::get_current_precision());get_points(); }; ptr_circpos=_circpos_matrix_GMP;}
    static const CGMP get_rat_integral(int element, const CGMP& s1, const CGMP& s2, const CGMP& s3);
#endif
    static box_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> BoxPoints;
    static bub_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> BubPoints;

	template <class T,class cutDbase> void compute_momenta(const eval_param<T>& ep,Cmom<T> (&momenta)[3],Cmom<T> (&momenta_m)[3],const cutDbase& cb);

    static void init();

private:
    static const bool get_points();
};

// The functions for the triangle rat integrals 
template <class T> std::complex<T> Rat_Int_Tri_1(const std::complex<T>& s1, const std::complex<T>& s2, const std::complex<T>& s3){
	return std::complex<T>(-1,0)/T(2);
}
	
template <class T> std::complex<T> Rat_Int_Tri_2(const std::complex<T>& s1, const std::complex<T>& s2, const std::complex<T>& s3){
	return T(8)*(s1+s2+s3)/T(-24);
}

template <class T> std::complex<T> Rat_Int_Tri_3(const std::complex<T>& s1, const std::complex<T>& s2, const std::complex<T>& s3){
	return (s1*s1+s2*s2+s3*s3+s1*s2+s1*s3+s2*s3)/T(-180);
}

template <class T> std::complex<T> Rat_Int_Tri_4(const std::complex<T>& s1, const std::complex<T>& s2, const std::complex<T>& s3){
	return std::complex<T>(0,0);
}
	
//template <class MomentumEvaluator> class General_RatTri_Combiner {
//public:
//	template <class T> std::complex<T> combine(
//			std::complex<T>* ampl,
//			std::complex<T>* ampl_err);
//};
//
//template <class MomentumEvaluator> class Normal_RatTri_Combiner {
//public:
//	template <class T> std::complex<T> combine(
//			std::complex<T>* ampl,
//			std::complex<T>* ampl_err);
//};


template <class MomentumEvaluator,class CutType> class Normal_RatTri_Corner_Tree_Strategy : public MomentumEvaluator {
public:
	template <class T> void fill_trees(
			const eval_param<T>&,
			std::vector<std::complex<T> > (&trees_result)[3],
			std::vector<size_t>& cut_mass,
			CutType& self);
};


template <class CutType> class Normal_RatTri_Specification {
public:
	enum { MUPOINTS=3 }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=2 }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=2 }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=7 }; // Number of C coefficients
	enum { CBUBPOINTS=3}; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=3}; // Number of points in y used to evaluate for the bubble

	typedef tri_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatTri_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
//	typedef Normal_RatTri_Combiner<MomentaEvaluator> Combiner;
};

template <class CutType, int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class General_RatTri_Specification {
public:
	enum { MUPOINTS=MUPOINTS_T }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=MUTRIPOINTS_T }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=MUBUBPOINTS_T }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=CPOINTS_T }; // Number of C coefficients
	enum { CBUBPOINTS=CBUBPOINTS_T }; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=YPOINTS_T}; // Number of points in y used to evaluate for the bubble

	typedef tri_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatTri_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
//	typedef General_RatTri_Combiner<MomentaEvaluator> Combiner;
};

template <class CutType> class Higgs_RatTri_Specification : public  General_RatTri_Specification<CutType,MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS> {};


// A class to store all the parameters to be passed between the bubble and the triangle
template <class T, int MUPOINTS_T, int YPOINTS_T> class triangle_param
{
public:
	momentum<std::complex<T> > K1; // The momenta of the corners
	Cmom<T> K1flatbc, chic; // The momenta of the basis vectors
	std::complex<T> S1, gammab; // Basis factors

	size_t tri_corner; // Stores the corner of the triangle that was opened to get the pentagon
	std::complex<T> bub_acc_pieces[YPOINTS_T]; // Stores the pieces required from the triangle for the error estimate

	std::complex<T> tri_sub_pieces[MUPOINTS_T-1]; // Stores the returned results from the triangle subtraction terms
//	std::complex<T> tri_sub_pieces_p[YPOINTS_T*MUPOINTS_T]; //Used for debugging
//	std::complex<T> tri_sub_pieces_m[YPOINTS_T*MUPOINTS_T]; //Used for debugging
//	std::complex<T> tri_sub_pieces_0p[YPOINTS_T*MUPOINTS_T]; //Used for debugging
//	std::complex<T> tri_sub_pieces_0m[YPOINTS_T*MUPOINTS_T]; //Used for debugging

	double accuracy; // Store the accuracy of our computation

	long int mcID; // Used in the mom_conf code
	size_t vertex_ref; // Used in the mom_conf code
};



//! Class for Darren's triangle cut diagrams
template <class BASE,  class RatTriSpecs/* = Normal_RatTri_Specification<BASE>*/ > class triangle_Rat : public BASE,  public RatTriSpecs::CornerTreeStrategy /*,  public RatTriSpecs::Combiner*/   {
	long int mcID;
	unsigned long int epID;
	std::vector<int> indID;

	C C0coeff;
	CHP C0coeff_HP;
	CVHP C0coeff_VHP;

	 // Stores the information we need to relate the c coeffs with different opened bubble corners
	momentum<std::complex<R> > K1fsmom, K2fsmom;
	momentum<std::complex<RHP> > K1fsmom_HP, K2fsmom_HP;
	momentum<std::complex<RVHP> > K1fsmom_VHP, K2fsmom_VHP;
	C coeffs[RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS];
	CHP coeffs_HP[RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS];
	CVHP coeffs_VHP[RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS];
	C gamma,gammap;
	CHP gamma_HP,gammap_HP;
	CVHP gamma_VHP,gammap_VHP;

	std::vector<C> coeffkeep[RatTriSpecs::MUTRIPOINTS], polekeep[RatTriSpecs::MUTRIPOINTS], denpole;
	std::vector<CHP> coeffkeep_HP[RatTriSpecs::MUTRIPOINTS], polekeep_HP[RatTriSpecs::MUTRIPOINTS], denpole_HP;
	std::vector<CVHP> coeffkeep_VHP[RatTriSpecs::MUTRIPOINTS], polekeep_VHP[RatTriSpecs::MUTRIPOINTS], denpole_VHP;
#if BH_USE_GMP
	CGMP C0coeff_GMP;
	momentum<std::complex<RGMP> > K1fsmom_GMP, K2fsmom_GMP;
	CGMP coeffs_GMP[RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS];
	CGMP gamma_GMP,gammap_GMP;
	std::vector<CGMP> coeffkeep_GMP[RatTriSpecs::MUTRIPOINTS], polekeep_GMP[RatTriSpecs::MUTRIPOINTS], denpole_GMP;

	#endif
protected:
	//This vector is used for storing the indexes that get passed to the cut trees
	std::vector<int> indlst[3];

	eval_param<R>* _ep[3];
	eval_param<RHP>* _ep_HP[3];
	eval_param<RVHP>* _ep_VHP[3];
#if BH_USE_GMP
	eval_param<RGMP>* _ep_GMP[3];
#endif
	std::vector<size_t> _cut_mass; // Store the unique mass_params of the cut legs


	bool _k1massive, _k2massive, _k3massive; // Keeps track of weather the K1, K2 and K3 legs are massive or not

public:
	mass_param_coll* _leg_masses; // Store the mass_params for each cut leg
	template <class T> triangle_Rat(T& t);
	template <class T> triangle_Rat(T* t);
	virtual ~triangle_Rat();

	virtual C eval(mom_conf&,const std::vector<int>&);
	virtual CHP eval(mom_conf_HP&,const std::vector<int>&);
	virtual CVHP eval(mom_conf_VHP&,const std::vector<int>&);

	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);

	//These versions should never be called
	virtual void get_sub_terms(mom_conf& mc, const std::vector<int>& ind, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;
	virtual void get_sub_terms(mom_conf_HP& mc, const std::vector<int>& ind, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;
	virtual void get_sub_terms(mom_conf_VHP& mc, const std::vector<int>& ind, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;

	virtual void get_sub_terms(const eval_param<R>& ep, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;
	virtual void get_sub_terms(const eval_param<RHP>& ep, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;
	virtual void get_sub_terms(const eval_param<RVHP>& ep, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;
#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>&,const std::vector<int>&);
	virtual CGMP eval(const eval_param<RGMP>& ep);
	virtual void get_sub_terms(momentum_configuration<RGMP>& mc, const std::vector<int>& ind, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;
	virtual void get_sub_terms(const eval_param<RGMP>& ep, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)=0;
#endif
	//Used for setting and checking weather this has already been computed
	inline void set_eval(long int mcset, const std::vector<int>& indset){mcID=mcset;indID=indset;};
	inline void set_eval(unsigned long int mcset){epID=mcset;};
	inline bool is_eval(long int mcset, const std::vector<int>& indset){if((mcset==mcID)&&(indset==indID)){return true;} return false;};
	inline bool is_eval(unsigned long int mcset){if(mcset==epID){return true;} return false;};

	// Add a cutD to the list
	void add_mass(const cutD& pcd);
	void add_mass(int m1, int m2,int m3);

	bool is_k1massive() const { return _k1massive;}
	bool is_k2massive() const { return _k2massive;}
	bool is_k3massive() const { return _k3massive;}

	// Some storage related items need to be set up after all other initialisation is done
	void final_initialisation(); // This is called after all the parent and daughter diagrams are set up but before eval is called

   	typedef typename RatTriSpecs::MomentaEvaluator MomentaEvaluatorType;
   	typedef typename daughter_info<triangle_Rat<BASE,RatTriSpecs> >::daughter_type* dau_type_p;
protected:
	// Computes the triangle coefficients
	template <class T> void get_coeffs_fn(momentum_configuration<T>& mc, const  std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref);
	void get_coeffs(momentum_configuration<R>& mc, const  std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref){get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};
	void get_coeffs(momentum_configuration<RHP>& mc, const  std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref){get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};
	void get_coeffs(momentum_configuration<RVHP>& mc, const  std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref){get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};

	template <class T> void get_coeffs_fn(const eval_param<T>& ep, double& accuracy);
	void get_coeffs(const eval_param<R>& ep, double& accuracy){get_coeffs_fn(ep,accuracy);};
	void get_coeffs(const eval_param<RHP>& ep, double& accuracy){get_coeffs_fn(ep,accuracy);};
	void get_coeffs(const eval_param<RVHP>& ep, double& accuracy){get_coeffs_fn(ep,accuracy);};

	// Gets the previously computed boxes for the subtraction of the boxes in the bubble
	inline void get_boxes(std::vector<C>*& coeffs, std::vector<C>*& poles, std::vector<C>*& denpoles){coeffs=coeffkeep;poles=polekeep;denpoles=&denpole;};
	inline void get_boxes(std::vector<CHP>*& coeffs, std::vector<CHP>*& poles, std::vector<CHP>*& denpoles){coeffs=coeffkeep_HP;poles=polekeep_HP;denpoles=&denpole_HP;};
	inline void get_boxes(std::vector<CVHP>*& coeffs, std::vector<CVHP>*& poles, std::vector<CVHP>*& denpoles){coeffs=coeffkeep_VHP;poles=polekeep_VHP;denpoles=&denpole_VHP;};

	//Used for setting the full coefficient
	void set_C0coeff(C value){C0coeff=value;};
	void set_C0coeff(CHP value){C0coeff_HP=value;C0coeff=to_double(value);};
	void set_C0coeff(CVHP value){C0coeff_VHP=value;C0coeff_HP=to_HP(value);C0coeff=to_double(value);};

	void get_C0coeff(C& value){value=C0coeff;};
	void get_C0coeff(CHP& value){value=C0coeff_HP;};
	void get_C0coeff(CVHP& value){value=C0coeff_VHP;};

	inline void boxcoeff_add(size_t in, C* coeff, C* pole, C& denval){for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
																								(coeffkeep[m0pole])[2*in]=coeff[2*m0pole];
																								(coeffkeep[m0pole])[1+2*in]=coeff[1+2*m0pole];
																								(polekeep[m0pole])[2*in]=pole[2*m0pole];
																								(polekeep[m0pole])[1+2*in]=pole[1+2*m0pole];
																							}
																							denpole[in]=denval;};
	inline void boxcoeff_add(size_t in, CHP* coeff, CHP* pole , CHP& denval){for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
																								(coeffkeep_HP[m0pole])[2*in]=coeff[2*m0pole];
																								(coeffkeep_HP[m0pole])[1+2*in]=coeff[1+2*m0pole];
																								(polekeep_HP[m0pole])[2*in]=pole[2*m0pole];
																								(polekeep_HP[m0pole])[1+2*in]=pole[1+2*m0pole];
																							}
																							denpole_HP[in]=denval;};
	inline void boxcoeff_add(size_t in, CVHP* coeff , CVHP* pole, CVHP& denval){for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
																								(coeffkeep_VHP[m0pole])[2*in]=coeff[2*m0pole];
																								(coeffkeep_VHP[m0pole])[1+2*in]=coeff[1+2*m0pole];
																								(polekeep_VHP[m0pole])[2*in]=pole[2*m0pole];
																								(polekeep_VHP[m0pole])[1+2*in]=pole[1+2*m0pole];
																							}
																							denpole_VHP[in]=denval;};

	void set_tri_param_basis_vectors(const momentum<std::complex<R> >& K1fmom, const momentum<std::complex<R> >& K2fmom){K1fsmom=K1fmom;K2fsmom=K2fmom;};
	void set_tri_param_basis_vectors(const momentum<std::complex<RHP> >& K1fmom, const momentum<std::complex<RHP> >& K2fmom){K1fsmom_HP=K1fmom;K2fsmom_HP=K2fmom;};
	void set_tri_param_basis_vectors(const momentum<std::complex<RVHP> >& K1fmom, const momentum<std::complex<RVHP> >& K2fmom){K1fsmom_VHP=K1fmom;K2fsmom_VHP=K2fmom;};

	void set_tri_param_gamma(const C& gam,const C& gamp){gamma=gam;gammap=gamp;};
	void set_tri_param_gamma(const CHP& gam,const CHP& gamp){gamma_HP=gam;gammap_HP=gamp;};
	void set_tri_param_gamma(const CVHP& gam,const CVHP& gamp){gamma_VHP=gam;gammap_VHP=gamp;};

	inline void coeffkeep_add(C *coeff){for(int icirc=0; icirc<RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS; icirc++){
															coeffs[icirc]=coeff[icirc];
														}};
	inline void coeffkeep_add(CHP *coeff){for(int icirc=0; icirc<RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS; icirc++){
															coeffs_HP[icirc]=coeff[icirc];
														}};
	inline void coeffkeep_add(CVHP *coeff){for(int icirc=0; icirc<RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS; icirc++){
															coeffs_VHP[icirc]=coeff[icirc];
														}};

	//Used for getting the vectors used to construct the triangle param this was computed using
	void get_tri_param_basis_vectors(momentum<std::complex<R> >& K1f, momentum<std::complex<R> >& K2f){K1f=K1fsmom;K2f=K2fsmom;};
	void get_tri_param_basis_vectors(momentum<std::complex<RHP> >& K1f, momentum<std::complex<RHP> >& K2f){K1f=K1fsmom_HP;K2f=K2fsmom_HP;};
	void get_tri_param_basis_vectors(momentum<std::complex<RVHP> >& K1f, momentum<std::complex<RVHP> >& K2f){K1f=K1fsmom_VHP;K2f=K2fsmom_VHP;};

	void get_tri_param_gamma(C& gam,C& gamp){gam=gamma;gamp=gammap;};
	void get_tri_param_gamma(CHP& gam,CHP& gamp){gam=gamma_HP;gamp=gammap_HP;};
	void get_tri_param_gamma(CVHP& gam,CVHP& gamp){gam=gamma_VHP;gamp=gammap_VHP;};

	inline void coeffkeep_get(C*& coeff){coeff=coeffs;};
	inline void coeffkeep_get(CHP*& coeff){coeff=coeffs_HP;};
	inline void coeffkeep_get(CVHP*& coeff){coeff=coeffs_VHP;};

	void get_ep(eval_param<R>**& ep){ep=_ep;}
	void get_ep(eval_param<RHP>**& ep){ep= _ep_HP;}
	void get_ep(eval_param<RVHP>**& ep){ep= _ep_VHP;}

#if BH_USE_GMP
	void get_coeffs(momentum_configuration<RGMP>& mc, const  std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref){get_coeffs_fn(mc,ind,accuracy,original_mcID,vertex_ref);};
	void get_coeffs(const eval_param<RGMP>& ep, double& accuracy){get_coeffs_fn(ep,accuracy);};
	inline void get_boxes(std::vector<CGMP>*& coeffs, std::vector<CGMP>*& poles, std::vector<CGMP>*& denpoles){coeffs=coeffkeep_GMP;poles=polekeep_GMP;denpoles=&denpole_GMP;};
	void set_C0coeff(CGMP value){C0coeff_GMP=value;C0coeff_HP=to_HP(value);C0coeff=to_double(value);};
	void get_C0coeff(CGMP& value){value=C0coeff_GMP;};
	inline void boxcoeff_add(size_t in, CGMP* coeff , CGMP* pole, CGMP& denval){for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
																								(coeffkeep_GMP[m0pole])[2*in]=coeff[2*m0pole];
																								(coeffkeep_GMP[m0pole])[1+2*in]=coeff[1+2*m0pole];
																								(polekeep_GMP[m0pole])[2*in]=pole[2*m0pole];
																								(polekeep_GMP[m0pole])[1+2*in]=pole[1+2*m0pole];
																							}
																							denpole_GMP[in]=denval;};

	void set_tri_param_basis_vectors(const momentum<std::complex<RGMP> >& K1fmom, const momentum<std::complex<RGMP> >& K2fmom){K1fsmom_GMP=K1fmom;K2fsmom_GMP=K2fmom;};
	void set_tri_param_gamma(const CGMP& gam,const CGMP& gamp){gamma_GMP=gam;gammap_GMP=gamp;};
	inline void coeffkeep_add(CGMP *coeff){for(int icirc=0; icirc<RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS; icirc++){
															coeffs_GMP[icirc]=coeff[icirc];
														}};

	//Used for getting the vectors used to construct the triangle param this was computed using
	void get_tri_param_basis_vectors(momentum<std::complex<RGMP> >& K1f, momentum<std::complex<RGMP> >& K2f){K1f=K1fsmom_GMP;K2f=K2fsmom_GMP;};
	void get_tri_param_gamma(CGMP& gam,CGMP& gamp){gam=gamma_GMP;gamp=gammap_GMP;};
	inline void coeffkeep_get(CGMP*& coeff){coeff=coeffs_GMP;};
	void get_ep(eval_param<RGMP>**& ep){if (_ep_GMP[0]->ref()->E().real().get_precision() < RGMP::get_current_precision() ){ init(); }; ep= _ep_GMP;}
#endif


	void init(); // We need to pass the cutD used for construction
};



template <class BASE, class RatTriSpecs/* = Normal_RatTri_Specification<BASE>*/ > class triangle_Rat_3mass : public triangle_Rat<BASE,RatTriSpecs> {
public:
	template <class T> triangle_Rat_3mass(T& t) : triangle_Rat<BASE,RatTriSpecs>(t) {};
	template <class T> triangle_Rat_3mass(T* t) : triangle_Rat<BASE,RatTriSpecs>(t) {};
	virtual ~triangle_Rat_3mass(){};

	virtual void get_sub_terms(mom_conf& mc, const std::vector<int>& ind, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(mc,ind,tp);};
	virtual void get_sub_terms(mom_conf_HP& mc, const std::vector<int>& ind, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(mc,ind,tp);};
	virtual void get_sub_terms(mom_conf_VHP& mc, const std::vector<int>& ind, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(mc,ind,tp);};

	virtual void get_sub_terms(const eval_param<R>& ep, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
	virtual void get_sub_terms(const eval_param<RHP>& ep, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
	virtual void get_sub_terms(const eval_param<RVHP>& ep, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
#if BH_USE_GMP
	virtual void get_sub_terms(momentum_configuration<RGMP>& mc, const std::vector<int>& ind, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(mc,ind,tp);};
	virtual void get_sub_terms(const eval_param<RGMP>& ep, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
#endif
private:
	template <class T> void get_sub_terms_work(momentum_configuration<T>& mc, const std::vector<int>& ind, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp);
	template <class T> void get_sub_terms_work(const eval_param<T>& ep, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp);
};

template <class BASE, class RatTriSpecs/* = Normal_RatTri_Specification<BASE>*/ > class triangle_Rat_plusminus : public triangle_Rat<BASE,RatTriSpecs> {
public:
	template <class T> triangle_Rat_plusminus(T& t) : triangle_Rat<BASE,RatTriSpecs>(t) {};
	template <class T> triangle_Rat_plusminus(T* t) : triangle_Rat<BASE,RatTriSpecs>(t) {};
	virtual ~triangle_Rat_plusminus(){};

	virtual void get_sub_terms(mom_conf& mc, const std::vector<int>& ind, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(mc,ind,tp);};
	virtual void get_sub_terms(mom_conf_HP& mc, const std::vector<int>& ind, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(mc,ind,tp);};
	virtual void get_sub_terms(mom_conf_VHP& mc, const std::vector<int>& ind, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(mc,ind,tp);};

	virtual void get_sub_terms(const eval_param<R>& ep, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
	virtual void get_sub_terms(const eval_param<RHP>& ep, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
	virtual void get_sub_terms(const eval_param<RVHP>& ep, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
#if BH_USE_GMP
	virtual void get_sub_terms(momentum_configuration<RGMP>& mc, const std::vector<int>& ind, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){get_sub_terms_work(mc,ind,tp);};
	virtual void get_sub_terms(const eval_param<RGMP>& ep, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){ get_sub_terms_work(ep,tp);};
#endif
	private:
	template <class T> void get_sub_terms_work(momentum_configuration<T>& mc, const std::vector<int>& ind, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp);
	template <class T> void get_sub_terms_work(const eval_param<T>& ep, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp);
};

template <class BASE, class RatTriSpecs/* = Normal_RatTri_Specification<BASE>*/ > class triangle_Rat_zero : public triangle_Rat<BASE,RatTriSpecs> {
public:
	template <class T> triangle_Rat_zero(T& t) : triangle_Rat<BASE,RatTriSpecs>(t) {};
	template <class T> triangle_Rat_zero(T* t) : triangle_Rat<BASE,RatTriSpecs>(t) {};
	virtual ~triangle_Rat_zero(){};

	virtual C eval(mom_conf&,const std::vector<int>&){return C(0.,0.);};
	virtual CHP eval(mom_conf_HP&,const std::vector<int>&){return CHP(0.,0.);};
	virtual CVHP eval(mom_conf_VHP&,const std::vector<int>&){return CVHP(0.,0.);};

	virtual C eval(const eval_param<R>& ep){return C(0.,0.);};
	virtual CHP eval(const eval_param<RHP>& ep){return CHP(0.,0.);};
	virtual CVHP eval(const eval_param<RVHP>& ep){return CVHP(0.,0.);};

	virtual void get_sub_terms(mom_conf& mc, const std::vector<int>& ind, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};
	virtual void get_sub_terms(mom_conf_HP& mc, const std::vector<int>& ind, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};
	virtual void get_sub_terms(mom_conf_VHP& mc, const std::vector<int>& ind, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};

	virtual void get_sub_terms(const eval_param<R>& ep, triangle_param<R,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};
	virtual void get_sub_terms(const eval_param<RHP>& ep, triangle_param<RHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};
	virtual void get_sub_terms(const eval_param<RVHP>& ep, triangle_param<RVHP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};

#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>&,const std::vector<int>&){return CGMP(0.,0.);};
	virtual CGMP eval(const eval_param<RGMP>& ep){return CGMP(0.,0.);};
	virtual void get_sub_terms(momentum_configuration<RGMP>& mc, const std::vector<int>& ind, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};
	virtual void get_sub_terms(const eval_param<RGMP>& ep, triangle_param<RGMP,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp){};
#endif
};


}
}

#endif /* TRIANGLE_RATEXT_H_ */
