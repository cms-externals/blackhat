/*
 * pentagon_ratext.h
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef PENTAGON_RATEXT_H_
#define PENTAGON_RATEXT_H_

#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "eval_param.h"
#include "momenta_evaluator_ratext.h"
#include "box_ratext.h"
#include "coeff_param.h"

namespace BH {

class cutD;

namespace ratext {

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> class pentagon_rat_eval_points : public MomentumEvaluator {
public:
	enum { MUPOINTS = MUPOINTS_T};
	enum { MUTRIPOINTS = MUTRIPOINTS_T};
	enum { MUBUBPOINTS = MUBUBPOINTS_T};
	enum { CPOINTS = CPOINTS_T};
	enum { CBUBPOINTS = CBUBPOINTS_T};
	enum { YPOINTS = YPOINTS_T};
	pentagon_rat_eval_points(){init();};
    virtual ~pentagon_rat_eval_points(){}

    static box_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> BoxPoints;

	template <class T,class cutDbase> void compute_momenta(const eval_param<T>& ep,Cmom<T> (&momenta)[5],Cmom<T> (&momenta_m)[5], std::complex<T>& pm0sol,const cutDbase& cb);
private:
    void init();
};

template <class MomentumEvaluator,class CutType> class Normal_RatPent_Corner_Tree_Strategy : public MomentumEvaluator {
public:
	template <class T> void fill_trees(
			const eval_param<T>&,
			std::vector<std::complex<T> > (&trees_result)[5],
			std::vector<size_t>& cut_mass,
			CutType& self);
};


template <class CutType> struct Normal_RatPent_Specification {
public:
	enum { MUPOINTS=3 }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=2 }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=2 }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=7 }; // Number of C coefficients
	enum { CBUBPOINTS=3}; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=3}; // Number of points in y used to evaluate for the bubble

	typedef pentagon_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatPent_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
};

template <class CutType, int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> struct General_RatPent_Specification {
public:
	enum { MUPOINTS=MUPOINTS_T }; // Number of powers of mu^2 in the box subtraction terms
	enum { MUTRIPOINTS=MUTRIPOINTS_T }; // Number of powers of mu^2 in the triangle subtraction terms
	enum { MUBUBPOINTS=MUBUBPOINTS_T }; // Number of powers of mu^2 in the bubble subtraction terms
	enum { CPOINTS=CPOINTS_T }; // Number of C coefficients
	enum { CBUBPOINTS=CBUBPOINTS_T }; // Number of points in t used to evaluate for the bubble
	enum { YPOINTS=YPOINTS_T}; // Number of points in y used to evaluate for the bubble

	typedef pentagon_rat_eval_points<MUPOINTS,MUTRIPOINTS,MUBUBPOINTS,CPOINTS,CBUBPOINTS,YPOINTS> MomentaEvaluator;
	typedef Normal_RatPent_Corner_Tree_Strategy<MomentaEvaluator,CutType>  CornerTreeStrategy;
};



template <class CutType> struct Higgs_RatPent_Specification : public  General_RatPent_Specification<CutType,MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS> {};


// A class to store all the parameters to be passed between the box and the pentagon
template <class T, int MUPOINTS_T> class pent_param
{
public:
	momentum<std::complex<T> > K1, K2, K4; // The momenta of the corners
	Cmom<T> K1flatc, K2flatc; // The momenta of the basis vectors
	momentum<std::complex<T> > vec1c, vec2c;
	std::complex<T> alp1, alp2; // Basis momenta factors

	std::complex<T> gamma; // Basis factor

	std::complex<T> boxpoles[2*MUPOINTS_T]; // The location of the poles in t for the box solution
	std::complex<T> pent_sub_ret[MUPOINTS_T*(MUPOINTS_T-2)]; // The pentagon subtraction term

	size_t pent_corner; // Stores the corner of the box that was opened to get the pentagon

	double accuracy; // Store the accuracy of our computation

	long int mcID; // Used in the mom_conf code
	size_t vertex_ref; // Used in the mom_conf code
};



//! Class for pentagon rat_ext diagrams
template <class BASE,  class RatPentSpecs/* = Normal_RatPent_Specification<BASE> */> class pentagon_Rat : public BASE,  public RatPentSpecs::CornerTreeStrategy {
	long int mcID;
	unsigned long int epID; // Used for tracking whether we have evaluated anything or not
	std::vector<int> indID;

	C E0coeff;
    CHP E0coeff_HP;
    CVHP E0coeff_VHP;
#if BH_USE_GMP
    CGMP E0coeff_GMP;
#endif
protected:
    //This vector is used for storing the indexes that get passed to the cut trees
 	std::vector<int> indlst[5];

	eval_param<R>* _ep[5];
	eval_param<RHP>* _ep_HP[5];
	eval_param<RVHP>* _ep_VHP[5];
#if BH_USE_GMP
	eval_param<RGMP>* _ep_GMP[5];
#endif

	std::vector<size_t> _cut_mass; // Store the unique mass_params of the cut legs


	bool _k1massive, _k2massive; // Keeps track of wether the K1 and K2 legs are massive or not

public:
	mass_param_coll* _leg_masses; // Store the mass_params for each cut leg
	template <class T> pentagon_Rat(T& t);
	template <class T> pentagon_Rat(T* t);
	virtual ~pentagon_Rat();

	//! Evaluation of box coefficients
	/** \param mc momentum configuration to be used \param ind indices of the momenta to be used \returns box coefficient*/
	virtual C eval(mom_conf& mc,const std::vector<int>& ind);
	virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind);
	virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind);

	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);
#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
	virtual CGMP eval(const eval_param<RGMP>& ep);
#endif
	// Computes the box subtraction terms for the triangle
	template <class T> void get_sub_terms(momentum_configuration<T>& mc,const std::vector<int>& ind, pent_param<T,RatPentSpecs::MUPOINTS>& pp);
	template <class T> void get_sub_terms(const eval_param<T>& ep, pent_param<T,RatPentSpecs::MUPOINTS>& pp);

	inline void set_eval(long int mcset, const std::vector<int>& indset){mcID=mcset;indID=indset;};
	inline void set_eval(unsigned long int mcset){epID=mcset;};
	inline bool is_eval(long int mcset, const std::vector<int>& indset){if((mcset==mcID)&&(indset==indID)){return true;} return false;};
	inline bool is_eval(unsigned long int mcset){if(mcset==epID){return true;} return false;};

	// Add a cutD to the list
	void add_mass(const cutD& pcd);
	void add_mass(int m1, int m2, int m3,int m4,int m5);

	void get_ep(eval_param<R>**& ep){ep=_ep;}
	void get_ep(eval_param<RHP>**& ep){ep=_ep_HP;}
	void get_ep(eval_param<RVHP>**& ep){ep=_ep_VHP;}
#if BH_USE_GMP
	void get_ep(eval_param<RGMP>**& ep){ep=_ep_GMP;}
#endif
	bool is_k1massive() const { return _k1massive;}
	bool is_k2massive() const { return _k2massive;}

   	typedef typename RatPentSpecs::MomentaEvaluator MomentaEvaluatorType;
protected:
	// Computes the pentagon coefficient
	template <class T> void get_coeffs(momentum_configuration<T>& mc, const std::vector<int>& ind, const long int original_mcID, const size_t vertex_ref);
	template <class T> void get_coeffs(const eval_param<T>& ep);

	void set_E0coeff(C value){E0coeff=value;};
	void set_E0coeff(CHP value){E0coeff_HP=value;E0coeff=to_double(value);};
	void set_E0coeff(CVHP value){E0coeff_VHP=value;E0coeff_HP=to_HP(value);E0coeff=to_double(value);};

	void get_E0coeff(C& value){value=E0coeff;};
	void get_E0coeff(CHP& value){value=E0coeff_HP;};
	void get_E0coeff(CVHP& value){value=E0coeff_VHP;};
#if BH_USE_GMP
	void set_E0coeff(CGMP value){E0coeff_GMP=value;E0coeff_VHP=to_VHP(value);E0coeff_HP=to_HP(value);E0coeff=to_double(value);};
	void get_E0coeff(CGMP& value){value=E0coeff_GMP;};
#endif
	void init();
};


}
}
#endif /* PENTAGON_RATEXT_H_ */
