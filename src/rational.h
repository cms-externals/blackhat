/*!\file rational.h
\brief Header file for the rational part
*/
#ifndef RATIONAL_H_
#define RATIONAL_H_

#include "amplitudes.h"
#include "rec_BB.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"
namespace BH {


typedef C (*Tree_Fn_Ptr)(momentum_configuration<R>&,const std::vector<int>&);
typedef CHP (*Tree_Fn_Ptr_HP)(momentum_configuration<RHP>&,const  std::vector<int>&);
typedef CVHP (*Tree_Fn_Ptr_VHP)(momentum_configuration<RVHP>&,const std::vector<int>&);
#ifdef BH_USE_GMP
typedef std::complex<RGMP> (*Tree_Fn_Ptr_GMP)(momentum_configuration<RGMP>&,const std::vector<int>&);
#endif

typedef C (*Tree_Fn_Ptr_eval)(const eval_param<R>&, const mass_param_coll&);
typedef CHP (*Tree_Fn_Ptr_eval_HP)(const eval_param<RHP>&, const mass_param_coll&);
typedef CVHP (*Tree_Fn_Ptr_eval_VHP)(const eval_param<RVHP>&, const mass_param_coll&);
#ifdef BH_USE_GMP
typedef std::complex<RGMP> (*Tree_Fn_Ptr_eval_GMP)(const eval_param<RGMP>&, const mass_param_coll&);
#endif

#ifndef SWIG

//! class for constant factorization functions
class Const_Fact_Fn : public Rec_BB{
	C _value;
	CHP _value_HP;
	CVHP _value_VHP;
#ifdef BH_USE_GMP
	CGMP _value_GMP;
	int d_den;
	int d_num;
	int d_GMP_precision;
#endif
public:
#ifdef BH_USE_GMP
	Const_Fact_Fn(C c,CHP chp ,CVHP cvhp,CGMP cgmp): _value(c),_value_HP(chp),_value_VHP(cvhp),_value_GMP(cgmp) {
		d_den=0;
		d_num=0;
		d_GMP_precision=RGMP::get_current_precision();
	};
	explicit Const_Fact_Fn(int i): _value(i),_value_HP(i),_value_VHP(i),_value_GMP(i) {
		d_den=0;
		d_num=i;
		d_GMP_precision=RGMP::get_current_precision();
	};
	Const_Fact_Fn(int num,int den): _value(double(num)/double(den)),_value_HP(RHP(num)/RHP(den)),_value_VHP(RVHP(num)/RVHP(den)),_value_GMP(RGMP(num)/RGMP(den)) {
		d_den=den;
		d_num=num;
		d_GMP_precision=RGMP::get_current_precision();
	};
#else
	Const_Fact_Fn(C c,CHP chp ,CVHP cvhp): _value(c),_value_HP(chp),_value_VHP(cvhp){};
	explicit Const_Fact_Fn(int i): _value(i),_value_HP(i),_value_VHP(i) {};
	Const_Fact_Fn(int num,int den): _value(double(num)/double(den)),_value_HP(RHP(num)/RHP(den)),_value_VHP(RVHP(num)/RVHP(den)) {};
#endif
	virtual C eval(mom_conf& mc, const std::vector<int>& ind){return _value;};
	virtual CHP eval(mom_conf_HP& mc, const std::vector<int>& ind){return _value_HP;};
	virtual CVHP eval(mom_conf_VHP& mc, const std::vector<int>& ind){return _value_VHP;};
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return eval_and_increase_precision_if_needed();};
#endif

	virtual C eval(const eval_param<R>& ep){return _value;};
	virtual CHP eval(const eval_param<RHP>& ep){return _value_HP;};
	virtual CVHP eval(const eval_param<RVHP>& ep){return _value_VHP;};
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(const eval_param<RGMP>& ep){return eval_and_increase_precision_if_needed();};
#endif
	virtual ~Const_Fact_Fn(){};//Destructor
	static Const_Fact_Fn* unity;
#ifdef BH_USE_GMP
private:
	std::complex<RGMP> eval_and_increase_precision_if_needed();
#endif

};

std::vector<particle_ID> possible_propagators(const process& pro);


//! class for pairs of elements in a recursive construction
/**Rec Pair contains two parts and a splitting function */
class Rec_Pair : public Rec_BB {
	part _part;
protected:
	int i;
	int j;

	// Store properties that do not change during evaluation
	size_t max,maxl,maxr,shifted_ind_j,shifted_ind_i;
	// The left and right hand indices that get passed to the daughters
	std::vector<int> indshiftl, indshiftr;
	// The left and right hand eval_params that get passed to the daughters
	eval_param<R> _epl, _epr;
	eval_param<RHP> _epl_HP, _epr_HP;
	eval_param<RVHP> _epl_VHP, _epr_VHP;
	
	Cmom<R> _PHat,_mPHat,_i_leg,_j_leg;
	Cmom<RHP> _PHat_HP,_mPHat_HP,_i_leg_HP,_j_leg_HP;
	Cmom<RVHP> _PHat_VHP,_mPHat_VHP,_i_leg_VHP,_j_leg_VHP;
#ifdef BH_USE_GMP
	eval_param<RGMP> _epl_GMP, _epr_GMP;
    Cmom<RGMP> _PHat_GMP,_mPHat_GMP,_i_leg_GMP,_j_leg_GMP;
#endif

public:
	//! returns the left part of the pair
	 Rec_BB* left(){return daughters[0];};
	//! returns the right part of the pair
	 Rec_BB* right(){return daughters[1];};
	//! returns the central "factorisation function" part of the pair
	 Rec_BB* fact(){return daughters[2];};
	 //! constructor
	 /** \param left left part \param right right part \param part partition of the external leg \param i [ part of the shift \param j > part of the shift \param Fact factorization function    */
#ifndef SWIG
	 Rec_Pair(Rec_BB* left,Rec_BB* right,part&,int i,int j, Const_Fact_Fn* fact=new Const_Fact_Fn(1));
	 //! return the partition of the external legs
#endif
	 part& get_part(){return _part;};
	 //! [ part of the shift
	 int get_i(){return i;};
	 //! > part of the shift
	 int get_j(){return j;};
	 int get_shifted_ind_i(){return shifted_ind_i;};
	 int get_shifted_ind_j(){return shifted_ind_j;};
	
	 template <class T> Cmom<T>& get_phat(){};
	 template <class T> Cmom<T>& get_mphat(){};
	 template <class T> Cmom<T>& get_i_leg(){};
	 template <class T> Cmom<T>& get_j_leg(){};


	 virtual C eval(mom_conf& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
#endif

	 virtual C eval(const eval_param<R>& ep){return Rec_Pair_eval(ep);};
	 virtual CHP eval(const eval_param<RHP>& ep){return Rec_Pair_eval(ep);};
	 virtual CVHP eval(const eval_param<RVHP>& ep){return Rec_Pair_eval(ep);};
#ifdef BH_USE_GMP
	 virtual CGMP eval(const eval_param<RGMP>& ep){return Rec_Pair_eval(ep);};
#endif

	 virtual ~Rec_Pair(){};

protected:
	template <class T> eval_param<T>& get_l_eval_param(){};
	template <class T> eval_param<T>& get_r_eval_param(){};	

private:
	template <class T> std::complex<T>  Rec_Pair_eval(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T>  Rec_Pair_eval(const eval_param<T>& ep);
};
	
template <> inline eval_param<R>& Rec_Pair::get_l_eval_param(){return _epl;};
template <> inline eval_param<RHP>& Rec_Pair::get_l_eval_param(){return _epl_HP;};
template <> inline eval_param<RVHP>& Rec_Pair::get_l_eval_param(){return _epl_VHP;};
	
template <> inline eval_param<R>& Rec_Pair::get_r_eval_param(){return _epr;};	
template <> inline eval_param<RHP>& Rec_Pair::get_r_eval_param(){return _epr_HP;};
template <> inline eval_param<RVHP>& Rec_Pair::get_r_eval_param(){return _epr_VHP;};
	
template <> inline Cmom<R>& Rec_Pair::get_phat(){return _PHat;};
template <> inline Cmom<RHP>& Rec_Pair::get_phat(){return _PHat_HP;};
template <> inline Cmom<RVHP>& Rec_Pair::get_phat(){return _PHat_VHP;};
	
template <> inline Cmom<R>& Rec_Pair::get_mphat(){return _mPHat;};
template <> inline Cmom<RHP>& Rec_Pair::get_mphat(){return _mPHat_HP;};
template <> inline Cmom<RVHP>& Rec_Pair::get_mphat(){return _mPHat_VHP;};

template <> inline Cmom<R>& Rec_Pair::get_i_leg(){return _i_leg;};
template <> inline Cmom<RHP>& Rec_Pair::get_i_leg(){return _i_leg_HP;};
template <> inline Cmom<RVHP>& Rec_Pair::get_i_leg(){return _i_leg_VHP;};

template <> inline Cmom<R>& Rec_Pair::get_j_leg(){return _j_leg;};
template <> inline Cmom<RHP>& Rec_Pair::get_j_leg(){return _j_leg_HP;};
template <> inline Cmom<RVHP>& Rec_Pair::get_j_leg(){return _j_leg_VHP;};

#ifdef BH_USE_GMP
template <> inline eval_param<RGMP>& Rec_Pair::get_l_eval_param(){return _epl_GMP;};
template <> inline eval_param<RGMP>& Rec_Pair::get_r_eval_param(){return _epr_GMP;};
template <> inline Cmom<RGMP>& Rec_Pair::get_phat(){return _PHat_GMP;};
template <> inline Cmom<RGMP>& Rec_Pair::get_mphat(){return _mPHat_GMP;};
template <> inline Cmom<RGMP>& Rec_Pair::get_i_leg(){return _i_leg_GMP;};
template <> inline Cmom<RGMP>& Rec_Pair::get_j_leg(){return _j_leg_GMP;};
#endif

/**Rec Pair contains two parts and a splitting function with a massive propagator*/
class Rec_Pair_massive : public Rec_Pair {
protected:
	// Stores the function for computing the ij shift
	size_t (*shift_ij)(momentum_configuration<R>&, std::vector<int>&, int, int, size_t, const std::complex<R>&);
	size_t (*shift_ij_HP)(momentum_configuration<RHP>&, std::vector<int>&, int, int, size_t, const std::complex<RHP>&);
	size_t (*shift_ij_VHP)(momentum_configuration<RVHP>&, std::vector<int>&, int, int, size_t, const std::complex<RVHP>&);

	// Stores the function for computing the ij shift
	void (*shift_ij_ep)(const eval_param<R>&, int, int, int, int, Cmom<R>&, Cmom<R>&, Cmom<R>&, const momentum<std::complex<R> >&, const std::complex<R>&, const Cmom<R>*&, const Cmom<R>*&);
	void (*shift_ij_ep_HP)(const eval_param<RHP>&, int, int, int, int, Cmom<RHP>&, Cmom<RHP>&, Cmom<RHP>&, const momentum<std::complex<RHP> >&, const std::complex<RHP>&, const Cmom<RHP>*&, const Cmom<RHP>*&);
	void (*shift_ij_ep_VHP)(const eval_param<RVHP>&, int, int, int, int, Cmom<RVHP>&, Cmom<RVHP>&, Cmom<RVHP>&, const momentum<std::complex<RVHP> >&, const std::complex<RVHP>&, const Cmom<RVHP>*&, const Cmom<RVHP>*&);
#ifdef BH_USE_GMP
	size_t (*shift_ij_GMP)(momentum_configuration<RGMP>&, std::vector<int>&, int, int, size_t, const std::complex<RGMP>&);
	void (*shift_ij_ep_GMP)(const eval_param<RGMP>&, int, int, int, int, Cmom<RGMP>&, Cmom<RGMP>&, Cmom<RGMP>&, const momentum<std::complex<RGMP> >&, const std::complex<RGMP>&, const Cmom<RGMP>*&, const Cmom<RGMP>*&);
#endif

private:
	int d_massive_shift;
	int _imass, _jmass;
public:
	 //! constructor
	 /** \param left left part \param right right part \param part partition of the external leg \param i [ part of the shift \param j > part of the shift \param Fact factorization function    */
	#ifndef SWIG
	Rec_Pair_massive(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, int imass=0, int jmass=0, int massive_shift=0, Const_Fact_Fn* fact=new Const_Fact_Fn(1));
	#endif
	 virtual ~Rec_Pair_massive(){};

	 size_t get_shifted_ij(momentum_configuration<R>& mc, std::vector<int>& ind, size_t P, const std::complex<R>& Psqr){return shift_ij(mc,ind,i,j,P,Psqr);};
	 size_t get_shifted_ij(momentum_configuration<RHP>& mc, std::vector<int>& ind, size_t P, const std::complex<RHP>& Psqr){return shift_ij_HP(mc,ind,i,j,P,Psqr);};
	 size_t get_shifted_ij(momentum_configuration<RVHP>& mc, std::vector<int>& ind, size_t P, const std::complex<RVHP>& Psqr){return shift_ij_VHP(mc,ind,i,j,P,Psqr);};

	 void get_shifted_ij(const eval_param<R>& ep, Cmom<R>& shifti, Cmom<R>& shiftj, Cmom<R>& Phat, const momentum<std::complex<R> >& P, const std::complex<R>& Psqr,const Cmom<R>*& ref_i,const Cmom<R>*& ref_j){shift_ij_ep(ep,i,j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};
	 void get_shifted_ij(const eval_param<RHP>& ep, Cmom<RHP>& shifti, Cmom<RHP>& shiftj, Cmom<RHP>& Phat,const momentum<std::complex<RHP> >& P, const std::complex<RHP>& Psqr,const Cmom<RHP>*& ref_i,const Cmom<RHP>*& ref_j){shift_ij_ep_HP(ep,i,j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};
	 void get_shifted_ij(const eval_param<RVHP>& ep, Cmom<RVHP>& shifti, Cmom<RVHP>& shiftj, Cmom<RVHP>& Phat, const momentum<std::complex<RVHP> >& P, const std::complex<RVHP>& Psqr,const Cmom<RVHP>*& ref_i,const Cmom<RVHP>*& ref_j){shift_ij_ep_VHP(ep,i,j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};

	 virtual C eval(mom_conf& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};

	 virtual C eval(const eval_param<R>& ep){return Rec_Pair_eval(ep);};
	 virtual CHP eval(const eval_param<RHP>& ep){return Rec_Pair_eval(ep);};
	 virtual CVHP eval(const eval_param<RVHP>& ep){return Rec_Pair_eval(ep);};
#ifdef BH_USE_GMP
	 size_t get_shifted_ij(momentum_configuration<RGMP>& mc, std::vector<int>& ind, size_t P, const std::complex<RGMP>& Psqr){return shift_ij_GMP(mc,ind,i,j,P,Psqr);};
	 void get_shifted_ij(const eval_param<RGMP>& ep, Cmom<RGMP>& shifti, Cmom<RGMP>& shiftj, Cmom<RGMP>& Phat, const momentum<std::complex<RGMP> >& P, const std::complex<RGMP>& Psqr,const Cmom<RGMP>*& ref_i,const Cmom<RGMP>*& ref_j){shift_ij_ep_GMP(ep,i,j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};
	 virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CGMP eval(const eval_param<RGMP>& ep){return Rec_Pair_eval(ep);};
#endif
	 int get_massive_shift(){return d_massive_shift;}
	 int get_imass() const {return _imass;}
	 int get_jmass() const {return _jmass;}
private:
	template <class T> std::complex<T>  Rec_Pair_eval(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T>  Rec_Pair_eval(const eval_param<T>& ep);
};

/**Rec Pair contains two parts and a splitting function with a massive propagator*/
class Rec_Pair_massive_prop : public Rec_Pair_massive {
	// Stores the mass_properties of the cut propagator leg
	size_t _mass_leg;

public:
	 //! constructor
	 /** \param left left part \param right right part \param part partition of the external leg \param i [ part of the shift \param j > part of the shift \param Fact factorization function    */
	#ifndef SWIG
	Rec_Pair_massive_prop(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, size_t in_mass_leg, int imass=0, int jmass=0, int massive_shift=0, Const_Fact_Fn* fact=new Const_Fact_Fn(1));
	#endif
	 virtual ~Rec_Pair_massive_prop(){};

	 size_t get_mass_leg(){return _mass_leg;};

	 virtual C eval(mom_conf& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};

	 virtual C eval(const eval_param<R>& ep){return Rec_Pair_eval(ep);};
	 virtual CHP eval(const eval_param<RHP>& ep){return Rec_Pair_eval(ep);};
	 virtual CVHP eval(const eval_param<RVHP>& ep){return Rec_Pair_eval(ep);};
#ifdef BH_USE_GMP
	 virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CGMP eval(const eval_param<RGMP>& ep){return Rec_Pair_eval(ep);};
#endif
private:
	template <class T> std::complex<T>  Rec_Pair_eval(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T>  Rec_Pair_eval(const eval_param<T>& ep);
};

/**Rec Pair contains two parts and a splitting function with a massive propagator*/
class Rec_Pair_massive_prop_massless_shift : public Rec_Pair {
	// Stores the mass_properties of the cut propagator leg
	size_t _mass_leg;

public:
	 //! constructor
	 /** \param left left part \param right right part \param part partition of the external leg \param i [ part of the shift \param j > part of the shift \param Fact factorization function    */
	#ifndef SWIG
	Rec_Pair_massive_prop_massless_shift(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, size_t in_mass_leg, Const_Fact_Fn* fact=new Const_Fact_Fn(1));
	#endif
	 virtual ~Rec_Pair_massive_prop_massless_shift(){};

	 size_t get_mass_leg(){return _mass_leg;};

	 virtual C eval(mom_conf& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};

	 virtual C eval(const eval_param<R>& ep){return Rec_Pair_eval(ep);};
	 virtual CHP eval(const eval_param<RHP>& ep){return Rec_Pair_eval(ep);};
	 virtual CVHP eval(const eval_param<RVHP>& ep){return Rec_Pair_eval(ep);};
#ifdef BH_USE_GMP
	 virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CGMP eval(const eval_param<RGMP>& ep){return Rec_Pair_eval(ep);};
#endif
private:
	template <class T> std::complex<T>  Rec_Pair_eval(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T>  Rec_Pair_eval(const eval_param<T>& ep);
};

/**Rec Pair contains two parts and a splitting function with a massive propagator*/
class Rec_Pair_massive_unshifted : public Rec_Pair {
public:
	 //! constructor
	 /** \param left left part \param right right part \param part partition of the external leg \param i [ part of the shift \param j > part of the shift \param Fact factorization function    */
	#ifndef SWIG
		Rec_Pair_massive_unshifted(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, Const_Fact_Fn* fact=new Const_Fact_Fn(1));
	#endif
		virtual ~Rec_Pair_massive_unshifted(){};

	 virtual C eval(mom_conf& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};

	 virtual C eval(const eval_param<R>& ep){return Rec_Pair_eval(ep);};
	 virtual CHP eval(const eval_param<RHP>& ep){return Rec_Pair_eval(ep);};
	 virtual CVHP eval(const eval_param<RVHP>& ep){return Rec_Pair_eval(ep);};
#ifdef BH_USE_GMP
	 virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return Rec_Pair_eval(mc,ind);};
	 virtual CGMP eval(const eval_param<RGMP>& ep){return Rec_Pair_eval(ep);};
#endif
private:
	template <class T> std::complex<T>  Rec_Pair_eval(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T>  Rec_Pair_eval(const eval_param<T>& ep);
};


std::ostream& operator<<(std::ostream& s, Rec_BB& RBB);
std::ostream& operator<<(std::ostream& s, Rec_Pair& RP);

void construct_partitions(const process& pro,const std::vector<particle_ID>& possible_props,std::vector<part>& leglist,int amp_i,int amp_j);

#endif

//! class for the recursive part of the rational part
class Rational_base : public Rec_BB {
protected:
	process d_process;
	
public:
	//! constructor
	Rational_base(const process& p): d_process(p) {};
	virtual ~Rational_base(){};
	virtual bool is_zero(){return false;};
	/* don't define
	virtual C eval(mom_conf&,const std::vector<int>&);
	virtual CHP eval(mom_conf_HP&,const std::vector<int>&);
	virtual CVHP eval(mom_conf_VHP&,const std::vector<int>&);
	*/
	const process& get_process() {return d_process;};
	virtual void dry_run(const std::vector<int>&){};
	
	// This gives the estimated accuracy of the computation
	virtual double get_accuracy()=0;
};





}
#endif /*RATIONAL_H_*/
