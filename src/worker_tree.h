/*
 * worker_tree.h
 *
 *  Created on: 10 Jul 2009
 *      Author: daniel
 */

#ifndef WORKER_TREE_H_
#define WORKER_TREE_H_

#include "BH_typedefs.h"
#include <vector>
#include "computable.h"
#include "spinor.h"
#include "eval_param.h"

namespace BH {
	template <class T> class momentum_configuration;
	template <class T> class eval_param;
	class mass_param_coll;

	class Rec_Tree;
	void write(Rec_Tree* RT,std::ostream& os,bool topLevel=true);

	class process;
}


namespace BH {

typedef C (*Tree_Fn_Ptr)(momentum_configuration<R>&,const std::vector<int>&);
typedef CHP (*Tree_Fn_Ptr_HP)(momentum_configuration<RHP>&,const  std::vector<int>&);
typedef CVHP (*Tree_Fn_Ptr_VHP)(momentum_configuration<RVHP>&,const std::vector<int>&);

typedef C (*Tree_Fn_Ptr_eval)(const eval_param<R>&, const mass_param_coll&);
typedef CHP (*Tree_Fn_Ptr_eval_HP)(const eval_param<RHP>&, const mass_param_coll&);
typedef CVHP (*Tree_Fn_Ptr_eval_VHP)(const eval_param<RVHP>&, const mass_param_coll&);

#ifdef BH_USE_GMP
typedef CGMP (*Tree_Fn_Ptr_GMP)(momentum_configuration<RGMP>&,const std::vector<int>&);
typedef CGMP (*Tree_Fn_Ptr_eval_GMP)(const eval_param<RGMP>&, const mass_param_coll&);
#endif

class Tree_Pair_base;

class worker_tree : public computable<std::complex> {
public:
/*      no need for this
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

*/
	virtual ~worker_tree(){};


};

class worker_tree_known  : public worker_tree {

	Tree_Fn_Ptr_eval d_eval_C_ep_ptr;
	Tree_Fn_Ptr_eval_HP d_eval_CHP_ep_ptr;
	Tree_Fn_Ptr_eval_VHP d_eval_CVHP_ep_ptr;

#ifdef BH_USE_GMP
	Tree_Fn_Ptr_eval_GMP d_eval_CGMP_ep_ptr;
#endif

	mass_param_coll* _masses;
public:
//	worker_tree_known(std::istream& is);
	worker_tree_known(const process& pro);
	virtual C eval(momentum_configuration<R>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
	virtual CHP eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
	virtual CVHP eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){  return eval_fn(mc,ind);};

	virtual C eval(const eval_param<R>& ep){ return (*d_eval_C_ep_ptr)(ep,*_masses);};
	virtual CHP eval(const eval_param<RHP>& ep){ return (*d_eval_CHP_ep_ptr)(ep,*_masses);};
	virtual CVHP eval(const eval_param<RVHP>& ep){ return (*d_eval_CVHP_ep_ptr)(ep,*_masses);};

#ifdef BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){  return eval_fn(mc,ind);};
	virtual CGMP eval(const eval_param<RGMP>& ep){ return (*d_eval_CGMP_ep_ptr)(ep,*_masses);};
#endif

private:
	void init(const process&);
	template <class T>  std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind){
		eval_param<T> ep(mc,ind);
		return eval(ep);
	};

};

class worker_tree_known_offset  : public worker_tree_known {
	int d_offset;
	int d_length;
public:
	worker_tree_known_offset(std::istream& is);
	worker_tree_known_offset(const process& pro,int offset);
	virtual C eval(momentum_configuration<R>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind); };
	virtual CHP eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind); };
	virtual CVHP eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind); };

	virtual C eval(const eval_param<R>& ep){ return eval_fn(ep); };
	virtual CHP eval(const eval_param<RHP>& ep){ return eval_fn(ep); };
	virtual CVHP eval(const eval_param<RVHP>& ep){ return eval_fn(ep); };
private:
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T> eval_fn(const eval_param<T>& ep);
};



class worker_tree_unknown : public worker_tree {
	int d_nbr_pairs;
	std::vector<Tree_Pair_base*> d_pairs;
public:
	worker_tree_unknown(std::istream& is);
	virtual C eval(momentum_configuration<R>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
	virtual CHP eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
	virtual CVHP eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};

	virtual C eval(const eval_param<R>& ep){ return eval_fn(ep);};
	virtual CHP eval(const eval_param<RHP>& ep){ return eval_fn(ep);};
	virtual CVHP eval(const eval_param<RVHP>& ep){ return eval_fn(ep);};

private:
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T> eval_fn(const eval_param<T>& ep);
//	not implemented
	worker_tree_unknown(const worker_tree_unknown&);
	worker_tree_unknown& operator=(const worker_tree_unknown&);

};
// need to define these structs because otherwise it is difficult to know the size
// of the eval_params before constructing them (which would be required because there is no default constructor)
// this way the indices_ are constructed first. If the eval params are not in a struct but as a member, they would
// need to be constructed first

struct indices_struct {
	int size;
	std::vector<int> indices;
	indices_struct(std::istream& is,bool prependMinusOne);
};



struct eval_params_struct {
	eval_param<R> _epl, _epr;
	eval_param<RHP> _epl_HP, _epr_HP;
	eval_param<RVHP> _epl_VHP, _epr_VHP;

#ifdef BH_USE_GMP
	eval_param<RGMP> _epl_GMP, _epr_GMP;
#endif
    eval_params_struct(int n_left,int n_right):
    	_epl(n_left), _epr(n_right),
    	_epl_HP(n_left), _epr_HP(n_right),
    	_epl_VHP(n_left), _epr_VHP(n_right)
#ifdef BH_USE_GMP
    	,_epl_GMP(n_left), _epr_GMP(n_right)
#endif
    {};
};

template <class T> struct momenta_struct {
	Cmom<T> PHat;
	Cmom<T> mPHat;
	Cmom<T> shifti;
	Cmom<T> shiftj;
};

class Tree_Pair_base {
	worker_tree* d_left;
	worker_tree* d_right;

	indices_struct d_left_indices;
	indices_struct d_right_indices;
	eval_params_struct d_eps;


public:
	std::vector<int> indshiftl, indshiftr;
	int maxl,maxr,max;
	int shifted_ind_j,shifted_ind_i;

	Tree_Pair_base(std::istream& is);
	int left_size(){return d_left_indices.size;};
	int right_size(){return d_right_indices.size;};
	int left_ind(int i){return d_left_indices.indices[i];};
	int right_ind(int i){return d_right_indices.indices[i];};
		virtual C eval(momentum_configuration<R>& mc,const std::vector<int>& ind)=0;
		virtual CHP eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind)=0;
		virtual CVHP eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind)=0;
	C eval_right(momentum_configuration<R>& mc,const std::vector<int>& ind){ return d_right->eval(mc,ind);};
	CHP eval_right(momentum_configuration<RHP>& mc,const std::vector<int>& ind){ return d_right->eval(mc,ind);};
	CVHP eval_right(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){ return d_right->eval(mc,ind);};
	C eval_left(momentum_configuration<R>& mc,const std::vector<int>& ind){ return d_left->eval(mc,ind);};
	CHP eval_left(momentum_configuration<RHP>& mc,const std::vector<int>& ind){ return d_left->eval(mc,ind);};
	CVHP eval_left(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){ return d_left->eval(mc,ind);};

	virtual C eval(const eval_param<R>& ep)=0;
	virtual CHP eval(const eval_param<RHP>& ep)=0;
	virtual CVHP eval(const eval_param<RVHP>& ep)=0;

	C eval_right(const eval_param<R>& ep){ return d_right->eval(ep);};
	CHP eval_right(const eval_param<RHP>& ep){ return d_right->eval(ep);};
	CVHP eval_right(const eval_param<RVHP>& ep){ return d_right->eval(ep);};

	C eval_left(const eval_param<R>& ep){ return d_left->eval(ep);};
	CHP eval_left(const eval_param<RHP>& ep){ return d_left->eval(ep);};
	CVHP eval_left(const eval_param<RVHP>& ep){ return d_left->eval(ep);};

	virtual ~Tree_Pair_base(){};

	 template <class T> Cmom<T>& get_phat(){};
	 template <class T> Cmom<T>& get_mphat(){};
	 template <class T> Cmom<T>& get_i_leg(){};
	 template <class T> Cmom<T>& get_j_leg(){};

#ifndef SWIG
	 template <class T> eval_param<T>& get_l_eval_param(){};
	 template <class T> eval_param<T>& get_r_eval_param(){};
#endif


private:
	Tree_Pair_base(const Tree_Pair_base&);
	Tree_Pair_base& operator=(const Tree_Pair_base&);


};


template <> inline eval_param<R>& Tree_Pair_base::get_l_eval_param(){return d_eps._epl;};
template <> inline eval_param<RHP>& Tree_Pair_base::get_l_eval_param(){return d_eps._epl_HP;};
template <> inline eval_param<RVHP>& Tree_Pair_base::get_l_eval_param(){return d_eps._epl_VHP;};

template <> inline eval_param<R>& Tree_Pair_base::get_r_eval_param(){return d_eps._epr;};
template <> inline eval_param<RHP>& Tree_Pair_base::get_r_eval_param(){return d_eps._epr_HP;};
template <> inline eval_param<RVHP>& Tree_Pair_base::get_r_eval_param(){return d_eps._epr_VHP;};


#ifdef BH_USE_GMP
template <> inline eval_param<RGMP>& Tree_Pair_base::get_l_eval_param(){return d_eps._epl_GMP;};
template <> inline eval_param<RGMP>& Tree_Pair_base::get_r_eval_param(){return d_eps._epr_GMP;};
#endif



class shift_base {
protected:
	int d_i,d_j;
public:
	shift_base(std::istream& is);
};


template <class Pair> class massless_shift : public shift_base {
public:
	massless_shift(std::istream& is): shift_base(is){};
	template <class T> std::complex<T> generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa);
};

template <class Pair> class massive_shift : public shift_base {
	size_t (*shift_ij)(momentum_configuration<R>&, std::vector<int>&, int, int, size_t, const std::complex<R>&);
	size_t (*shift_ij_HP)(momentum_configuration<RHP>&, std::vector<int>&, int, int, size_t, const std::complex<RHP>&);
	size_t (*shift_ij_VHP)(momentum_configuration<RVHP>&, std::vector<int>&, int, int, size_t, const std::complex<RVHP>&);
	void (*shift_ij_ep)(const eval_param<R>&, int, int, int, int, Cmom<R>&, Cmom<R>&, Cmom<R>&, const momentum<std::complex<R> >&, const std::complex<R>&, const Cmom<R>*&, const Cmom<R>*&);
	void (*shift_ij_ep_HP)(const eval_param<RHP>&, int, int, int, int, Cmom<RHP>&, Cmom<RHP>&, Cmom<RHP>&, const momentum<std::complex<RHP> >&, const std::complex<RHP>&, const Cmom<RHP>*&, const Cmom<RHP>*&);
	void (*shift_ij_ep_VHP)(const eval_param<RVHP>&, int, int, int, int, Cmom<RVHP>&, Cmom<RVHP>&, Cmom<RVHP>&, const momentum<std::complex<RVHP> >&, const std::complex<RVHP>&, const Cmom<RVHP>*&, const Cmom<RVHP>*&);
#ifdef BH_USE_GMP
	size_t (*shift_ij_GMP)(momentum_configuration<RGMP>&, std::vector<int>&, int, int, size_t, const std::complex<RGMP>&);
	void (*shift_ij_ep_GMP)(const eval_param<RGMP>&, int, int, int, int, Cmom<RGMP>&, Cmom<RGMP>&, Cmom<RGMP>&, const momentum<std::complex<RGMP> >&, const std::complex<RGMP>&, const Cmom<RGMP>*&, const Cmom<RGMP>*&);
#endif
	int _imass, _jmass;

public:
	 size_t get_shifted_ij(momentum_configuration<R>& mc, std::vector<int>& ind, size_t P, const std::complex<R>& Psqr){return shift_ij(mc,ind,d_i,d_j,P,Psqr);};
	 size_t get_shifted_ij(momentum_configuration<RHP>& mc, std::vector<int>& ind, size_t P, const std::complex<RHP>& Psqr){return shift_ij_HP(mc,ind,d_i,d_j,P,Psqr);};
	 size_t get_shifted_ij(momentum_configuration<RVHP>& mc, std::vector<int>& ind, size_t P, const std::complex<RVHP>& Psqr){return shift_ij_VHP(mc,ind,d_i,d_j,P,Psqr);};
	 void get_shifted_ij(const eval_param<R>& ep, Cmom<R>& shifti, Cmom<R>& shiftj, Cmom<R>& Phat, const momentum<std::complex<R> >& P, const std::complex<R>& Psqr,const Cmom<R>*& ref_i,const Cmom<R>*& ref_j){shift_ij_ep(ep,d_i,d_j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};
	 void get_shifted_ij(const eval_param<RHP>& ep, Cmom<RHP>& shifti, Cmom<RHP>& shiftj, Cmom<RHP>& Phat,const momentum<std::complex<RHP> >& P, const std::complex<RHP>& Psqr,const Cmom<RHP>*& ref_i,const Cmom<RHP>*& ref_j){shift_ij_ep_HP(ep,d_i,d_j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};
	 void get_shifted_ij(const eval_param<RVHP>& ep, Cmom<RVHP>& shifti, Cmom<RVHP>& shiftj, Cmom<RVHP>& Phat, const momentum<std::complex<RVHP> >& P, const std::complex<RVHP>& Psqr,const Cmom<RVHP>*& ref_i,const Cmom<RVHP>*& ref_j){shift_ij_ep_VHP(ep,d_i,d_j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};

	 massive_shift(std::istream& is);
	template <class T> std::complex<T> generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa);
#ifdef BH_USE_GMP
	 size_t get_shifted_ij(momentum_configuration<RGMP>& mc, std::vector<int>& ind, size_t P, const std::complex<RGMP>& Psqr){return shift_ij_GMP(mc,ind,d_i,d_j,P,Psqr);};
	 void get_shifted_ij(const eval_param<RGMP>& ep, Cmom<RGMP>& shifti, Cmom<RGMP>& shiftj, Cmom<RGMP>& Phat, const momentum<std::complex<RGMP> >& P, const std::complex<RGMP>& Psqr,const Cmom<RGMP>*& ref_i,const Cmom<RGMP>*& ref_j){shift_ij_ep_GMP(ep,d_i,d_j,_imass,_jmass,shifti,shiftj,Phat,P,Psqr,ref_i,ref_j);};
#endif

};

template <class Pair> class massive_prop_shift : public massive_shift<Pair> {
	size_t _mass_leg;
public:
	massive_prop_shift(std::istream& is);
	template <class T> std::complex<T> generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa);
};

template <class Pair> class massive_prop_massless_shift : public shift_base {
//	size_t (*shift_ij)(momentum_configuration<R>&, std::vector<int>&, int, int, size_t, const std::complex<R>&);
//	size_t (*shift_ij_HP)(momentum_configuration<RHP>&, std::vector<int>&, int, int, size_t, const std::complex<RHP>&);
//	size_t (*shift_ij_VHP)(momentum_configuration<RVHP>&, std::vector<int>&, int, int, size_t, const std::complex<RVHP>&);
int d_mass_leg;

public:
//	 size_t get_shifted_ij(momentum_configuration<R>& mc, std::vector<int>& ind, size_t P, const std::complex<R>& Psqr){return shift_ij(mc,ind,d_i,d_j,P,Psqr);};
//	 size_t get_shifted_ij(momentum_configuration<RHP>& mc, std::vector<int>& ind, size_t P, const std::complex<RHP>& Psqr){return shift_ij_HP(mc,ind,d_i,d_j,P,Psqr);};
//	 size_t get_shifted_ij(momentum_configuration<RVHP>& mc, std::vector<int>& ind, size_t P, const std::complex<RVHP>& Psqr){return shift_ij_VHP(mc,ind,d_i,d_j,P,Psqr);};
	massive_prop_massless_shift(std::istream& is): shift_base(is){ is >> d_mass_leg; };
	template <class T> std::complex<T> generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa);
};

template <class Pair> class massive_unshifted_shift : public shift_base {
public:
	massive_unshifted_shift(std::istream& is): shift_base(is){};
	template <class T> std::complex<T> generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa);
};

template < template <class> class ShiftType> class Tree_Pair: public Tree_Pair_base, public ShiftType<Tree_Pair_base> {
public:
	Tree_Pair(std::istream& is): Tree_Pair_base(is), ShiftType<Tree_Pair_base>(is) {};
	virtual C eval(momentum_configuration<R>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
	virtual CHP eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
	virtual CVHP eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};

	virtual C eval(const eval_param<R>& ep){ return eval_fn(ep);};
	virtual CHP eval(const eval_param<RHP>& ep){ return eval_fn(ep);};
	virtual CVHP eval(const eval_param<RVHP>& ep){ return eval_fn(ep);};

private:
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T> eval_fn(const eval_param<T>& ep);

};


typedef Tree_Pair<massless_shift> massless_pair;
typedef Tree_Pair<massive_shift> massive_pair;
typedef Tree_Pair<massive_prop_shift> massive_prop_pair;
typedef Tree_Pair<massive_prop_massless_shift> massive_prop_massless_pair;
typedef Tree_Pair<massive_unshifted_shift> massive_unshifted_pair;

worker_tree* create_worker_tree(std::istream& is);

class worker_tree_factory  {
public:
	worker_tree* new_tree(const process& pro);
};



}
#endif /* WORKER_TREE_H_ */
