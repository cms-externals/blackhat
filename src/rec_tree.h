/*
 * rec_tree.h
 *
 *  Created on: 24-Nov-2008
 *      Author: daniel
 */

#ifndef REC_TREE_H_
#define REC_TREE_H_

#include "rational.h"  // needed for Rec_BB
#include "counted.h"
#ifdef BH_USE_GMP
#include "gmp_r.h"
#endif

namespace BH {

class coupling_summary;

template <class T> std::complex<T>  ZeroF_eval(const eval_param<T>& ep,const mass_param_coll& masses);

//! class for trees included in a recursive construction
class Rec_Tree : public HelAmpl, public Rec_BB  {
public:
	//! constructor
	/** \param pro process for the amplitude */
	Rec_Tree(const process& pro): HelAmpl(pro) {};
	//! process
	/** \return process of the amplitude */
	const process& get_process() const {return d_process;};
	virtual ~Rec_Tree(){};
};

//! class for recursive trees for whom an analytic evaluation is possible
class Known_Rec_Tree_base : public Rec_Tree {
#ifndef SWIG
	friend std::ostream& operator<<(std::ostream& ,const Known_Rec_Tree_base& );
#endif
	protected:

	Tree_Fn_Ptr_eval _eval_C_ep_ptr;
	Tree_Fn_Ptr_eval_HP _eval_CHP_ep_ptr;
	Tree_Fn_Ptr_eval_VHP _eval_CVHP_ep_ptr;

#ifdef BH_USE_GMP
	Tree_Fn_Ptr_eval_GMP _eval_CGMP_ep_ptr;
#endif

	mass_param_coll* _masses;
public:
	//! constructor
	Known_Rec_Tree_base(const process& pro): Rec_Tree(pro) {_masses=new mass_param_coll(pro);};

	virtual C eval(mom_conf& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};

#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
#endif

	virtual C eval(const eval_param<R>& ep){return (*_eval_C_ep_ptr)(ep,*_masses);};
	virtual CHP eval(const eval_param<RHP>& ep){return (*_eval_CHP_ep_ptr)(ep,*_masses);};
	virtual CVHP eval(const eval_param<RVHP>& ep){return (*_eval_CVHP_ep_ptr)(ep,*_masses);};
#ifdef BH_USE_GMP
	virtual CGMP eval(const eval_param<RGMP>& ep){return (*_eval_CGMP_ep_ptr)(ep,*_masses);};
#endif

	virtual ~Known_Rec_Tree_base(){delete _masses;};

private:
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind){eval_param<T> ep(mc,ind); return eval(ep);};
};

//! class for recursive trees for whom an analytic evaluation is possible
class Known_Rec_Tree : public Known_Rec_Tree_base {
public:
	//! constructor
	Known_Rec_Tree(const process& pro);
	virtual ~Known_Rec_Tree(){};
};

class Known_Rec_Tree_offset : public Known_Rec_Tree_base {
	size_t _offset;
	size_t _length;
public:
	Known_Rec_Tree_offset(const process& pro,size_t offset=0);

	// eval(mc,ind) is taken over by the bse class, relaying it to the ep eval

	//! evaluation of the rational term in simple precision
	virtual C eval(const eval_param<R>& ep);
	//! evaluation of the rational term in double precision
	virtual CHP eval(const eval_param<RHP>& ep);
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(const eval_param<RVHP>& ep);

#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(const eval_param<RGMP>& ep);
#endif


	int get_offset(){return _offset;};
	virtual ~Known_Rec_Tree_offset(){};//Destructor to delete all the new rational terms we have computed

};


class Known_Rec_Tree_permutation : public Known_Rec_Tree_base {
	std::vector<int> d_perm_ind;
public:
	Known_Rec_Tree_permutation(const process& pro,const std::vector<int>& per);

	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);

#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(const eval_param<RGMP>& ep);
#endif


	virtual ~Known_Rec_Tree_permutation(){};

	const std::vector<int>& get_perm_ind(){return d_perm_ind;}
};


//! class for recursive trees for whom an analytic evaluation is not known but has to be done using recursion
class Unknown_Rec_Tree : public Rec_Tree {
protected:
	size_t i, j; //The two shifted legs

public:
	//! constructor
	Unknown_Rec_Tree(const process& amp_pro,const std::vector<particle_ID>&); //Default constructor to set up HelAmpl
	Unknown_Rec_Tree(const process& amp_pro,const std::vector<particle_ID>&,int OneB_i,int OneB_j); //Default constructor to set up HelAmpl
	Unknown_Rec_Tree(const process& amp_pro); //Default constructor to set up HelAmpl
	virtual ~Unknown_Rec_Tree(){};//Destructor to delete all the new rational terms we have computed

	virtual C eval(momentum_configuration<double>&,const std::vector<int>&);
	virtual CHP eval(momentum_configuration<dd_real>&,const std::vector<int>&);
	virtual CVHP eval(momentum_configuration<qd_real>&,const std::vector<int>&);
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
#endif

	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> eval(const eval_param<RGMP>& ep);
#endif

	//! Retrieve which (one-based) legs were shifted
	void get_shifted_legs(size_t* amp_i, size_t* amp_j){*amp_i=i+1;*amp_j=j+1;};
	//! [ part of the shift
	int get_i(){return i+1;};
	//! > part of the shift
	int get_j(){return j+1;};

private:
	//Computes the recursive rational terms
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T> eval_tree(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> std::complex<T> eval_tree(eval_param<T>& ep);
	void clean_up(const coupling_summary&);
	void generate_recursion(const std::vector<particle_ID>&);
	void generate_recursion(const std::vector<particle_ID>&,int OneB_i,int OneB_j);
};

//! factory class for trees in a recursive construction
class Tree_factory {
public :
	//! returns a new tree object
	Rec_Tree* new_tree(const process& pro);
	Rec_Tree* new_tree(const process& pro, const std::vector<int>& labels){
		throw BHerror("Not implemented");
	}; // Here we ignore any labels as on-shell recursion does not use them
};


//! a class for storing Rec_Trees
template <class FAC, class T> class Rec_Tree_cache;
template <class FAC, class T> std::ostream& operator<<(std::ostream&, const Rec_Tree_cache<FAC,T> &);

template <class FAC, class T> class Rec_Tree_cache {
private:
	std::map<process,T*> _RTmap;
#ifndef SWIG
	friend std::ostream& operator<< <FAC,T> (std::ostream&, const Rec_Tree_cache<FAC,T>&);
#endif
public:
	Rec_Tree_cache(){};
	~Rec_Tree_cache(){flush();}; // We must delete all objects here because we do not know if they are being used elsewhere in the code and so only when we destroy this do we know to get rid of them

	T* get_Rec_Tree(const process& pro); // Use to add a new element or to get a pointer to the previously added one
	T* get_Rec_Tree(const process& pro, const std::vector<int>& labels);
	void flush(); // Delete all the elements in the map and then empty the map

	int size() {return _RTmap.size();};
};

//! factory class for trees in a recursive construction with caching of the structures to reduce memory consumption
template <class FAC, class T> class Tree_factory_WC;
template <class FAC, class T> std::ostream& operator<<(std::ostream&, const Tree_factory_WC<FAC,T>&);

template <class FAC, class T> class Tree_factory_WC {
private:
	static Rec_Tree_cache<FAC,T> _storage; // We want to cache trees over the whole evaluation
#ifndef SWIG
	friend std::ostream& operator<< <FAC,T> (std::ostream&, const Tree_factory_WC<FAC,T>&);
#endif
public :
	Tree_factory_WC() {};
	~Tree_factory_WC() {};

	//! returns a new tree object
	T* new_tree(const process& pro) {return _storage.get_Rec_Tree(pro);};
	T* new_tree(const process& pro, const std::vector<int>& labels); // Construct tree assuming some legs can be reused, those that cannot are marked as -1 in the labels
	void delete_tree(T* del) {}; // We do not want to delete trees here, this will be done when the cache is flushed

	static int cache_size() {return _storage.size();};
	static void flush() { _storage.flush();}
};


std::ostream& operator<<(std::ostream& s, Rec_Tree& RT);
std::ostream& operator<<(std::ostream& s, Known_Rec_Tree& RT);
std::ostream& operator<<(std::ostream& s, Unknown_Rec_Tree& RT);

    
/*
 *
 *
 * A cached rec_tree class
 *
 *
 */

class Rec_Tree_WC : public Rec_Tree {
private:
    Rec_Tree* _rtn;
public:
    Rec_Tree_WC(Rec_Tree* rtn) : Rec_Tree(rtn->get_process()), _rtn(rtn) {};
    ~Rec_Tree_WC() {};
    
    virtual C eval(const eval_param<R>& ep){return _rtn->eval(ep);};
    virtual CHP eval(const eval_param<RHP>& ep){return _rtn->eval(ep);};
    virtual CVHP eval(const eval_param<RVHP>& ep){return _rtn->eval(ep);};
    
    virtual C eval(mom_conf& mc,const std::vector<int>& ind){return _rtn->eval(mc,ind);};
    virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return _rtn->eval(mc,ind);};
    virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return _rtn->eval(mc,ind);};
    
    virtual C get_value(const eval_param<R>& ep){return _rtn->get_value(ep);};
    virtual CHP get_value(const eval_param<RHP>& ep){return _rtn->get_value(ep);};
    virtual CVHP get_value(const eval_param<RVHP>& ep){return _rtn->get_value(ep);};
    
    virtual C get_value(mom_conf& mc,const std::vector<int>& ind){return _rtn->get_value(mc,ind);};
    virtual CHP get_value(mom_conf_HP& mc,const std::vector<int>& ind){return _rtn->get_value(mc,ind);};
    virtual CVHP get_value(mom_conf_VHP& mc,const std::vector<int>& ind){return _rtn->get_value(mc,ind);};
};

/*
 *
 *
 * Factory class for the Rec_Tree_WC functions
 *
 *
 */

typedef Tree_factory_WC<Tree_factory,Rec_Tree> REC_TREE_FACTORY_WC;

//! factory class for trees in a recursive construction
class Rec_Tree_WC_factory {
    public :
    //! returns a new tree object
    Rec_Tree_WC* new_tree(const process& pro) {REC_TREE_FACTORY_WC TF;Rec_Tree* prtn=TF.new_tree(pro); return new Rec_Tree_WC(prtn);};
    Rec_Tree_WC* new_tree(const process& pro, const std::vector<int>& labels){return new_tree(pro);}; // Here we ignore any labels as on-shell recursion does not use them
};

}

#endif /* REC_TREE_H_ */
