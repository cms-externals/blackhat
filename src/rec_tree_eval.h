/*
 * rec_tree_eval.h
 *
 *  Created on: Apr 7, 2009
 *      Author: dforde
 */

#ifndef REC_TREE_EVAL_H_
#define REC_TREE_EVAL_H_

#include "rational_eval.h"  // needed for Rec_BB
#include "eval_param.h"

namespace BH {

template <class T> complex<T>  ZeroF_eval(const eval_param<T>& ep);

//! class for trees included in a recursive construction
class Rec_Tree_eval : public HelAmpl, public Rec_BB_eval {
public:
	//! constructor
	/** \param pro process for the amplitude */
	Rec_Tree_eval(const process& pro): HelAmpl(pro) {};
	//! process
	/** \return process of the amplitude */
	const process& get_process() const {return d_process;};
	virtual ~Rec_Tree_eval(){};
};

//! class for recursive trees for whom an analytic evaluation is possible
class Known_Rec_Tree_base_eval : public Rec_Tree_eval {
protected:
	Tree_Fn_Ptr_eval _eval_C_ptr;
	Tree_Fn_Ptr_eval_HP _eval_CHP_ptr;
	Tree_Fn_Ptr_eval_VHP _eval_CVHP_ptr;
public:
	//! constructor
	Known_Rec_Tree_base_eval(const process& pro): Rec_Tree_eval(pro){};
	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);
	virtual ~Known_Rec_Tree_base_eval(){};
};

//! class for recursive trees for whom an analytic evaluation is possible
class Known_Rec_Tree_eval : public Known_Rec_Tree_base_eval {
public:
	//! constructor
	Known_Rec_Tree_eval(const process& pro);
	virtual ~Known_Rec_Tree_eval(){};
};

class Known_Rec_Tree_offset_eval : public Known_Rec_Tree_base_eval {
	size_t _offset;
	size_t _length;
public:
	Known_Rec_Tree_offset_eval(const process& pro,size_t offset=0);
	//! evaluation of the rational term in simple precision
	virtual C eval(const eval_param<R>& ep);
	//! evaluation of the rational term in double precision
	virtual CHP eval(const eval_param<RHP>& ep);
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(const eval_param<RVHP>& ep);
	virtual ~Known_Rec_Tree_offset_eval(){};//Destructor to delete all the new rational terms we have computed

};


class Known_Rec_Tree_permutation_eval : public Known_Rec_Tree_eval {
	std::vector<int> _perm_ind;
public:
	Known_Rec_Tree_permutation_eval(const process& pro,const vector<int>& per);
	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);
	virtual ~Known_Rec_Tree_permutation_eval(){};

};


//! class for recursive trees for whom an analytic evaluation is not known but has to be done using recursion
class Unknown_Rec_Tree_eval : public Rec_Tree_eval {
protected:
	size_t i, j; //The two shifted legs

public:
	//! constructor
	Unknown_Rec_Tree_eval(const process& amp_pro,const std::vector<particle_ID>&); //Default constructor to set up HelAmpl
	Unknown_Rec_Tree_eval(const process& amp_pro); //Default constructor to set up HelAmpl
	virtual ~Unknown_Rec_Tree_eval(){};//Destructor to delete all the new rational terms we have computed

	virtual C eval(const eval_param<R>& ep);
	virtual CHP eval(const eval_param<RHP>& ep);
	virtual CVHP eval(const eval_param<RVHP>& ep);

	//! Retrieve which (one-based) legs were shifted
	void get_shifted_legs(size_t* amp_i, size_t* amp_j){*amp_i=i+1;*amp_j=j+1;};
	//! [ part of the shift
	int get_i(){return i+1;};
	//! > part of the shift
	int get_j(){return j+1;};

private:
	//Computes the recursive rational terms
	template <class T> complex<T> eval_tree(eval_param<T>& ep);
	void clean_up(const coupling_summary_eval&);
	void generate_recursion(const std::vector<particle_ID>&);
};

//! factory class for trees in a recursive construction
class Tree_factory_eval {
public :
	//! returns a new tree object
	Rec_Tree_eval* new_tree(const process& pro);
};

std::ostream& operator<<(std::ostream& s, Rec_Tree_eval& RT);
std::ostream& operator<<(std::ostream& s, Known_Rec_Tree_eval& RT);
std::ostream& operator<<(std::ostream& s, Unknown_Rec_Tree_eval& RT);

}



#endif /* REC_TREE_EVAL_H_ */
