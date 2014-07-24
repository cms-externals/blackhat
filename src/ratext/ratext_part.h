/*
 * ratext_part.h
 *
 *  Created on: 17-May-2009
 *      Author: daniel
 */

#ifndef RATEXT_PART_H_
#define RATEXT_PART_H_


#include <vector>
#include "partitions.h"



#include "eval_param.h"
#include "ratext/rat_ext.h"
#include "ratext/bubble_ratext.h"
#include "ratext/triangle_ratext.h"
#include "ratext/box_ratext.h"
#include "ratext/pentagon_ratext.h"
#include "process_utils.h"

#ifdef BH_PUBLIC
#include "worker_tree.h"
#endif


namespace BH {

namespace ratext {


template <class PENT,class BOX,class TRI,class BUB>  class ratext_part {
private:
	double _accuracy;
	
protected:
	std::vector<PENT*> d_pentagons;
	std::vector<BOX*> d_boxes;
    std::vector<TRI*> d_triangles;
    std::vector<BUB*> d_bubbles;
	symmetry_factor d_colour_struct_factor;

	TREE_TYPE* d_tree_ptr;
public:
	ratext_part (const process& p) : _accuracy(0) { TREE_FACTORY_TYPE TF; d_tree_ptr=TF.new_tree(fix_flavors(p));};

	//! box cut diagram
	/** \param i (1-based) index \return boxD object representing the ith box diagram */
	PENT* pentagon(size_t i) {return d_pentagons[i-1];};
	//! box cut diagram
	/** \param i (1-based) index \return boxD object representing the ith box diagram */
	BOX* box(size_t i) {return d_boxes[i-1];};
	//! triangle cut diagram
	/** \param i (1-based) index \return triangleD object representing the ith triangle diagram */
	TRI* triangle(size_t i) {return d_triangles[i-1];};
	//! bubble cut diagram
	/** \param i (1-based) index \return bubbleD object representing the ith bubble diagram */
	BUB* bubble(size_t i) {return d_bubbles[i-1];};
	/** \return number of pentagon cut diagrams */
	size_t nbr_pentagons() const {return d_pentagons.size();};
	/** \return number of box cut diagrams */
	size_t nbr_boxes() const {return d_boxes.size();};
	/** \return number of triangle cut diagrams */
	size_t nbr_triangles() const {return d_triangles.size();};
	/** \return number of bubble cut diagrams */
	size_t nbr_bubbles() const {return d_bubbles.size();};
	void  set_colour_fac(const symmetry_factor& sy){d_colour_struct_factor=sy;};
	void  set_colour_fac(int num,int den){d_colour_struct_factor=symmetry_factor(num,den);};
	template <class T> T get_colour_fac(){return d_colour_struct_factor.eval<T>();};
	template <class T> void get_colour_fac(T& col){col=d_colour_struct_factor.eval<T>();};

	//! returns the tree
	template <class T> std::complex<T> get_tree(momentum_configuration<T>& mc, const std::vector<int>& ind){return d_tree_ptr->get_value(mc,ind);}
	// computable not implemented for eval param, so we need to recompute
	template <class T> std::complex<T> get_tree(const eval_param<T>& ep){return d_tree_ptr->eval(ep);}

	// Define the type of objects stored in the vectors
	typedef PENT pent_type;
	typedef BOX box_type;
	typedef TRI tri_type;
	typedef BUB bub_type;

	virtual ~ratext_part();
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return _accuracy;};

protected:
	// Does the actual evaluation for the rational terms
	template <class T> std::complex<T> eval_rat(momentum_configuration<T>& mc, const std::vector<int>& ind);
	template <class T> std::complex<T> eval_rat_ep(const eval_param<T>& ep);

	// Does the actual evaluation for the rational terms
	virtual std::complex<R> eval_rat_IR_checked(momentum_configuration<R>& mc, const std::vector<int>& ind);
	virtual std::complex<RHP> eval_rat_IR_checked(momentum_configuration<RHP>& mc, const std::vector<int>& ind);
	virtual std::complex<RVHP> eval_rat_IR_checked(momentum_configuration<RVHP>& mc, const std::vector<int>& ind){return eval_rat(mc,ind);};

	virtual std::complex<R> eval_rat_IR_checked_ep(const eval_param<R>& ep);
	virtual std::complex<RHP> eval_rat_IR_checked_ep(const eval_param<RHP>& ep);
	virtual std::complex<RVHP> eval_rat_IR_checked_ep(const eval_param<RVHP>& ep){return eval_rat_ep(ep);};

#if BH_USE_GMP
	virtual std::complex<RGMP> eval_rat_IR_checked(momentum_configuration<RGMP>& mc, const std::vector<int>& ind);
	virtual std::complex<RGMP> eval_rat_IR_checked_ep(const eval_param<RGMP>& ep){return eval_rat_ep(ep);};
#endif

};




}

}

#endif /* RATEXT_PART_H_ */
