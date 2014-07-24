/*
 * standard_cut_part.h
 *
 *  Created on: Jun 29, 2009
 *      Author: daniel
 */

#ifndef STANDARD_CUT_PART_H_
#define STANDARD_CUT_PART_H_


#include <vector>
#include <map>
#include <complex>
#include <iosfwd>
#include "cut_part.h"
#include "cut_Darren.h"
#include "cached_integral.h"

namespace BH {

namespace cut {

template <class BOX,class TRI,class BUB>  class standard_cut_part : public Cut_Part_base {
private:
	double _accuracy;
    SeriesC<R> _conjugate_cut_part;
    SeriesC<RHP> _conjugate_cut_part_HP;
    SeriesC<RVHP> _conjugate_cut_part_VHP;
#ifdef BH_USE_GMP
    SeriesC<RGMP> _conjugate_cut_part_GMP;
#endif
	std::vector<CachedIntegral::Cached_Box_Integral_User*> d_box_integrals;
	std::vector<CachedIntegral::Cached_Triangle_Integral_User*> d_triangle_integrals;
	std::vector<CachedIntegral::Cached_Bubble_Integral_User*> d_bubble_integrals;
protected:
	std::vector<BOX*> d_boxes;
    std::vector<TRI*> d_triangles;
    std::vector<BUB*> d_bubbles;
public:
    typedef BOX BoxType;
    typedef TRI TriType;
    typedef BUB BubType;
	standard_cut_part (const process& PRO):Cut_Part_base(PRO){};
#ifndef BH_PUBLIC
	template <class OrigType> standard_cut_part (OrigType& orig, option* opt,cutD_factory* cf);
#endif
	//	virtual C eval(mom_conf);
//	virtual bool is_zero();
	//! box cut diagram
	/** \param i (1-based) index \return boxD object representing the ith box diagram */
	BOX* box(size_t i) {return d_boxes[i-1];};
	//! triangle cut diagram
	/** \param i (1-based) index \return triangleD object representing the ith triangle diagram */
	TRI* triangle(size_t i) {return d_triangles[i-1];};
	//! bubble cut diagram
	/** \param i (1-based) index \return bubbleD object representing the ith bubble diagram */
	BUB* bubble(size_t i) {return d_bubbles[i-1];};
	/** \return number of box cut diagrams */
	size_t nbr_boxes() const {return d_boxes.size();};
	/** \return number of triangle cut diagrams */
	size_t nbr_triangles() const {return d_triangles.size();};
	/** \return number of bubble cut diagrams */
	size_t nbr_bubbles() const {return d_bubbles.size();};
	virtual ~standard_cut_part();
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return _accuracy;};
    	SeriesC<R> get_conjugate_cut_part(){return _conjugate_cut_part;};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){return _conjugate_cut_part_HP;};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){return _conjugate_cut_part_VHP;};
#ifdef BH_USE_GMP
    	SeriesC<RGMP> get_conjugate_cut_part_GMP(){return _conjugate_cut_part_GMP;};
#endif

	
protected:
//    template <class T> SeriesC<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind);
      SeriesC<R> eval_fn(momentum_configuration<R>& mc,const std::vector<int>& ind);
      SeriesC<RHP> eval_fn(momentum_configuration<RHP>& mc,const std::vector<int>& ind);
      SeriesC<RVHP> eval_fn(momentum_configuration<RVHP>& mc,const std::vector<int>& ind);
#ifdef BH_USE_GMP
      SeriesC<RGMP> eval_fn(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
#endif
      template <class T> SeriesC<T> eval_without_check(momentum_configuration<T>& mc,const std::vector<int>& ind);
      SeriesC<R> eval_with_check(momentum_configuration<R>& mc,const std::vector<int>& ind);
      SeriesC<R> eval_with_check_wCI(momentum_configuration<R>& mc,const std::vector<int>& ind);
      SeriesC<RHP> eval_with_check(momentum_configuration<RHP>& mc,const std::vector<int>& ind);
      SeriesC<RHP> eval_with_check_wCI(momentum_configuration<RHP>& mc,const std::vector<int>& ind);
      SeriesC<RVHP> eval_with_check(momentum_configuration<RVHP>& mc,const std::vector<int>& ind);
      SeriesC<RVHP> eval_with_check_wCI(momentum_configuration<RVHP>& mc,const std::vector<int>& ind);

    template <class T> SeriesC<T> eval_fn(const eval_param<T>& ep);
    template <class T> SeriesC<T> eval_without_check(const eval_param<T>& ep);
	SeriesC<R> eval_with_check(const eval_param<R>& ep);
	SeriesC<RHP> eval_with_check(const eval_param<RHP>& ep);
	SeriesC<RVHP> eval_with_check(const eval_param<RVHP>& ep);

#ifdef BH_USE_GMP
      SeriesC<RGMP> eval_with_check(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
      SeriesC<RGMP> eval_with_check_wCI(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
  	SeriesC<RGMP> eval_with_check(const eval_param<RGMP>& ep);
#endif

protected:
	void construct_integral_table();


};



}
}



#endif /* STANDARD_CUT_PART_H_ */
