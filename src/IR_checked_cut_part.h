/*
 * IR_checked_cut_part.h
 *
 *  Created on: 2 Jul 2009
 *      Author: daniel
 */

#ifndef IR_CHECKED_CUT_PART_H_
#define IR_CHECKED_CUT_PART_H_

#include "IR_checked.h"

namespace BH {

//! class for Cut_Part that has an IR test (big IR)
class IR_checked_Cut_Part : public Cut_Part_base, public IR_checked {
	Cut_Part_base* d_cut_part_p;
	
	//double _accuracy;
    //SeriesC<R> _conjugate_cut_part;
public :
	//! evaluation of the cut part
	virtual SeriesC<R> eval(mom_conf&,const std::vector<int>&);
	//! evaluation of the cut part
	virtual SeriesC<RHP> eval(mom_conf_HP&,const std::vector<int>&);
	//! evaluation of the cut part
	virtual SeriesC<RVHP> eval(mom_conf_VHP&,const std::vector<int>&);
	//! evaluation of the cut part
	virtual SeriesC<R> eval(const eval_param<R>&);
	//! evaluation of the cut part
	virtual SeriesC<RHP> eval(const eval_param<RHP>&);
	//! evaluation of the cut part
	virtual SeriesC<RVHP> eval(const eval_param<RVHP>&);
	//! IR tolerances
	/** \retun the tolerance set for this particular object. */

#ifdef BH_USE_GMP
	virtual SeriesC<RGMP> eval(momentum_configuration<RGMP>&,const std::vector<int>&);
	virtual SeriesC<RGMP> eval(const eval_param<RGMP>&);
#endif

	template <class T> inline T IR_tolerance();

	/* return the Cut_Part_base */
	Cut_Part_base* cut_part(){return d_cut_part_p;};
//! constructor
	IR_checked_Cut_Part(Cut_Part_base* ptr): Cut_Part_base(ptr->get_process()), IR_checked(1e-4,1e-9,1e-30) , d_cut_part_p(ptr)/*, _accuracy(0) ,_conjugate_cut_part(SeriesC<R>(-2,0))*/ {};
	virtual ~IR_checked_Cut_Part(){ delete d_cut_part_p;};
//	void set_tolerances(R t_R,R t_RHP,R t_RVHP){_IR_tolerance_R= t_R;_IR_tolerance_RHP= t_RHP;_IR_tolerance_RVHP = t_RVHP; };

	// This gives the estimated accuracy of the computation
	double get_accuracy(){ return (d_cut_part_p->get_accuracy());};
    	SeriesC<R> get_conjugate_cut_part(){ return (d_cut_part_p->get_conjugate_cut_part());};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){ return (d_cut_part_p->get_conjugate_cut_part_HP());};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){ return (d_cut_part_p->get_conjugate_cut_part_VHP());};
#ifdef BH_USE_GMP
    	SeriesC<RGMP> get_conjugate_cut_part_GMP(){ return (d_cut_part_p->get_conjugate_cut_part_GMP());};
#endif

};

}
#endif /* IR_CHECKED_CUT_PART_H_ */
