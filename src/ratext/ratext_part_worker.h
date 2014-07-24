/*
 * ratext_part_worker.h
 *
 *  Created on: 7 Jul 2009
 *      Author: daniel
 */

#ifndef RATEXT_PART_WORKER_H_
#define RATEXT_PART_WORKER_H_

#include "rational_factory.h"
#include "ratext_part.h"
#include "ratext/rat_worker.h"
#include "settings.h"

namespace BH {

namespace ratext {

template <class RatBubSpecs, class RatTriSpecs, class RatBoxSpecs, class RatPentSpecs> class general_worker_ratext : public Rational_base, private ratext_part<pentagon_Rat<rat_worker,RatPentSpecs>,box_Rat<rat_worker,RatBoxSpecs>,triangle_Rat<rat_worker,RatTriSpecs>,bubble_Rat<rat_worker,RatBubSpecs> > {
public:
	typedef pentagon_Rat<rat_worker,RatPentSpecs> pent_type;
	typedef box_Rat<rat_worker,RatBoxSpecs> box_type;
	typedef triangle_Rat<rat_worker,RatTriSpecs> tri_type;
	typedef bubble_Rat<rat_worker,RatBubSpecs> bub_type;
	typedef ratext_part<pentagon_Rat<rat_worker,RatPentSpecs>,box_Rat<rat_worker,RatBoxSpecs>,triangle_Rat<rat_worker,RatTriSpecs>,bubble_Rat<rat_worker,RatBubSpecs> > RP;
	general_worker_ratext(std::istream& is);
	virtual C eval(mom_conf& mc, const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual CHP eval(mom_conf_HP& mc, const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual CVHP eval(mom_conf_VHP& mc, const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual C eval(const eval_param<R>& ep){return RP::eval_rat_IR_checked_ep(ep);};
	virtual CHP eval(const eval_param<RHP>& ep){return RP::eval_rat_IR_checked_ep(ep);};
	virtual CVHP eval(const eval_param<RVHP>& ep){return RP::eval_rat_IR_checked_ep(ep);};
#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc, const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual CGMP eval(const eval_param<RGMP>& ep){return RP::eval_rat_IR_checked_ep(ep);};
#endif
	virtual double get_accuracy(){return RP::get_accuracy();};

private :
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc, const std::vector<int>& ind){
		if(settings::general::s_use_ep_only){
			eval_param<T> ep(mc,ind);
			return eval(ep);
		} else {
			if(settings::rational_settings::s_use_IR_in_ratext){
				return RP::eval_rat_IR_checked(mc,ind);
			} else {
				return RP::eval_rat(mc,ind);
			}
		}
	};

};

typedef general_worker_ratext<Normal_RatBub_Specification<rat_worker>, Normal_RatTri_Specification<rat_worker>, Normal_RatBox_Specification<rat_worker>, Normal_RatPent_Specification<rat_worker> > worker_ratext;
typedef general_worker_ratext<Normal_RatBub_Specification<rat_worker>, Normal_RatTri_Specification<rat_worker>, Normal_RatBox_Specification<rat_worker>, Normal_RatPent_Specification<rat_worker> > higgs_worker_ratext;

class worker_rational_factory : public Rational_factory<Rational_base> {
public:
	virtual Rational_base* new_rational(const process& pro,color_structure cs);
	virtual Rational_base* new_rational(const process& pro,const std::vector<particle_ID>& possible_props,option* opt){throw BHerror("Not implemented");};
	static Rational_factory<Rational_base>* s_default_WRF;
};


} /* ratext */

} /* BH */


#endif /* RATEXT_PART_WORKER_H_ */
