/*!\file IR_checked.h
\brief header for IR checked objects (One Loop Helicity Amplitudes, Spurious poles)
 */
#ifndef IR_CHECKED_H_
#define IR_CHECKED_H_

#include "OneLoopHelAmpl.h"
#include "cut_part_factory.h"

namespace BH {


//! class for IR checked objects
/** the IR_checked objects contain three tolerances, one for R, RHP and RVHP */
class IR_checked {
protected:
	R _IR_tolerance_R;
	RHP _IR_tolerance_RHP;
	RVHP _IR_tolerance_RVHP;
#ifdef BH_USE_GMP
	RGMP _IR_tolerance_RGMP;
#endif
public:
	template <class T> inline T IR_tolerance();
	IR_checked(R tol_r=1e-4,RHP tol_rhp=1e-18,RVHP tol_rvhp=1e-40): _IR_tolerance_R(tol_r), _IR_tolerance_RHP(tol_rhp), _IR_tolerance_RVHP(tol_rvhp) {} ;
	void set_tolerances(R t_R,R t_RHP,R t_RVHP){_IR_tolerance_R= t_R;_IR_tolerance_RHP= t_RHP;_IR_tolerance_RVHP = t_RVHP; };

};

//! Class for IR checked One loop helicity amplitudes
/** IR checked One loop helicity amplitudes have an IR_checked cut part and the rational part is build up from IR_checked spurious poles.
 	A IR_checked_OLHA has two different types of tolerances. The one for the  "Big IR" test and the one for the "little IR" test. They are set by two different member functions
   \sa IR_checked_Cut_Part IR_checked_Spurious_Pole

  */

class IR_checked_cut_part_factory : public cut_part_factory<Cut_Part_base> {
	public:
		virtual Cut_Part_base* new_cut_part(const process&,color_structure);
		static IR_checked_cut_part_factory* s_default_IR_checked_cut_part_factory;
		virtual ~IR_checked_cut_part_factory(){};
};


class IR_checked_OLHA : public One_Loop_Helicity_Amplitude, public IR_checked {
public:
	//! constructor
//	IR_checked_OLHA(const process& p,const std::vector<particle_ID>& possible_props, cutD_factory* cf= cutD_factory::default_CF, option* op= option::always_true);
	//! constructor
	IR_checked_OLHA(const process& p, color_structure cs );
	IR_checked_OLHA(const process& p, color_structure cs, Rational_factory<Rational_base>* RRF , cut_part_factory<IR_checked_Cut_Part>* CPF) ;
	//! evaluation of the amplitude
	SeriesC<R> eval(mom_conf&,const std::vector<int>&);
	//! evaluation of the amplitude
	SeriesC<RHP> eval(mom_conf_HP&,const std::vector<int>&);
	//! evaluation of the amplitude
	SeriesC<RVHP> eval(mom_conf_VHP&,const std::vector<int>&);

	//! evaluation of the amplitude
	SeriesC<R> eval(const eval_param<R>&);
	//! evaluation of the amplitude
	SeriesC<RHP> eval(const eval_param<RHP>&);
	//! evaluation of the amplitude
	SeriesC<RVHP> eval(const eval_param<RVHP>&);

#ifdef BH_USE_GMP
	SeriesC<RGMP> eval(momentum_configuration<RGMP>&,const std::vector<int>&);
	SeriesC<RGMP> eval(const eval_param<RGMP>&);
#endif

	//! sets the tolerances for the big IR test
	void set_bigIR_tolerances(R t_R,R t_RHP,R t_RVHP);
	//! sets the tolerances for the big IR test
	template <class T> inline T IR_tolerance();
	//! sets the tolerances for the big IR test
	void set_littleIR_tolerances(R t_R,R t_RHP,R t_RVHP);

	virtual ~IR_checked_OLHA(){};
private:
	void construct();
};



}
#endif /*IR_CHECKED_H_*/
