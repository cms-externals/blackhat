/*!\file cut_part.h
\brief header for the cut part
*/
#ifndef CUT_PART_H_
#define CUT_PART_H_

#include "cut_part_base.h"
#include "amplitudes.h"
#include "BH_utilities.h"
#include "integrals.h"
#include "mode_dependent_typedefs.h"
#include "multi_precision_constant.h"
namespace BH {




typedef SeriesC<R> (*Cut_Fn_Ptr)(momentum_configuration<R>&,const std::vector<int>&,int);
typedef SeriesC<RHP> (*Cut_Fn_Ptr_HP)(momentum_configuration<RHP>&,const  std::vector<int>&,int);
typedef SeriesC<RVHP> (*Cut_Fn_Ptr_VHP)(momentum_configuration<RVHP>&,const std::vector<int>&,int);

typedef SeriesC<R> (*Cut_Fn_eval_Ptr)(const eval_param<R>&,const R&);
typedef SeriesC<RHP> (*Cut_Fn_eval_Ptr_HP)(const eval_param<RHP>&,const RHP&);
typedef SeriesC<RVHP> (*Cut_Fn_eval_Ptr_VHP)(const eval_param<RVHP>&,const RVHP&);



#ifdef BH_USE_OLD_KNOWN_CUT_PART
//! Cut_Part represent the cut-part of an OneLoopHelAmpl
class Known_Cut_Part : public Cut_Part_base {
		Cut_Fn_Ptr d_eval_C_ptr;
		Cut_Fn_Ptr_HP d_eval_CHP_ptr;
		Cut_Fn_Ptr_VHP d_eval_CVHP_ptr;

		Cut_Fn_eval_Ptr d_eval_C_ep_ptr;
		Cut_Fn_eval_Ptr_HP d_eval_CHP_ep_ptr;
		Cut_Fn_eval_Ptr_VHP d_eval_CVHP_ep_ptr;
	public :
		SeriesC<R> eval(mom_conf& mc,const std::vector<int>& ind){if (d_mu_index==0){d_mu_index = DefineMu<R>(mc,m_MU);}; return (*d_eval_C_ptr)(mc,ind,d_mu_index);};
		SeriesC<RHP> eval(mom_conf_HP& mc,const std::vector<int>& ind){if (d_mu_index_HP==0){d_mu_index_HP = DefineMu<RHP>(mc,m_MU);};return (*d_eval_CHP_ptr)(mc,ind,d_mu_index_HP);};
		SeriesC<RVHP> eval(mom_conf_VHP& mc,const std::vector<int>& ind){if (d_mu_index_VHP==0){d_mu_index_VHP = DefineMu<RVHP>(mc,m_MU);};return (*d_eval_CVHP_ptr)(mc,ind,d_mu_index_VHP);};

		SeriesC<R> eval(const eval_param<R>& ep){multi_precision_constant mmu(get_mu()); return (*d_eval_C_ep_ptr)(ep,mmu*mmu);};
		SeriesC<RHP> eval(const eval_param<RHP>& ep){multi_precision_constant mmu(get_mu()); return (*d_eval_CHP_ep_ptr)(ep,mmu*mmu);};
		SeriesC<RVHP> eval(const eval_param<RVHP>& ep){multi_precision_constant mmu(get_mu()); return (*d_eval_CVHP_ep_ptr)(ep,mmu*mmu);};

		//! Constructor
		Known_Cut_Part(const process& pro,color_structure cs);
		virtual ~Known_Cut_Part(){};
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
    	SeriesC<R> get_conjugate_cut_part(){return SeriesC<R>(-2,0);};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){return SeriesC<RHP>(-2,0);};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){return SeriesC<RVHP>(-2,0);};
    	//template <class T> SeriesC<T> get_conjugate_cut_part(){return SeriesC<T>(-2,0);};

};

class Known_Cut_Part_offset : public Cut_Part_base {
		size_t d_offset;
		Cut_Fn_Ptr d_eval_C_ptr;
		Cut_Fn_Ptr_HP d_eval_CHP_ptr;
		Cut_Fn_Ptr_VHP d_eval_CVHP_ptr;

		Cut_Fn_eval_Ptr d_eval_C_ep_ptr;
		Cut_Fn_eval_Ptr_HP d_eval_CHP_ep_ptr;
		Cut_Fn_eval_Ptr_VHP d_eval_CVHP_ep_ptr;
	public :
		SeriesC<R> eval(mom_conf& mc,const std::vector<int>& ind);
		SeriesC<RHP> eval(mom_conf_HP& mc,const std::vector<int>& ind);
		SeriesC<RVHP> eval(mom_conf_VHP& mc,const std::vector<int>& ind);
		SeriesC<R> eval(const eval_param<R>& ep);
		SeriesC<RHP> eval(const eval_param<RHP>& ep);
		SeriesC<RVHP> eval(const eval_param<RVHP>& ep);

		size_t get_offset() const {return d_offset;};
		//! Constructor
		Known_Cut_Part_offset(const process& pro,color_structure cs,size_t offset);
		virtual ~Known_Cut_Part_offset(){};
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
    	//template <class T> SeriesC<T> get_conjugate_cut_part(){return SeriesC<T>(-2,0);};
    	SeriesC<R> get_conjugate_cut_part(){return SeriesC<R>(-2,0);};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){return SeriesC<RHP>(-2,0);};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){return SeriesC<RVHP>(-2,0);};
};

#endif /*BH_USE_OLD_KNOWN_CUT_PART*/

#ifndef BH_PUBLIC

//! Cut_Part represent the cut-part of an OneLoopHelAmpl
class Cut_Part : public Cut_Part_base {
protected:
	std::vector<computable<std::complex>*> eval_bubbles;
	std::vector<computable<std::complex>*> eval_triangles;
	std::vector<computable<std::complex>*> eval_boxes;
	template <class T> SeriesC<T> Cut_Part_eval(momentum_configuration<T>& mc,const std::vector<int>& ind);
	template <class T> SeriesC<T> Cut_Part_eval(const eval_param<T>& ep);

	SeriesC<R> eval_with_check(momentum_configuration<R>& mc,const std::vector<int>& ind);
	SeriesC<RHP> eval_with_check(momentum_configuration<RHP>& mc,const std::vector<int>& ind);
	SeriesC<RVHP> eval_with_check(momentum_configuration<RVHP>& mc,const std::vector<int>& ind);
	SeriesC<R> eval_with_check(const eval_param<R>& ep);
	SeriesC<RHP> eval_with_check(const eval_param<RHP>& ep);
	SeriesC<RVHP> eval_with_check(const eval_param<RVHP>& ep);

	std::vector<bubbleD*> bubbles;
	std::vector<triangleD*> triangles;
	std::vector<boxD*> boxes;
public :
	virtual SeriesC<R> eval(mom_conf&,const std::vector<int>&);
	virtual SeriesC<RHP> eval(mom_conf_HP&,const std::vector<int>&);
	virtual SeriesC<RVHP> eval(mom_conf_VHP&,const std::vector<int>&);

	virtual SeriesC<R> eval(const eval_param<R>& ep);
	virtual SeriesC<RHP> eval(const eval_param<RHP>& ep);
	virtual SeriesC<RVHP> eval(const eval_param<RVHP>& ep);
#ifdef BH_USE_GMP
	virtual SeriesC<RGMP> eval(momentum_configuration<RGMP>&,const std::vector<int>&);
#endif

	//	virtual C eval(mom_conf);
//	virtual bool is_zero();
	//! box cut diagram
	/** \param i (1-based) index \return boxD object representing the ith box diagram */
	boxD* box(size_t i) {return boxes[i-1];};
	//! triangle cut diagram
	/** \param i (1-based) index \return triangleD object representing the ith triangle diagram */
	triangleD* triangle(size_t i) {return triangles[i-1];};
	//! bubble cut diagram
	/** \param i (1-based) index \return bubbleD object representing the ith bubble diagram */
	bubbleD* bubble(size_t i) {return bubbles[i-1];};
	/** \return number of box cut diagrams */
	size_t nbr_boxes() const {return boxes.size();};
	/** \return number of triangle cut diagrams */
	size_t nbr_triangles() const {return triangles.size();};
	/** \return number of bubble cut diagrams */
	size_t nbr_bubbles() const {return bubbles.size();};
//	//! Constructor
//	Cut_Part() {} ;
	//! Constructor
	Cut_Part(Cut_Part* cp, cutD_factory* = cutD_factory::default_CF,option* opt=option::always_true);
	//! Constructor
	Cut_Part(const process& p,const std::vector<particle_ID>& possible_props, cutD_factory* = cutD_factory::default_CF, option* = option::always_true) ;
	//! Constructor
	Cut_Part(const process& p,const std::vector<particle_ID>& possible_props, const ordering_constraint& ordered, cutD_factory* = cutD_factory::default_CF, option* = option::always_true) ;
	//! Constructor
	explicit Cut_Part(const process& p) : Cut_Part_base(p) {};
	//! Constructor needed for Cut_part_D_Dims constructor
	/** \param A pointer to the OnleLoopAmplitude for which the helicity amplitude is constructed
		\param hel vector containing the helicities of the external legs.
	 */
//	OneLoopHelAmpl(OneLoopAmplitude* A,vector<short> hel) ;
	virtual ~Cut_Part();
	void set_bubble_eval_daughter(size_t n,computable<std::complex>* cp){eval_bubbles[n-1]=cp;};
	computable<std::complex>* eval_bubble(size_t i) {return eval_bubbles[i-1];};
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
    	//template <class T> SeriesC<T> get_conjugate_cut_part(){return SeriesC<T>(-2,0);};
    	SeriesC<R> get_conjugate_cut_part(){return SeriesC<R>(-2,0);};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){return SeriesC<RHP>(-2,0);};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){return SeriesC<RVHP>(-2,0);};
#ifdef BH_USE_GMP
    	SeriesC<RGMP> get_conjugate_cut_part_GMP(){return SeriesC<RGMP>(-2,0);};
#endif

};

#endif


}

#endif /*CUT_PART_H_*/
