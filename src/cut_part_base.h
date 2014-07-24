/*
 * cut_part_base.h
 *
 *  Created on: 05.07.2012
 *      Author: daniel
 */

#ifndef CUT_PART_BASE_H_
#define CUT_PART_BASE_H_

#include "amplitudes.h"
#include "BH_utilities.h"
#include "mode_dependent_typedefs.h"
#include "multi_precision_constant.h"

namespace BH {




//! Cut_Part represent the cut-part of an OneLoopHelAmpl
class Cut_Part_base : public HelAmpl , public computable<SeriesC> {
	TREE_TYPE* d_tree_ptr;
protected:
		multi_precision_constant m_MU;
		int d_mu_index;
		int d_mu_index_HP;
		int d_mu_index_VHP;
#ifdef BH_USE_GMP
		int d_mu_index_GMP;
#endif

	    std::vector<int> _masses; // Stores the mass_labels of the external particles
	public :
//		virtual SeriesC<R> eval(mom_conf& mc,const std::vector<int>& ind){	int david_mu=DefineMu<R>(mc,m_MU);return (*_eval_C_ptr)(mc,ind,david_mu);};
//		virtual SeriesC<RHP> eval(mom_conf_HP& mc,const std::vector<int>& ind){int david_mu=DefineMu<RHP>(mc,m_MU);return (*_eval_CHP_ptr)(mc,ind,david_mu);};
//		virtual SeriesC<RVHP> eval(mom_conf_VHP& mc,const std::vector<int>& ind){int david_mu=DefineMu<RVHP>(mc,m_MU);return (*_eval_CVHP_ptr)(mc,ind,david_mu);};
		const process& get_process() const {return d_process;};
		//! Constructor
		Cut_Part_base(const process& pro);
		virtual ~Cut_Part_base();
		//! returns the tree
		template <class T> std::complex<T> get_tree(momentum_configuration<T>& mc, const std::vector<int>& ind);
		template <class T> std::complex<T> get_tree(const eval_param<T>& ep);
		void set_mu(int index){d_mu_index=index;};
		void set_mu_HP(int index){d_mu_index_HP=index;};
		void set_mu_VHP(int index){d_mu_index_VHP=index;};
#ifdef BH_USE_GMP
		void set_mu_GMP(int index){d_mu_index_GMP=index;};
		void set_mu(R val,RHP val_HP,RVHP val_VHP,const RGMP& val_GMP){m_MU.set(val,val_HP,val_VHP,val_GMP);d_mu_index=0;d_mu_index_HP=0;d_mu_index_VHP=0;};
#else
		void set_mu(R val,RHP val_HP,RVHP val_VHP){m_MU.set(val,val_HP,val_VHP);d_mu_index=0;d_mu_index_HP=0;d_mu_index_VHP=0;};
#endif
		void set_mu(const multi_precision_constant& multi){m_MU=multi;d_mu_index=0;d_mu_index_HP=0;d_mu_index_VHP=0;};
		const multi_precision_constant& get_mu(){return m_MU;};
		template <class T> inline int get_mu_index();
		template <class T> inline void set_mu_index(int index);
		virtual void dry_run(const std::vector<int>&){};

		virtual double get_accuracy()=0;
        virtual SeriesC<R> get_conjugate_cut_part()=0;
        virtual SeriesC<RHP> get_conjugate_cut_part_HP()=0;
        virtual SeriesC<RVHP> get_conjugate_cut_part_VHP()=0;
#if BH_USE_GMP
        virtual SeriesC<RGMP> get_conjugate_cut_part_GMP()=0;
#endif

};

template <> inline void Cut_Part_base::set_mu_index<R>(int ind){
	d_mu_index=ind;
}
template <> inline void Cut_Part_base::set_mu_index<RHP>(int ind){
	d_mu_index_HP=ind;
}
template <> inline void Cut_Part_base::set_mu_index<RVHP>(int ind){
	d_mu_index_VHP=ind;
}
template <> inline int Cut_Part_base::get_mu_index<R>(){
	return d_mu_index;
}
template <> inline int Cut_Part_base::get_mu_index<RHP>(){
	return d_mu_index_HP;
}
template <> inline int Cut_Part_base::get_mu_index<RVHP>(){
	return d_mu_index_VHP;
}
#ifdef BH_USE_GMP
template <> inline int Cut_Part_base::get_mu_index<RGMP>(){
	return d_mu_index_GMP;
}
template <> inline void Cut_Part_base::set_mu_index<RGMP>(int ind){
	d_mu_index_GMP=ind;
}
#endif

}

#endif /* CUT_PART_BASE_H_ */
