
#include "rec_rational_factory.h"
#include "particles.h"
#include "amplitudes.h"
#include "cut_part.h"
#include "rational.h"
#include "tree_amp.h"
#include <vector>
#include <memory>
#include <algorithm>
#include "cut_part_factory.h"
#include "process_utils.h"
#ifndef BH_PUBLIC
#include "shifts.h"
#endif
#include "IR_checked.h"
#include "IR_checked_cut_part.h"
//#include "cut/amplitudes_cut.h"
#include "iterators.h"
#include "cached_integral.h"
#include "ratext/ratext_part.h"
#ifndef BH_PUBLIC
#include "cut_part_normal.h"
#endif
#include "cut_part_worker.h"
#include "settings.h"
#include "BH_debug.h"

using namespace std;

namespace BH {


Cut_Part_base* Known_cut_part_factory::new_cut_part(const process& pro,color_structure cs){

	if ( settings::general::s_use_cached_integrals  ) {
		vector<int> ind; for (int i=1;i<=pro.n();i++){ind.push_back(i);};
		CachedIntegral::Cut_Part_wCI* Cptr = CachedIntegral::CwCI_Ptr(pro,cs,ind);
		if ( Cptr ) {
			BH_DEBUG_MESSAGE2("using known cut part wCI for ",pro);
			return new CachedIntegral::Known_Cut_wCI(pro,Cptr);
		}
	}
#ifdef NOTDEFINED /*  we only wnat to use cached integrals */

	if ( C_Ptr<R>(pro,cs) ) {
		    BH_DEBUG_MESSAGE2("using known cut part for ",pro);
		return new Known_Cut_Part(pro,cs);
	} else {
		return 0;
	}
#endif
	return 0;
}

Known_cut_part_factory global_KCPF;
Known_cut_part_factory* Known_cut_part_factory::s_default_KCPF=&global_KCPF;

#ifdef BH_USE_OLD_CPF
template <class CUTD_FACTORY> Cut_Part_base* custom_cut_part_factory<CUTD_FACTORY>::new_cut_part(const process& pro,color_structure cs){

	if ( settings::general::s_use_known_formulae  ) {

		if ( settings::general::s_use_cached_integrals  ) {
			vector<int> ind; for (int i=1;i<=pro.n();i++){ind.push_back(i);};
			CachedIntegral::Cut_Part_wCI* Cptr = CachedIntegral::CwCI_Ptr(pro,cs,ind);
			if ( Cptr ) {
				BH_DEBUG_MESSAGE2("using known cut part wCI for ",pro);
				return new CachedIntegral::Known_Cut_wCI(pro,Cptr);
			} else {
				for (int offset=1;offset<= pro.n();offset++) {
				vector<particle_ID> pp;
				for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((offset+(i-2))%pro.n()+1)) ; }
				process new_pro(pp);
					CachedIntegral::Cut_Part_wCI* Cptr = CachedIntegral::CwCI_Ptr(new_pro,cs,ind);
					if( Cptr ){
				        BH_DEBUG_MESSAGE2("using known cut part wCI for ",pro);
						return new CachedIntegral::Known_Cut_wCI_offset(pro,Cptr,offset);
					}
				 }
			}
		}
/*  we only wnat to use cached integrals */
#ifdef BH_USE_OLD_CUT_PART
		if ( C_Ptr<R>(pro,cs) ) {
		    BH_DEBUG_MESSAGE2("using known cut part for ",pro);
			return new Known_Cut_Part(pro,cs);
		}
		else {
			for (int offset=1;offset<= pro.n();offset++) {
			vector<particle_ID> pp;
			for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((offset+(i-2))%pro.n()+1)) ; }
				process new_pro(pp);
				if( C_Ptr<R>(new_pro,cs) != 0){
				    BH_DEBUG_MESSAGE2("using known cut part for ",pro);
					return new Known_Cut_Part_offset(pro,cs,offset);
				}
			 }
		}
#else
		throw BHerror("should always use cached integrals");
#endif
	}



//	vector<int> flavors;
//	size_t max_internal_flavor=1;
//	for(int i=1;i<=pro.n();i++){
//		  if( (*pro.p(i).type()) == quark ) {
//			  vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
//			  if ( pos == flavors.end() ){
//				  flavors.push_back(pro.p(i).flavor());
//				  if (pro.p(i).flavor() > max_internal_flavor) {max_internal_flavor = pro.p(i).flavor();}
//			  }
//
//		  }
//	}
	CUTD_FACTORY  Factory;
//	Darren_factory  Darren_F;


		vector<ph_type> pp;
		process new_pro=pro;
		option* the_opt=get_possible_props_and_options(pro,cs,pp,new_pro);
		std::auto_ptr<option> ap_op(the_opt);
		ordering_constraint oc=get_ordering_constraint(pro);
//		_PRINT(oc);
		if ( (! oc.weak.empty()) || oc.strong.size()!=1 || oc.strong[0].size() != pro.n()  ){
			if ( settings::general::s_use_cached_integrals  ) {
				return new CachedIntegral::Unknown_Cut_Part_wCI(new_pro,pp,oc,&Factory,the_opt);
			} else {
				return new Cut_Part(new_pro,pp,oc,&Factory,the_opt);
			}
		}
		else
			if ( settings::general::s_use_cached_integrals  ) {
				return new CachedIntegral::Unknown_Cut_Part_wCI(new_pro,pp,&Factory,the_opt);
			} else {
				return new Cut_Part(new_pro,pp,&Factory,the_opt);
			}
}


template <class CUTD_FACTORY> Cut_Part_base* custom_IRC_cut_part_factory<CUTD_FACTORY>::new_cut_part(const process& pro,color_structure cs){

	custom_cut_part_factory<CUTD_FACTORY> CPF;

		return new IR_checked_Cut_Part(CPF.new_cut_part(pro,cs));



}

#endif

#ifndef BH_PUBLIC
cut_part_factory_Darren global_CPFD;
cut_part_factory_FHZ global_CPFF;
cut_part_factory_higgs global_CPFH;
#endif
cut_part_factory_worker global_CPFW;



template <> cut_part_factory<Cut_Part_base>* cut_part_factory<Cut_Part_base>::s_default_cut_part_factory(const process& pro){

	process::const_iterator it = find_if(pro.begin(),pro.end(),is_of_type(higgs));

	bool has_a_higgs = ( it != pro.end() ) ;
//	bool has_a_higgs = true;
#ifndef BH_PUBLIC
	if ( has_a_higgs ) {
			switch ( settings::general::s_cut_type ){
			case  settings::general::Darren:
				return &global_CPFH;
			case  settings::general::FHZ:
				_WARNING("cut_part_factory<Cut_Part_base>::s_default_higgs_cut_part_factory() : No FHZ factory exists using Darren extended factory");
				return &global_CPFH;
			case  settings::general::worker:
				_WARNING("cut_part_factory<Cut_Part_base>::s_default_higgs_cut_part_factory() : No worker factory exists using Darren extended factory");
				return &global_CPFH;
			}
	} else {
		switch ( settings::general::s_cut_type ){
		case  settings::general::Darren:
			return &global_CPFD;
		case  settings::general::FHZ:
			return &global_CPFF;
		case  settings::general::worker:
			return &global_CPFW;
		}
	}
#else
	return &global_CPFW;
#endif
};


#if _USE_GCC
#ifdef BH_USE_OLD_CPF


template class custom_cut_part_factory<Darren_factory_new>;
template class custom_cut_part_factory<FHZ_factory>;

template class custom_IRC_cut_part_factory<Darren_factory_new>;
template class custom_IRC_cut_part_factory<FHZ_factory>;
#endif

#endif

#if _USE_PGCC

template class custom_cut_part_factory<Darren_factory_new>;
template class custom_cut_part_factory<FHZ_factory>;

template class custom_IRC_cut_part_factory<Darren_factory_new>;
template class custom_IRC_cut_part_factory<FHZ_factory>;

template <class T> bool custom_cut_part_factory<T>::_use_known_formulae=true;
template <class T> bool custom_cut_part_factory<T>::s_use_cached_integrals=true;


#endif


}
