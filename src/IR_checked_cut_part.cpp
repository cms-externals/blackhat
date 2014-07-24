/*
 * IR_checked_cut_part.cpp
 *
 *  Created on: 2 Jul 2009
 *      Author: daniel
 */
#include "IR_checked_cut_part.h"
#include "polylog.h"  //  for pi
#include "BH_debug.h"
#include <limits>
using namespace std;

#define _VERBOSE 0
namespace BH {

/* obsolete old way. new function can also deal with analytics.
// return true if test is passed
template <class T> bool IR_check_Cut_Part(Cut_Part_base* cpb,momentum_configuration<T>& mc,const vector<int>& ind,T tolerance){
	Cut_Part* CP=dynamic_cast<Cut_Part*>(cpb);
	if (CP){
		complex<T> Atree=CP->get_tree(mc,ind);
		complex<T> res;
#if _VERBOSE && 0
_PRINT(mc.p(1));
_PRINT(mc.p(2));
_PRINT(mc.p(3));
_PRINT(mc.p(4));
cout<<setprecision(16);
for(int kk=1;kk<=7;kk++){
	cout<<real(mc.p(kk).E())<<" "<<real(mc.p(kk).X())<<" "<<real(mc.p(kk).Y())<<" "<<real(mc.p(kk).Z())<<" ";
}
cout<<endl;
//_PRINT(mc.m2(7));
//_PRINT(*(CP->bubble(1)));
cout<<"{ "<<ind[0]<<" , "<<ind[1]<<" , "<<ind[2]<<" , "<<ind[3]<<" , "<<ind[4]<<" , "<<ind[5]<<" , "<<ind[6]<<" }\n";
_PRINT(CP->get_process());
_PRINT(CP->get_tree(mc,ind));
_PRINT(CP->eval(mc,ind)/CP->get_tree(mc,ind));
#endif
		for (int i=1;i<=CP->nbr_bubbles();i++){
			res+=CP->bubble(i)->get_value(mc,ind);
		};
		res/=Atree;
		if ( abs(res.imag()) > tolerance){
#if _VERBOSE
		_MESSAGE4("IR in  IR_check_Cut_Part failed: ",abs(res.imag())," > " ,tolerance )	;
#endif
		return false;
	}
	else{
#if _VERBOSE
		_MESSAGE4("IR in  IR_check_Cut_Part passed: ",abs(res.imag())," < " ,tolerance )	;
#endif
	return true;
	};
	} else return true ; // this happens if the cut part is known, in which case the big IR test will be satisfied

}

// return true if test is passed
template <class T> bool IR_check_Cut_Part(Cut_Part_base* cpb,const eval_param<T>&ep ,T tolerance){
	Cut_Part* CP=dynamic_cast<Cut_Part*>(cpb);
	if (CP){
		complex<T> Atree=CP->get_tree(ep);
		complex<T> res;
		for (int i=1;i<=CP->nbr_bubbles();i++){
			res+=CP->bubble(i)->get_value(ep);
		};
		res/=Atree;
		if ( abs(res.imag()) > tolerance){
#if _VERBOSE
		_MESSAGE4("IR in  IR_check_Cut_Part failed: ",abs(res.imag())," > " ,tolerance )	;
#endif
		return false;
	}
	else{
#if _VERBOSE
		_MESSAGE4("IR in  IR_check_Cut_Part passed: ",abs(res.imag())," < " ,tolerance )	;
#endif
	return true;
	};
	} else return true ; // this happens if the cut part is known, in which case the big IR test will be satisfied

}
*/

// return true if test is passed
template <class T> bool IR_check_Cut_Part(complex<T> Atree,complex<T> single_pole, T tolerance){

  BH_DEBUG_MESSAGE2("Atree: ",Atree);  
  T imag_ratio=abs((single_pole/Atree/complex<T>(pi<T>(),0)).imag())+tolerance;
	// we assume here that the the momentum configuration has been rescaled to order one momenta
	// this variable will then be used to ensure that the tree did not blow up
	T inverse_tree=1./abs(Atree);
	if ( abs(to_double(imag_ratio) )> std::numeric_limits<int>::max() ){
	  return false;
	}
	int int_part=int(to_double(imag_ratio));

	return ((0<imag_ratio-int_part)&&(imag_ratio-int_part<T(2)*tolerance)&&(inverse_tree>tolerance/T(1000)));
}



Cut_Part_base* IR_checked_cut_part_factory::new_cut_part(const process& PRO,color_structure cs){
	Cut_Part_base* ptr = cut_part_factory<Cut_Part_base>::s_default_cut_part_factory(PRO)->new_cut_part(PRO, cs);
	return new IR_checked_Cut_Part(ptr);
};

IR_checked_cut_part_factory global_IRCCPF;
IR_checked_cut_part_factory* IR_checked_cut_part_factory::s_default_IR_checked_cut_part_factory=&global_IRCCPF;



template<> inline R IR_checked_Cut_Part::IR_tolerance(){return _IR_tolerance_R;}
template<> inline RHP IR_checked_Cut_Part::IR_tolerance(){return _IR_tolerance_RHP;}
template<> inline RVHP IR_checked_Cut_Part::IR_tolerance(){return _IR_tolerance_RVHP;}

SeriesC<R> IR_checked_Cut_Part::eval(mom_conf& mc,const vector<int>& ind){
	d_cut_part_p->set_mu_index<R>(this->get_mu_index<R>());
	if (this->get_mu_index<R>() == 0){
		d_cut_part_p->set_mu_index<R>( DefineMu<R>(mc,m_MU));
	}
	complex<R> Atree=d_cut_part_p->get_tree(mc,ind);
	SeriesC<R> cut_result=d_cut_part_p->eval(mc,ind);
	
	multi_precision_reader* mpr=dynamic_cast<multi_precision_reader*>(&mc);

	if (!IR_check_Cut_Part(Atree,cut_result[-1],_IR_tolerance_R) ) {
		if ( mpr ==0 ){
//			_WARNING("pointer is not a multi_precision_reader, we improvise");
//			throw BHerror("IR_checked_Cut_Part::eval called with an argument that is not a multi_precision_reader");
			SeriesC<RHP> HP_result;
			vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
			mom_conf_HP mcHP=mc.extend<RHP>(ind);
			// if d_index_mu is not set, we use the values from the multi_precision_constant,
			// otherwise, we construct a mu vector in the extended mc corresponding to the energy component of the vector d_mu+index in the NP mc.
			int new_mu_index_HP;
			int old_mu_index_HP=this->get_mu_index<RHP>();
//_PRINT(this->get_mu_index<RHP>());
//_PRINT(mc.mom(this->get_mu_index<R>()).E().real());
//_PRINT(this->get_mu());
//			if (this->get_mu_index<RHP>() == 0) {
//				new_mu_index_HP = DefineMu<RHP>(mcHP,this->get_mu());
//			} else {
//				new_mu_index_HP = DefineMu<RHP>(mcHP,RHP(mc.mom(this->get_mu_index<R>()).E().real()));
//			}
			new_mu_index_HP = DefineMu<RHP>(mcHP,RHP(mc.mom(this->get_mu_index<R>()).E().real()));
			this->set_mu_HP(new_mu_index_HP);
			HP_result=eval(mcHP,new_ind);
			cut_result= to_double(HP_result);
			this->set_mu_HP(old_mu_index_HP);
		} else {
			_MESSAGE("IR for cut failed, increasing precision to double double...!");
			SeriesC<RHP> cut_result_HP =eval(mpr->mc_HP(),ind);
			cut_result= to_double(cut_result_HP);
		}
//		_PRINT(cut_result);
	}
	
	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result;
}

SeriesC<RHP> IR_checked_Cut_Part::eval(mom_conf_HP& mc,const vector<int>& ind){
	d_cut_part_p->set_mu_index<RHP>(this->get_mu_index<RHP>());
	if (this->get_mu_index<RHP>() == 0){
		d_cut_part_p->set_mu_index<RHP>( DefineMu<RHP>(mc,m_MU));
	}
	complex<RHP> Atree_HP=d_cut_part_p->get_tree(mc,ind);
	SeriesC<RHP> cut_result_HP=d_cut_part_p->eval(mc,ind);

	if (!IR_check_Cut_Part(Atree_HP,cut_result_HP[-1],_IR_tolerance_RHP) ) {
#if _VERBOSE
		_MESSAGE("IR for cut failed even with HP, increasing precision to VHP...!");
#endif
		multi_precision_reader_HP* mpr=dynamic_cast<multi_precision_reader_HP*>(&mc);
		if ( mpr ==0 ){
//			_WARNING("pointer is not a multi_precision_reader_HP, improvise");
			// throw BHerror("IR_checked_Cut_Part::eval called with an argument that is not a multi_precision_reader");
			SeriesC<RVHP> VHP_result;
			vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
			mom_conf_VHP mcVHP=mc.extend<RVHP>(ind);
			// if d_index_mu is not set, we use the values from the multi_precision_constant,
			// otherwise, we construct a mu vector in the extended mc corresponding to the energy component of the vector d_mu+index in the NP mc.
			int new_mu_index_VHP;
			int old_mu_index_VHP=this->get_mu_index<RVHP>();
//			if (this->get_mu_index<RVHP>() == 0) {
//				new_mu_index_VHP = DefineMu<RVHP>(mcVHP,this->get_mu());
//			} else {
//				new_mu_index_VHP = DefineMu<RVHP>(mcVHP,RVHP(mc.mom(this->get_mu_index<RHP>()).E().real()));
//			}
			new_mu_index_VHP = DefineMu<RVHP>(mcVHP,RVHP(mc.mom(this->get_mu_index<RHP>()).E().real()));
			this->set_mu_VHP(new_mu_index_VHP);
			VHP_result=eval(mcVHP,new_ind);
			cut_result_HP = to_HP(VHP_result);
			this->set_mu_VHP(old_mu_index_VHP);
		} else {
			SeriesC<RVHP> cut_result_VHP =eval(mpr->mc_VHP(),ind);
			cut_result_HP= to_HP(cut_result_VHP);
		}
	}
	
	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result_HP;
}

SeriesC<RVHP> IR_checked_Cut_Part::eval(mom_conf_VHP& mc,const vector<int>& ind){


		d_cut_part_p->set_mu_index<RVHP>(this->get_mu_index<RVHP>());
		if (this->get_mu_index<RVHP>() == 0){
			d_cut_part_p->set_mu_index<RVHP>( DefineMu<RVHP>(mc,m_MU));
		}
		complex<RVHP> Atree_VHP=d_cut_part_p->get_tree(mc,ind);
		SeriesC<RVHP> cut_result_VHP =d_cut_part_p->eval(mc,ind);
		if(!IR_check_Cut_Part(Atree_VHP,cut_result_VHP[-1],_IR_tolerance_RVHP) ) {
#if _VERBOSE
			_MESSAGE("IR for cut failed even with VHP, no way to do better, sorry...!");
#endif
			}
//		_MESSAGE("done");
	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result_VHP;
}

#ifdef BH_USE_GMP
SeriesC<RGMP> IR_checked_Cut_Part::eval(momentum_configuration<RGMP>& mc,const vector<int>& ind){


		d_cut_part_p->set_mu_index<RGMP>(this->get_mu_index<RGMP>());
		if (this->get_mu_index<RGMP>() == 0){
			d_cut_part_p->set_mu_index<RGMP>( DefineMu<RGMP>(mc,m_MU));
		}
		complex<RGMP> Atree_GMP=d_cut_part_p->get_tree(mc,ind);
		SeriesC<RGMP> cut_result_GMP =d_cut_part_p->eval(mc,ind);
		if(!IR_check_Cut_Part(Atree_GMP,cut_result_GMP[-1],_IR_tolerance_RGMP) ) {
#if _VERBOSE
			_MESSAGE("IR for cut failed even with VHP, no way to do better, sorry...!");
#endif
			}
//		_MESSAGE("done");
	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result_GMP;
}

#endif

SeriesC<R> IR_checked_Cut_Part::eval(const eval_param<R>& ep){
	//We need to pass the chosen scale onto the cut_part that does the computation
	multi_precision_constant mus(get_mu());
	d_cut_part_p->set_mu(mus);
	complex<R> Atree=d_cut_part_p->get_tree(ep);
	SeriesC<R> cut_result=d_cut_part_p->eval(ep);

	if (!IR_check_Cut_Part(Atree,cut_result[-1],_IR_tolerance_R) ) {
		SeriesC<RHP> HP_result;
		vector<Cmom<RHP> > *new_mom=extend_momenta<R,RHP>(ep);
		eval_param<RHP> epHP(*new_mom);
		HP_result=eval(epHP);
		cut_result= to_double(HP_result);
		delete new_mom;
	}
	
	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result;
}

SeriesC<RHP> IR_checked_Cut_Part::eval(const eval_param<RHP>& ep){
	//We need to pass the chosen scale onto the cut_part that does the computation
	multi_precision_constant mus(get_mu());
	d_cut_part_p->set_mu(mus);
	complex<RHP> Atree_HP=d_cut_part_p->get_tree(ep);
	SeriesC<RHP> cut_result_HP=d_cut_part_p->eval(ep);

	if (!IR_check_Cut_Part(Atree_HP,cut_result_HP[-1],_IR_tolerance_RHP) ) {
		SeriesC<RVHP> VHP_result;
		vector<Cmom<RVHP> > *new_mom=extend_momenta<RHP,RVHP>(ep);
		eval_param<RVHP> epVHP(*new_mom);
		VHP_result=eval(epVHP);
		cut_result_HP= to_HP(VHP_result);
		delete new_mom;
	}
	
	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result_HP;
}

SeriesC<RVHP> IR_checked_Cut_Part::eval(const eval_param<RVHP>& ep){
	//We need to pass the chosen scale onto the cut_part that does the computation
	multi_precision_constant mus(get_mu());
	d_cut_part_p->set_mu(mus);
	complex<RVHP> Atree_VHP=d_cut_part_p->get_tree(ep);
	SeriesC<RVHP> cut_result_VHP =d_cut_part_p->eval(ep);
	if (!IR_check_Cut_Part(Atree_VHP,cut_result_VHP[-1],_IR_tolerance_RVHP) ) {
#if _VERBOSE
			_MESSAGE("IR for cut failed even with VHP, no way to do better, sorry...!");
#endif
	}
	
	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result_VHP;
}

#ifdef BH_USE_GMP
SeriesC<RGMP> IR_checked_Cut_Part::eval(const eval_param<RGMP>& ep){
	//We need to pass the chosen scale onto the cut_part that does the computation
	multi_precision_constant mus(get_mu());
	d_cut_part_p->set_mu(mus);
	complex<RGMP> Atree_GMP=d_cut_part_p->get_tree(ep);
	SeriesC<RGMP> cut_result_GMP =d_cut_part_p->eval(ep);
	if (!IR_check_Cut_Part(Atree_GMP,cut_result_GMP[-1],_IR_tolerance_RGMP) ) {
#if _VERBOSE
			_MESSAGE("IR for cut failed even with GMP, no way to do better, sorry...!");
#endif
	}

	//_accuracy=d_cut_part_p->get_accuracy();
	return cut_result_GMP;
}

#endif

}
