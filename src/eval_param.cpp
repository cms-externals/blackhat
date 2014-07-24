/*
 * eval_param.cpp
 *
 *  Created on: Apr 9, 2009
 *      Author: Darren Forde
 */

#include "eval_param.h"
#include "BH_debug.h"
#if BH_USE_GMP
#include "gmp_r.h"
#endif

#define _VERBOSE 0

using std::complex;

namespace BH {

// Initialize static variables

// We choose a momenta that does not require that qd is initialized and so we do not loose precision when the static variable is initialised
//template <class T> const Cmom<T> eval_param<T>::_ep_quark_ref(momentum<complex<T> >(complex<T>(1,1),complex<T>(0,1),complex<T>(1,1),complex<T>(1,0)),_mt_unknown);
template <class T> const Cmom<T> eval_param<T>::_ep_quark_ref(momentum<complex<T> >(sqrt(complex<T>(-2,2)),complex<T>(0,1),complex<T>(1,1),complex<T>(0,1)),_mt_massless);
//template <class T> const Cmom<T> eval_param<T>::_ep_quark_ref(momentum<complex<T> >(complex<T>(1,0),complex<T>(0,0),complex<T>(0,0),complex<T>(1,0)),_mt_massless);
//template <class T> const Cmom<T> eval_param<T>::_ep_quark_ref(momentum<complex<T> >(sqrt(complex<T>(3,0)),complex<T>(1,0),complex<T>(1,0),complex<T>(1,0)),_mt_massless);
#ifdef BH_USE_GMP
template <> eval_param<RGMP>::eval_param(int n){
	_size=n; _em=new const Cmom<RGMP>*[n];  _epstate=new eval_param_state(_size);
	static Cmom<RGMP> better_ep_quark_ref;
	better_ep_quark_ref=Cmom<RGMP>(momentum<complex<RGMP> >(sqrt(complex<RGMP>(-2,2)),complex<RGMP>(0,1),complex<RGMP>(1,1),complex<RGMP>(0,1)),_mt_massless);
	_ref_vector=&better_ep_quark_ref;
	BH_DEBUG_MESSAGE("Using new constructor");
}

template <> eval_param<RGMP>::eval_param(const vector<Cmom<RGMP> >& vec) {
	_size=vec.size();
	_em=new const Cmom<RGMP>*[vec.size()];
	for(int i=0;i<_size;i++){_em[i]=&vec[i];} ;
	static Cmom<RGMP> better_ep_quark_ref;
	better_ep_quark_ref=Cmom<RGMP>(momentum<complex<RGMP> >(sqrt(complex<RGMP>(-2,2)),complex<RGMP>(0,1),complex<RGMP>(1,1),complex<RGMP>(0,1)),_mt_massless);
	_ref_vector=&better_ep_quark_ref;
	_epstate=new eval_param_state(_size);
};
template <> eval_param<RGMP>::eval_param(momentum_configuration<RGMP>& mc, std::vector<int>& ind) {
	_size=ind.size();
	_em=new const Cmom<RGMP>*[ind.size()];
	for(int i=0;i<_size;i++){_em[i]=(&mc[ind[i]]);};
	static Cmom<RGMP> better_ep_quark_ref;
	better_ep_quark_ref=Cmom<RGMP>(momentum<complex<RGMP> >(sqrt(complex<RGMP>(-2,2)),complex<RGMP>(0,1),complex<RGMP>(1,1),complex<RGMP>(0,1)),_mt_massless);
	_ref_vector=&better_ep_quark_ref;
	_epstate=new eval_param_state(_size);
};
template <> eval_param<RGMP>::eval_param(momentum_configuration<RGMP>& mc, const std::vector<int>& ind) {
	_size=ind.size();
	_em=new const Cmom<RGMP>*[ind.size()];
	for(int i=0;i<_size;i++){_em[i]=(&mc[ind[i]]);};
	static Cmom<RGMP> better_ep_quark_ref;
	better_ep_quark_ref=Cmom<RGMP>(momentum<complex<RGMP> >(sqrt(complex<RGMP>(-2,2)),complex<RGMP>(0,1),complex<RGMP>(1,1),complex<RGMP>(0,1)),_mt_massless);
	_ref_vector=&better_ep_quark_ref;
	_epstate=new eval_param_state(_size);
};

#endif

template <class T> eval_param<T>::eval_param(int n) {
	_size=n;
	_em=new const Cmom<T>*[n];
	_ref_vector=&_ep_quark_ref;
	_epstate=new eval_param_state(_size);
};

template <class T> eval_param<T>::eval_param(const vector<Cmom<T> >& vec) {
	_size=vec.size();
	_em=new const Cmom<T>*[vec.size()];
	for(int i=0;i<_size;i++){_em[i]=&vec[i];} ;
	_ref_vector=&_ep_quark_ref; _epstate=new eval_param_state(_size);
};
template <class T> eval_param<T>::eval_param(momentum_configuration<T>& mc, std::vector<int>& ind) {
	_size=ind.size();
	_em=new const Cmom<T>*[ind.size()];
	for(int i=0;i<_size;i++){_em[i]=(&mc[ind[i]]);};
	_ref_vector=&_ep_quark_ref;
	_epstate=new eval_param_state(_size);
};
template <class T> eval_param<T>::eval_param(momentum_configuration<T>& mc, const std::vector<int>& ind) {
	_size=ind.size();
	_em=new const Cmom<T>*[ind.size()];
	for(int i=0;i<_size;i++){_em[i]=(&mc[ind[i]]);};
	_ref_vector=&_ep_quark_ref;
	_epstate=new eval_param_state(_size);
};


template <class T> void eval_param<T>::update(momentum_configuration<T>& mc, const std::vector<int>& ind){
        delete[] _em; delete _epstate;
        _size=ind.size(); 
        _em=new const Cmom<T>*[ind.size()]; 
        for(int i=0;i<_size;i++){_em[i]=(&mc[ind[i]]);} 
        _ref_vector=&_ep_quark_ref; 
        _epstate=new eval_param_state(_size);};


unsigned long int  eval_param_state::eval_param_next_ID(1);

int mass_param::_next_mass_label(0);

mass_param::mass_param()
{
	_label=0;
}

mass_param::mass_param(const C mass, const CHP mass_HP, const CVHP mass_VHP) :
	_original(mass),
	_original_sqr(mass*mass),
	_actual(mass),
	_actual_sqr(mass*mass),
	_actual_HP(mass_HP),
	_actual_sqr_HP(mass_HP*mass_HP),
	_original_HP(mass_HP),
	_original_sqr_HP(mass_HP*mass_HP),
	_actual_VHP(mass_VHP),
	_actual_sqr_VHP(mass_VHP*mass_VHP),
	_original_VHP(mass_VHP),
	_original_sqr_VHP(mass_VHP*mass_VHP)
{
	_label=++mass_param::_next_mass_label;
}

mass_param::mass_param(const mass_param& mp) :
	_original(mp._original),
	_original_sqr(mp._original_sqr),
	_actual(mp._actual),
	_actual_sqr(mp._actual_sqr),
	_actual_HP(mp._actual_HP),
	_actual_sqr_HP(mp._actual_sqr_HP),
	_original_HP(mp._original_HP),
	_original_sqr_HP(mp._original_sqr_HP),
	_actual_VHP(mp._actual_VHP),
	_actual_sqr_VHP(mp._actual_sqr_HP),
	_original_VHP(mp._original_VHP),
	_original_sqr_VHP(mp._original_sqr_HP)
{
	_label=mp._label;
}

mass_param& mass_param::operator=(const mass_param& mp){
	_original=mp._original;
	_original_sqr=mp._original_sqr;
	_actual=mp._actual;
	_actual_sqr=mp._actual_sqr;
	_actual_HP=mp._actual_HP;
	_actual_sqr_HP=mp._actual_sqr_HP;
	_original_HP=mp._original_HP;
	_original_sqr_HP=mp._original_sqr_HP;
	_actual_VHP=mp._actual_VHP;
	_actual_sqr_VHP=mp._actual_sqr_HP;
	_original_VHP=mp._original_VHP;
	_original_sqr_VHP=mp._original_sqr_HP;
	_label=mp._label;

	return *this;
}

mass_param_library::mass_param_library(mass_param& mp){
	// If the label we are asked for is actually not the next in order we pad the list
	// out until we get to it.
	for(int i=_mass_params.size();i<mp.label();i++){
		mass_param empty;
		_mass_params.push_back(empty);
	}
	_mass_params.push_back(mp);
}


void mass_param_library::add(mass_param& mp){
	//We insert the mass at the location stored in the mass label of the particle_ID
	// If this is below the current max we just insert it
	if(_mass_params.size()>mp.label()){
		_mass_params[mp.label()-1]=mp;
	}else{
		// If the label we are asked for is actually not the next in order we pad the list
		// out until we get to it.
		for(int i=_mass_params.size();i<mp.label();i++){
			mass_param empty;
			_mass_params.push_back(empty);
		}
		_mass_params.push_back(mp);
	}
}

#ifdef BH_USE_GMP
const CGMP& mass_param::getGMPvalue(CGMP& value) const {
	if (value.real().get_precision()<RGMP::get_current_precision()){
		value=CGMP(RGMP(value.real(),RGMP::get_current_precision()),RGMP(value.imag(),RGMP::get_current_precision()));
	}
	return value;
}

#endif
// Adapted from Fernando and harald's extend for the mom_conf in extend.cpp
template <class T_low,class T_high> std::vector<Cmom<T_high> >* extend_momenta(const eval_param<T_low>& ep)
{
	std::vector<Cmom<T_high> > *momenta=new std::vector<Cmom<T_high> >;
	momentum<complex<T_high> > sum;
	size_t n=ep.size();
	for (size_t i=0;i<n-2;i++){
		momenta->push_back(Cmom<T_high>(lambdat<T_high>(ep.p(i)->Lt()),lambda<T_high>(ep.p(i)->L())));
		sum+=momenta->back().P();
	}

	lambdat<T_high> Z1(complex<T_high>(1,0),complex<T_high>(0,0));
	lambdat<T_high> Z2(complex<T_high>(0,0),complex<T_high>(1,0));
	lambda<T_high> Z3(complex<T_high>(1,0),complex<T_high>(0,0));
	lambda<T_high> Z4(complex<T_high>(0,0),complex<T_high>(1,0));

	// completed k_{n-1}
	Cmom<T_high> knm1c(lambdat<T_high>(ep.p(n-2)->Lt()),lambda<T_high>(ep.p(n-2)->L()));
	// notice that the sqrt is always taken in the vicinity of 1, away of its branch cut
	complex<T_high> sqrtt=sqrt(-sum*sum/(complex<T_high>(2,0)*sum*knm1c.P()));
	// final extended k_{n-1}
	Cmom<T_high> knm1f(sqrtt*lambdat<T_high>(ep.p(n-2)->Lt()),sqrtt*lambda<T_high>(ep.p(n-2)->L()));

	sum+=knm1f.P();
	// extended k_{n}, just by momentum conservation
	Cmom<T_high> knf=(-sum);
	// to make sure we didn't hit a branch cut
	// completed k_{n}
	Cmom<T_high> knc(lambdat<T_high>(ep.p(n-1)->Lt()),lambda<T_high>(ep.p(n-1)->L()));

	complex<T_high> phase(1,0);
	complex<T_high> phaseinv(1,0);

	if(abs((knf.L()-knc.L())*Z1)>T_high(1)/T_high(100000) ||
		abs((knf.L()-knc.L())*Z2)>T_high(1)/T_high(100000) ||
		abs(Z3*(knf.Lt()-knc.Lt()))>T_high(1)/T_high(100000) ||
		abs(Z4*(knf.Lt()-knc.Lt()))>T_high(1)/T_high(100000)
	){
#if _VERBOSE
		std::cout<<"\nfixing phases in extended momentum --- extend of complex momenta\n";
#endif
		complex<T_high> knca=(knc.L()*Z1+knc.L()*Z2)/T_high(2);
		complex<T_high> knfa=(knf.L()*Z1+knf.L()*Z2)/T_high(2);
		phase=knca/knfa;
		phaseinv=knfa/knca;

#if _VERBOSE
		_PRINT(knc.L());
		_PRINT(knc.Lt());
		_PRINT(knf.L());
		_PRINT(knf.Lt());
		_PRINT(phase*knf.L());
		_PRINT(phaseinv*knf.Lt());
		std::cout<<std::endl;
#endif
	}

	momenta->push_back(Cmom<T_high>(sqrtt*lambdat<T_high>(ep.p(n-2)->Lt()),sqrtt*lambda<T_high>(ep.p(n-2)->L())));
	momenta->push_back(Cmom<T_high>(phaseinv*knf.Lt(),phase*knf.L()));


	return momenta;
}


/*
 *
 *
 * cout functions
 *
 *
 */

template <class T> std::ostream& operator<<(std::ostream& s,const eval_param<T>& ep){
	s << "eval_param size: " << ep.size() << " : {";
	if(ep.size()!=0){
		s << *(ep.p(0));
	}
	for(int i=1;i<ep._size;i++){
		s << "," << std::endl << *(ep.p(i));
	}
	s << "}" << std::endl;
//	s << "ref_mom:" << *ep.ref() << std::endl;
//	s << "masses:" << ep._masses;
	return s;
}

std::ostream& operator<<(std::ostream& s,const mass_param_coll& mpc){
	s << mpc._size << " mass_params in collection:{";
	for(int i=0;i<mpc._size;i++){
		s << mpc._em[i] << ",";
	}
	s << "}";
	return s;
}

std::ostream& operator<<(std::ostream& s,const mass_param_library& mpl){
	s << mpl._mass_params.size() << " mass_params in library :{" << std::endl;
	for(std::vector<mass_param>::const_iterator it=mpl._mass_params.begin();it!=mpl._mass_params.end();it++){
		s << *it << std::endl;
	}
	s << "}";
	return s;
}

std::ostream& operator<<(std::ostream& s,const mass_param& mp){
	s << "mass_param label:" << mp._label << " mass:" << mp._actual << " mass^2:" << mp._actual_sqr;
	return s;
}

/*
 *
 *
 * eval_param_state eval functions
 *
 *
 */

bool operator==(const eval_param_state& eps1, const eval_param_state& eps2){
	if(eps1._ep_ID!=eps2._ep_ID){return false;}
	if(eps1._size!=eps2._size){return false;}

	for(int i=0;i<eps1._size;i++){
		if(eps1._eps_legs[i]!=eps2._eps_legs[i]) {return false;};
	}

	return true;
}

bool operator!=(const eval_param_state& eps1, const eval_param_state& eps2){
	if(eps1._ep_ID!=eps2._ep_ID){return true;}
	if(eps1._size!=eps2._size){return true;}

	for(int i=0;i<eps1._size;i++){
		if(eps1._eps_legs[i]!=eps2._eps_legs[i]) {return true;};
	}

	return false;
}

bool eval_param_state::is_match(eval_param_state* eps){
	return eps->get_full_state()==_full_state;
}

void eval_param_state::toggle_state() {
	_state^=1; // Use _state xor 1 to flip the bit to 1 or to 0

	// Now that we have changed the state we need to update the total state check number
	_full_state=0;
	for(int i=0;i<_size;i++){
		_full_state+=(_eps_legs[i]->get_state()<<i);
	}
}

// Explicit Instantiation

template std::vector<Cmom<RHP> >* extend_momenta<R,RHP>(const eval_param<R>& ep);
template std::vector<Cmom<RVHP> >* extend_momenta<R,RVHP>(const eval_param<R>& ep);
template std::vector<Cmom<RVHP> >* extend_momenta<RHP,RVHP>(const eval_param<RHP>& ep);

template const Cmom<R> eval_param<R>::_ep_quark_ref;
template const Cmom<RHP> eval_param<RHP>::_ep_quark_ref;
template const Cmom<RVHP> eval_param<RVHP>::_ep_quark_ref;


template void eval_param<R>::update(momentum_configuration<R>& mc, const std::vector<int>& ind);
template void eval_param<RHP>::update(momentum_configuration<RHP>& mc, const std::vector<int>& ind);
template void eval_param<RVHP>::update(momentum_configuration<RVHP>& mc, const std::vector<int>& ind);



template std::ostream& operator<<(std::ostream& s,const eval_param<R>& ep);
template std::ostream& operator<<(std::ostream& s,const eval_param<RHP>& ep);
template std::ostream& operator<<(std::ostream& s,const eval_param<RVHP>& ep);


template class eval_param<double>;
template class eval_param<dd_real>;
template class eval_param<qd_real>;

#if BH_USE_GMP

template const Cmom<RGMP> eval_param<RGMP>::_ep_quark_ref;
template std::ostream& operator<<(std::ostream& s,const eval_param<RGMP>& ep);

#endif



}
