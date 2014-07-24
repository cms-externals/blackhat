/*
 * standard_cut_part.cpp
 *
 *  Created on: Jun 29, 2009
 *      Author: daniel
 */


#define _VERBOSE 0 // 0 is silent, 1 gives some information



#include "standard_cut_part.hpp"

namespace BH {
namespace cut {

#if BH_USE_GMP

void extendPrecisionComplex(CGMP& c){
	c.real().set_prec(RGMP::get_current_precision());
	c.imag().set_prec(RGMP::get_current_precision());
};


class givemeaccess {
	public:
	CGMP a1;
	CGMP a2;
};

template <class SP> SP extendPrecision(const SP& sp){
	SP newsp(sp);
	givemeaccess& s=reinterpret_cast<givemeaccess&>(newsp);
	extendPrecisionComplex(s.a1);
	extendPrecisionComplex(s.a2);
	return newsp;
}

std::vector<Cmom<RGMP> >* extend_momentaGMP(const eval_param<RGMP>& ep)
{
	std::vector<Cmom<RGMP> > *momenta=new std::vector<Cmom<RGMP> >;
	momentum<complex<RGMP> > sum;
	size_t n=ep.size();
	for (size_t i=0;i<n-2;i++){
		lambdat<RGMP> Lt=extendPrecision(ep.p(i)->Lt());
		lambda<RGMP> L=extendPrecision(ep.p(i)->L());
		momenta->push_back(Cmom<RGMP>(Lt,L));
		sum+=momenta->back().P();
	}

	lambdat<RGMP> Z1(complex<RGMP>(1,0),complex<RGMP>(0,0));
	lambdat<RGMP> Z2(complex<RGMP>(0,0),complex<RGMP>(1,0));
	lambda<RGMP> Z3(complex<RGMP>(1,0),complex<RGMP>(0,0));
	lambda<RGMP> Z4(complex<RGMP>(0,0),complex<RGMP>(1,0));

	// completed k_{n-1}
	Cmom<RGMP> knm1c(extendPrecision(ep.p(n-2)->Lt()),extendPrecision(ep.p(n-2)->L()));
	// notice that the sqrt is always taken in the vicinity of 1, away of its branch cut
	complex<RGMP> sqrtt=sqrt(-sum*sum/(complex<RGMP>(2,0)*sum*knm1c.P()));
	// final extended k_{n-1}
	Cmom<RGMP> knm1f(sqrtt*extendPrecision(ep.p(n-2)->Lt()),sqrtt*extendPrecision(ep.p(n-2)->L()));

	sum+=knm1f.P();
	// extended k_{n}, just by momentum conservation
	Cmom<RGMP> knf=(-sum);
	// to make sure we didn't hit a branch cut
	// completed k_{n}
	Cmom<RGMP> knc(extendPrecision(ep.p(n-1)->Lt()),extendPrecision(ep.p(n-1)->L()));

	complex<RGMP> phase(1,0);
	complex<RGMP> phaseinv(1,0);

	if(abs((knf.L()-knc.L())*Z1)>RGMP(1)/RGMP(100000) ||
		abs((knf.L()-knc.L())*Z2)>RGMP(1)/RGMP(100000) ||
		abs(Z3*(knf.Lt()-knc.Lt()))>RGMP(1)/RGMP(100000) ||
		abs(Z4*(knf.Lt()-knc.Lt()))>RGMP(1)/RGMP(100000)
	){
#if _VERBOSE
		std::cout<<"\nfixing phases in extended momentum --- extend of complex momenta\n";
#endif
		complex<RGMP> knca=(knc.L()*Z1+knc.L()*Z2)/RGMP(2);
		complex<RGMP> knfa=(knf.L()*Z1+knf.L()*Z2)/RGMP(2);
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

	momenta->push_back(Cmom<RGMP>(sqrtt*extendPrecision(ep.p(n-2)->Lt()),sqrtt*extendPrecision(ep.p(n-2)->L())));
	momenta->push_back(Cmom<RGMP>(phaseinv*knf.Lt(),phase*knf.L()));


	return momenta;
}
#endif
/*
spinor<R> to_double(const spinor<RGMP>& sp){
	givemeaccess& s=reinterpret_cast<givemeaccess&>(sp);
	C a1=to_double(s.a1);
	C a2=to_double(s.a2);
	return spinor<R>(a1,a2);
}


void GMPtoDouble(momentum_configuration<RGMP>& mc,momentum_configuration<R>& mcR)
{
	for (int i=0;i<mc.n();i++){
		lambdat<R> Lt(to_double(mc.Lt(i)));
		lambda<R> L=to_double(mc.L(i)));
		mcR.insert(Lt,L);
	}
}
*/



template class standard_cut_part<worker::worker_boxDarren,worker::worker_triangleDarren,worker::worker_bubbleDarren>;
#ifndef BH_PUBLIC
template class standard_cut_part<normal_boxDarren,normal_triangleDarren,normal_bubbleDarren>;
template class standard_cut_part<boxFHZ,triangleFHZ,bubbleFHZ>;
template class standard_cut_part<boxD,triangleD,bubbleD>;
template class standard_cut_part<higgs_boxDarren,higgs_triangleDarren,higgs_bubbleDarren>;
#endif

typedef standard_cut_part<worker::worker_boxDarren,worker::worker_triangleDarren,worker::worker_bubbleDarren> SCP_worker ;
#ifndef BH_PUBLIC
typedef standard_cut_part<normal_boxDarren,normal_triangleDarren,normal_bubbleDarren> SCP_Darren ;
typedef standard_cut_part<boxFHZ,triangleFHZ,bubbleFHZ> SCP_FHZ ;
typedef standard_cut_part<higgs_boxDarren,higgs_triangleDarren,higgs_bubbleDarren> SCP_higgs ;
#endif
template SeriesC<R> SCP_worker::eval_fn(const eval_param<R>& ep);
template SeriesC<RHP> SCP_worker::eval_fn(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_worker::eval_fn(const eval_param<RVHP>& ep);

//template SeriesC<R>  SCP_worker::eval_fn(momentum_configuration<R>& mc,const vector<int>& ind);
//template SeriesC<RHP>  SCP_worker::eval_fn(momentum_configuration<RHP>& mc,const vector<int>& ind);
//template SeriesC<RVHP>  SCP_worker::eval_fn(momentum_configuration<RVHP>& mc,const vector<int>& ind);
#ifndef BH_PUBLIC

template SeriesC<R> SCP_Darren::eval_fn(const eval_param<R>& ep);
template SeriesC<RHP> SCP_Darren::eval_fn(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_Darren::eval_fn(const eval_param<RVHP>& ep);

//template SeriesC<R>  SCP_Darren::eval_fn(momentum_configuration<R>& mc,const vector<int>& ind);
//template SeriesC<RHP>  SCP_Darren::eval_fn(momentum_configuration<RHP>& mc,const vector<int>& ind);
//template SeriesC<RVHP>  SCP_Darren::eval_fn(momentum_configuration<RVHP>& mc,const vector<int>& ind);

template SeriesC<R> SCP_FHZ::eval_fn(const eval_param<R>& ep);
template SeriesC<RHP> SCP_FHZ::eval_fn(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_FHZ::eval_fn(const eval_param<RVHP>& ep);

//template SeriesC<R>  SCP_FHZ::eval_fn(momentum_configuration<R>& mc,const vector<int>& ind);
//template SeriesC<RHP>  SCP_FHZ::eval_fn(momentum_configuration<RHP>& mc,const vector<int>& ind);
//template SeriesC<RVHP>  SCP_FHZ::eval_fn(momentum_configuration<RVHP>& mc,const vector<int>& ind);

template SeriesC<R> SCP_higgs::eval_fn(const eval_param<R>& ep);
template SeriesC<RHP> SCP_higgs::eval_fn(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_higgs::eval_fn(const eval_param<RVHP>& ep);

//template SeriesC<R>  SCP_higgs::eval_fn(momentum_configuration<R>& mc,const vector<int>& ind);
//template SeriesC<RHP>  SCP_higgs::eval_fn(momentum_configuration<RHP>& mc,const vector<int>& ind);
//template SeriesC<RVHP>  SCP_higgs::eval_fn(momentum_configuration<RVHP>& mc,const vector<int>& ind);

#endif /*BH_PUBLIC*/

#if BH_USE_GMP

template SeriesC<RGMP> SCP_worker::eval_fn(const eval_param<RGMP>& ep);
#ifndef BH_PUBLIC
template SeriesC<RGMP> SCP_Darren::eval_fn(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_FHZ::eval_fn(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_higgs::eval_fn(const eval_param<RGMP>& ep);
#endif
#endif

}
}
