#include "BH_error.h"
#include "particles.h"
#include "BH_Ampl.h"
#include "BH_Ampl_processes.h"
#include "mom_conf.h"
#include "integrals.h"
#include "OneLoopHelAmpl.h"
#include "cached_OLHA.h"
#include "constants.h"

/*
	const double Nc=3.;
	const double Nf=5.;
	const double Ns=0.;
*/

using namespace std;

using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {

#if 0
int BH_interface::d_which_study=0;
#endif

#define _VERBOSE 0
#define SPA(i,j) mc.spa(ind[i-1],ind[j-1])
#define SPB(i,j) mc.spb(ind[i-1],ind[j-1])
#define S(i,j) mc.s(ind[i-1],ind[j-1])

/*
class BH_interface_concrete {
	momentum_configuration<double>* d_mc_p;
	settings_table* d_settings_p;
public:
	BH_interface_concrete();
	~BH_interface_concrete();

	//! Return a pointer to a new BH_Ampl object.
	BH_amplitude* new_ampl(const std::vector<int>&);
	//! Sets the named setting to the value provided.
	// The type is deduced from the argument. It is the responsibility of the caller to ensure that the type matches the type of the setting
	template <class T> void set(const std::string& name,T value);
	//! Prints a table of the settings and their current values.
	void print_settings();
	virtual double operator()(BHinput& in);
	momentum_configuration<double>* get_mc(){return d_mc_p;}
};
*/

/*
template <class T> complex<T> AqppqmCut(int eporder, momentum_configuration<T>& mc, const std::vector<int>& ind, int mu)
{
//define corner vectors
 std::vector<int> c1;  c1.push_back(ind[1-1]);
 std::vector<int> c2;  c2.push_back(ind[2-1]);
 std::vector<int> c3;  c3.push_back(ind[3-1]);
 std::vector<int> c4;  c4.push_back(ind[4-1]);
 std::vector<int> c5;  c5.push_back(ind[5-1]);

 std::vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 std::vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
 std::vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
 std::vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
 std::vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);

 std::vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind[i-1]);}
 std::vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind[i-1]);}
 std::vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind[i-1]);}
 std::vector<int> c451;  for(int i = 4; i<=5; i++) {c451.push_back(ind[i-1]);}
                      c451.push_back(ind[1-1]);
 std::vector<int> c145 = c451;
                        // WARNING: Triangle sign fudge to account for
                        // David's backwards sign in the triangles.
                        // Later this will be +1.
 complex<T> trianglesign = complex<T>(1.,0.);
 return( (complex<T>(0.,-1./2.)*(-(pow(SPA(1,3),2)*pow(SPA(4,5),2)*
  pow(SPB(5,1),2)*S(2,3))+pow(SPA(1,3),2)*pow(SPA(4,5),2)*
 pow(SPB(5,1),2)*S(4,5)))/(pow(S(2,3)-S(4,5),2)*
 pow(SPA(4,5),2)*SPA(1,2)*SPA(2,3)*SPB(5,4)) +
 (complex<T>(0.,1./2.)*(trianglesign*complex<T>(-2.,0.)*
 Int(eporder,mc,mu,c4,c5,c123)*
 pow(SPA(2,3),2)*pow(SPA(3,4),2)*pow(SPA(4,5),2)*
 pow(SPB(3,2),2)*pow(SPB(5,4),2)+trianglesign*complex<T>(-2.,0.)*
 Int(eporder,mc,mu,c4,c5,c123)*pow(SPA(3,4),2)*
 pow(SPA(4,5),4)*pow(SPB(5,4),4)-pow(SPA(1,3),2)*
 pow(SPA(4,5),2)*pow(SPB(5,1),2)*S(2,3)+
 trianglesign*complex<T>(4.,0.)*
 Int(eporder,mc,mu,c4,c5,c123)*pow(SPA(3,4),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),3)*S(2,3) +
pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(5,1),2)*S(4,5) +
Int(eporder,mc,mu,c45,c123)*pow(SPA(1,3),2)*pow(SPA(4,5),2)*
 pow(SPB(5,1),2)*S(4,5)+Int(eporder,mc,mu,c3,c2,c1,c45)*
 pow(S(2,3)-S(4,5),2)*pow(SPA(3,4),2)*S(1,2)*S(2,3)*
 S(4,5)+complex<T>(2.,0.)*Int(eporder,mc,mu,c45,c123)*pow(SPA(4,5),3)*
 pow(SPB(5,4),2)*SPA(1,3)*SPA(3,4)*SPB(5,1) +
complex<T>(-2.,0.)*Int(eporder,mc,mu,c45,c123)*S(2,3)*S(4,5)*SPA(1,3)*
 SPA(3,4)*SPA(4,5)*SPB(5,1)+Int(eporder,mc,mu,c23,c145)*
 S(4,5)*(complex<T>(3.,0.)*pow(S(2,3)-S(4,5),2)*pow(SPA(3,4),2) -
  pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(5,1),2) +
  complex<T>(-2.,0.)*(-S(2,3)+S(4,5))*SPA(1,3)*SPA(3,4)*SPA(4,5)*
   SPB(5,1))))/(pow(S(2,3)-S(4,5),2)*pow(SPA(4,5),2)*
 SPA(1,2)*SPA(2,3)*SPB(5,4))

       );
}

template <class T> complex<T> AqppqmCut2(int eporder, momentum_configuration<T>& mc, const std::vector<int>& ind, int mu)
{
//define corner vectors
 std::vector<int> c1;  c1.push_back(ind[1-1]);
 std::vector<int> c2;  c2.push_back(ind[2-1]);
 std::vector<int> c3;  c3.push_back(ind[3-1]);
 std::vector<int> c4;  c4.push_back(ind[4-1]);
 std::vector<int> c5;  c5.push_back(ind[5-1]);

 std::vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 std::vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
 std::vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
 std::vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
 std::vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);

 std::vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind[i-1]);}
 std::vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind[i-1]);}
 std::vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind[i-1]);}
 std::vector<int> c451;  for(int i = 4; i<=5; i++) {c451.push_back(ind[i-1]);}
                      c451.push_back(ind[1-1]);
 std::vector<int> c145 = c451;
                        // WARNING: Triangle sign fudge to account for
                        // David's backwards sign in the triangles.
                        // Later this will be +1.
 complex<T> trianglesign = complex<T>(1.,0.);
 return( (complex<T>(0,0.5)*(Int(eporder,mc,mu,c12,c3,c45))*(S(1,2) - S(4,5))*
	      pow(SPA(3,4),2))/(SPA(1,2)*SPA(2,3)*SPA(4,5)) +
	   (complex<T>(0,0.5)*(Int(eporder,mc,mu,c1,c23,c45))*(S(2,3) - S(4,5))*
	      pow(SPA(3,4),2))/(SPA(1,2)*SPA(2,3)*SPA(4,5)) -
	   (complex<T>(0,0.5)*(Int(eporder,mc,mu,c1,c2,c345))*pow(SPA(3,4),2)*
	      SPB(2,1))/(SPA(2,3)*SPA(4,5)) -
	   (complex<T>(0,0.5)*(Int(eporder,mc,mu,c2,c3,c145))*pow(SPA(3,4),2)*
	      SPB(3,2))/(SPA(1,2)*SPA(4,5)) +
	   (complex<T>(0,0.5)*(Int(eporder,mc,mu,c3,c2,c1,c45))*pow(SPA(3,4),2)*
	      SPB(2,1)*SPB(3,2))/SPA(4,5) +
	   (complex<T>(0,0.5)*(Int(eporder,mc,mu,c23,c145))*
	      (complex<T>(3,0)*pow(S(2,3),2)*pow(SPA(3,4),2) +
	    		  complex<T>(3,0)*pow(S(4,5),2)*pow(SPA(3,4),2) -
	        complex<T>(2,0)*S(4,5)*SPA(1,3)*SPA(3,4)*SPA(4,5)*SPB(5,1) -
	        pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(5,1),2) +
	        complex<T>(2,0)*S(2,3)*SPA(3,4)*(-complex<T>(3,0)*S(4,5)*SPA(3,4) +
	           SPA(1,3)*SPA(4,5)*SPB(5,1))))/
	    (pow(S(2,3) - S(4,5),2)*SPA(1,2)*SPA(2,3)*SPA(4,5)) +
	   (complex<T>(0,0.5)*(Int(eporder,mc,mu,c45,c123))*SPA(1,3)*SPB(5,1)*
	      (-complex<T>(2,0)*S(2,3)*SPA(3,4) + SPA(4,5)*
	         (SPA(1,3)*SPB(5,1) + complex<T>(2,0)*SPA(3,4)*SPB(5,4))))/
	    (pow(S(2,3) - S(4,5),2)*SPA(1,2)*SPA(2,3))

       );
}
// this seems to give the right real part for the cut and is faster than AqppqmCut2
template <class T> complex<T> AqppqmCut3(int eporder, momentum_configuration<T>& mc, const std::vector<int>& ind, int mu)
{
//define corner vectors
 std::vector<int> c1;  c1.push_back(ind[1-1]);
 std::vector<int> c2;  c2.push_back(ind[2-1]);
 std::vector<int> c3;  c3.push_back(ind[3-1]);
 std::vector<int> c4;  c4.push_back(ind[4-1]);
 std::vector<int> c5;  c5.push_back(ind[5-1]);

 std::vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 std::vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
 std::vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
 std::vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
 std::vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);

 std::vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind[i-1]);}
 std::vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind[i-1]);}
 std::vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind[i-1]);}
 std::vector<int> c451;  for(int i = 4; i<=5; i++) {c451.push_back(ind[i-1]);}
                      c451.push_back(ind[1-1]);
 std::vector<int> c145 = c451;
                        // WARNING: Triangle sign fudge to account for
                        // David's backwards sign in the triangles.
                        // Later this will be +1.
 complex<T> trianglesign = complex<T>(1.,0.);
 return(
		 (complex<T>(0,0.5)*(Int(eporder,mc,mu,c3,c2,c1,c45)*
	        pow(S(2,3) - S(4,5),2)*SPA(1,2)*SPA(2,3)*pow(SPA(3,4),2)*
	        SPB(2,1)*SPB(3,2) + Int(eporder,mc,mu,c23,c145)*
	        (complex<T>(3,0)*pow(S(2,3),2)*pow(SPA(3,4),2) +
	          complex<T>(3,0)*pow(S(4,5),2)*pow(SPA(3,4),2) -
	          complex<T>(2,0)*S(4,5)*SPA(1,3)*SPA(3,4)*SPA(4,5)*SPB(5,1) -
	          pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(5,1),2) +
	          complex<T>(2,0)*S(2,3)*SPA(3,4)*(-complex<T>(3,0)*S(4,5)*SPA(3,4) +
	             SPA(1,3)*SPA(4,5)*SPB(5,1))) +
	       SPA(4,5)*(-complex<T>(2,0)*Int(eporder,mc,mu,c4,c5,c123)*
	           pow(S(2,3) - S(4,5),2)*pow(SPA(3,4),2)*SPB(5,4) +
	          Int(eporder,mc,mu,c45,c123)*SPA(1,3)*SPB(5,1)*
	           (-complex<T>(2,0)*S(2,3)*SPA(3,4) +
	             SPA(4,5)*(SPA(1,3)*SPB(5,1) + complex<T>(2,0)*SPA(3,4)*SPB(5,4))))))/
	   (pow(S(2,3) - S(4,5),2)*SPA(1,2)*SPA(2,3)*SPA(4,5))

       );
}



template <class T> complex<T> AqppqmRat(momentum_configuration<T> mc,const std::vector<int> ind){


return( (complex<T>(0.,1./2.)*(-(pow(SPA(1,3),2)*pow(SPA(4,5),2)*
 pow(SPB(5,1),2)*S(2,3))+pow(SPA(1,3),2)*pow(SPA(4,5),2)*
pow(SPB(5,1),2)*S(4,5)))/(pow(S(2,3)-S(4,5),2)*
pow(SPA(4,5),2)*SPA(1,2)*SPA(2,3)*SPB(5,4))
	   );
}



template <class T> complex<T> AqpqmpRat(momentum_configuration<T> mc,const std::vector<int>& ind){
 return( complex<T>(0.,1)*((complex<T>(1./2.,0.)*pow(SPA(2,4),2))/
(SPA(1,3)*SPA(2,3)*
SPA(4,5))+(complex<T>(-1./2.,0.)*pow(SPA(1,4),2)*pow(SPB(3,1),2))/
 ((S(2,3)-S(4,5))*SPA(1,3)*SPA(4,5)*SPB(3,2)) -
(SPA(1,2)*SPA(3,4)*SPB(3,1)*SPB(5,3))/
(S(1,2)*(S(1,2)-S(4,5))*
SPA(1,3))+(complex<T>(-1./2.,0.)*(SPB(3,2)*SPB(5,1)+
SPB(3,1)*SPB(5,2))*
SPB(5,3))/(SPA(1,3)*SPB(2,1)*SPB(3,2)*SPB(5,4)))
	);
}

template <class T> complex<T> Aqppqm(int eporder, momentum_configuration<T>& mc, const std::vector<int>& ind, int mu){
if (eporder==0){
	return AqppqmCut(eporder,mc, ind, mu)+AqppqmRat(mc, ind);
} else {
	return AqppqmCut(eporder,mc, ind, mu);
}

}

template <class T> complex<T> AqpqmpCut(int eporder, momentum_configuration<T> mc, const std::vector<int>& ind, int mu)
{
//define corner vectors
 std::vector<int> c1;  c1.push_back(ind[0]);
 std::vector<int> c2;  c2.push_back(ind[1]);
 std::vector<int> c3;  c3.push_back(ind[2]);
 std::vector<int> c4;  c4.push_back(ind[3]);
 std::vector<int> c5;  c5.push_back(ind[4]);

 std::vector<int> c12;  c12.push_back(ind[0]); c12.push_back(ind[1]);
 std::vector<int> c23;  c23.push_back(ind[1]); c23.push_back(ind[2]);
 std::vector<int> c34;  c34.push_back(ind[2]); c34.push_back(ind[3]);
 std::vector<int> c45;  c45.push_back(ind[3]); c45.push_back(ind[4]);
 std::vector<int> c51;  c51.push_back(ind[4]); c51.push_back(ind[0]);

 std::vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind[i-1]);}
 std::vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind[i-1]);}
 std::vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind[i-1]);}
 std::vector<int> c451;  for(int i = 4; i<=5; i++) {c451.push_back(ind[i-1]);}
                      c451.push_back(ind[0]);
 std::vector<int> c145 = c451;
 std::vector<int> c245;  c245.push_back(ind[1]);
                    for(int i = 4; i<=5; i++) {c245.push_back(ind[i-1]);}
                        // WARNING: Triangle sign fudge to account for
                        // David's backwards sign in the triangles.
                        // Later this will be +1.
 complex<T> trianglesign = complex<T>(1.,0.);
 return(
complex<T>(0.,-1.)*((complex<T>(1./2.,0.)*pow(SPA(2,4),2))/
(SPA(1,3)*SPA(2,3)*
 SPA(4,5))+(complex<T>(-1./2.,0.)*pow(SPA(1,4),2)*pow(SPB(3,1),2))/
((S(2,3)-S(4,5))*SPA(1,3)*SPA(4,5)*SPB(3,2)) -
 (SPA(1,2)*SPA(3,4)*SPB(3,1)*SPB(5,3))/(S(1,2)*(S(1,2)-S(4,5))*
 SPA(1,3))+(complex<T>(-1./2.,0.)*(SPB(3,2)*SPB(5,1)+SPB(3,1)*SPB(5,2))*
 SPB(5,3))/(SPA(1,3)*SPB(2,1)*SPB(3,2)*SPB(5,4))) +
complex<T>(0.,1)*((complex<T>(-1./2.,0.)*
(complex<T>(-1.,0.)+Int(eporder,mc,mu,c45,c123))*
 pow(SPA(2,4),2))/(SPA(1,3)*SPA(2,3)*SPA(4,5)) +
 (pow(SPA(2,4),2)*(complex<T>(2.,0.)*Int(eporder,mc,mu,c45,c123) -
  trianglesign*Int(eporder,mc,mu,c1,c2,c345)*S(1,2)))/
(SPA(1,3)*SPA(2,3)*SPA(4,5)) +
 (pow(SPA(2,4),2)*(trianglesign*Int(eporder,mc,mu,c1,c2,c345)*
   S(1,2)+trianglesign*Int(eporder,mc,mu,c1,c3,c245)*
   S(1,3)+complex<T>(1./2.,0.)*Int(eporder,mc,mu,c3,c1,c2,c45)*S(1,2)*
   S(1,3)-trianglesign*Int(eporder,mc,mu,c4,c5,c123)*
   S(4,5)))/(SPA(1,3)*SPA(2,3)*SPA(4,5)) +
 (pow(SPA(1,4),2)*(trianglesign*Int(eporder,mc,mu,c1,c2,c345)*
   S(1,2)+trianglesign*Int(eporder,mc,mu,c2,c3,c145)*
   S(2,3)+complex<T>(1./2.,0.)*Int(eporder,mc,mu,c3,c2,c1,c45)*S(1,2)*
   S(2,3)-trianglesign*Int(eporder,mc,mu,c4,c5,c123)*S(4,5))*
 SPA(2,3))/(pow(SPA(1,3),3)*SPA(4,5)) +
 ((trianglesign*Int(eporder,mc,mu,c1,c2,c345)*S(1,2) +
  trianglesign*Int(eporder,mc,mu,c2,c3,c145)*S(2,3) +
  complex<T>(1./2.,0.)*Int(eporder,mc,mu,c3,c2,c1,c45)*S(1,2)*S(2,3) -
  trianglesign*Int(eporder,mc,mu,c4,c5,c123)*S(4,5))*SPA(2,4)*
 (-(SPA(1,4)*SPA(2,3))+SPA(1,2)*SPA(3,4)))/
(pow(SPA(1,3),2)*SPA(2,3)*SPA(4,5)) +
 ((Int(eporder,mc,mu,c23,c145)-Int(eporder,mc,mu,c45,c123))*
 pow(SPA(1,4),2)*SPA(2,3)*SPB(3,1))/(pow(SPA(1,3),2)*
 (S(2,3)-S(4,5))*SPA(4,5)) +
 (complex<T>(-2.,0.)*(Int(eporder,mc,mu,c23,c145)-Int(eporder,mc,mu,c45,
   c123))*SPA(1,4)*SPA(2,4)*SPB(3,1))/((S(2,3)-S(4,5))*
 SPA(1,3)*SPA(4,5))-((Int(eporder,mc,mu,c12,c345) -
  Int(eporder,mc,mu,c45,c123))*SPA(1,2)*SPA(1,4)*SPA(3,4)*
 SPB(3,1))/(pow(SPA(1,3),2)*(S(1,2)-S(4,5))*SPA(4,5)) +
 (complex<T>(-1./2.,0.)*pow(SPA(1,4),2)*pow(SPB(3,1),2)*
 (S(2,3)+Int(eporder,mc,mu,c23,c145)*S(2,3) -
  Int(eporder,mc,mu,c45,c123)*S(2,3)-S(4,5)))/
(pow(S(2,3)-S(4,5),2)*SPA(1,3)*SPA(4,5)*SPB(3,2)) -
 ((S(1,2)+Int(eporder,mc,mu,c12,c345)*S(1,2) -
  Int(eporder,mc,mu,c45,c123)*S(1,2)-S(4,5))*SPA(1,2)*
 SPA(3,4)*SPB(3,1)*SPB(5,3))/(pow(S(1,2)-S(4,5),2)*S(1,2)*
 SPA(1,3))+
(complex<T>(-1./2.,0.)*(SPB(3,2)*SPB(5,1)+SPB(3,1)*SPB(5,2))*
 SPB(5,3))/(SPA(1,3)*SPB(2,1)*SPB(3,2)*SPB(5,4)))
);
}
*/
/*
class BH_Ampl_concrete : public BH_Ampl {
private:
	process _process;
	TreeHelAmpl* _tree;
	std::vector<OneLoopHelAmpl*> _Ampls;
	std::complex<double> _factor;
	SeriesC<R> _last_result;
	C _last_Atree;
public:
	BH_Ampl_concrete(size_t pc);
	virtual double operator()(BHinput& in);
	virtual double get_single_pole();
	virtual double get_double_pole();
};


BH_Ampl_concrete::BH_Ampl_concrete(size_t pc){
	switch (pc){
	case 220: {
		_process=process(lm,lbp,qm,qbp);
		_tree=new TreeHelAmpl(_process);
		_Ampls.push_back(new OneLoopHelAmpl(_process,leading_color));
		_factor=C(2.,0);break;
	}

	}
	for (int j=0;j<_Ampls.size();j++){
		_Ampls[j]->set_scheme(HV);
	}
}


double BH_Ampl_concrete::get_single_pole(){

	return (_factor*_last_result[-1]/_last_Atree).real();
}

double BH_Ampl_concrete::get_double_pole(){

	return (_factor*_last_result[-2]/_last_Atree).real();
}

double BH_Ampl_concrete::operator()(BHinput& in){
	mom_conf mc;
// two first momenta are incoming
	mc.insert(Cmom<double>(-in.m_momenta[0][0],-in.m_momenta[0][1],-in.m_momenta[0][2],-in.m_momenta[0][3]));
	mc.insert(Cmom<double>(-in.m_momenta[1][0],-in.m_momenta[1][1],-in.m_momenta[1][2],-in.m_momenta[1][3]));
// others are outgoing



	std::vector<int> ind;
	ind.push_back(1);
	ind.push_back(2);

	for (int i=2;i < in.m_momenta.size();i++){
		mc.insert(Cmom<double>(in.m_momenta[i][0],in.m_momenta[i][1],in.m_momenta[i][2],in.m_momenta[i][3]));
		ind.push_back(i+1);
	}

	SeriesC<R> res(-2,0);
	for (int i=0;i<_Ampls.size();i++ ){
		_Ampls[i]->set_mu(in.m_mu,in.m_mu,in.m_mu);
		res+=_Ampls[i]->get_value(mc,ind);
	}
	_last_result=res;
	_last_Atree=_tree->eval(mc,ind);
	return (_factor*res[0]/_last_Atree).real();
}
*/
/*
class BH_Ampl_ee3jet : public BH_Ampl {
protected:
	std::vector<TreeHelAmpl*> _trees;
	std::vector<OneLoopHelAmpl*> _Ampls;
	std::complex<double> _factor;
	C _last_result_2;
	C _last_result_1;
	C _last_result_0;
	std::vector<int> new_ind;
	std::vector<int> new_ind2;
	std::vector<int> new_ind3;
	std::vector<int> new_ind4;
	std::vector<int> new_ind_flip45;
	std::vector<int> new_ind2_flip45;
	std::vector<int> ind;
	std::vector<int> ind_flip13;
	std::vector<int> ind_flip45;
	std::vector<int> ind_flip1345;
public:
	BH_Ampl_ee3jet();
	BH_Ampl_ee3jet(std::vector<int> mom_assg);
	// storing the momenta in the BH canonical order q g g qb photon
	std::vector<int> momenta_assignment;
	virtual double operator()(BHinput& in);
	virtual double get_single_pole();
	virtual double get_double_pole();
};

class BH_Ampl_ee3jet_LC : public BH_Ampl_ee3jet {
public:
	BH_Ampl_ee3jet_LC() : BH_Ampl_ee3jet() {};
	virtual double operator()(BHinput& in);
};

class BH_Ampl_ee3jet_SLC : public BH_Ampl_ee3jet {
public:
	BH_Ampl_ee3jet_SLC() : BH_Ampl_ee3jet() {};
	virtual double operator()(BHinput& in);
};



BH_Ampl_ee3jet::BH_Ampl_ee3jet(){
		momenta_assignment.push_back(0);
		momenta_assignment.push_back(1);
		momenta_assignment.push_back(2);
		momenta_assignment.push_back(3);
		momenta_assignment.push_back(4);

	process PRO1(qm,m,qbp,lm,lbp);
	process PRO2(qm,p,qbp,lm,lbp);
	process PRO3(qp,m,qbm,lm,lbp);
	process PRO4(qp,p,qbm,lm,lbp);

	process PRO9(qp,qbm,p,lm,lbp);
	process PRO10(qp,qbm,m,lm,lbp);
	process PRO11(p,qm,qbp,lm,lbp);
	process PRO12(m,qm,qbp,lm,lbp);

	_Ampls.push_back(new OneLoopHelAmpl(PRO1,leading_color));
	_Ampls.push_back(new OneLoopHelAmpl(PRO2,leading_color));
	_Ampls.push_back(new OneLoopHelAmpl(PRO3,leading_color));
	_Ampls.push_back(new OneLoopHelAmpl(PRO4,leading_color));

	_trees.push_back(new TreeHelAmpl(PRO1));
	_trees.push_back(new TreeHelAmpl(PRO2));
	_trees.push_back(new TreeHelAmpl(PRO3));
	_trees.push_back(new TreeHelAmpl(PRO4));

	_trees.push_back(new TreeHelAmpl(PRO9));
	_trees.push_back(new TreeHelAmpl(PRO10));
	_trees.push_back(new TreeHelAmpl(PRO11));
	_trees.push_back(new TreeHelAmpl(PRO12));



	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);

	ind_flip13.push_back(3);
	ind_flip13.push_back(2);
	ind_flip13.push_back(1);
	ind_flip13.push_back(4);
	ind_flip13.push_back(5);

	ind_flip1345.push_back(3);
	ind_flip1345.push_back(2);
	ind_flip1345.push_back(1);
	ind_flip1345.push_back(5);
	ind_flip1345.push_back(4);

	ind_flip45.push_back(1);
	ind_flip45.push_back(2);
	ind_flip45.push_back(3);
	ind_flip45.push_back(5);
	ind_flip45.push_back(4);

	new_ind.push_back(1);
	new_ind.push_back(3);
	new_ind.push_back(2);
	new_ind.push_back(4);
	new_ind.push_back(5);

	new_ind2.push_back(3);
	new_ind2.push_back(1);
	new_ind2.push_back(2);
	new_ind2.push_back(5);
	new_ind2.push_back(4);

	new_ind3.push_back(2);
	new_ind3.push_back(1);
	new_ind3.push_back(3);
	new_ind3.push_back(4);
	new_ind3.push_back(5);

	new_ind4.push_back(2);
	new_ind4.push_back(3);
	new_ind4.push_back(1);
	new_ind4.push_back(5);
	new_ind4.push_back(4);

	new_ind_flip45.push_back(1);
	new_ind_flip45.push_back(3);
	new_ind_flip45.push_back(2);
	new_ind_flip45.push_back(5);
	new_ind_flip45.push_back(4);

	new_ind2_flip45.push_back(3);
	new_ind2_flip45.push_back(1);
	new_ind2_flip45.push_back(2);
	new_ind2_flip45.push_back(4);
	new_ind2_flip45.push_back(5);


	_factor=C(1.,0);

	for (int j=0;j<_Ampls.size();j++){
		_Ampls[j]->set_scheme(HV);
	}



}

BH_Ampl_ee3jet::BH_Ampl_ee3jet(std::vector<int> mom_assg): momenta_assignment(mom_assg){

	process PRO1(qm,m,qbp,lm,lbp);
	process PRO2(qm,p,qbp,lm,lbp);
	process PRO3(qp,m,qbm,lm,lbp);
	process PRO4(qp,p,qbm,lm,lbp);

	process PRO9(qp,qbm,p,lm,lbp);
	process PRO10(qp,qbm,m,lm,lbp);
	process PRO11(p,qm,qbp,lm,lbp);
	process PRO12(m,qm,qbp,lm,lbp);

	_Ampls.push_back(new OneLoopHelAmpl(PRO1,leading_color));
	_Ampls.push_back(new OneLoopHelAmpl(PRO2,leading_color));
	_Ampls.push_back(new OneLoopHelAmpl(PRO3,leading_color));
	_Ampls.push_back(new OneLoopHelAmpl(PRO4,leading_color));

	_trees.push_back(new TreeHelAmpl(PRO1));
	_trees.push_back(new TreeHelAmpl(PRO2));
	_trees.push_back(new TreeHelAmpl(PRO3));
	_trees.push_back(new TreeHelAmpl(PRO4));

	_trees.push_back(new TreeHelAmpl(PRO9));
	_trees.push_back(new TreeHelAmpl(PRO10));
	_trees.push_back(new TreeHelAmpl(PRO11));
	_trees.push_back(new TreeHelAmpl(PRO12));



	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);

	ind_flip13.push_back(3);
	ind_flip13.push_back(2);
	ind_flip13.push_back(1);
	ind_flip13.push_back(4);
	ind_flip13.push_back(5);

	ind_flip1345.push_back(3);
	ind_flip1345.push_back(2);
	ind_flip1345.push_back(1);
	ind_flip1345.push_back(5);
	ind_flip1345.push_back(4);

	ind_flip45.push_back(1);
	ind_flip45.push_back(2);
	ind_flip45.push_back(3);
	ind_flip45.push_back(5);
	ind_flip45.push_back(4);

	new_ind.push_back(1);
	new_ind.push_back(3);
	new_ind.push_back(2);
	new_ind.push_back(4);
	new_ind.push_back(5);

	new_ind2.push_back(3);
	new_ind2.push_back(1);
	new_ind2.push_back(2);
	new_ind2.push_back(5);
	new_ind2.push_back(4);

	new_ind3.push_back(2);
	new_ind3.push_back(1);
	new_ind3.push_back(3);
	new_ind3.push_back(4);
	new_ind3.push_back(5);

	new_ind4.push_back(2);
	new_ind4.push_back(3);
	new_ind4.push_back(1);
	new_ind4.push_back(5);
	new_ind4.push_back(4);

	new_ind_flip45.push_back(1);
	new_ind_flip45.push_back(3);
	new_ind_flip45.push_back(2);
	new_ind_flip45.push_back(5);
	new_ind_flip45.push_back(4);

	new_ind2_flip45.push_back(3);
	new_ind2_flip45.push_back(1);
	new_ind2_flip45.push_back(2);
	new_ind2_flip45.push_back(4);
	new_ind2_flip45.push_back(5);


	_factor=C(1.,0);

	for (int j=0;j<_Ampls.size();j++){
		_Ampls[j]->set_scheme(HV);
	}




}

double BH_Ampl_ee3jet_LC::operator()(BHinput& in){
// the momenta are entered in the mom_conf so that the quark is firsst, the gluon second and that anti-quark third

	mom_conf mc;

	// q
	if(momenta_assignment[0]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[0]][0],-in.m_momenta[momenta_assignment[0]][1],
			-in.m_momenta[momenta_assignment[0]][2],-in.m_momenta[momenta_assignment[0]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[0]][0],in.m_momenta[momenta_assignment[0]][1],
			in.m_momenta[momenta_assignment[0]][2],in.m_momenta[momenta_assignment[0]][3]));
	}
	// gluon
	if(momenta_assignment[1]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[1]][0],-in.m_momenta[momenta_assignment[1]][1],
			-in.m_momenta[momenta_assignment[1]][2],-in.m_momenta[momenta_assignment[1]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[1]][0],in.m_momenta[momenta_assignment[1]][1],
			in.m_momenta[momenta_assignment[1]][2],in.m_momenta[momenta_assignment[1]][3]));
	}
	// qb
	if(momenta_assignment[2]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[2]][0],-in.m_momenta[momenta_assignment[2]][1],
			-in.m_momenta[momenta_assignment[2]][2],-in.m_momenta[momenta_assignment[2]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[2]][0],in.m_momenta[momenta_assignment[2]][1],
			in.m_momenta[momenta_assignment[2]][2],in.m_momenta[momenta_assignment[2]][3]));
	}
	// l
	if(momenta_assignment[3]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[3]][0],-in.m_momenta[momenta_assignment[3]][1],
			-in.m_momenta[momenta_assignment[3]][2],-in.m_momenta[momenta_assignment[3]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[3]][0],in.m_momenta[momenta_assignment[3]][1],
			in.m_momenta[momenta_assignment[3]][2],in.m_momenta[momenta_assignment[3]][3]));
	}
	// lb
	if(momenta_assignment[4]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[4]][0],-in.m_momenta[momenta_assignment[4]][1],
			-in.m_momenta[momenta_assignment[4]][2],-in.m_momenta[momenta_assignment[4]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[4]][0],in.m_momenta[momenta_assignment[4]][1],
			in.m_momenta[momenta_assignment[4]][2],in.m_momenta[momenta_assignment[4]][3]));
	}
//	// quark
//	mc.insert(Cmom<double>(in.m_momenta[3][0],in.m_momenta[3][1],in.m_momenta[3][2],in.m_momenta[3][3]));
//	// gluon
//	mc.insert(Cmom<double>(in.m_momenta[2][0],in.m_momenta[2][1],in.m_momenta[2][2],in.m_momenta[2][3]));
//	// antiquark
//	mc.insert(Cmom<double>(in.m_momenta[4][0],in.m_momenta[4][1],in.m_momenta[4][2],in.m_momenta[4][3]));
//
//	// two first momenta are incoming
//	mc.insert(Cmom<double>(-in.m_momenta[0][0],-in.m_momenta[0][1],-in.m_momenta[0][2],-in.m_momenta[0][3]));
//	mc.insert(Cmom<double>(-in.m_momenta[1][0],-in.m_momenta[1][1],-in.m_momenta[1][2],-in.m_momenta[1][3]));
//
	std::vector<int> ind;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);

	int the_david_mu = DefineMu(mc,in.m_mu);
	SeriesC<R> res(-2,0);
	C treesum(0,0);

	for (int i=0;i<_Ampls.size();i++ ){
				_Ampls[i]->set_mu(in.m_mu,in.m_mu,in.m_mu);
		res += (_Ampls[i]->get_value(mc,ind)*conj(_trees[i]->eval(mc,ind)));

	}



	C TSL1=_trees[4]->eval(mc,new_ind);
	C TSL2=_trees[4]->eval(mc,new_ind2);
	C TSL3=_trees[4]->eval(mc,new_ind2_flip45);
	C TSL4=_trees[4]->eval(mc,new_ind_flip45);

	C trees_SL =TSL1*conj(TSL1)+TSL2*conj(TSL2)+TSL3*conj(TSL3)+TSL4*conj(TSL4);

	double renorm = - (Nc*(11./3. ))/2.;

	_last_result_2 = Nc * (res[-2] )/trees_SL;
	_last_result_1 = Nc * (res[-1] )/trees_SL  + renorm ;

	_last_result_0 = Nc * ( res[ 0])    ;
	return (_factor*_last_result_0).real();
}

double BH_Ampl_ee3jet_SLC::operator()(BHinput& in){
// the momenta are entered in the mom_conf so that the quark is firsst, the gluon second and that anti-quark third

	mom_conf mc;

	// q
	if(momenta_assignment[0]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[0]][0],-in.m_momenta[momenta_assignment[0]][1],
			-in.m_momenta[momenta_assignment[0]][2],-in.m_momenta[momenta_assignment[0]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[0]][0],in.m_momenta[momenta_assignment[0]][1],
			in.m_momenta[momenta_assignment[0]][2],in.m_momenta[momenta_assignment[0]][3]));
	}
	// gluon
	if(momenta_assignment[1]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[1]][0],-in.m_momenta[momenta_assignment[1]][1],
			-in.m_momenta[momenta_assignment[1]][2],-in.m_momenta[momenta_assignment[1]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[1]][0],in.m_momenta[momenta_assignment[1]][1],
			in.m_momenta[momenta_assignment[1]][2],in.m_momenta[momenta_assignment[1]][3]));
	}
	// qb
	if(momenta_assignment[2]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[2]][0],-in.m_momenta[momenta_assignment[2]][1],
			-in.m_momenta[momenta_assignment[2]][2],-in.m_momenta[momenta_assignment[2]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[2]][0],in.m_momenta[momenta_assignment[2]][1],
			in.m_momenta[momenta_assignment[2]][2],in.m_momenta[momenta_assignment[2]][3]));
	}
	// l
	if(momenta_assignment[3]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[3]][0],-in.m_momenta[momenta_assignment[3]][1],
			-in.m_momenta[momenta_assignment[3]][2],-in.m_momenta[momenta_assignment[3]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[3]][0],in.m_momenta[momenta_assignment[3]][1],
			in.m_momenta[momenta_assignment[3]][2],in.m_momenta[momenta_assignment[3]][3]));
	}
	// lb
	if(momenta_assignment[4]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[4]][0],-in.m_momenta[momenta_assignment[4]][1],
			-in.m_momenta[momenta_assignment[4]][2],-in.m_momenta[momenta_assignment[4]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[4]][0],in.m_momenta[momenta_assignment[4]][1],
			in.m_momenta[momenta_assignment[4]][2],in.m_momenta[momenta_assignment[4]][3]));
	}
//	// quark
//	mc.insert(Cmom<double>(in.m_momenta[3][0],in.m_momenta[3][1],in.m_momenta[3][2],in.m_momenta[3][3]));
//	// gluon
//	mc.insert(Cmom<double>(in.m_momenta[2][0],in.m_momenta[2][1],in.m_momenta[2][2],in.m_momenta[2][3]));
//	// antiquark
//	mc.insert(Cmom<double>(in.m_momenta[4][0],in.m_momenta[4][1],in.m_momenta[4][2],in.m_momenta[4][3]));
//
//	// two first momenta are incoming
//	mc.insert(Cmom<double>(-in.m_momenta[0][0],-in.m_momenta[0][1],-in.m_momenta[0][2],-in.m_momenta[0][3]));
//	mc.insert(Cmom<double>(-in.m_momenta[1][0],-in.m_momenta[1][1],-in.m_momenta[1][2],-in.m_momenta[1][3]));

	int the_david_mu = DefineMu(mc,in.m_mu);

	C TSL1=_trees[4]->eval(mc,new_ind);
	C TSL2=_trees[4]->eval(mc,new_ind2);
	C TSL3=_trees[4]->eval(mc,new_ind2_flip45);
	C TSL4=_trees[4]->eval(mc,new_ind_flip45);


	double scheme_shift=-1./2.;
	C LSL2a=(AqpqmpCut(-2,mc,new_ind,the_david_mu)*conj(TSL1));
	C LSL1a=(AqpqmpCut(-1,mc,new_ind,the_david_mu)*conj(TSL1));
	C LSL0a=((AqpqmpCut( 0,mc,new_ind,the_david_mu)+AqpqmpRat(mc,new_ind)+scheme_shift*TSL1)*conj(TSL1));

	C LSL2b=(AqpqmpCut(-2,mc,new_ind2,the_david_mu)*conj(TSL2));
	C LSL1b=(AqpqmpCut(-1,mc,new_ind2,the_david_mu)*conj(TSL2));
	C LSL0b=((AqpqmpCut( 0,mc,new_ind2,the_david_mu)+AqpqmpRat(mc,new_ind2)+scheme_shift*TSL2)*conj(TSL2));

	C LSL2c=(AqpqmpCut(-2,mc,new_ind2_flip45,the_david_mu)*conj(TSL3));
	C LSL1c=(AqpqmpCut(-1,mc,new_ind2_flip45,the_david_mu)*conj(TSL3));
	C LSL0c=((AqpqmpCut( 0,mc,new_ind2_flip45,the_david_mu)+AqpqmpRat(mc,new_ind2_flip45)+scheme_shift*TSL3)*conj(TSL3));

	C LSL2d=(AqpqmpCut(-2,mc,new_ind_flip45,the_david_mu)*conj(TSL4));
	C LSL1d=(AqpqmpCut(-1,mc,new_ind_flip45,the_david_mu)*conj(TSL4));
	C LSL0d=((AqpqmpCut( 0,mc,new_ind_flip45,the_david_mu)+AqpqmpRat(mc,new_ind_flip45)+scheme_shift*TSL4)*conj(TSL4));



	C trees_SL =TSL1*conj(TSL1)+TSL2*conj(TSL2)+TSL3*conj(TSL3)+TSL4*conj(TSL4);



	double inv=1./(Nc*Nc);
	double renorm = - (Nc*( - (2*Nf)/(3.*Nc) ))/2.;

	_last_result_2 = Nc * ( - inv *(LSL2a+LSL2b+LSL2c+LSL2d))/trees_SL;
	_last_result_1 = Nc * ( - inv *(LSL1a+LSL1b+LSL1c+LSL1d))/trees_SL  + renorm ;

	_last_result_0 = Nc * ( - inv *(LSL0a+LSL0b+LSL0c+LSL0d));
	return (_factor*_last_result_0).real();
}


double BH_Ampl_ee3jet::operator()(BHinput& in){
// the momenta are entered in the mom_conf so that the quark is firsst, the gluon second and that anti-quark third

	mom_conf mc;

	// q
	if(momenta_assignment[0]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[0]][0],-in.m_momenta[momenta_assignment[0]][1],
			-in.m_momenta[momenta_assignment[0]][2],-in.m_momenta[momenta_assignment[0]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[0]][0],in.m_momenta[momenta_assignment[0]][1],
			in.m_momenta[momenta_assignment[0]][2],in.m_momenta[momenta_assignment[0]][3]));
	}
	// gluon
	if(momenta_assignment[1]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[1]][0],-in.m_momenta[momenta_assignment[1]][1],
			-in.m_momenta[momenta_assignment[1]][2],-in.m_momenta[momenta_assignment[1]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[1]][0],in.m_momenta[momenta_assignment[1]][1],
			in.m_momenta[momenta_assignment[1]][2],in.m_momenta[momenta_assignment[1]][3]));
	}
	// qb
	if(momenta_assignment[2]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[2]][0],-in.m_momenta[momenta_assignment[2]][1],
			-in.m_momenta[momenta_assignment[2]][2],-in.m_momenta[momenta_assignment[2]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[2]][0],in.m_momenta[momenta_assignment[2]][1],
			in.m_momenta[momenta_assignment[2]][2],in.m_momenta[momenta_assignment[2]][3]));
	}
	// l
	if(momenta_assignment[3]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[3]][0],-in.m_momenta[momenta_assignment[3]][1],
			-in.m_momenta[momenta_assignment[3]][2],-in.m_momenta[momenta_assignment[3]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[3]][0],in.m_momenta[momenta_assignment[3]][1],
			in.m_momenta[momenta_assignment[3]][2],in.m_momenta[momenta_assignment[3]][3]));
	}
	// lb
	if(momenta_assignment[4]<2){
		mc.insert(Cmom<double>(-in.m_momenta[momenta_assignment[4]][0],-in.m_momenta[momenta_assignment[4]][1],
			-in.m_momenta[momenta_assignment[4]][2],-in.m_momenta[momenta_assignment[4]][3]));
	}
	else{
		mc.insert(Cmom<double>(in.m_momenta[momenta_assignment[4]][0],in.m_momenta[momenta_assignment[4]][1],
			in.m_momenta[momenta_assignment[4]][2],in.m_momenta[momenta_assignment[4]][3]));
	}
//	// quark
//	mc.insert(Cmom<double>(in.m_momenta[3][0],in.m_momenta[3][1],in.m_momenta[3][2],in.m_momenta[3][3]));
//	// gluon
//	mc.insert(Cmom<double>(in.m_momenta[2][0],in.m_momenta[2][1],in.m_momenta[2][2],in.m_momenta[2][3]));
//	// antiquark
//	mc.insert(Cmom<double>(in.m_momenta[4][0],in.m_momenta[4][1],in.m_momenta[4][2],in.m_momenta[4][3]));
//
//	// two first momenta are incoming
//	mc.insert(Cmom<double>(-in.m_momenta[0][0],-in.m_momenta[0][1],-in.m_momenta[0][2],-in.m_momenta[0][3]));
//	mc.insert(Cmom<double>(-in.m_momenta[1][0],-in.m_momenta[1][1],-in.m_momenta[1][2],-in.m_momenta[1][3]));
// others are outgoing


	int the_david_mu = DefineMu(mc,in.m_mu);
//	int the_david_mu = DefineMu(mc,1.);
//	_PRINT(mc.s(1,2));
//	_PRINT(mc.s(1,3));
//	_PRINT(mc.s(1,4));
//	_PRINT(mc.s(1,5));
//	_PRINT(mc.s(3,4));
//	_PRINT(mc.s(3,5));
//	_PRINT(mc.s(4,5));
	SeriesC<R> res(-2,0);
	C treesum(0,0);

//	_PRINT(in.m_mu);
//	_PRINT(log(mc.s(1,2))+log(mc.s(2,3)));
//	_PRINT(log(mc.s(1,3)));
#if 0
	for (int i=0;i<_Ampls.size();i++ ){
				_Ampls[i]->set_mu(in.m_mu,in.m_mu,in.m_mu);
//				_Ampls[i]->set_mu(1.,1.,1.);
//				_PRINT(_Ampls[i]->get_value(mc,ind));
//				_PRINT(_Ampls[i]->get_cut(mc,ind)*conj(_trees[i]->eval(mc,ind)));
//				_PRINT(_Ampls[i]->get_rational(mc,ind)*conj(_trees[i]->eval(mc,ind)));
//				_PRINT(1./2.*_trees[i]->eval(mc,ind)*conj(_trees[i]->eval(mc,ind)));
//				_PRINT(_Ampls[i]->get_value(mc,ind)*conj(_trees[i]->eval(mc,ind)));
		res += (_Ampls[i]->get_value(mc,ind)*conj(_trees[i]->eval(mc,ind)));
//	_PRINT(	_trees[i]->eval(mc,ind)*conj(_trees[i]->eval(mc,ind)));

	}
#endif
	C res0(0,0),res1(0,0),res2(0,0);

	res2+=(AqppqmCut2(-2,mc,ind,the_david_mu))*conj(_trees[3]->eval(mc,ind));
	res2+=(AqppqmCut2(-2,mc,ind_flip13,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip13));
	res2+=(AqppqmCut2(-2,mc,ind_flip45,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip45));
	res2+=(AqppqmCut2(-2,mc,ind_flip1345,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip1345));


	res1+=(AqppqmCut2(-1,mc,ind,the_david_mu))*conj(_trees[3]->eval(mc,ind));
	res1+=(AqppqmCut2(-1,mc,ind_flip13,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip13));
	res1+=(AqppqmCut2(-1,mc,ind_flip45,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip45));
	res1+=(AqppqmCut2(-1,mc,ind_flip1345,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip1345));


	res0+=(AqppqmCut3(0,mc,ind,the_david_mu)+AqppqmRat(mc,ind)-1./2.*_trees[3]->eval(mc,ind))*conj(_trees[3]->eval(mc,ind));
	res0+=(AqppqmCut3(0,mc,ind_flip13,the_david_mu)+AqppqmRat(mc,ind_flip13)-1./2.*_trees[3]->eval(mc,ind_flip13))*conj(_trees[3]->eval(mc,ind_flip13));
	res0+=(AqppqmCut3(0,mc,ind_flip45,the_david_mu)+AqppqmRat(mc,ind_flip45)-1./2.*_trees[3]->eval(mc,ind_flip45))*conj(_trees[3]->eval(mc,ind_flip45));
	res0+=(AqppqmCut3(0,mc,ind_flip1345,the_david_mu)+AqppqmRat(mc,ind_flip1345)-1./2.*_trees[3]->eval(mc,ind_flip1345))*conj(_trees[3]->eval(mc,ind_flip1345));

//	_MESSAGE2("cut: ",(AqppqmCut(0,mc,ind,the_david_mu))*conj(_trees[3]->eval(mc,ind)));
//	_MESSAGE2("rat: ",AqppqmRat(mc,ind)*conj(_trees[3]->eval(mc,ind)));
//	_MESSAGE2("scheme shift:",(-1./2.*_trees[3]->eval(mc,ind))*conj(_trees[3]->eval(mc,ind)));
//	_PRINT((AqppqmCut(0,mc,ind,the_david_mu)+AqppqmRat(mc,ind)-1./2.*_trees[3]->eval(mc,ind))*conj(_trees[3]->eval(mc,ind)));
//	C cut=(AqppqmCut(0,mc,ind_flip13,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip13));
//	C rat= AqppqmRat(mc,ind_flip13)*conj(_trees[3]->eval(mc,ind_flip13));
//	C sch=(-1./2.*_trees[3]->eval(mc,ind_flip13))*conj(_trees[3]->eval(mc,ind_flip13));
//	_MESSAGE2("cut: ",cut);
//	_MESSAGE2("rat: ",rat);
//	_MESSAGE2("scheme shift: ",sch);
//	_PRINT(sch+rat+cut);
//	_PRINT((AqppqmCut(0,mc,ind_flip13,the_david_mu)+AqppqmRat(mc,ind_flip13)-1./2.*_trees[3]->eval(mc,ind_flip13))*conj(_trees[3]->eval(mc,ind_flip13)));
//	_MESSAGE2("cut: ",(AqppqmCut(0,mc,ind_flip45,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip45)));
//	_MESSAGE2("rat: ",AqppqmRat(mc,ind_flip45)*conj(_trees[3]->eval(mc,ind_flip45)));
//	_MESSAGE2("scheme shift:",(-1./2.*_trees[3]->eval(mc,ind_flip45))*conj(_trees[3]->eval(mc,ind_flip45)));
//	_PRINT((AqppqmCut(0,mc,ind_flip45,the_david_mu)+AqppqmRat(mc,ind_flip45)-1./2.*_trees[3]->eval(mc,ind_flip45))*conj(_trees[3]->eval(mc,ind_flip45)));
//	_MESSAGE2("cut: ",(AqppqmCut(0,mc,ind_flip1345,the_david_mu))*conj(_trees[3]->eval(mc,ind_flip1345)));
//	_MESSAGE2("rat: ",AqppqmRat(mc,ind_flip1345)*conj(_trees[3]->eval(mc,ind_flip1345)));
//	_MESSAGE2("scheme shift:",(-1./2.*_trees[3]->eval(mc,ind_flip1345))*conj(_trees[3]->eval(mc,ind_flip1345)));
//	_PRINT((AqppqmCut(0,mc,ind_flip1345,the_david_mu)+AqppqmRat(mc,ind_flip1345)-1./2.*_trees[3]->eval(mc,ind_flip1345))*conj(_trees[3]->eval(mc,ind_flip1345)));

//	_PRINT(res2);
//	_PRINT(res1);
//	_PRINT(res0);
//	_PRINT(res[0]-res0);

//	_PRINT(_trees[0]->eval(mc,ind));
//	_PRINT(_trees[1]->eval(mc,ind));
//	_PRINT(_trees[2]->eval(mc,ind));
//	_PRINT(_trees[3]->eval(mc,ind));
//
//	_PRINT(_trees[3]->eval(mc,ind_flip13));
//	_PRINT(_trees[3]->eval(mc,ind_flip45));
//	_PRINT(_trees[3]->eval(mc,ind_flip1345));
//
//	_PRINT(_Ampls[0]->get_cut(mc,ind)/_trees[0]->eval(mc,ind));
//	_PRINT(_Ampls[1]->get_cut(mc,ind)/_trees[1]->eval(mc,ind));
//	_PRINT(_Ampls[2]->get_cut(mc,ind)/_trees[2]->eval(mc,ind));
//	_PRINT(_Ampls[3]->get_cut(mc,ind)/_trees[3]->eval(mc,ind));

//	_PRINT(_Ampls[0]->get_value(mc,ind)*_trees[0]->eval(mc,ind));
//	_PRINT(_Ampls[1]->get_value(mc,ind)*_trees[1]->eval(mc,ind));
//	_PRINT(_Ampls[2]->get_value(mc,ind)*_trees[2]->eval(mc,ind));
//	_PRINT(_Ampls[3]->get_value(mc,ind)*_trees[3]->eval(mc,ind));
//
//	_PRINT(AqppqmCut(0,mc,ind_flip13,the_david_mu)/_trees[3]->eval(mc,ind_flip13));
//	_PRINT(AqppqmCut(0,mc,ind_flip45,the_david_mu)/_trees[3]->eval(mc,ind_flip45));
//	_PRINT(AqppqmCut(0,mc,ind_flip1345,the_david_mu)/_trees[3]->eval(mc,ind_flip1345));
//
//	_PRINT((AqppqmCut2(0,mc,ind_flip13,the_david_mu))/_trees[3]->eval(mc,ind_flip13));
//	_PRINT((AqppqmCut2(0,mc,ind_flip45,the_david_mu))/_trees[3]->eval(mc,ind_flip45));
//	_PRINT((AqppqmCut2(0,mc,ind_flip1345,the_david_mu))/_trees[3]->eval(mc,ind_flip1345));
//	_PRINT((AqppqmCut2(0,mc,ind,the_david_mu))/_trees[3]->eval(mc,ind));
//
//	_PRINT((AqppqmCut3(0,mc,ind_flip13,the_david_mu))/_trees[3]->eval(mc,ind_flip13));
//	_PRINT((AqppqmCut3(0,mc,ind_flip45,the_david_mu))/_trees[3]->eval(mc,ind_flip45));
//	_PRINT((AqppqmCut3(0,mc,ind_flip1345,the_david_mu))/_trees[3]->eval(mc,ind_flip1345));
//
//	_PRINT((AqppqmCut2(0,mc,ind_flip13,the_david_mu)-AqppqmCut3(0,mc,ind_flip13,the_david_mu))/_trees[3]->eval(mc,ind_flip13));
//	_PRINT((AqppqmCut2(0,mc,ind_flip45,the_david_mu)-AqppqmCut3(0,mc,ind_flip45,the_david_mu))/_trees[3]->eval(mc,ind_flip45));
//	_PRINT((AqppqmCut2(0,mc,ind_flip1345,the_david_mu)-AqppqmCut3(0,mc,ind_flip1345,the_david_mu))/_trees[3]->eval(mc,ind_flip1345));
//
//	_PRINT(_Ampls[0]->get_rational(mc,ind)/_trees[0]->eval(mc,ind));
//	_PRINT(_Ampls[1]->get_rational(mc,ind)/_trees[1]->eval(mc,ind));
//	_PRINT(_Ampls[2]->get_rational(mc,ind)/_trees[2]->eval(mc,ind));
//	_PRINT(_Ampls[3]->get_rational(mc,ind)/_trees[3]->eval(mc,ind));
//
//	_PRINT(AqppqmRat(mc,ind)/_trees[3]->eval(mc,ind));
//	_PRINT(AqppqmRat(mc,ind_flip13)/_trees[3]->eval(mc,ind_flip13));
//	_PRINT(AqppqmRat(mc,ind_flip45)/_trees[3]->eval(mc,ind_flip45));
//	_PRINT(AqppqmRat(mc,ind_flip1345)/_trees[3]->eval(mc,ind_flip1345));
//	_PRINT(Aqppqm(0,mc,ind,the_david_mu)/_trees[3]->eval(mc,ind));


//	_PRINT(mc[1]);
//	_PRINT(mc[2]);
//	_PRINT(mc[3]);
//	_PRINT(mc[4]);
//	_PRINT(mc[5]);
//	_PRINT(mc[1]*mc[1]);
//	_PRINT(mc[2]*mc[2]);
//	_PRINT(mc[3]*mc[3]);
//	_PRINT(mc[4]*mc[4]);
//	_PRINT(mc[5]*mc[5]);
//	_PRINT(mc[1]+mc[2]+mc[3]+mc[4]+mc[5]);

	C TSL1=_trees[4]->eval(mc,new_ind);
	C TSL2=_trees[4]->eval(mc,new_ind2);
	C TSL3=_trees[4]->eval(mc,new_ind2_flip45);
	C TSL4=_trees[4]->eval(mc,new_ind_flip45);

//	_PRINT(_trees[0]->eval(mc,ind));
//	_PRINT(_trees[1]->eval(mc,ind));
//	_PRINT(_trees[2]->eval(mc,ind));
//	_PRINT(_trees[3]->eval(mc,ind));
//	_PRINT(_trees[4]->eval(mc,new_ind));
//	_PRINT(_trees[4]->eval(mc,new_ind2));
//	_PRINT(_trees[4]->eval(mc,new_ind_flip45));
//	_PRINT(_trees[4]->eval(mc,new_ind2_flip45));
//
//	_PRINT(AqpqmpCut(-2,mc,new_ind,the_david_mu));
//	_PRINT(AqpqmpCut(-2,mc,new_ind2,the_david_mu));
//	_PRINT(AqpqmpCut(-2,mc,new_ind2_flip45,the_david_mu));
//	_PRINT(AqpqmpCut(-2,mc,new_ind_flip45,the_david_mu));
//
//	_PRINT(TSL1);
//	_PRINT(TSL2);
//	_PRINT(TSL3);
//	_PRINT(TSL4);

	double scheme_shift=-1./2.;
	C LSL2a=(AqpqmpCut(-2,mc,new_ind,the_david_mu)*conj(TSL1));
	C LSL1a=(AqpqmpCut(-1,mc,new_ind,the_david_mu)*conj(TSL1));
	C LSL0a=((AqpqmpCut( 0,mc,new_ind,the_david_mu)+AqpqmpRat(mc,new_ind)+scheme_shift*TSL1)*conj(TSL1));

	C LSL2b=(AqpqmpCut(-2,mc,new_ind2,the_david_mu)*conj(TSL2));
	C LSL1b=(AqpqmpCut(-1,mc,new_ind2,the_david_mu)*conj(TSL2));
	C LSL0b=((AqpqmpCut( 0,mc,new_ind2,the_david_mu)+AqpqmpRat(mc,new_ind2)+scheme_shift*TSL2)*conj(TSL2));

	C LSL2c=(AqpqmpCut(-2,mc,new_ind2_flip45,the_david_mu)*conj(TSL3));
	C LSL1c=(AqpqmpCut(-1,mc,new_ind2_flip45,the_david_mu)*conj(TSL3));
	C LSL0c=((AqpqmpCut( 0,mc,new_ind2_flip45,the_david_mu)+AqpqmpRat(mc,new_ind2_flip45)+scheme_shift*TSL3)*conj(TSL3));

	C LSL2d=(AqpqmpCut(-2,mc,new_ind_flip45,the_david_mu)*conj(TSL4));
	C LSL1d=(AqpqmpCut(-1,mc,new_ind_flip45,the_david_mu)*conj(TSL4));
	C LSL0d=((AqpqmpCut( 0,mc,new_ind_flip45,the_david_mu)+AqpqmpRat(mc,new_ind_flip45)+scheme_shift*TSL4)*conj(TSL4));

//	_PRINT((LSL2a));
//	_PRINT((LSL2b));
//	_PRINT((LSL2c));
//	_PRINT((LSL2d));
//
//	_PRINT((LSL2a)/(TSL1*conj(TSL1)));
//	_PRINT((LSL2b)/(TSL2*conj(TSL2)));
//	_PRINT((LSL2c)/(TSL3*conj(TSL3)));
//	_PRINT((LSL2d)/(TSL4*conj(TSL4)));


	C trees_SL =TSL1*conj(TSL1)+TSL2*conj(TSL2)+TSL3*conj(TSL3)+TSL4*conj(TSL4);

//	_PRINT((LSL2a+LSL2b+LSL2c+LSL2d)/trees_SL);
//	_PRINT((LSL1a+LSL1b+LSL1c+LSL1d)/trees_SL);
//	_PRINT((LSL0a+LSL0b+LSL0c+LSL0d)/trees_SL);


	double inv=1./(Nc*Nc);
//	double inv=0./(Nc*Nc);
	double renorm = - (Nc*(11./3. - (2.*Nf)/(3.*Nc) - Ns/(3.*Nc)))/2.;
//	double renorm = 0;

	_last_result_2 = Nc * (res2 - inv *(LSL2a+LSL2b+LSL2c+LSL2d))/trees_SL;
	_last_result_1 = Nc * (res1 - inv *(LSL1a+LSL1b+LSL1c+LSL1d))/trees_SL  + renorm ;
//	_PRINT(Nc*res[0]/trees_SL);
//	_PRINT((LSL0a+LSL0b+LSL0c+LSL0d)/trees_SL/Nc);

	_last_result_0 = Nc * ( res0 - inv *(LSL0a+LSL0b+LSL0c+LSL0d))/trees_SL;
//	_last_Atree=_tree->eval(mc,ind);
	return (_factor*_last_result_0).real();
}

double BH_Ampl_ee3jet::get_single_pole(){

	return (_factor*_last_result_1).real();
}

double BH_Ampl_ee3jet::get_double_pole(){

	return (_factor*_last_result_2).real();
}
*/
#if 0
BH_Ampl* BH_factory::new_ampl(const std::vector<int>& labels){
	// gluon 21
	// photon 22
	// d 1
	// u 2
	//s c b t
	//e- 11
	//ne 12
	//mu- 13
	//nmu 14
	//tau- 15
	//ntau- 16

	// updates the values of the constant
	constants::update_constants(d_settings_p);

	std::vector<particle> types;

size_t nbr_gluons=0,nbr_quarks=0,nbr_leptons=0,nbr_photon=0;

for (int i=0;i<labels.size();i++){
	switch (std::abs(labels[i])){
		case 1: case 2: case 3: case 4: case 5: types.push_back(quark); ++nbr_quarks; break;
		case 11: case 13: case 15: types.push_back(lepton); ++nbr_leptons;break;
		case 12: case 14: case 16: types.push_back(lepton); ++nbr_leptons;break;
		case 21: types.push_back(gluon);++nbr_gluons; break;
		case 22: types.push_back(photon);++nbr_photon; break;
	}
}

size_t pc=nbr_gluons+10*nbr_quarks+100*nbr_leptons+100000*nbr_photon;

	int case4q=-1;
	bool up_down_quark=0;
	// default zero, which is for only off-shell photon
	int photonZW=0;
	// labels_flipped corrects for the all outgoing BH conventions if initial fermions are present
	std::vector<int> labels_flipped;
		if((std::abs(labels[0])>0 && std::abs(labels[0])<6) || (std::abs(labels[0])>10 && std::abs(labels[0])<17))
			labels_flipped.push_back(-labels[0]);
		else
			labels_flipped.push_back(labels[0]);
		if((std::abs(labels[1])>0 && std::abs(labels[1])<6) || (std::abs(labels[1])>10 && std::abs(labels[1])<17))
			labels_flipped.push_back(-labels[1]);
		else
			labels_flipped.push_back(labels[1]);
		for(int k=2;k<labels.size();k++)
			labels_flipped.push_back(labels[k]);

	switch(pc){
		  case 41: case 42:{
			// needs modification for W --- so far only drops the run
			int quark1=0,quark2=0,quark3=0,quark4=0;
			for(int i=0;i<labels_flipped.size();i++){
				if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark1==0)
					quark1=labels_flipped[i];
				else if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark2==0)
					quark2=labels_flipped[i];
				else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark3==0)
					quark3=labels_flipped[i];
				else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark4==0)
					quark4=labels_flipped[i];
			}
			if(quark1==0 || quark2==0 || quark3==0 || quark4==0 ){
				std::cout<<"wrong assignmnet of particle labels in 2q2Q1g amps\n";
				throw BHerror("Sorry, can't do that yet. (from 2q2Q1g)");
			}
			else if((quark1==-quark3 && quark2==-quark4) || (quark1==-quark4 && quark2==-quark3)){
				if(quark1==quark2 && quark1%2==0)
					// uu
					case4q=0;
				else if(quark1==quark2)
					// dd
					case4q=5;
				else if(std::abs(quark1-quark2)%2==0 && quark1%2==0)
					// uup
					case4q=1;
				else if(std::abs(quark1-quark2)%2==0)
					// ddp
					case4q=4;
				else if(quark1%2==0)
					// ud
					case4q=2;
				else
					// du
					case4q=3;
			}
			else{
				std::cout<<"wrong assignmnet of particle labels in 2q2Q1g amps!\n";
				throw BHerror("Sorry, can't do that yet.");
			}
		  } break;
		case 220: case 221: case 222: case 223: case 224:{
		int lepton1=0,lepton2=0;
		int quark1=0,quark2=0;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17 && lepton1==0)
				lepton1=labels_flipped[i];
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)
				lepton2=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark1==0)
				quark1=labels_flipped[i];
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark2==0)
				quark2=labels_flipped[i];
		}
		if(quark1==0 || quark2==0 || lepton1==0 || lepton2==0){
			std::cout<<"wrong assignmnet of particle labels in 2q2g2l amps\n";
			throw BHerror("Sorry, can't do that yet.");
		}
		else if(lepton1==-lepton2 && (std::abs(lepton1)==11 || std::abs(lepton1)==13 /*|| std::abs(lepton1)==15*/)){
			// gamma+Z case
			photonZW=1;
		}
		else if(lepton1==-lepton2  && /* abuse of tauon for this switch */ (std::abs(lepton1)==15)){
			// forced gamma decaying into leptons
			photonZW=0;
		}
		else if(lepton1==-lepton2  && (std::abs(lepton1)==12 || std::abs(lepton1)==14 || std::abs(lepton1)==16)){
			// pure Z case decaying into neutrinos
			photonZW=2;
		}
		else if((((lepton1==11 || lepton1==13 || lepton1==15) && (lepton1+lepton2)==-1) ||
			((lepton1==12 || lepton1==14 || lepton1==16) && (lepton1+lepton2)==1)) &&
			std::abs(quark1+quark2)%2==1){
			// W case
			photonZW=3;
		}
		else{
			std::cout<<"wrong assignmnet of particle labels in 2q2g2l amps\n";
			throw BHerror("Sorry, can't do that yet.");
		}
		// W case keeps default value up_down_quark=0
		if(std::abs(quark1)==std::abs(quark2)){
			switch (std::abs(quark1)){
				// up quarks
				case 2: case 4:	up_down_quark=0; break;
				// down quarks
				case 1: case 3: case 5:up_down_quark=1; break;
			}
		}
	  } break;
	  case 240: case 241:{
		// needs modification for W --- so far only drops the run
		int lepton1=0,lepton2=0;
		int quark1=0,quark2=0,quark3=0,quark4=0;
		for(int i=0;i<labels_flipped.size();i++){
			if(std::abs(labels_flipped[i])>10 && std::abs(labels_flipped[i])<17 && lepton1==0)
				lepton1=labels_flipped[i];
			else if(std::abs(labels_flipped[i])>10 && std::abs(labels_flipped[i])<17)
				lepton2=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark1==0)
				quark1=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark2==0)
				quark2=labels_flipped[i];
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark3==0)
				quark3=labels_flipped[i];
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark4==0)
				quark4=labels_flipped[i];
		}
		if(quark1==0 || quark2==0 || quark3==0 || quark4==0 || lepton1==0 || lepton2==0){
			std::cout<<"wrong assignmnet of particle labels in 2q2Q1g2l amps\n";
			throw BHerror("Sorry, can't do that yet. (from 2q2Q1g2l)");
		}
		// Z & gamma cases
		else if((quark1==-quark3 && quark2==-quark4) || (quark1==-quark4 && quark2==-quark3)){
			if(quark1==quark2 && quark1%2==0)
				// uu
				case4q=0;
			else if(quark1==quark2)
				// dd
				case4q=5;
			else if(std::abs(quark1-quark2)%2==0 && quark1%2==0)
				// uup
				case4q=1;
			else if(std::abs(quark1-quark2)%2==0)
				// ddp
				case4q=4;
			else if(quark1%2==0)
				// ud
				case4q=2;
			else
				// du
				case4q=3;

			if(lepton1==-lepton2 && (std::abs(lepton1)==11 || std::abs(lepton1)==13 /* used to switch to photon only || std::abs(lepton1)==15*/)){
				// gamma+Z case
				photonZW=1;
			}
			else if(lepton1==-lepton2  && /* abuse of tauon for this switch */ (std::abs(lepton1)==15)){
				// forced gamma decaying into leptons
				photonZW=0;
			}
			else if(lepton1==-lepton2){
				// pure Z case decaying into neutrinos
				photonZW=2;
			}
		}
		// W cases
		else if(quark1==-quark3 && std::abs(quark2+quark4)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}
		else if(quark1==-quark4 && std::abs(quark2+quark3)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}
		else if(quark2==-quark3 && std::abs(quark1+quark4)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}
		else if(quark2==-quark4 && std::abs(quark1+quark3)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}

		else{
			std::cout<<"wrong assignmnet of particle labels in 2q2Q1g2l amps!\n";
			throw BHerror("Sorry, can't do that yet.");
		}
	  } break;
	  case 100040:{
		int quark1=-1,quark2=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark1<0)
				quark1=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6)
				quark2=labels_flipped[i];
		}
		if(quark1==quark2)
			case4q=0;
		else if(std::abs(quark1-quark2)%2==0)
			case4q=1;
		else
			case4q=2;
	  } break;
	  default: break;
	}
#if _VERBOSE
_PRINT(case4q);
_PRINT(photonZW);
_PRINT(up_down_quark);
#endif

//_PRINT(pc);

// with this we get the correct assigment of momenta --- notice that first to momenta are assumed to be initial
	std::vector<int> mom_assg;
	switch(pc){
		 case 41:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,g1pos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]==21&& g1pos==-1) {g1pos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos==-1){qpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos==-1){q2pos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qbpos==-1){qbpos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qb2pos==-1){qb2pos=i;}
		}
		mom_assg.push_back(qpos);
		mom_assg.push_back(qb2pos);
		mom_assg.push_back(q2pos);
		mom_assg.push_back(qbpos);
		mom_assg.push_back(g1pos);
	  } break;
       		 case 42:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,g1pos=-1,g2pos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]==21&& g1pos==-1) {g1pos=i;}
			if(labels_flipped[i]==21&& g2pos==-1) {g2pos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos==-1){qpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos==-1){q2pos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qbpos==-1){qbpos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qb2pos==-1){qb2pos=i;}
		}
		mom_assg.push_back(qpos);
		mom_assg.push_back(qb2pos);
		mom_assg.push_back(q2pos);
		mom_assg.push_back(qbpos);
		mom_assg.push_back(g1pos);
		mom_assg.push_back(g2pos);
	  } break;
		case 220: case 221: case 222: case 223: case 224:{
		int qpos=-1,qbpos=-1,lpos=-1,lbpos=-1;
		std::vector<int> gpos;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17)	{lpos=i;}
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)	{lbpos=i;}
			else if(labels_flipped[i]==21)	{gpos.push_back(i);}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6)	{qpos=i;}
			else	{qbpos=i;}
		}
		mom_assg.push_back(qpos);
		for(int k=0;k<gpos.size();k++)
			mom_assg.push_back(gpos[k]);
		mom_assg.push_back(qbpos);
		mom_assg.push_back(lpos);
		mom_assg.push_back(lbpos);
	  } break;
	  case 240:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,lpos=-1,lbpos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17)	{lpos=i;}
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)	{lbpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos==-1){qpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos==-1){q2pos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qbpos==-1){qbpos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qb2pos==-1){qb2pos=i;}
		}
		// Z/gamma case
		if(photonZW!=3){
			switch (case4q){
				case 0: case 5: break;
				default:{
					if(labels_flipped[qpos]!=-labels_flipped[qbpos]){
						int qbp(qbpos),qb2p(qb2pos);
						qbpos=qb2p;
						qb2pos=qbp;
					}
				} break;
			}
		}
		// W case
		else{
			if(labels_flipped[qpos]==-labels_flipped[qbpos]){
				int qp(qpos),q2p(q2pos);
				int qbp(qbpos),qb2p(qb2pos);
				qpos=q2p;
				q2pos=qp;
				qbpos=qb2p;
				qb2pos=qbp;
			}
			else if(labels_flipped[qpos]==-labels_flipped[qb2pos]){
				int qp(qpos),q2p(q2pos);
				qpos=q2p;
				q2pos=qp;
			}
			else if(labels_flipped[q2pos]==-labels_flipped[qbpos]){
				int qbp(qbpos),qb2p(qb2pos);
				qbpos=qb2p;
				qb2pos=qbp;
			}
		}
		mom_assg.push_back(qpos);
		mom_assg.push_back(qb2pos);
		mom_assg.push_back(q2pos);
		mom_assg.push_back(qbpos);
		mom_assg.push_back(lpos);
		mom_assg.push_back(lbpos);
	  } break;
	  case 241:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,g1pos=-1,lpos=-1,lbpos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17)	{lpos=i;}
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)	{lbpos=i;}
			else if(labels_flipped[i]==21) {g1pos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos==-1){qpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos==-1){q2pos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qbpos==-1){qbpos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qb2pos==-1){qb2pos=i;}
		}
		// Z/gamma case
		if(photonZW!=3){
			switch (case4q){
				case 0: case 5: break;
				default:{
					if(labels_flipped[qpos]!=-labels_flipped[qbpos]){
						int qbp(qbpos),qb2p(qb2pos);
						qbpos=qb2p;
						qb2pos=qbp;
					}
				} break;
			}
		}
		// W case
		else{
			if(labels_flipped[qpos]==-labels_flipped[qbpos]){
				int qp(qpos),q2p(q2pos);
				int qbp(qbpos),qb2p(qb2pos);
				qpos=q2p;
				q2pos=qp;
				qbpos=qb2p;
				qb2pos=qbp;
			}
			else if(labels_flipped[qpos]==-labels_flipped[qb2pos]){
				int qp(qpos),q2p(q2pos);
				qpos=q2p;
				q2pos=qp;
			}
			else if(labels_flipped[q2pos]==-labels_flipped[qbpos]){
				int qbp(qbpos),qb2p(qb2pos);
				qbpos=qb2p;
				qb2pos=qbp;
			}
		}
		mom_assg.push_back(qpos);
		mom_assg.push_back(qb2pos);
		mom_assg.push_back(q2pos);
		mom_assg.push_back(qbpos);
		mom_assg.push_back(g1pos);
		mom_assg.push_back(lpos);
		mom_assg.push_back(lbpos);
	  } break;
       	 case 100021: case 100022:{
		int qpos=-1,qbpos=-1,ypos=-1;
		std::vector<int> gpos;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]==22)	{ypos=i;}
			else if(labels_flipped[i]==21)	{gpos.push_back(i);}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6)	{qpos=i;}
			else	{qbpos=i;}
		}
		mom_assg.push_back(qpos);
		for(int k=0;k<gpos.size();k++)
			mom_assg.push_back(gpos[k]);
		mom_assg.push_back(qbpos);
		mom_assg.push_back(ypos);
	  } break;
	  case 100040:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,ypos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]==22)	{ypos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos<0){
				qpos=i;
				for(int j=0;j<labels_flipped.size();j++){
					if(labels_flipped[i]==-labels_flipped[j]){
						qbpos=j;
						break;
					}
				}
			}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos<0){
				q2pos=i;
				for(int j=0;j<labels_flipped.size();j++){
					if(labels_flipped[i]==-labels_flipped[j] && j!=qbpos){
						qb2pos=j;
						break;
					}
				}
			}
		}
		mom_assg.push_back(qpos);
		mom_assg.push_back(qb2pos);
		mom_assg.push_back(q2pos);
		mom_assg.push_back(qbpos);
		mom_assg.push_back(ypos);
	  } break;
	  default:{
		mom_assg.push_back(0);
		mom_assg.push_back(1);
		mom_assg.push_back(2);
		mom_assg.push_back(3);
		mom_assg.push_back(4);
		mom_assg.push_back(5);
		mom_assg.push_back(6);
		mom_assg.push_back(7);
	  } break;
	}
#if _VERBOSE
for(int kk=0;kk<mom_assg.size();kk++)
	_PRINT(mom_assg[kk]);
#endif


int color,tree_color;

switch (settings::BH_interface_settings::s_BH_color_mode)
	{
        case settings::BH_interface_settings::full_color: {color=0;tree_color=0;}; break;
        case settings::BH_interface_settings::leading_color : {color=1;
		switch (settings::BH_interface_settings::s_BH_tree_color_mode){
        		case settings::BH_interface_settings::full_color : {tree_color=0;}; break;
        		case settings::BH_interface_settings::leading_color : {tree_color=1;}; break;
			};
		}; break;
        case settings::BH_interface_settings::full_minus_leading_color : {color=2;tree_color=0;}; break;
	}



switch (pc){
case 41:return new BH_Ampl_2q2Q1g(case4q,color,tree_color,mom_assg);
case 42:return new BH_Ampl_2q2Q2g(case4q,color,tree_color,mom_assg);
//case 220:return new BH_Ampl_concrete(pc);
//case 221:return new BH_Ampl_ee3jet(mom_assg);
case 220:return new BH_Ampl_2q2l(up_down_quark,photonZW,color,tree_color,mom_assg);
case 221:return new BH_Ampl_2q1g2l(up_down_quark,photonZW,color,tree_color,mom_assg);
case 222:return new BH_Ampl_2q2g2l(up_down_quark,photonZW,color,tree_color,mom_assg);
case 223:return new BH_Ampl_2q3g2l(up_down_quark,photonZW,color,tree_color,mom_assg);
case 224:return new BH_Ampl_2q4g2l(up_down_quark,photonZW,color,tree_color,mom_assg);
case 240:return new BH_Ampl_2q2Q2l(photonZW,case4q,color,tree_color,mom_assg);
case 241:return new BH_Ampl_2q2Q1g2l(photonZW,case4q,color,tree_color,mom_assg);
case 100021: return new BH_Ampl_2q1g1y(up_down_quark,color,mom_assg);
case 100022: return new BH_Ampl_2q2g1y(up_down_quark,color,tree_color,mom_assg);
case 100040:	return new BH_Ampl_2q2Q1y(up_down_quark,case4q,color,tree_color,mom_assg);
default : throw BHerror("Sorry, can't do that yet.");
}

}

#endif

#if 0
BH_Ampl* BH_interface::new_ampl(const std::vector<int>& labels){
	// gluon 21
	// photon 22
	// d 1
	// u 2
	//s c b t
	//e- 11
	//ne 12
	//mu- 13
	//nmu 14
	//tau- 15
	//ntau- 16
//_MESSAGE("BH_interface::new_ampl CALLED");
//_MESSAGE("library in ~ffebres/workspace/W3jdists/blackhat-lib-0.6.3/");

  d_ampl_index++;
  std::cout<<"\n// process "<<d_ampl_index<<":"<<endl;
  std::cout<<"	int array"<<d_ampl_index<<"[]={";
	for(int jj=0;jj<(labels.size()-1);jj++)
		std::cout<<labels[jj]<<",";
		std::cout<<labels[labels.size()-1]<<"};\n";
  std::cout<<"	std::vector<int> prop_n"<<d_ampl_index<<"(array"<<d_ampl_index<<",array"<<d_ampl_index<<"+sizeof(array"<<d_ampl_index<<")/sizeof(int));\n";


	// updates the values of the constant
	prop_hel_fn::update_constants(d_settings_p);

	std::vector<particle> types;

size_t nbr_gluons=0,nbr_quarks=0,nbr_leptons=0,nbr_photon=0;

for (int i=0;i<labels.size();i++){
	switch (std::abs(labels[i])){
		case 1: case 2: case 3: case 4: case 5: types.push_back(quark); ++nbr_quarks; break;
		case 11: case 13: case 15: types.push_back(lepton); ++nbr_leptons;break;
		case 12: case 14: case 16: types.push_back(lepton); ++nbr_leptons;break;
		case 21: types.push_back(gluon);++nbr_gluons; break;
		case 22: types.push_back(photon);++nbr_photon; break;
	}
}

size_t pc=nbr_gluons+10*nbr_quarks+100*nbr_leptons+100000*nbr_photon;

	int case4q=-1;
	bool up_down_quark=0;
	// default zero, which is for only off-shell photon
	int photonZW=0;
	// labels_flipped corrects for the all outgoing BH conventions if initial fermions are present
	std::vector<int> labels_flipped;
		if((std::abs(labels[0])>0 && std::abs(labels[0])<6) || (std::abs(labels[0])>10 && std::abs(labels[0])<17))
			labels_flipped.push_back(-labels[0]);
		else
			labels_flipped.push_back(labels[0]);
		if((std::abs(labels[1])>0 && std::abs(labels[1])<6) || (std::abs(labels[1])>10 && std::abs(labels[1])<17))
			labels_flipped.push_back(-labels[1]);
		else
			labels_flipped.push_back(labels[1]);
		for(int k=2;k<labels.size();k++)
			labels_flipped.push_back(labels[k]);
	switch(pc){
	  case 220: case 221: case 222: case 223: case 224:{
		int lepton1=0,lepton2=0;
		int quark1=0,quark2=0;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17 && lepton1==0)
				lepton1=labels_flipped[i];
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)
				lepton2=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark1==0)
				quark1=labels_flipped[i];
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark2==0)
				quark2=labels_flipped[i];
		}
		if(quark1==0 || quark2==0 || lepton1==0 || lepton2==0){
			std::cout<<"wrong assignmnet of particle labels in 2q2g2l amps\n";
			throw BHerror("Sorry, can't do that yet.");
		}
		else if(lepton1==-lepton2 && (std::abs(lepton1)==11 || std::abs(lepton1)==13 /*|| std::abs(lepton1)==15*/)){
			// gamma+Z case
			photonZW=1;
		}
		else if(lepton1==-lepton2  && /* abuse of tauon for this switch */ (std::abs(lepton1)==15)){
			// forced gamma decaying into leptons
			photonZW=0;
		}
		else if(lepton1==-lepton2  && (std::abs(lepton1)==12 || std::abs(lepton1)==14 || std::abs(lepton1)==16)){
			// pure Z case decaying into neutrinos
			photonZW=2;
		}
		else if((((lepton1==11 || lepton1==13 || lepton1==15) && (lepton1+lepton2)==-1) ||
			((lepton1==12 || lepton1==14 || lepton1==16) && (lepton1+lepton2)==1)) &&
			std::abs(quark1+quark2)%2==1){
			// W case
			photonZW=3;
		}
		else{
			std::cout<<"wrong assignmnet of particle labels in 2q2g2l amps\n";
			throw BHerror("Sorry, can't do that yet.");
		}
		// W case keeps default value up_down_quark=0
		if(std::abs(quark1)==std::abs(quark2)){
			switch (std::abs(quark1)){
				// up quarks
				case 2: case 4:	up_down_quark=0; break;
				// down quarks
				case 1: case 3: case 5:up_down_quark=1; break;
			}
		}
	  } break;
	  case 240: case 241:{
		// needs modification for W --- so far only drops the run
		int lepton1=0,lepton2=0;
		int quark1=0,quark2=0,quark3=0,quark4=0;
		for(int i=0;i<labels_flipped.size();i++){
			if(std::abs(labels_flipped[i])>10 && std::abs(labels_flipped[i])<17 && lepton1==0)
				lepton1=labels_flipped[i];
			else if(std::abs(labels_flipped[i])>10 && std::abs(labels_flipped[i])<17)
				lepton2=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark1==0)
				quark1=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark2==0)
				quark2=labels_flipped[i];
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark3==0)
				quark3=labels_flipped[i];
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && quark4==0)
				quark4=labels_flipped[i];
		}
		if(quark1==0 || quark2==0 || quark3==0 || quark4==0 || lepton1==0 || lepton2==0){
			std::cout<<"wrong assignmnet of particle labels in 2q2Q2l amps\n";
			throw BHerror("Sorry, can't do that yet.");
		}
		// Z & gamma cases
		else if((quark1==-quark3 && quark2==-quark4) || (quark1==-quark4 && quark2==-quark3)){
			if(quark1==quark2 && quark1%2==0)
				// uu
				case4q=0;
			else if(quark1==quark2)
				// dd
				case4q=5;
			else if(std::abs(quark1-quark2)%2==0 && quark1%2==0)
				// uup
				case4q=1;
			else if(std::abs(quark1-quark2)%2==0)
				// ddp
				case4q=4;
			else if(quark1%2==0)
				// ud
				case4q=2;
			else
				// du
				case4q=3;

			if(lepton1==-lepton2 && (std::abs(lepton1)==11 || std::abs(lepton1)==13 /* used to switch to photon only || std::abs(lepton1)==15*/)){
				// gamma+Z case
				photonZW=1;
			}
			else if(lepton1==-lepton2  && /* abuse of tauon for this switch */ (std::abs(lepton1)==15)){
				// forced gamma decaying into leptons
				photonZW=0;
			}
			else if(lepton1==-lepton2){
				// pure Z case decaying into neutrinos
				photonZW=2;
			}
		}
		// W cases
		else if(quark1==-quark3 && std::abs(quark2+quark4)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}
		else if(quark1==-quark4 && std::abs(quark2+quark3)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}
		else if(quark2==-quark3 && std::abs(quark1+quark4)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}
		else if(quark2==-quark4 && std::abs(quark1+quark3)%2==1){
			photonZW=3;
			if(quark1==quark2)	case4q=0;
			else if(quark3==quark4)	case4q=1;
			else	case4q=2;
		}

		else{
			std::cout<<"wrong assignmnet of particle labels in 2q2Q2l amps!\n";
			throw BHerror("Sorry, can't do that yet.");
		}
	  } break;
	  case 100040:{
		int quark1=-1,quark2=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>0 && labels_flipped[i]<6 && quark1<0)
				quark1=labels_flipped[i];
			else if(labels_flipped[i]>0 && labels_flipped[i]<6)
				quark2=labels_flipped[i];
		}
		if(quark1==quark2)
			case4q=0;
		else if(std::abs(quark1-quark2)%2==0)
			case4q=1;
		else
			case4q=2;
	  } break;
	  default: break;
	}
#if _VERBOSE
_PRINT(case4q);
_PRINT(photonZW);
_PRINT(up_down_quark);
#endif

//_PRINT(pc);

// with this we get the correct assigment of momenta --- notice that first to momenta are assumed to be initial
	std::vector<int> mom_assg;
	switch(pc){
		case 220: case 221: case 222: case 223: case 224:{
		int qpos=-1,qbpos=-1,lpos=-1,lbpos=-1;
		std::vector<int> gpos;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17)	{lpos=i;}
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)	{lbpos=i;}
			else if(labels_flipped[i]==21)	{gpos.push_back(i);}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6)	{qpos=i;}
			else	{qbpos=i;}
		}
		mom_assg.push_back(qpos+1);
		for(int k=0;k<gpos.size();k++)
			mom_assg.push_back(gpos[k]+1);
		mom_assg.push_back(qbpos+1);
		mom_assg.push_back(lpos+1);
		mom_assg.push_back(lbpos+1);
	  } break;
	  case 240:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,lpos=-1,lbpos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17)	{lpos=i;}
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)	{lbpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos==-1){qpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos==-1){q2pos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qbpos==-1){qbpos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qb2pos==-1){qb2pos=i;}
		}
		// Z/gamma case
		if(photonZW!=3){
			switch (case4q){
				case 0: case 5: break;
				default:{
					if(labels_flipped[qpos]!=-labels_flipped[qbpos]){
						int qbp(qbpos),qb2p(qb2pos);
						qbpos=qb2p;
						qb2pos=qbp;
					}
				} break;
			}
		}
		// W case
		else{
			if(labels_flipped[qpos]==-labels_flipped[qbpos]){
				int qp(qpos),q2p(q2pos);
				int qbp(qbpos),qb2p(qb2pos);
				qpos=q2p;
				q2pos=qp;
				qbpos=qb2p;
				qb2pos=qbp;
			}
			else if(labels_flipped[qpos]==-labels_flipped[qb2pos]){
				int qp(qpos),q2p(q2pos);
				qpos=q2p;
				q2pos=qp;
			}
			else if(labels_flipped[q2pos]==-labels_flipped[qbpos]){
				int qbp(qbpos),qb2p(qb2pos);
				qbpos=qb2p;
				qb2pos=qbp;
			}
		}
		// NOTICE that for BH_interface mom_assg is base 1, as it'll be passed as input ind
		mom_assg.push_back(qpos+1);
		mom_assg.push_back(qb2pos+1);
		mom_assg.push_back(q2pos+1);
		mom_assg.push_back(qbpos+1);
		mom_assg.push_back(lpos+1);
		mom_assg.push_back(lbpos+1);
	  } break;
	  case 241:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,g1pos=-1,lpos=-1,lbpos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]>10 && labels_flipped[i]<17)	{lpos=i;}
			else if(labels_flipped[i]<-10 && labels_flipped[i]>-17)	{lbpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos==-1){qpos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos==-1){q2pos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qbpos==-1){qbpos=i;}
			else if(labels_flipped[i]<0 && labels_flipped[i]>-6 && qb2pos==-1){qb2pos=i;}
			else if(labels_flipped[i]==21 && g1pos==-1){g1pos=i;}
		}
		// Z/gamma case
		if(photonZW!=3){
			switch (case4q){
				case 0: case 5: break;
				default:{
					if(labels_flipped[qpos]!=-labels_flipped[qbpos]){
						int qbp(qbpos),qb2p(qb2pos);
						qbpos=qb2p;
						qb2pos=qbp;
					}
				} break;
			}
		}
		// W case
		else{
			if(labels_flipped[qpos]==-labels_flipped[qbpos]){
				int qp(qpos),q2p(q2pos);
				int qbp(qbpos),qb2p(qb2pos);
				qpos=q2p;
				q2pos=qp;
				qbpos=qb2p;
				qb2pos=qbp;
			}
			else if(labels_flipped[qpos]==-labels_flipped[qb2pos]){
				int qp(qpos),q2p(q2pos);
				qpos=q2p;
				q2pos=qp;
			}
			else if(labels_flipped[q2pos]==-labels_flipped[qbpos]){
				int qbp(qbpos),qb2p(qb2pos);
				qbpos=qb2p;
				qb2pos=qbp;
			}
		}
		// NOTICE that for BH_interface mom_assg is base 1, as it'll be passed as input ind
		mom_assg.push_back(qpos+1);
		mom_assg.push_back(qb2pos+1);
		mom_assg.push_back(q2pos+1);
		mom_assg.push_back(qbpos+1);
		mom_assg.push_back(g1pos+1);
		mom_assg.push_back(lpos+1);
		mom_assg.push_back(lbpos+1);
	  } break;
       	 case 100021: case 100022:{
		int qpos=-1,qbpos=-1,ypos=-1;
		std::vector<int> gpos;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]==22)	{ypos=i;}
			else if(labels_flipped[i]==21)	{gpos.push_back(i);}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6)	{qpos=i;}
			else	{qbpos=i;}
		}
		mom_assg.push_back(qpos+1);
		for(int k=0;k<gpos.size();k++)
			mom_assg.push_back(gpos[k]+1);
		mom_assg.push_back(qbpos+1);
		mom_assg.push_back(ypos+1);
	  } break;
	  case 100040:{
		int qpos=-1,qbpos=-1,q2pos=-1,qb2pos=-1,ypos=-1;
		for(int i=0;i<labels_flipped.size();i++){
			if(labels_flipped[i]==22)	{ypos=i;}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && qpos<0){
				qpos=i;
				for(int j=0;j<labels_flipped.size();j++){
					if(labels_flipped[i]==-labels_flipped[j]){
						qbpos=j;
						break;
					}
				}
			}
			else if(labels_flipped[i]>0 && labels_flipped[i]<6 && q2pos<0){
				q2pos=i;
				for(int j=0;j<labels_flipped.size();j++){
					if(labels_flipped[i]==-labels_flipped[j] && j!=qbpos){
						qb2pos=j;
						break;
					}
				}
			}
		}
		// NOTICE that for BH_interface mom_assg is base 1, as it'll be passed as input ind
		mom_assg.push_back(qpos+1);
		mom_assg.push_back(qb2pos+1);
		mom_assg.push_back(q2pos+1);
		mom_assg.push_back(qbpos+1);
		mom_assg.push_back(ypos+1);
	  } break;
	  default:{
		// NOTICE that for BH_interface mom_assg is base 1, as it'll be passed as input ind
		mom_assg.push_back(0+1);
		mom_assg.push_back(1+1);
		mom_assg.push_back(2+1);
		mom_assg.push_back(3+1);
		mom_assg.push_back(4+1);
		mom_assg.push_back(5+1);
		mom_assg.push_back(6+1);
		mom_assg.push_back(7+1);
	  } break;
	}

#if _VERBOSE
for(int kk=0;kk<mom_assg.size();kk++)
	_PRINT(mom_assg[kk]);
#endif

int color,tree_color;

switch (settings::BH_interface_settings::s_BH_color_mode)
	{
        case settings::BH_interface_settings::full_color: {color=0;tree_color=0;}; break;
        case settings::BH_interface_settings::leading_color : {color=1;
		switch (settings::BH_interface_settings::s_BH_tree_color_mode){
        		case settings::BH_interface_settings::full_color : {tree_color=0;}; break;
        		case settings::BH_interface_settings::leading_color : {tree_color=1;}; break;
			};
		}; break;
        case settings::BH_interface_settings::full_minus_leading_color : {color=2;tree_color=0;}; break;
	}



switch (pc){
case 41:{
	// to use IR_checked_OLHA's instead of OneLoopHelAmpl's
	CachedOLHA::use_IR_checked();
	return new BH_Ampl_2q2Q1g(case4q,color,tree_color,mom_assg,this,d_ampl_index);
	} break;
//case 220:return new BH_Ampl_concrete(pc);
//case 221:return new BH_Ampl_ee3jet(mom_assg);
case 220:return new BH_Ampl_2q2l(up_down_quark,photonZW,color,tree_color,mom_assg,this);
case 221:return new BH_Ampl_2q1g2l(up_down_quark,photonZW,color,tree_color,mom_assg,this);
case 222:return new BH_Ampl_2q2g2l(up_down_quark,photonZW,color,tree_color,mom_assg,this);
case 223:{
	// to use IR_checked_OLHA's instead of OneLoopHelAmpl's
	CachedOLHA::use_IR_checked();
	return new BH_Ampl_2q3g2l(up_down_quark,photonZW,color,tree_color,mom_assg,this,d_ampl_index);
	} break;
case 224:{
	// to use IR_checked_OLHA's instead of OneLoopHelAmpl's
	CachedOLHA::use_IR_checked();
	return new BH_Ampl_2q4g2l(up_down_quark,photonZW,color,tree_color,mom_assg,this,d_ampl_index);
	} break;
case 240:return new BH_Ampl_2q2Q2l(photonZW,case4q,color,tree_color,mom_assg,this);
case 241:{
	// to use IR_checked_OLHA's instead of OneLoopHelAmpl's
	CachedOLHA::use_IR_checked();
	return new BH_Ampl_2q2Q1g2l(photonZW,case4q,color,tree_color,mom_assg,this,d_ampl_index);
	} break;
case 100021: return new BH_Ampl_2q1g1y(up_down_quark,color,mom_assg,this);
case 100022: return new BH_Ampl_2q2g1y(up_down_quark,color,tree_color,mom_assg,this);
case 100040:	return new BH_Ampl_2q2Q1y(up_down_quark,case4q,color,tree_color,mom_assg,this);
default : throw BHerror("Sorry, can't do that yet.");
}

}

template <class T> void BH_interface::set(const std::string& name,T value){
	d_settings_p->set(name,value);
}

BH_interface::BH_interface() : d_settings_p(new settings_table()),d_phase_space_point(0),d_ampl_index(0),
		d_file_base(200), d_file_extension(0) {
	d_settings_p->add("Z_mass",91.188);
	d_settings_p->add("Z_width",2.49);
	d_settings_p->add("W_mass",80.419);
	d_settings_p->add("W_width",2.06);
	d_settings_p->add("sin_th_2",0.23);
	d_settings_p->add("sin_2th",0.8416650165000326);
	d_settings_p->add("alpha_S",0.118);
	d_settings_p->add("alpha_QED",1./128.802);
	//
    d_settings_p->add("G3_Lambda2",10e5);
	d_mc_p= new momentum_configuration<R>();

	// to employ cached integrals
	QCD_cut_part_factory::s_use_cached_integrals=true;

	///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	// for runs when producing/using stored ME's
	///////////////////////////////////////////////////////

	// this integer defines wheter the run is for:
	//		0:	regular BH+Sherpa run
	//		1:	warming up grids, i.e. ME's (normalized by tree) returned are set to 1
	//		2:	writing PS points
	//		3:	returning stored matrix elements
//	BH_interface::d_which_study=3;

	switch (d_which_study){
		case 2:{
			d_ps_output_file.open("PSpointsOUTW3jvirt_test.dat",ios::app);
			d_ps_output_file<<setprecision(16);
			d_ps_output_file<<"testing open"<<endl;
		} break;
		case 3:{
			d_file_base=200;
			d_file_extension=0;
		} break;
		default: break;
	}

}


BH_interface::~BH_interface(){

	switch (d_which_study){
		case 2:{
			d_ps_output_file<<"testing close"<<endl;
			if(d_ps_output_file.is_open())
				d_ps_output_file.close();
		} break;
		case 3:{
			if(d_matrix_elements.is_open())
				d_matrix_elements.close();
		} break;
		default: break;
	}

	delete d_settings_p;
}



void BH_interface::print_settings(){
	_MESSAGE("-- settings list --------");
	d_settings_p->display();
	_MESSAGE("-------------------------");
}

void BH_interface::operator()(BHinput& in){
	if(d_which_study==0 || d_which_study==1){
	d_mc_p->clear();
	mom_conf& mc=*d_mc_p;

///////////////////////////////////////////////////////////////////
// rescaling momenta and dimensionful objects to O(1) units
	R _GeV;	_GeV=double(in.m_momenta.size())/abs(in.m_momenta[0][0]+in.m_momenta[1][0]);
//_PRINT(_GeV);
	prop_hel_fn::set_GeV(_GeV);
///////////////////////////////////////////////////////////////////

		// flip sign of incoming momenta (first two passed)
		mc.insert(Cmom<double>(-in.m_momenta[0][0]*_GeV,
				-in.m_momenta[0][2]*_GeV,
				-in.m_momenta[0][3]*_GeV,
				-in.m_momenta[0][1]*_GeV));
		mc.insert(Cmom<double>(-in.m_momenta[1][0]*_GeV,
				-in.m_momenta[1][2]*_GeV,
				-in.m_momenta[1][3]*_GeV,
				-in.m_momenta[1][1]*_GeV));
		for(int k=2;k<in.m_momenta.size();k++){
			mc.insert(Cmom<double>(in.m_momenta[k][0]*_GeV,
					in.m_momenta[k][2]*_GeV,
					in.m_momenta[k][3]*_GeV,
					in.m_momenta[k][1]*_GeV));
		}

	d_mu=in.m_mu*_GeV;
  }   else if(d_which_study==2){
		d_mc_p->clear();
		mom_conf& mc=*d_mc_p;

			// flip sign of incoming momenta (first two passed)
			mc.insert(Cmom<double>(-in.m_momenta[0][0],-in.m_momenta[0][1],-in.m_momenta[0][2],-in.m_momenta[0][3]));
			mc.insert(Cmom<double>(-in.m_momenta[1][0],-in.m_momenta[1][1],-in.m_momenta[1][2],-in.m_momenta[1][3]));
			for(int k=2;k<in.m_momenta.size();k++){
				mc.insert(Cmom<double>(in.m_momenta[k][0],in.m_momenta[k][1],in.m_momenta[k][2],in.m_momenta[k][3]));
			}

		d_ps_output_file<<setprecision(16);
		d_ps_output_file<<in.m_momenta[0][0]<<" "<<in.m_momenta[0][2]<<" "<<in.m_momenta[0][3]<<" "<<in.m_momenta[0][1]<<" ";
		d_ps_output_file<<in.m_momenta[1][0]<<" "<<in.m_momenta[1][2]<<" "<<in.m_momenta[1][3]<<" "<<in.m_momenta[1][1]<<" ";
		for(int kk=2;kk<in.m_momenta.size();kk++){
			d_ps_output_file<<in.m_momenta[kk][0]<<" "<<in.m_momenta[kk][2]<<" "<<in.m_momenta[kk][3]<<" "<<in.m_momenta[kk][1]<<" ";
		}
		d_ps_output_file<<endl;

		d_mu=in.m_mu;
	  }
	  else if(d_which_study==3){
		// notice that his works only for 7-pt, although generalization is straightforward
		const char * filename;
		if(d_phase_space_point%d_file_base==0){
			std::string sfilename="data/W3jvirt.";
			std::string s;
			std::stringstream out;
			out << d_file_extension;
			sfilename+=out.str();
			filename=sfilename.c_str();

			if(d_matrix_elements.is_open())
				d_matrix_elements.close();
			d_matrix_elements.open(filename,ios::in);

			d_file_extension++;
		}
		d_phase_space_point++;
		d_mc_p->clear();
		for(int ii=0;ii<7;ii++){
			double e,x,y,z;
			d_matrix_elements>>e;
			d_matrix_elements>>x;
			d_matrix_elements>>y;
			d_matrix_elements>>z;
			d_mc_p->insert(Cmom<R>(e,x,y,z));
		}

		d_mu=in.m_mu;

	// SANITY check
	if((d_phase_space_point-1)%10000==0){
		std::cout<<setprecision(16);
		_PRINT(filename);
		_PRINT(d_phase_space_point);
		_PRINT(in.m_mu);
		_PRINT(d_mc_p->p(1));
	// getting the PS points
		std::cout<<in.m_momenta[0][0]<<" "<<in.m_momenta[0][2]<<" "<<in.m_momenta[0][3]<<" "<<in.m_momenta[0][1]<<" ";
		std::cout<<in.m_momenta[1][0]<<" "<<in.m_momenta[1][2]<<" "<<in.m_momenta[1][3]<<" "<<in.m_momenta[1][1]<<" ";
		for(int kk=2;kk<d_mc_p->n();kk++){
			std::cout<<in.m_momenta[kk][0]<<" "<<in.m_momenta[kk][2]<<" "<<in.m_momenta[kk][3]<<" "<<in.m_momenta[kk][1]<<" ";
		}
		std::cout<<endl;
	}
	// end of SANITIY check
	  }
}


#endif

double BH_Ampl::operator ()(BHinput& in){
	return 0.;
}

double BH_Ampl::get_single_pole(){return 0.;}
double BH_Ampl::get_double_pole(){return 0.;}

#if 0

template <class T> void BH_factory::set(const std::string& name,T value){
	d_settings_p->set(name,value);
}


BH_factory::BH_factory() : d_settings_p(new settings_table()) {
	d_settings_p->add("Z_mass",91.188);
	d_settings_p->add("Z_width",2.49);
	d_settings_p->add("W_mass",80.419);
	d_settings_p->add("W_width",2.06);
	d_settings_p->add("sin_th_2",0.23);
	d_settings_p->add("sin_2th",0.8416650165000326);
	d_settings_p->add("alpha_S",0.118);
	d_settings_p->add("alpha_QED",1./128.802);
}



BH_factory::~BH_factory(){
	delete d_settings_p;
}



void BH_factory::print_settings(){
	_MESSAGE("-- settings list --------");
	d_settings_p->display();
	_MESSAGE("-------------------------");
}


template void BH_factory::set(const std::string& name,double value);
template void BH_factory::set(const std::string& name,int value);
template void BH_factory::set(const std::string& name,bool value);

#endif

#if 0
template void BH_interface::set(const std::string& name,double value);
template void BH_interface::set(const std::string& name,int value);
template void BH_interface::set(const std::string& name,bool value);
#endif

}





