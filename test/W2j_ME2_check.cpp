/*
 *
 * program to check: 
 * 	(1) ME2 of W+2j process
 * 
 * PS-point: gkm6
 * 
*/

#include <mom_conf.h>
#include "assembly.h"
#include "settings.h"
#include "settings_reader.h"

#include "cached_OLHA.h"
#include "Interface/BH_interface.h"
#include "Interface/BH_Ampl.h"
#include "polylog.h"  // for pi

#include <ctime>


#define GET_BORN 1 	//[0,1] for [don't use,use] get_born() to evaluate morn ME2.

#define write_TEX 0 	//[0,1] for [no,yes] LaTeX output.
#define write_C 0 	//[0,1] for [no,yes] c++ usable output.
#define do_CHECKS 1	//[0,1] for [no,yes] performing checks.


// to choose which partial amplitudes and ME2 should we computed
#define CHECK_ME2_4q_dist_quark 1
#define CHECK_ME2_4q_ident_quark 1
#define CHECK_ME2_4q_ident_antiq 1

#define CHECK_ME2_2q2g 1



#define EPS 0.00000000000001
#define _CHECK(X,Y,tol,mess,err) if ( abs( 1. - Y/(X+EPS)) < tol ){ _MESSAGE2(mess," OK");	} else { _MESSAGE5(mess," NOT OK:",X," != ", Y );err ++;	}

using namespace BH;
using namespace std;

string isgn(C a){if(imag(a)>0.){return(" + i ");} else {return(" - i ");};}


int main(void) {
        	
		unsigned int old_cw;
	        fpu_fix_start(&old_cw);

int err(0);
double tol(0.00001);
std::cout<<setprecision(10);

#include "std_momenta.hpp"
settings::use_setting("COLOR_MODE full_color");

  ind6.clear();
  ind6.push_back(1);
  ind6.push_back(2);
  ind6.push_back(3);
  ind6.push_back(4);
  ind6.push_back(5);
  ind6.push_back(6);

	vector<double> p1;
     	vector<double> p2;
        vector<double> p3;
        vector<double> p4;
        vector<double> p5;
        vector<double> p6;

	vector<int> particles;
	
	double mu_momenta6pt(6.);
	double c(1.);
p1.push_back(-gkm6.p(1).E().real()*c);p1.push_back(-gkm6.p(1).X().real()*c);p1.push_back(-gkm6.p(1).Y().real()*c);p1.push_back(-gkm6.p(1).Z().real()*c);
p2.push_back(-gkm6.p(2).E().real()*c);p2.push_back(-gkm6.p(2).X().real()*c);p2.push_back(-gkm6.p(2).Y().real()*c);p2.push_back(-gkm6.p(2).Z().real()*c);
p3.push_back(gkm6.p(3).E().real()*c);p3.push_back(gkm6.p(3).X().real()*c);p3.push_back(gkm6.p(3).Y().real()*c);p3.push_back(gkm6.p(3).Z().real()*c);
p4.push_back(gkm6.p(4).E().real()*c);p4.push_back(gkm6.p(4).X().real()*c);p4.push_back(gkm6.p(4).Y().real()*c);p4.push_back(gkm6.p(4).Z().real()*c);
p5.push_back(gkm6.p(5).E().real()*c);p5.push_back(gkm6.p(5).X().real()*c);p5.push_back(gkm6.p(5).Y().real()*c);p5.push_back(gkm6.p(5).Z().real()*c);
p6.push_back(gkm6.p(6).E().real()*c);p6.push_back(gkm6.p(6).X().real()*c);p6.push_back(gkm6.p(6).Y().real()*c);p6.push_back(gkm6.p(6).Z().real()*c);


        vector<vector <double> > momenta6pt;
        momenta6pt.push_back(p1);
        momenta6pt.push_back(p2);
        momenta6pt.push_back(p3);
        momenta6pt.push_back(p4);
        momenta6pt.push_back(p5);
        momenta6pt.push_back(p6);

	C me2_finite;
	C me2_single_pole;
	C me2_double_pole;
	C me2_born;

	me2_finite=C(0.,0.);
	me2_single_pole=C(0.,0.);
	me2_double_pole=C(0.,0.);
	me2_born=C(0.,0.);

	particles.clear();

	BH_interface bhi;
	BH_Ampl* me2;
#include "theconstants.h"

vector< vector<C> > res_me2(3);
vector<C> res_me2_born;

// ME2 distinct quarks 4q2l 
res_me2[0].push_back(C(32.82387766,0));
res_me2[1].push_back(C(-13.07271972,0));
res_me2[2].push_back(C(-5.333333333,0));
res_me2_born.push_back(C(3.679728142e-06,0));
// ME2 ident quarks 4q2l 
res_me2[0].push_back(C(32.21883588,0));
res_me2[1].push_back(C(-13.12542597,0));
res_me2[2].push_back(C(-5.333333333,0));
res_me2_born.push_back(C(3.78041248e-06,0));
// ME2 identical anti-quarks 4q2l 
res_me2[0].push_back(C(26.53713426,0));
res_me2[1].push_back(C(-14.01005141,0));
res_me2[2].push_back(C(-5.333333333,0));
res_me2_born.push_back(C(4.878681245e-06,0));
// ME2 2q2g2l 
res_me2[0].push_back(C(3.128007358,0));
res_me2[1].push_back(C(-25.26180849,0));
res_me2[2].push_back(C(-8.666666667,0));
res_me2_born.push_back(C(2.704287436e-05,0));



#if CHECK_ME2_4q_dist_quark 

/*****************************************************************/
/*****************************************************************/
  std::cout<< "**********************************************" <<std::endl;
   	cout << " ME2 distinct quarks 4q2l " << endl;
{
	
                particles.push_back(-2);
		particles.push_back(3);
                particles.push_back(3);
                particles.push_back(-1);
                particles.push_back(11);
                particles.push_back(-12);
        
	BHinput input(momenta6pt,mu_momenta6pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();

#if GET_BORN
	me2_born=me2->get_born();
#endif

_PRINT(me2_finite);
_PRINT(me2_single_pole);
_PRINT(me2_double_pole);

#if GET_BORN
_PRINT(me2_born);
#endif

  std::cout<< "-----------------------------------------------" <<std::endl;

  
#if write_TEX 
  cout<<fixed<< "$(\\bar u c -> c \\bar d e^- \\abr\\nu)$"  <<
	  "  & $ "<<me2_double_pole << 
	  " $ & $ "<<me2_single_pole<< 
	  " $ & $ "<<me2_finite<< 
	  " $\\\\ \\hline "<<endl;
  cout<<
	  "  & $ "<< 
	  " $ & $ "<< 
	  " $ & $ "<<me2_born<< 
	  " $\\\\ \\hline "<<endl;
#endif

#if write_C
	cout<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][0],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][0],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][0],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[0],me2_born,tol,": born",err) 
#endif
  
  std::cout<< "-----------------------------------------------" <<std::endl;
}
  std::cout<< "**********************************************" <<std::endl;
/*****************************************************************/
#endif


#if CHECK_ME2_4q_ident_quark 

/*****************************************************************/
/*****************************************************************/
  std::cout<< "**********************************************" <<std::endl;
   	cout << " ME2 ident quarks 4q2l " << endl;
{
	particles.clear();	
                particles.push_back(-2);
		particles.push_back(2);
                particles.push_back(2);
                particles.push_back(-1);
                particles.push_back(11);
                particles.push_back(-12);

        BHinput input(momenta6pt,mu_momenta6pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();

#if GET_BORN
	me2_born=me2->get_born();
#endif

_PRINT(me2_finite);
_PRINT(me2_single_pole);
_PRINT(me2_double_pole);
#if GET_BORN
_PRINT(me2_born);
#endif

  std::cout<< "-----------------------------------------------" <<std::endl;

  
#if write_TEX 
  cout<<fixed<< "$(\\bar u u -> u \\bar d e^- \\abr\\nu)$"  <<
	  "  & $ "<<me2_double_pole << 
	  " $ & $ "<<me2_single_pole<< 
	  " $ & $ "<<me2_finite<< 
	  " $\\\\ \\hline "<<endl;
  cout<<
	  "  & $ "<< 
	  " $ & $ "<< 
	  " $ & $ "<<me2_born<< 
	  " $\\\\ \\hline "<<endl;
#endif

#if write_C
	cout<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][1],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][1],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][1],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[1],me2_born,tol,": born",err) 
#endif
  
  std::cout<< "-----------------------------------------------" <<std::endl;
}
  std::cout<< "**********************************************" <<std::endl;
/*****************************************************************/
#endif


#if CHECK_ME2_4q_ident_antiq 
/*****************************************************************/
/*****************************************************************/
  std::cout<< "**********************************************" <<std::endl;
   	cout << " ME2 identical anti-quarks 4q2l " << endl;
{
		
	particles.clear();	
                particles.push_back(-2);
                particles.push_back(1);
		particles.push_back(1);
                particles.push_back(-1);
                particles.push_back(11);
                particles.push_back(-12);

        BHinput input(momenta6pt,mu_momenta6pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();

#if GET_BORN
	me2_born=me2->get_born();
#endif

_PRINT(me2_finite);
_PRINT(me2_single_pole);
_PRINT(me2_double_pole);
#if GET_BORN
_PRINT(me2_born);
#endif

  std::cout<< "-----------------------------------------------" <<std::endl;

  
#if write_TEX 
  cout<<fixed<< "$(\\bar u d -> d \\bar b e^- \\bar\\nu)$"  <<
	  "  & $ "<<me2_double_pole << 
	  " $ & $ "<<me2_single_pole<< 
	  " $ & $ "<<me2_finite<< 
	  " $\\\\ \\hline "<<endl;
  cout<<
	  "  & $ "<< 
	  " $ & $ "<< 
	  " $ & $ "<<me2_born<< 
	  " $\\\\ \\hline "<<endl;
#endif

#if write_C
	cout<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][2],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][2],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][2],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[2],me2_born,tol,": born",err) 
#endif
  
  std::cout<< "-----------------------------------------------" <<std::endl;
}
  std::cout<< "**********************************************" <<std::endl;
/*****************************************************************/
#endif


#if CHECK_ME2_2q2g
/*****************************************************************/
/*****************************************************************/
  std::cout<< "**********************************************" <<std::endl;
   	cout << " ME2 2q2g2l " << endl;
{
		
	particles.clear();	
                particles.push_back(-2);
                particles.push_back(21);
                particles.push_back(21);
                particles.push_back(-1);
                particles.push_back(11);
                particles.push_back(-12);

        BHinput input(momenta6pt,mu_momenta6pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();
#if GET_BORN
 	me2_born=me2->get_born();
#endif

_PRINT(me2_finite);
_PRINT(me2_single_pole);
_PRINT(me2_double_pole);
#if GET_BORN
_PRINT(me2_born);
#endif

  std::cout<< "-----------------------------------------------" <<std::endl;

  
#if write_TEX 
  cout<<fixed<< "$(\\bar u g -> g \\bar d e^- \\bar\\nu)$"  <<
	  "  & $ "<<me2_double_pole << 
	  " $ & $ "<<me2_single_pole<< 
	  " $ & $ "<<me2_finite<< 
	  " $\\\\ \\hline "<<endl;
  cout<<
	  "  & $ "<< 
	  " $ & $ "<< 
	  " $ & $ "<<me2_born<< 
	  " $\\\\ \\hline "<<endl;
#endif

#if write_C
	cout<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][3],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][3],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][3],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[3],me2_born,tol,": born",err) 
#endif
  
  std::cout<< "-----------------------------------------------" <<std::endl;
}
  std::cout<< "**********************************************" <<std::endl;
/*****************************************************************/
#endif





#if do_CHECKS
cout<< "*********************************************"<<endl;
cout<< "*                                           *"<<endl;
cout<< "* total number of errors: "<<err<<"                *"  << endl;
cout<< "*                                           *"<<endl;
cout<< "*********************************************"<<endl;
#endif

//BH::CachedOLHA::Cached_OLHA_factory::default_COLHA->print_state();

 
fpu_fix_end(&old_cw);


return err;
}
