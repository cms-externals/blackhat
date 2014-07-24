/*
 *
 * program to check: 
 * 	(1) partial amplitudes of the W+3jets process
 * 	(2) ME2 of W+3j process
 * 
 * PS-point: gkm4
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
#define CHECK_ME2_2q1g 1



#define _CHECK(X,Y,tol,mess,err) if ( abs( X - Y) < tol ){ _MESSAGE2(mess," OK");	} else { _MESSAGE5(mess," NOT OK:",X," != ", Y );err ++;	}

using namespace BH;
using namespace std;

string isgn(C a){if(imag(a)>0.){return(" + i ");} else {return(" - i ");};}


int main(void) {
        	
		unsigned int old_cw;
	        fpu_fix_start(&old_cw);

   	settings::use_setting("COLOR_MODE full_color");

#include "std_momenta.hpp"

  cout<<setprecision(8);
  double tol=0.00001;

  int err=0;



  ind4.clear();
  ind4.push_back(1);
  ind4.push_back(2);
  ind4.push_back(3);
  ind4.push_back(4);

	vector<double> p1;
     	vector<double> p2;
        vector<double> p3;
        vector<double> p4;

	vector<int> particles;
	
	double mu_momenta4pt(4.);
	double c(1.);
p1.push_back(-gkm4.p(1).E().real()*c);p1.push_back(-gkm4.p(1).X().real()*c);p1.push_back(-gkm4.p(1).Y().real()*c);p1.push_back(-gkm4.p(1).Z().real()*c);
p2.push_back(-gkm4.p(2).E().real()*c);p2.push_back(-gkm4.p(2).X().real()*c);p2.push_back(-gkm4.p(2).Y().real()*c);p2.push_back(-gkm4.p(2).Z().real()*c);
p3.push_back(gkm4.p(3).E().real()*c);p3.push_back(gkm4.p(3).X().real()*c);p3.push_back(gkm4.p(3).Y().real()*c);p3.push_back(gkm4.p(3).Z().real()*c);
p4.push_back(gkm4.p(4).E().real()*c);p4.push_back(gkm4.p(4).X().real()*c);p4.push_back(gkm4.p(4).Y().real()*c);p4.push_back(gkm4.p(4).Z().real()*c);


        vector<vector <double> > momenta4pt;
        momenta4pt.push_back(p1);
        momenta4pt.push_back(p2);
        momenta4pt.push_back(p3);
        momenta4pt.push_back(p4);

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

vector<C> res_me2;
C res_me2_born;


res_me2.push_back(C(-17.52675772,0.00000000));
res_me2.push_back(C(-19.75199044,0.00000000));
res_me2.push_back(C(-5.66666667,0.00000000));
res_me2_born=C(8.7445108,0);


#if CHECK_ME2_2q1g
/*****************************************************************/
/*****************************************************************/
  std::cout<< "**********************************************" <<std::endl;
   	cout << " ME2 2q1g1y " << endl;
{
		
	particles.clear();	
                particles.push_back(-2);
                particles.push_back(21);
                particles.push_back(-2);
                particles.push_back(22);
                

        BHinput input(momenta4pt,mu_momenta4pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();
#if GET_BORN
 	me2_born=me2->get_born();
#endif

_PRINT(particles);
_PRINT(me2_finite);
_PRINT(me2_single_pole);
_PRINT(me2_double_pole);
#if GET_BORN
_PRINT(me2_born);
#endif

  std::cout<< "-----------------------------------------------" <<std::endl;

  
#if write_TEX 
  cout<<fixed<< "$(\\bar u g -> \\bar d e^- \\bar\\nu)$"  <<
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
  _CHECK(res_me2[2],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born,me2_born,tol,": born",err) 
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


 
fpu_fix_end(&old_cw);


return err;
}
