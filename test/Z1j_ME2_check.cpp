/*
 *
 * program to check: 
 * 	(1) ME2 of Z+1j process
 * 
 * PS-point: gkm5
 * 
*/

#include <mom_conf.h>
#include "assembly.h"
#include "settings.h"
#include "settings_reader.h"
#include "polylog.h"  // for pi


#include "Interface/BH_interface.h"
#include "Interface/BH_Ampl.h"

#include <ctime>


#define write_TEX 0 	//[0,1] for [no,yes] LaTeX output.
#define write_C 0	//[0,1] for [no,yes] c++ usable output.
#define do_CHECKS 1	//[0,1] for [no,yes] performing checks.


// to choose which partial amplitudes and ME2 should we computed
#define CHECK_ME2_2q1g 1

#define _CHECK(X,Y,tol,mess,err) if ( ((abs(X) == 0.) && (abs(Y) > tol)) || ((abs(X) != 0.) && abs( X - Y)/abs(X) > tol) || std::isnan(abs(X)) ){  _MESSAGE8(mess," NOT OK:",X," != ", Y ," (only " ,-log(abs( (X - Y)/Y))/log(10.)," digits)");err ++;	}



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
#if !oldBH
settings::use_setting("COLOR_MODE full_color");
#endif

  ind5.clear();
  ind5.push_back(1);
  ind5.push_back(2);
  ind5.push_back(3);
  ind5.push_back(4);
  ind5.push_back(5);

	vector<double> p1;
     	vector<double> p2;
        vector<double> p3;
        vector<double> p4;
        vector<double> p5;

	vector<int> particles;
	
	double mu_momenta5pt(5.);
	double c(1.);
p1.push_back(-gkm5.p(1).E().real()*c);p1.push_back(-gkm5.p(1).X().real()*c);p1.push_back(-gkm5.p(1).Y().real()*c);p1.push_back(-gkm5.p(1).Z().real()*c);
p2.push_back(-gkm5.p(2).E().real()*c);p2.push_back(-gkm5.p(2).X().real()*c);p2.push_back(-gkm5.p(2).Y().real()*c);p2.push_back(-gkm5.p(2).Z().real()*c);
p3.push_back(gkm5.p(3).E().real()*c);p3.push_back(gkm5.p(3).X().real()*c);p3.push_back(gkm5.p(3).Y().real()*c);p3.push_back(gkm5.p(3).Z().real()*c);
p4.push_back(gkm5.p(4).E().real()*c);p4.push_back(gkm5.p(4).X().real()*c);p4.push_back(gkm5.p(4).Y().real()*c);p4.push_back(gkm5.p(4).Z().real()*c);
p5.push_back(gkm5.p(5).E().real()*c);p5.push_back(gkm5.p(5).X().real()*c);p5.push_back(gkm5.p(5).Y().real()*c);p5.push_back(gkm5.p(5).Z().real()*c);


        vector<vector <double> > momenta5pt;
        momenta5pt.push_back(p1);
        momenta5pt.push_back(p2);
        momenta5pt.push_back(p3);
        momenta5pt.push_back(p4);
        momenta5pt.push_back(p5);

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
#if oldBH
	BH2sherpa* me2;
#endif
#if !oldBH
	BH_Ampl* me2;
#endif
#include "theconstants.h"

vector<C> res_me2;
C res_me2_born;

res_me2.push_back(C(-22.4260007896,0.00000));
res_me2.push_back(C(-22.010516,0.000000));
res_me2.push_back(C(-5.666667,0.000000));
res_me2_born=C(0.4005217108,0.000000);

#if CHECK_ME2_2q1g
  std::cout<< "**********************************************" <<std::endl;
   	cout << " ME2 2q1g2l " << endl;
{
		
	particles.clear();	
                particles.push_back(-2);
                particles.push_back(21);
                particles.push_back(-2);
                particles.push_back(11);
                particles.push_back(-11);

        BHinput input(momenta5pt,mu_momenta5pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();
 	me2_born=me2->get_born();


  
#if write_TEX 
  cout<<fixed<< "$(\\bar u g -> g \\bar u e^- e^+)$"  <<
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
	cout<<"res_me2_born=C"<<me2_born<<";"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born,me2_born,tol,": born",err) 
#endif
  
}
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
