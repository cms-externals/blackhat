/*
 *
 * program to check: 
 * 	(1) ME2 of photon+2j process
 *  (2) need to use automated assembly
 *
 * PS-point: gkm5
 * 
*/


#include <mom_conf.h>
#include "assembly.h"
#include "settings.h"
#include "settings_reader.h"
#include "polylog.h"  // for pi

#include "cached_OLHA.h"
#include "Interface/BH_interface.h"
#include "Interface/BH_Ampl.h"


#include <ctime>


#define GET_BORN 1 	//[0,1] for [don't use,use] get_born() to evaluate morn ME2.

#define write_TEX 0 	//[0,1] for [no,yes] LaTeX output.
#define write_C 0 	//[0,1] for [no,yes] c++ usable output.
#define do_CHECKS 1	//[0,1] for [no,yes] performing checks.


// to choose which partial amplitudes and ME2 should we computed
#define CHECK_ME2_4q_dist_quark 1
#define CHECK_ME2_4q_ident_quark 1
#define CHECK_ME2_2q2g 1



#define _CHECK(X,Y,tol,mess,err) if ( abs( 1. - Y/X) < tol ){ /*_MESSAGE2(mess," OK");*/	} else { _MESSAGE5(mess," NOT OK:",X," != ", Y );err ++;	}

using namespace BH;
using namespace std;

string isgn(C a){if(imag(a)>0.){return(" + i ");} else {return(" - i ");};}


int main(void) {
        	
	unsigned int old_cw;
    fpu_fix_start(&old_cw);

	settings::use_setting("USE_HIGHER_PRECISION no");
	#include "std_momenta.hpp"

	cout<<setprecision(10);
	double tol=0.00001;

	int err=0;



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
	BH_Ampl* me2;
#include "theconstants.h"

vector< vector<C> > res_me2(3);
vector<C> res_me2_born;

//res_me2[0].push_back(C(4.31801662,0.00000000));
res_me2[0].push_back(C(5.529166141,0.00000000));
res_me2[1].push_back(C(-23.81423441,0.00000000));
res_me2[2].push_back(C(-5.33333333,0.00000000));
res_me2_born.push_back(C(24.72500015,0));

//res_me2[0].push_back(C(2.38082130,0.00000000));
res_me2[0].push_back(C(1.548813934,0.00000000));
res_me2[1].push_back(C(-23.81423441,0.00000000));
res_me2[2].push_back(C(-5.33333333,0.00000000));
res_me2_born.push_back(C(74.52206889,0));

//res_me2[0].push_back(C(-5.88002546,0.00000000));
res_me2[0].push_back(C(-5.630878751,0.00000000));
res_me2[1].push_back(C(-24.05351261,0.00000000));
res_me2[2].push_back(C(-5.33333333,0.00000000));
res_me2_born.push_back(C(26.04682837,0));

//old target
//res_me2[0].push_back(C(-55.05121687,0.00000000));
//per cent level difference by including full helicity sum 
res_me2[0].push_back(C(-55.64223004,0.00000000));
res_me2[1].push_back(C(-42.05300015,0.00000000));
res_me2[2].push_back(C(-8.66666667,0.00000000));
res_me2_born.push_back(C(596.3508444,0));


#if CHECK_ME2_4q_dist_quark 
cout<<endl;
/*****************************************************************/
/*****************************************************************/
{
	
                particles.push_back(-2);
	        	particles.push_back(3);
                particles.push_back(3);
                particles.push_back(-2);
                particles.push_back(22);
                
        
	BHinput input(momenta5pt,mu_momenta5pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();

#if GET_BORN
	me2_born=me2->get_born();
#endif




  
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
 
  
}
#endif


#if CHECK_ME2_4q_dist_quark 
cout<<endl;

{
            	particles.clear();	
                particles.push_back(-2);
        		particles.push_back(4);
                particles.push_back(4);
                particles.push_back(-2);
                particles.push_back(22);
                

        BHinput input(momenta5pt,mu_momenta5pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();

#if GET_BORN
	me2_born=me2->get_born();
#endif


  
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
  
}
#endif


#if CHECK_ME2_4q_ident_quark 
cout<<endl;
{
		
	particles.clear();	
                particles.push_back(-1);
                particles.push_back(1);
		particles.push_back(1);
                particles.push_back(-1);
                particles.push_back(22);
                

        BHinput input(momenta5pt,mu_momenta5pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();

#if GET_BORN
	me2_born=me2->get_born();
#endif



  
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
 
  
}
#endif


#if CHECK_ME2_2q2g
cout<<endl;
{
		
	particles.clear();	
                particles.push_back(-2);
                particles.push_back(21);
                particles.push_back(21);
                particles.push_back(-2);
                particles.push_back(22);
                

        BHinput input(momenta5pt,mu_momenta5pt);
	me2=bhi.new_ampl(particles);
	bhi.operator()(input);
	

 	me2_finite=me2->get_finite();
 	me2_single_pole=me2->get_single_pole();
 	me2_double_pole=me2->get_double_pole();
#if GET_BORN
 	me2_born=me2->get_born();
#endif



  
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
  
}
#endif





#if do_CHECKS
cout<<endl;
cout<< "*********************************************"<<endl;
cout<< "*                                           *"<<endl;
cout<< "* total number of errors: "<<err<<"                *"  << endl;
cout<< "*                                           *"<<endl;
cout<< "*********************************************"<<endl;
#endif

#if 0 
BH::CachedOLHA::Cached_OLHA_factory::default_COLHA->print_state();
cout<< "*********************************************"<<endl;
BH::CachedTHA::Cached_THA_factory::default_CTHA->print_state();
#endif
 
fpu_fix_end(&old_cw);


return err;
}
