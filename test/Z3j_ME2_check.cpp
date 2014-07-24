/*
 *
 * program to check: 
 * 	(1) ME2 of Z+3j process
 * 
 * PS-point: bbdfk7pt 
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

#define HigherPrecision 0 //[0,1] for [double,HP] aritm.

#define write_TEX 0 	//[0,1] for [no,yes] LaTeX output.
#define write_C 0 	//[0,1] for [no,yes] c++ usable output.
#define do_CHECKS 1	//[0,1] for [no,yes] performing checks.


// to choose which partial amplitudes and ME2 should we computed
#define CHECK_ME2_4q_dist_up_down_quark 1
#define CHECK_ME2_4q_dist_up_quark 1
#define CHECK_ME2_4q_ident_up_quark 1

#define CHECK_ME2_up_2q3g 1



#define EPS 0.00000000000001
#define _CHECK(X,Y,tol,mess,err) if ( ((abs(X) == 0.) && (abs(Y) > tol)) || ((abs(X) != 0.) && abs( X - Y)/abs(X) > tol) || std::isnan(abs(X)) ){  _MESSAGE8(mess," NOT OK:",X," != ", Y ," (only " ,-log(abs( (X - Y)/Y))/log(10.)," digits)");err ++; }



using namespace BH;
using namespace std;

string isgn(C a){if(imag(a)>0.){return(" + i ");} else {return(" - i ");};}


int main(void) {
        	
int err(0);
//old accuracy aim: 
//double tol(0.00000015);
  double tol(0.0000002);

cout<<endl;
_PRINT(tol);
cout<<endl;
std::cout<<setprecision(10);

#include "std_momenta.hpp"


settings::use_setting("USE_HIGHER_PRECISION no");
settings::use_setting("INTERFACE_MODE normal");
settings::use_setting("USE_AUTOMATED_ASSEMBLY yes");

  ind7.clear();
  ind7.push_back(1);
  ind7.push_back(2);
  ind7.push_back(3);
  ind7.push_back(4);
  ind7.push_back(5);
  ind7.push_back(6);
  ind7.push_back(7);

	vector<double> p1;
     	vector<double> p2;
        vector<double> p3;
        vector<double> p4;
        vector<double> p5;
        vector<double> p6;
        vector<double> p7;

	vector<int> particles;
	
	double mu_momenta7pt(7.);
	double c(1.);
p1.push_back(-bbdfk7pt.p(1).E().real()*c);p1.push_back(-bbdfk7pt.p(1).X().real()*c);p1.push_back(-bbdfk7pt.p(1).Y().real()*c);p1.push_back(-bbdfk7pt.p(1).Z().real()*c);
p2.push_back(-bbdfk7pt.p(2).E().real()*c);p2.push_back(-bbdfk7pt.p(2).X().real()*c);p2.push_back(-bbdfk7pt.p(2).Y().real()*c);p2.push_back(-bbdfk7pt.p(2).Z().real()*c);
p3.push_back(bbdfk7pt.p(3).E().real()*c);p3.push_back(bbdfk7pt.p(3).X().real()*c);p3.push_back(bbdfk7pt.p(3).Y().real()*c);p3.push_back(bbdfk7pt.p(3).Z().real()*c);
p4.push_back(bbdfk7pt.p(4).E().real()*c);p4.push_back(bbdfk7pt.p(4).X().real()*c);p4.push_back(bbdfk7pt.p(4).Y().real()*c);p4.push_back(bbdfk7pt.p(4).Z().real()*c);
p5.push_back(bbdfk7pt.p(5).E().real()*c);p5.push_back(bbdfk7pt.p(5).X().real()*c);p5.push_back(bbdfk7pt.p(5).Y().real()*c);p5.push_back(bbdfk7pt.p(5).Z().real()*c);
p6.push_back(bbdfk7pt.p(6).E().real()*c);p6.push_back(bbdfk7pt.p(6).X().real()*c);p6.push_back(bbdfk7pt.p(6).Y().real()*c);p6.push_back(bbdfk7pt.p(6).Z().real()*c);
p7.push_back(bbdfk7pt.p(7).E().real()*c);p7.push_back(bbdfk7pt.p(7).X().real()*c);p7.push_back(bbdfk7pt.p(7).Y().real()*c);p7.push_back(bbdfk7pt.p(7).Z().real()*c);


        vector<vector <double> > momenta7pt;
        momenta7pt.push_back(p1);
        momenta7pt.push_back(p2);
        momenta7pt.push_back(p3);
        momenta7pt.push_back(p4);
        momenta7pt.push_back(p5);
        momenta7pt.push_back(p6);
        momenta7pt.push_back(p7);

    clock_t before, after;
	
    C me2_finite;
	C me2_single_pole;
	C me2_double_pole;
	C me2_born;

	me2_finite=C(0.,0.);
	me2_single_pole=C(0.,0.);
	me2_double_pole=C(0.,0.);
	me2_born=C(0.,0.);

	particles.clear();

	BH_Ampl* me2;
	BHinput input(momenta7pt,mu_momenta7pt);

vector< vector<C> > res_me2(3);
vector<C> res_me2_born;

//ME2 distinct up/down-type quarks 4q1g2l 
res_me2[0].push_back(C(5.237471628,0));
res_me2[1].push_back(C(-32.37456065,0));
res_me2[2].push_back(C(-8.333333333,0));
res_me2_born.push_back(C(129.655116,0));
//ME2 distinct up-type quarks 4q1g2l
res_me2[0].push_back(C(2.246735884,0));
res_me2[1].push_back(C(-32.38497697,0));
res_me2[2].push_back(C(-8.333333333,0));
res_me2_born.push_back(C(179.494319,0));
//ME2 ident up-quarks 4q1g2l 
res_me2[0].push_back(C(0.4387412994,0));
res_me2[1].push_back(C(-32.5180903,0));
res_me2[2].push_back(C(-8.333333333,0));
res_me2_born.push_back(C(201.4076989,0));
// ME2 2q3g2l up-quarks
res_me2[0].push_back(C(-15.23260829,0));
res_me2[1].push_back(C(-42.32792665,0));
res_me2[2].push_back(C(-11.66666667,0));
res_me2_born.push_back(C(289.6247791,0));

cout<<endl;
cout<<endl;

#if CHECK_ME2_4q_dist_up_down_quark 

{
	
	BH_interface bhi;
#include "theconstants.h"

	bhi.operator()(input);
	
        particles.push_back(-2);
      	particles.push_back(1);
        particles.push_back(1);
        particles.push_back(-2);
        particles.push_back(21);
        particles.push_back(11);
        particles.push_back(-11);
        
        before=clock();
        me2=bhi.new_ampl(particles);
        after=clock();
        _MESSAGE3("\n\tconstr. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s");
	
        before=clock();
#if HigherPrecision==0
        me2_finite=me2->get_finite();
        me2_single_pole=me2->get_single_pole();
        me2_double_pole=me2->get_double_pole();
#endif


#if HigherPrecision==1
        me2_finite=me2->get_finite_HP();
        me2_single_pole=me2->get_single_pole_HP();
        me2_double_pole=me2->get_double_pole_HP();
#endif

	me2_born=me2->get_born();

        after=clock();
        _MESSAGE3("\teval. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s\n");

  
#if write_TEX 
  cout<<fixed<< "$(\\bar u d -> d \\bar u  g e^- e^+)$"  <<
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



#if CHECK_ME2_4q_dist_up_quark 

{
	
	BH_interface bhi;
#include "theconstants.h"
	bhi.operator()(input);
       
	particles.clear();	
        particles.push_back(-2);
        particles.push_back(4);
        particles.push_back(4);
        particles.push_back(-2);
        particles.push_back(21);
        particles.push_back(11);
        particles.push_back(-11);
        
        before=clock();
        me2=bhi.new_ampl(particles);
        after=clock();
        _MESSAGE3("\n\tconstr. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s");

        before=clock();
#if HigherPrecision==0
        me2_finite=me2->get_finite();
        me2_single_pole=me2->get_single_pole();
        me2_double_pole=me2->get_double_pole();
#endif


#if HigherPrecision==1
        me2_finite=me2->get_finite_HP();
        me2_single_pole=me2->get_single_pole_HP();
        me2_double_pole=me2->get_double_pole_HP();
#endif

    	me2_born=me2->get_born();
        after=clock();
        _MESSAGE3("\teval. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s\n");



  
#if write_TEX 
  cout<<fixed<< "$(\\bar u c -> c \\bar u  g e^- e^+)$"  <<
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







#if CHECK_ME2_4q_ident_up_quark 

{
	
	BH_interface bhi;
#include "theconstants.h"
	bhi.operator()(input);

	particles.clear();	
        particles.push_back(-2);
	particles.push_back(2);
        particles.push_back(2);
        particles.push_back(-2);
        particles.push_back(21);
        particles.push_back(11);
        particles.push_back(-11);

        before=clock();
        me2=bhi.new_ampl(particles);
        after=clock();
        _MESSAGE3("\n\tconstr. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s");
	
        before=clock();
#if HigherPrecision==0
        me2_finite=me2->get_finite();
        me2_single_pole=me2->get_single_pole();
        me2_double_pole=me2->get_double_pole();
#endif


#if HigherPrecision==1
        me2_finite=me2->get_finite_HP();
        me2_single_pole=me2->get_single_pole_HP();
        me2_double_pole=me2->get_double_pole_HP();
#endif

	me2_born=me2->get_born();
        after=clock();
        _MESSAGE3("\teval. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s\n");



  
#if write_TEX 
  cout<<fixed<< "$(\\bar u u -> u \\bar du g e^- e^+)$"  <<
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
  
//	delete me2;  
}
#endif


#if CHECK_ME2_up_2q3g
{
	BH_interface bhi;
#include "theconstants.h"
	bhi.operator()(input);
		
	particles.clear();	
        particles.push_back(-2);
        particles.push_back(21);
        particles.push_back(21);
        particles.push_back(21);
        particles.push_back(-2);
        particles.push_back(11);
        particles.push_back(-11);

        before=clock();
        me2=bhi.new_ampl(particles);
        after=clock();
        _MESSAGE3("\n\tconstr. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s");
	
	
        before=clock();
	
#if HigherPrecision==0
        me2_finite=me2->get_finite();
        me2_single_pole=me2->get_single_pole();
        me2_double_pole=me2->get_double_pole();
#endif


#if HigherPrecision==1
        me2_finite=me2->get_finite_HP();
        me2_single_pole=me2->get_single_pole_HP();
        me2_double_pole=me2->get_double_pole_HP();
#endif

 	me2_born=me2->get_born();


        after=clock();
        _MESSAGE3("\teval. time:\t",double(after-before)/double(CLOCKS_PER_SEC)," s\n");

  
#if write_TEX 
  cout<<fixed<< "$(\\bar u g -> g g g \\bar u e^- e^+)$"  <<
	  "  & $ "<<me2_double_pole << 
	  " $ & $ "<<me2_single_pole<< 
	  " $ & $ "<<me2_finite<< 
                particles.push_back(21);
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
cout<< "*********************************************"<<endl;
cout<< "*                                           *"<<endl;
cout<< "* total number of errors: "<<err<<"                *"  << endl;
cout<< "*                                           *"<<endl;
cout<< "*********************************************"<<endl;
#endif


 


return err;
}
