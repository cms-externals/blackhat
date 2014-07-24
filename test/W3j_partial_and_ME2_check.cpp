/*
 *
 * program to check: 
 * 	(1) partial amplitudes of the W+3jets process
 * 	(2) ME2 of W+3j process
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
#
#include <ctime>

#define HigherPrecision 0 //[0,1] for [double,HP] aritm.

#define write_TEX 0 	//[0,1] for [no,yes] LaTeX output.
#define write_C 0	//[0,1] for [no,yes] c++ usable output.
#define do_CHECKS 1	//[0,1] for [no,yes] performing checks.


// to choose which partial amplitudes and ME2 should we computed
#define CHECK_ME2_4q1g_dist_quark 1
#define CHECK_ME2_4q1g_ident_quark 1
#define CHECK_ME2_4q1g_ident_antiq 1

#define CHECK_ME2_2q3g 1

#define print_cache 0

#define EPS 0.0000000000001
#define _CHECK(X,Y,tol,mess,err) if ( ((abs(X) == 0.) && (abs(Y) > tol)) || ((abs(X) != 0.) && abs( X - Y)/abs(X) > tol) || std::isnan(abs(X)) ){  _MESSAGE8(mess," NOT OK:",X," != ", Y ," (only " ,-log(abs( (X - Y)/Y))/log(10.)," digits)");err ++;	}

using namespace BH;
using namespace std;

string isgn(C a){if(imag(a)>0.){return(" + i ");} else {return(" - i ");};}


int main(void) {
        	
		unsigned int old_cw;
	    fpu_fix_start(&old_cw);

	settings::use_setting("INTERFACE_MODE normal");

#include "std_momenta.hpp"

  int mu_index = DefineMu<double>(bbdfk7pt,7); //use mu=7 for partial amplitudes
  int mu_index_HP = DefineMu<RHP>(bbdfk7pt_HP,7); //use mu=7 for partial amplitudes
 
  int err(0);
  cout<<setprecision(10);
#if HigherPrecision==1
  cout<<setprecision(10);
#endif
  double tol=0.000001;

    clock_t before, after;

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
	BHinput input(momenta7pt,mu_momenta7pt);
	bhi.operator()(input);


vector< vector<C> > res_me2(3);
vector<C> res_me2_born;

// ME2 distinct quarks 4q1g2l 
res_me2[0].push_back(C(1.7780613303,0.0000000000));
res_me2[1].push_back(C(-32.3767720981,0.0000000000));
res_me2[2].push_back(C(-8.3333333333,0.0000000000));
res_me2_born.push_back(C(0.0002632647,0.0000000000));
// ME2 identical quarks 4q1g2l 
res_me2[0].push_back(C(1.0350002564,0.0000000000));
res_me2[1].push_back(C(-32.4080716525,0.0000000000));
res_me2[2].push_back(C(-8.3333333333,0.0000000000));
res_me2_born.push_back(C(0.0002675986,0.0000000000));
// ME2 identical anti-quarks 4q1g2l 
res_me2[0].push_back(C(0.4788030624,0.0000000000));
res_me2[1].push_back(C(-32.5075013637,0.0000000000));
res_me2[2].push_back(C(-8.3333333333,0.0000000000));
res_me2_born.push_back(C(0.0002967181,0.0000000000));
// ME2 2q3g2l 
res_me2[0].push_back(C(-13.9799122531,0.0000000000));
res_me2[1].push_back(C(-42.3430362824,0.0000000000));
res_me2[2].push_back(C(-11.6666666667,0.0000000000));
res_me2_born.push_back(C(0.0005316697,0.0000000000));

cout<<endl;
cout<<endl;
#if CHECK_ME2_4q1g_dist_quark 

{
	
                particles.push_back(-2);
		particles.push_back(3);
                particles.push_back(3);
                particles.push_back(-1);
                particles.push_back(21);
                particles.push_back(11);
                particles.push_back(-12);
        
        before=clock();
	me2=bhi.new_ampl(particles);
    after=clock();

    _MESSAGE3("\tconstr. time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");

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


    _MESSAGE3("\teval time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");
    cout<<endl;
  
#if write_TEX 
  cout<<fixed<< "$(\\bar u c -> c \\bar d g e^- \\abr\\nu)$"  <<
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
	cout<<fixed<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<fixed<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<fixed<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<fixed<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][0],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][0],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][0],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[0],me2_born,tol,": born",err) 
#endif
  
}
#endif


#if CHECK_ME2_4q1g_ident_quark 

{
	particles.clear();	
                particles.push_back(-2);
		particles.push_back(2);
                particles.push_back(2);
                particles.push_back(-1);
                particles.push_back(21);
                particles.push_back(11);
                particles.push_back(-12);

        before=clock();
	me2=bhi.new_ampl(particles);
	    after=clock();

    _MESSAGE3("\tconstr. time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");

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


    _MESSAGE3("\teval time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");
    cout<<endl;
 
#if write_TEX 
  cout<<fixed<< "$(\\bar u u -> u \\bar d g e^- \\abr\\nu)$"  <<
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
	cout<<fixed<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<fixed<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<fixed<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<fixed<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][1],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][1],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][1],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[1],me2_born,tol,": born",err) 
#endif
  
}
#endif


#if CHECK_ME2_4q1g_ident_antiq 
{
		
	particles.clear();	
                particles.push_back(-2);
                particles.push_back(1);
		particles.push_back(1);
                particles.push_back(-1);
                particles.push_back(21);
                particles.push_back(11);
                particles.push_back(-12);

        before=clock();
	me2=bhi.new_ampl(particles);
		    after=clock();

    _MESSAGE3("\tconstr. time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");

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


    _MESSAGE3("\teval time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");
    cout<<endl;
 



  
#if write_TEX 
  cout<<fixed<< "$(\\bar u d -> d \\bar b g e^- \\bar\\nu)$"  <<
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
	cout<<fixed<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<fixed<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<fixed<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<fixed<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][2],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][2],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][2],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[2],me2_born,tol,": born",err) 
#endif
  
}
#endif


#if CHECK_ME2_2q3g
{
		
	particles.clear();	
                particles.push_back(-2);
                particles.push_back(21);
                particles.push_back(21);
                particles.push_back(21);
                particles.push_back(-1);
                particles.push_back(11);
                particles.push_back(-12);

        before=clock();
	me2=bhi.new_ampl(particles);
		    after=clock();

    _MESSAGE3("\tconstr. time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");

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


    _MESSAGE3("\teval time:    \t", double(after-before)/double(CLOCKS_PER_SEC)," s");
    cout<<endl;
 



  
#if write_TEX 
  cout<<fixed<< "$(\\bar u g -> g g \\bar d e^- \\bar\\nu)$"  <<
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
	cout<<fixed<<"res_me2[0].push_back(C"<<me2_finite<<");"<<endl;
	cout<<fixed<<"res_me2[1].push_back(C"<<me2_single_pole<<");"<<endl;
	cout<<fixed<<"res_me2[2].push_back(C"<<me2_double_pole<<");"<<endl;
	cout<<fixed<<"res_me2_born.push_back(C"<<me2_born<<");"<<endl;
#endif


#if do_CHECKS
  _CHECK(res_me2[2][3],me2_double_pole,tol,":-2",err) 
  _CHECK(res_me2[1][3],me2_single_pole,tol,":-1",err)
  _CHECK(res_me2[0][3],me2_finite,tol,": 0",err) 
  _CHECK(res_me2_born[3],me2_born,tol,": born",err) 
#endif
  
}
#endif


#if print_cache
BH::CachedOLHA::Cached_OLHA_factory::default_COLHA->print_state();
#endif




#if do_CHECKS
cout<< "*********************************************"<<endl;
cout<< "*                                           *"<<endl;
cout<< "* total number of errors: "<<err<<"                *"  << endl;
cout<< "*                                           *"<<endl;
cout<< "*********************************************"<<endl;
#endif

BH::CachedOLHA::Cached_OLHA_factory::default_COLHA->print_state();

 
fpu_fix_end(&old_cw);


return err ;
}
