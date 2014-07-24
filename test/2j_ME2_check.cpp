/*
 *
 * program to check: 
 * ME2 of 2jets process
 * QED order 0
 * QCD LO
 * values checked against madgraph5
*/

#include <mom_conf.h>
#include "assembly.h"
#include "settings.h"
#include "settings_reader.h"

#include "cached_OLHA.h"
#include "Interface/BH_interface.h"
#include "Interface/BH_Ampl.h"
#include <ctime>
#include "polylog.h"  // for pi

#define EPS 0.00000000000001
#define _CHECK(X,Y,tol,mess,err) if ( abs( 1. - Y/(X+Y*EPS)) < tol ){ _MESSAGE2(mess,"OK");	} else { _MESSAGE5(mess,"NOT OK:",X," != ", Y );err ++;	}
using namespace BH;
using namespace std;

int main(void) {
        	
    unsigned int old_cw;
    fpu_fix_start(&old_cw);

//settings::use_setting("USE_AUTOMATED_ASSEMBLY yes");
//settings::use_setting("REAL_ASSEMBLY_ONLY yes");


string process_names[]={
   "ME2_indentical_quarks_4q",
   "ME2_distinct_quarks_4q",
   "ME2_2q2g",
   "ME2_4g"};

int process[]={
   2,-2, 2,-2,
   2,-2, 1,-1,
   2,-2, 21,21,
   21,21,21,21};

bool check[]={
  true,
  true,
  true,
  true
  };




    int err(0);
    double tol(0.00001);
    int target_number;

    std::cout<<setprecision(14);

	double mu_momenta4pt(4.);
    vector<vector <double> > momenta4pt;

	BH_interface bhi;
	BH_Ampl* me2;

	bhi.set("W_mass",80.4190000);
 	bhi.set("W_width",2.04759951);
	bhi.set("alpha_S",0.118000000);
	bhi.set("alpha_QED",1./132.506980);
	double sin_th_2=1.-pow(80.419/91.188,2);
	bhi.set("sin_th_2",sin_th_2);
	bhi.set("sin_2th",sin(2.*asin(sqrt(sin_th_2))));
	double GF=1.16639*pow(10.,-5);
	double alpha_QED=sqrt(2.)*GF*pow(80.419,2)*sin_th_2/3.141592653589793;
	
  double p_a1[]={0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03};//  0.0000000E+00
  double p_a2[]={0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03};//  0.0000000E+00
  double p_a3[]={0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03};//  0.9725608E-05
  double p_a4[]={0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03};//  0.1009274E-04


    vector<double> p1(p_a1,p_a1+4), p2(p_a2,p_a2+4), p3(p_a3,p_a3+4), p4(p_a4,p_a4+4);
     momenta4pt.push_back(p1);
     momenta4pt.push_back(p2);
     momenta4pt.push_back(p3);
     momenta4pt.push_back(p4);

    C me2_born=C(0,0);
    vector<C> res_me2_born;
    vector<C> me2_loop;
    vector<vector<C> > res_me2_loop(4);


    double normalization; 
    res_me2_born.clear();

    normalization=6.*6.*(1); 
    res_me2_born.push_back(C(normalization*2.82769286,0));
    res_me2_loop[0].push_back(C(-279.27597016122,0));
    res_me2_loop[0].push_back(C(48.958066856582,0));
    res_me2_loop[0].push_back(C(-5.3333333310582,0));

    res_me2_born.push_back(C(normalization*0.566450061,0));
    res_me2_loop[1].push_back(C(-288.86414417921,0));
    res_me2_loop[1].push_back(C(47.62480706301,0));
    res_me2_loop[1].push_back(C(-5.3333333370724,0));

    normalization=6.*6.*(2); 
    res_me2_born.push_back(C(normalization*1.89410151,0));
    res_me2_loop[2].push_back(C(-422.20699951351,0));
    res_me2_loop[2].push_back(C(78.3091660889,0));
    res_me2_loop[2].push_back(C(-8.6666666666666666,0));
    
    normalization=16.*16.*(2); 
    res_me2_born.push_back(C(normalization*55.1792506,0));
    res_me2_loop[3].push_back(C(-556.449693072,0));
    res_me2_loop[3].push_back(C(108.8994149369,0));
    res_me2_loop[3].push_back(C(-12,0));



    clock_t before, after, before_c, after_c;
    int nbr_points(1);

        BHinput input(momenta4pt,mu_momenta4pt);

for (int i=0; i<4; i++){ 
    if(check[i]){
        _MESSAGE("====");
        vector<int> particles(process+4*i,process+4*(i+1));
       

        before_c=clock();
            me2=bhi.new_ampl(particles);
        after_c=clock();
        _MESSAGE("-----");
        _MESSAGE3(process_names[i],":\t",particles);
    	_MESSAGE3("\tconstr. time:\t",double(after_c-before_c)/double(CLOCKS_PER_SEC)," s");

       
        before=clock();
       for(int j=0;j<nbr_points;j++){ 
        //    bhi.operator()(input);
           me2_loop.clear();
        bhi.operator()(input);

            me2_born=me2->get_born();
        	me2_loop.push_back(me2->get_finite());
        	me2_loop.push_back(me2->get_single_pole());
        	me2_loop.push_back(me2->get_double_pole());

        }
        after=clock();
    }

    if(check[i]){
        _MESSAGE3("\teval time:    \t", double(after-before)/double(CLOCKS_PER_SEC)/double(nbr_points)," s");
        _CHECK(res_me2_loop[i][2],me2_loop[2],tol,"\tdouble:\t",err) 
        _CHECK(res_me2_loop[i][1],me2_loop[1],tol,"\tsingle:\t",err) 
        _CHECK(res_me2_loop[i][0],me2_loop[0],tol,"\tfinite:\t",err) 
        _CHECK(res_me2_born[i],me2_born,tol,"\tmatch target:\t",err) 
    }
}
    cout<< "*********************************************"<<endl;
    cout<< "*                                           *"<<endl;
    cout<< "* total number of errors: "<<err<<"                *"  << endl;
    cout<< "*                                           *"<<endl;
    cout<< "*********************************************"<<endl;


//BH::CachedOLHA::Cached_OLHA_factory::default_COLHA->print_state();
//BH::CachedTHA::Cached_THA_factory::default_CTHA->print_state();


    fpu_fix_end(&old_cw);


return err;
}
