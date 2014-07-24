/*
 *
 * program to check: 
 * ME2 of 3jets process
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
#define _CHECK(X,Y,tol,mess,err) if ( abs( 1. - Y/(X+Y*EPS)) < tol ){ _MESSAGE2(mess,"OK");	} else { _MESSAGE7(mess,"NOT OK:",X," != ", Y ," rel err " ,1.-Y/(X+Y*EPS));err ++;	}
using namespace BH;
using namespace std;

int main(void) {
        	
    unsigned int old_cw;
    fpu_fix_start(&old_cw);

//settings::use_setting("USE_AUTOMATED_ASSEMBLY yes");
//settings::use_setting("REAL_ASSEMBLY_ONLY yes");


string process_names[]={
   "virt ME2",
   "virt ME2",
   "virt ME2",
   "virt ME2",
   "virt ME2",
   "virt ME2"};

int process[]={
    1, 1, 21, 1, 1,
    1,-1, 21, 1,-1,
    1, -1,21, 2, -2,
    1,-1,21,21,21,
   21,21,21,1,-1,
    21,21,21,21,21};

bool check[]={
    true,
    true,
    true,
    true,
    true,
    true};
/*
bool check[]={
  false,
  false,
  false,
  false,
  true,
  false};
*/


    int err(0);
    double tol(0.00001);
    int target_number;

    std::cout<<setprecision(10);

	double mu_momenta5pt(90.);
	//double mu_momenta5pt(1.);
    vector<vector <double> > momenta5pt;

	BH_interface bhi;
	BH_Ampl* me2;
#if 0
	bhi.set("W_mass",80.4190000);
 	bhi.set("W_width",2.04759951);
	bhi.set("alpha_S",0.118000000);
	bhi.set("alpha_QED",1./132.506980);
	double sin_th_2=1.-pow(80.419/91.188,2);
	bhi.set("sin_th_2",sin_th_2);
	bhi.set("sin_2th",sin(2.*asin(sqrt(sin_th_2))));
	double GF=1.16639*pow(10.,-5);
	double alpha_QED=sqrt(2.)*GF*pow(80.419,2)*sin_th_2/3.141592653589793;
#endif
#include "constants.h"

/* sherpa 10 digits right for p^2*/
   double p_a1[]={849.734687032927,0,0,849.734687032927};
   double p_a2[]={849.734687032927,0,0,-849.734687032927};
   double p_a3[]={669.684727662989,26.8464568574744,657.500757397685,-124.296646136262};
   double p_a4[]={253.791480540209,-96.6527838445971,27.7879086518811,-233.015422456616};
   double p_a5[]={775.993165862656,69.8063269871227,-685.288666049566,357.312068592878};
/**/
    /* mc5 
double p_a1[]={-1.86058205404180962,-0.22322341633486098,-1.29381396176358134,-1.31832557381242467};
double p_a2[]={-2.22702886899219324,1.97843090981622047,-0.18090821896703733,-1.00635030417771722};
double p_a3[]={1.77552105445334353,-1.61888957780145734,0.46152052119123634,-0.56450895317284531};
double p_a4[]={-4.59937300385497873,4.17419448896838902,-0.96173464232271022,-1.67493249852413707};
double p_a5[]={-1.26375897363236766,-0.80009741768557219,-0.97450805959914479,-0.08523442629315951};
*//*
#include "std_momenta.hpp"
double p_a1[]={gkm5.p(1).E().real(),gkm5.p(1).X().real(),gkm5.p(1).Z().real(),gkm5.p(1).Y().real()};
double p_a2[]={gkm5.p(2).E().real(),gkm5.p(2).X().real(),gkm5.p(2).Z().real(),gkm5.p(2).Y().real()};
double p_a3[]={gkm5.p(3).E().real(),gkm5.p(3).X().real(),gkm5.p(3).Z().real(),gkm5.p(3).Y().real()};
double p_a4[]={gkm5.p(4).E().real(),gkm5.p(4).X().real(),gkm5.p(4).Z().real(),gkm5.p(4).Y().real()};
double p_a5[]={gkm5.p(5).E().real(),gkm5.p(5).X().real(),gkm5.p(5).Z().real(),gkm5.p(5).Y().real()};
*/

    vector<double> p1(p_a1,p_a1+4), p2(p_a2,p_a2+4), p3(p_a3,p_a3+4), p4(p_a4,p_a4+4), p5(p_a5,p_a5+4);
     momenta5pt.push_back(p1);
     momenta5pt.push_back(p2);
     momenta5pt.push_back(p3);
     momenta5pt.push_back(p4);
     momenta5pt.push_back(p5);

    C me2_born=C(0,0);
    vector<C> me2_loop;
    vector<C> res_me2_born;
    vector<vector<C> > res_me2_loop(6);

    //cut part only
    double normalization; 
    res_me2_born.clear();
    res_me2_born.push_back(C(0.2753607004,0));
    //cut part res_me2_loop[0].push_back(C(-61.969248,0));
    res_me2_loop[0].push_back(C(-58.46640046,0));
    res_me2_loop[0].push_back(C(29.02548571629 ,0));
    res_me2_loop[0].push_back(C(-8.3333333333333,0));

    res_me2_born.push_back(C(0.003242323504,0));
    //cut part res_me2_loop[1].push_back(C(-64.344137,0));
    res_me2_loop[1].push_back(C(-62.35432842,0));
    res_me2_loop[1].push_back(C(26.640840581622,0));
    res_me2_loop[1].push_back(C(-8.3333333333333,0));

    res_me2_born.push_back(C(0.001042280806,0));
    //cut part res_me2_loop[2].push_back(C(-87.839009,0));
    res_me2_loop[2].push_back(C(-83.87072177,0));
    res_me2_loop[2].push_back(C(26.631664665353,0));
    res_me2_loop[2].push_back(C(-8.3333333333333,0));

    res_me2_born.push_back(C(0.1226909135,0));
    //cut part res_me2_loop[3].push_back(C(-77.888732,0));
    res_me2_loop[3].push_back(C(-78.32950033,0));
    res_me2_loop[3].push_back(C(29.961263232251,0));
    res_me2_loop[3].push_back(C(-11.66666666666,0));

    res_me2_born.push_back(C(0.07665952372,0));
    //cut part res_me2_loop[4].push_back(C(-93.47800989,0));
    res_me2_loop[4].push_back(C(-76.83884167,0));
    //res_me2_loop[4].push_back(C(32.212763914293,0));
    res_me2_loop[4].push_back(C(32.10740152483,0));
    res_me2_loop[4].push_back(C(-11.66666666666,0));

    res_me2_born.push_back(C(55.64624626,0));
    //cut part res_me2_loop[4].push_back(C(-93.47800989,0));
    res_me2_loop[5].push_back(C(-93.257644823867));
    res_me2_loop[5].push_back(C(42.712311640779,0));
    res_me2_loop[5].push_back(C(-15,0));



    clock_t before, after, before_c, after_c;
    int nbr_points(1);

cout<<setprecision(15);
cout<<"alpha_S: "<<constants::alpha_S<<endl;

        BHinput input(momenta5pt,mu_momenta5pt);

for (int i=0; i<6; i++){ 
    if(check[i]){
        _MESSAGE("====");
        vector<int> particles(process+5*i,process+5*(i+1));
       

        before_c=clock();
            me2=bhi.new_ampl(particles);
        after_c=clock();
        _MESSAGE("-----");
        _MESSAGE3(process_names[i],":\t",particles);
    	_MESSAGE3("\tconstr. time:\t",double(after_c-before_c)/double(CLOCKS_PER_SEC)," s");

    {
        // quick fix for avoiding HP-computation
/*        for(int j=0;j<1;j++){ 
            bhi.operator()(input);
           me2_loop.clear();

        	me2_loop.push_back(me2->get_finite());
        	me2_loop.push_back(me2->get_single_pole());
        	me2_loop.push_back(me2->get_double_pole());

        	me2_born=me2->get_born();
        }
*/
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
    }
    if(check[i]){
        _MESSAGE3("\teval time:    \t", double(after-before)/double(CLOCKS_PER_SEC)/double(nbr_points)," s");
        _CHECK(res_me2_loop[i][2],me2_loop[2],tol,"\tdouble:\t",err) 
        _CHECK(res_me2_loop[i][1],me2_loop[1],tol,"\tsingle:\t",err) 
        _CHECK(res_me2_loop[i][0],me2_loop[0],tol,"\tfinite:\t",err) 
        _CHECK(res_me2_born[i],me2_born,tol,"\tborn:\t",err) 
    }
}
    cout<< "*********************************************"<<endl;
    cout<< "*                                           *"<<endl;
    cout<< "* total number of errors: "<<err<<"                *"  << endl;
    cout<< "*                                           *"<<endl;
    cout<< "*********************************************"<<endl;

    //BH::CachedOLHA::Cached_OLHA_factory::default_COLHA->print_state();
 
    fpu_fix_end(&old_cw);


return err;
}
