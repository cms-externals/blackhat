/*
 * partial_amplitudes.cpp
 *
 *  Created on: Aug 21, 2008
 *      Author: daniel
 */

//#define _CHECK(X,Y,tol,mess,err) if ( abs( X - Y)/abs(X) < tol ){ _MESSAGE4 (mess," OK (",-log(abs( (X - Y)/Y))/log(10.)," digits)");} else { _MESSAGE5(mess," NOT OK:",X," != ", Y );err ++;	}
#define _CHECK(X,Y,tol,mess,err) if ( abs( X - Y)/abs(X) > tol ){ _MESSAGE5(mess," NOT OK:",X," != ", Y );err ++;	}

#include <mom_conf.h>
#include "assembly.h"
#include "settings.h"
#include "settings_reader.h"
#include "path.h"
#include "cached_OLHA.h"
#include "Interface/BH_interface.h"
#include "Interface/BH_Ampl.h"
#include <ctime>
#include "polylog.h"  // for pi

#define WRITE_TARGETS 1    //[0,1] ...  [write to, compare with] targets
// has to be run in study dir.

using namespace std;
using namespace BH;


int main(){

	int err=0;
	double tol(0.0001);
	int pspoint=1;

	BH_interface bhi;
	std::vector<BH_Ampl*> me2;

	bhi.set("W_mass",80.4190000);
 	bhi.set("W_width",2.04759951);
	bhi.set("Z_mass",91.188);
 	bhi.set("Z_width",2.49);
	bhi.set("alpha_S",0.118000000);
	bhi.set("alpha_QED",1./132.506980);
	double sin_th_2=1.-pow(80.419/91.188,2);
	bhi.set("sin_th_2",sin_th_2);
	//bhi.set("sin_2th",sin(2.*asin(sqrt(sin_th_2))));
	double GF=1.16639*pow(10.,-5);
	double alpha_QED=sqrt(2.)*GF*pow(80.419,2)*sin_th_2/3.141592653589793;

	// building the virtual process
	std::vector<std::vector<int> > processes;
#include THEPATHTOCOLLECTPS
	cout <<"======================="<<endl;
	cout <<"ME construction starts"<<endl;
	for(size_t i=0; i<processes.size(); i++){
		me2.push_back(bhi.new_ampl(processes[i]));
	}
	cout<<"nbr of subprocesses: "<<me2.size()<<endl;
	cout <<"======================="<<endl;


	// running on PS points
        cout<<setprecision(16);
        string MOMENTAfile;
	MOMENTAfile=GetSrcPath()+string("/test/PStest/")+ string(THEPROCESSTOTEST)+string("/PSpoints.dat");
	size_t multiplicity(processes[0].size());

        string TARGETfile;
	TARGETfile=GetSrcPath()+string("/test/PStest/")+ string(THEPROCESSTOTEST)+string("/target.dat");


	ifstream momenta_file;
	momenta_file.open(MOMENTAfile.c_str());

#if WRITE_TARGETS==1
	ifstream targets_file;
	targets_file.open(TARGETfile.c_str());
#else
	ofstream targets_file;
	targets_file.open(TARGETfile.c_str());
	targets_file<<setprecision(16);
#endif

    	while ( momenta_file.good() && pspoint<10001 ) {

		string line;
		getline(momenta_file,line);
		if(line.size() < 10 ) break;
		stringstream momenta_str(line);

		int pro, ps;
		momenta_str>>pro;


		vector<vector <double> > momenta;
		for(size_t i=0;i<multiplicity;i++){
			vector<double> p1;
			for(size_t j=0;j<4;j++){
				double val;
				momenta_str>>val;
				p1.push_back(val);
			}
			momenta.push_back(p1);
		}

		double vals[3];

		double mu_momenta(100.);
		BHinput input(momenta,mu_momenta);
		bhi.operator()(input);

		vals[0]=me2[pro-1]->get_finite();
		vals[1]=me2[pro-1]->get_single_pole();
		vals[2]=me2[pro-1]->get_double_pole();

		int err_before(err);
#if WRITE_TARGETS==1
		//compare to targets
		targets_file>>ps;
		double targets[3];
		targets_file>>targets[0];
		targets_file>>targets[1];
		targets_file>>targets[2];

		_CHECK(targets[0],vals[0],tol,"finite:  ",err);
		_CHECK(targets[1],vals[1],tol,"1/eps:   ",err);
		_CHECK(targets[2],vals[2],tol,"1/eps^2: ",err);

		if(err!=err_before){
			err=err_before+1;
			cout<<"--------------------------------------"<<endl;
			cout <<"#: "<< pspoint;
			cout <<": particles: "<< processes[pro-1] <<endl;
			cout<<"--------------------------------------"<<endl;
		}
#else
		//write targets
		targets_file<<pspoint<<"\t";
		targets_file<< vals[0]<<"\t";
		targets_file<< vals[1]<<"\t";
		targets_file<< vals[2]<<endl;
#endif
	pspoint++;
	};

   	momenta_file.close();
	targets_file.close();

	cout<<"***************************"<<endl;
	cout<<"* total errors: "<<err<<"\t  *"<<endl;
	cout<<"* PS points:    "<<(pspoint-1)<<"\t  *"<<endl;
	cout<<"* tolerance:    "<<tol<<"\t  *"<<endl;
	cout<<"***************************"<<endl;

	return err;
}
