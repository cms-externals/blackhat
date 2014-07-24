/*
 * LH_interface.cpp
 *
 *  Created on: Jun 16, 2009
 *      Author: daniel
 */



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <string>
#include <iterator>
#include <iomanip>
#include <algorithm>


#include "Interface/BH_interface.h"
#include "Interface/BH_interface_impl.h"
#include "Interface/BH_Ampl.h"
#include "mom_conf.h"
#include "settings.h"

#include "Interface/OptionsHandler.h"
#include "settings_reader.h"

using namespace std;
using namespace BH::Tools;
namespace BH {

namespace LesHouches {

template <class T> struct do_delete : public std::unary_function<T*,void> {
	void operator()(T* ptr) { delete ptr;}
};

class syntax_error : std::exception {
	std::string d_line;
public:
	syntax_error(const std::string& line)  {
		d_line= std::string("Syntax error in line: \'") + line + string("\'");
	};
	virtual const char* what() const throw() {
		return d_line.c_str();
	};
	virtual ~syntax_error() throw() {};

};


#define SINGLE_VALUE_OPTION(NAME,VALUE)		case NAME ## _set : {\
			string atype;                                              \
			settingline >> atype;                                      \
			if  ( atype == #VALUE ){                                   \
				if ( generate_contract_file ){                         \
					os << " \t | OK " << endl;                         \
				} else {                                               \
				}                                                      \
			}                                                          \
			else {                                                     \
				os << " \t | ERROR: unsupported " << #NAME << " value \'" << atype << "\'. The only supported type is " << #VALUE << endl; \
				cerr << " ERROR: unsupported " << #NAME << " value \'" << atype << "\'. The only supported type is " << #VALUE << endl; \
				return false;    \
			}                    \
		} break;

class BHOptionsHandler;
class LHOptionsHandler;

class LH_interface {
public:
	bool read_LH_order(const char* filename,const char* contract_filename,bool generate_contract_file);
	static LH_interface* Instance(){return &UniqueInstance;};
	int nbr_particles(int i){ return d_nbr_particles[i-1];}
	void eval(int Label,BHinput& bhi,double (&result)[4]);
	void printHelp(std::ostream& os) const;
	BHOptionsHandler* getBHOptionsHandler() { return d_BHOptionsHandler_p; }
	LHOptionsHandler* getLHOptionsHandler() { return d_LHOptionsHandler_p; }
	template <class T> bool apply_setting(const std::string& name, T& param);
	virtual ~LH_interface();
private:
	vector<BH_Ampl*> d_amplitudes;
	BH_interface* d_bhi;
	static LH_interface UniqueInstance;
	LH_interface();
	std::vector<int> d_nbr_particles;
	settings_table d_options;
	BHOptionsHandler* d_BHOptionsHandler_p;
	LHOptionsHandler* d_LHOptionsHandler_p;

};

//
//


class  LHOptionsHandler : public OptionsHandler {
public:
	int d_alphasPower;
	int d_alphaPower;
	LHOptionsHandler(setable* ST);
	virtual ~LHOptionsHandler(){}
};


LHOptionsHandler::LHOptionsHandler(setable* ST){
	add(new singleValueOption("MassiveParticleScheme","OnShell","Specifies the Scheme for the massive particle."));
	add(new singleValueOption("CorrectionType","QCD","Specifies the type of correction."));
	add(new singleValueOption("MatrixElementSquareType","CHsummed","Specifies the type of matrix element squared."));
	add(new singleValueOption("IRsubtractionMethod","None","Specifies whether the result has an IR subtraction, an which one."));
	add(new singleValueOption("SubdivideProcess","No","Specifies whether to subdivide the result in several contributions when possible."));
	add(new singleValueOption("OperationMode","NormalizedByBorn","Specifies the operation mode."));
	add(new singleValueOption("ResonanceTreatment","ComplexMassScheme","How to treat resonances."));
	add(new singleValueOption("EWRenormalisationScheme","alphaMZ","Sets the EW scheme. This option will be ignored."));
	add(new AlwaysErrorOption("alpha_S","Couplings are all set to 1. This option will be ignored.","Sets the value of the strong coupling constant."));
	add(new AlwaysErrorOption("alpha_QED","Couplings are all set to 1. This option will be ignored.","Sets the value of the QED coupling constant."));
	add(new AlwaysErrorOption("ModelFile","Reading of LH model files is not yet implemented, sorry.","Sets the filename of the model file to use."));
//	add(new AlwaysErrorOption("mode","Only one mode is supported, sorry.","Sets the operating mode."));
	add(new multipleValueOptionWithTableUpdate("IRregularisation","tHV","CDR","Sets the IR regularisation scheme.",ST));
	add(new ValueSettingOption<int>("AlphasPower",d_alphasPower,"Sets alphas power."));
	add(new ValueSettingOption<int>("AlphaPower",d_alphaPower,"Sets alpha power."));
	//enableDebug();

	d_alphasPower=-1;
	d_alphaPower=-1;
}






class  BHOptionsHandler : public OptionsHandler {
public:
	BHOptionsHandler(BH_interface* ST);
	virtual ~BHOptionsHandler(){}
};

BHOptionsHandler::BHOptionsHandler(BH_interface* bhi){
	setable* ST=bhi;
	add(new SettingsTableOption<double>("Z_mass",bhi,"Sets the mass of the Z boson (value should be given in GeV)"));
	add(new SettingsTableOption<double>("W_mass",bhi,"Sets the mass of the W boson (value should be given in GeV)"));
	add(new SettingsTableOption<double>("H_mass",bhi,"Sets the mass of the Higgs boson (value should be given in GeV)"));
	add(new SettingsTableOption<double>("top_mass",bhi,"Sets the mass of the top quark (value should be given in GeV)"));
	add(new SettingsTableOption<double>("bottom_mass",bhi,"Sets the mass of the bottom quark (value should be given in GeV)"));
	add(new SettingsTableOption<double>("Z_width",bhi,"Sets the width of the Z boson (value should be given in GeV)"));
	add(new SettingsTableOption<double>("W_width",bhi,"Sets the width of the W boson (value should be given in GeV)"));
	add(new SettingsTableOption<double>("H_width",bhi,"Sets the width of the Higgs boson (value should be given in GeV)"));
	add(new SettingsTableOption<double>("top_width",bhi,"Sets the width of the top quark (value should be given in GeV)"));
	add(new SettingsTableOption<double>("bottom_width",bhi,"Sets the width of the bottom quark (value should be given in GeV)"));
	add(new SettingsTableOption<double>("sin_th_2",bhi,"Sets the square of sin(\\Theta_W), NOTE: sin_2th should be set to its new value too."));
	add(new SettingsTableOption<double>("sin_2th",bhi,"Sets sin(2\\Theta_W) NOTE: sin_th_2 should be set to its new value too."));
	//enableDebug();
}

template <class T> bool LH_interface::apply_setting(const std::string& name, T& param){
	if (d_options.apply_setting(name,param)){
		return true;
	} else {
		if (d_bhi->apply_setting(name,param)){
				return true;
			} else {
				return false;
			}
	}

};


bool get_setting(char* buffer,ostream& os,bool generate_contract_file,BH_interface* bhi,settings_table* ST=0){

	string message;
	stringstream settingline(buffer);
	stringstream settingline2(buffer);
	stringstream settingline3(buffer);
	bool success;
	BHOptionsHandler* BHOH=LH_interface::Instance()->getBHOptionsHandler();
	LHOptionsHandler* LHOH=LH_interface::Instance()->getLHOptionsHandler();
	success=BHOH->process(settingline,message);
	if ( success ) {
		os << " \t | OK " << endl;
		return true;
	} else {
		if ( BHOH->d_state == OptionsHandler::failed ){
			cerr << message << endl;
			os << " \t | ERROR: " << message <<  endl;
			return false;
		} else {
			success=LHOH->process(settingline2,message);
			if ( success ) {
				os << " \t | OK " << endl;
				return true;
			} else {
					if (BH::settings::read_from_stream(settingline3)) {
										os << " \t | OK " << endl;
										return true;
									} else {
											cerr << message << endl;
											os << " \t | ERROR: " << message <<  endl;
											return false;
					}
					cerr << message << endl;
					os << " \t | ERROR: " << message <<  endl;
					return false;
			}
		}
	}
	// if we get there, the everything should be fine
	return true;
}

void LH_interface::printHelp(ostream& os) const {
	os << "\nBinoth Les Houches interface options: \n\n";
	d_LHOptionsHandler_p->printHelp(os);
	os << "\nBlackHat options: \n\n";
	d_BHOptionsHandler_p->printHelp(os);
}

std::vector<int> read_subprocess_old(char* buffer,int& nbr_tot){
	std::istringstream command(buffer);
	int in_number;
	int out_number;
	string separator;
	command >> in_number;
	command >> separator;
	command >> out_number;

	nbr_tot=in_number+out_number;
	vector<int> particles;

	int nbr_particles=in_number+out_number;
	int newpart;
	for (int i=1; i<= nbr_particles ; i++){
		command >> newpart;
		particles.push_back(newpart);
	}
	return particles;

}
std::vector<int> read_subprocess(char* buffer,int& nbr_tot){
	std::istringstream command(buffer);
	int pdg_1;
	int pdg_2;
	string separator;
	command >> pdg_1;
	command >> pdg_2;
	command >> separator;

	if ( separator!="->"){
		throw syntax_error(buffer);
	}



	vector<int> particles;

	particles.push_back(pdg_1);
	particles.push_back(pdg_2);

	std::istream_iterator<int> it(command);
	std::istream_iterator<int> end;

	copy(it,end,back_inserter(particles));

	nbr_tot=particles.size();
	return particles;

}




struct BLHA_options {
	std::string d_MatrixElementSquareType;
	std::string d_CorrectionType;
	std::string d_IRregularisation;
	std::string d_MassiveParticleScheme;
	std::string d_IRsubtractionMethod;
	std::string d_OperationMode;
	std::string d_SubdivideProcess;
	std::string d_CouplingStrippedOff;
	BLHA_options():
		d_MatrixElementSquareType("CHsummed"),
		d_CorrectionType("QCD"),
		d_IRregularisation("tHV"),
		d_MassiveParticleScheme("OnShell"),
		d_IRsubtractionMethod("None"),
		d_OperationMode("FullColor"),
		d_SubdivideProcess("No"),
		d_CouplingStrippedOff("Yes")
		{}
} Default_BLHA_options;


LH_interface::LH_interface() : d_bhi(new BH_interface){
	d_options.add("MatrixElementSquareType",Default_BLHA_options.d_MatrixElementSquareType);
	d_options.add("CorrectionType",Default_BLHA_options.d_CorrectionType);
	d_options.add("IRregularisation",Default_BLHA_options.d_IRregularisation);
	d_options.add("IRsubtractionMethod",Default_BLHA_options.d_IRsubtractionMethod);
	d_options.add("CouplingStrippedOff",Default_BLHA_options.d_CouplingStrippedOff);
	d_options.add("MassiveParticleScheme",Default_BLHA_options.d_MassiveParticleScheme);
	d_options.add("OperationMode",Default_BLHA_options.d_OperationMode);
	d_options.add("SubdivideProcess",Default_BLHA_options.d_SubdivideProcess);

	// there
	d_BHOptionsHandler_p=new BHOptionsHandler(d_bhi);
	d_LHOptionsHandler_p=new LHOptionsHandler(&d_options);

}

LH_interface LH_interface::UniqueInstance;

LH_interface::~LH_interface(){
	delete d_BHOptionsHandler_p;
	delete d_LHOptionsHandler_p;
	delete d_bhi;
};


bool LH_interface::read_LH_order(const char* order_filename,const char* contract_filename,bool generate_contract_file){
	ifstream order_file;
	ofstream contract_file;
	order_file.open(order_filename);
	if ( generate_contract_file ) { contract_file.open(contract_filename); };
	char buffer[256];
	while (order_file.getline(buffer,256)){
		;
		int pos=-1;
		while ( buffer[++pos] == ' ');

		switch ( buffer[pos]){
		case '#' : {
			if ( generate_contract_file ){
				contract_file << buffer << endl;
			}
		} break;
		case '0': case '1' : case '2':case '3': case '4': case '5': case '6': case '7': case '8': case '9': case '-': {
//			cout << "subprocess: " << buffer << endl;
			int nbr_tot;
			vector<int> particles;
			try {
				particles = read_subprocess(buffer,nbr_tot);
				//			copy(particles.begin(),particles.end(),ostream_iterator<int>(cout," ")); cout << endl;
				try {
					BH_Ampl* new_a=d_bhi->new_ampl(particles);
					int interfacePower=d_LHOptionsHandler_p->d_alphasPower;
					if (interfacePower> 0 && new_a->get_order_qcd()!= interfacePower){
						if ( generate_contract_file ){
							contract_file << buffer <<   "\t| ERROR: mismatch in the power of alphas, process has power "<<new_a->get_order_qcd()<<" but alphas power " << d_LHOptionsHandler_p->d_alphasPower << " is requested." << endl;
						}
					}
					interfacePower=d_LHOptionsHandler_p->d_alphaPower;
					if (interfacePower> 0 && new_a->get_order_qed()!= interfacePower){
						if ( generate_contract_file ){
							contract_file << buffer <<   "\t| ERROR: mismatch in the power of alpha, process has power "<<new_a->get_order_qed()<<" but alpha power " << d_LHOptionsHandler_p->d_alphasPower << " is requested." << endl;
						}
					}

					d_amplitudes.push_back(new_a);
					d_nbr_particles.push_back(nbr_tot);
					if ( generate_contract_file ){
						contract_file << buffer <<  "\t| 1 " << d_amplitudes.size() << endl ;
					}
				} catch (...) {
					if ( generate_contract_file ){
						contract_file << buffer <<   "\t| ERROR: sorry, cannot do that. " << endl;
					}
				}

			}
			catch (syntax_error& se) {
				cerr << se.what() << endl;
				if ( generate_contract_file ){
					contract_file << buffer <<   "\t| ERROR: syntax error. " << endl;
				}

			}

		} break;
		case ' ': case '\0': {
//			cout << "blank line " << endl;
			if ( generate_contract_file ){
				contract_file <<  "\n" ;
			}
		};
		break;
		default: {
//			cout << "setting: " << buffer << endl;
			if ( generate_contract_file ){
				contract_file << buffer << "\t" ;
			}
			bool success=get_setting(buffer,contract_file,generate_contract_file,d_bhi,&d_options);
			if ( !success ){
				cerr << "Option \'" << buffer << "\' lead to the termination of the initialization phase. Please fix the problem and run again." << endl;
				if ( generate_contract_file ){
					contract_file << "Option \'" << buffer << "\' lead to the termination of the initialization phase. Please fix the problem and run again." << endl;
				}
				return false;
			}
		} break;
		}
	}

	if ( generate_contract_file ){
		contract_file << "\n#\n# options\n#\n" ;
		d_options.display(contract_file,"# ");
		contract_file << "#\n";

	}
	return true;
}

void LH_interface::eval(int Label,BHinput& bhi,double (&result)[4]){
	UniqueInstance.d_bhi->operator ()(bhi);
	result[2]=d_amplitudes[Label-1]->get_finite();

	result[0]=d_amplitudes[Label-1]->get_double_pole();
	result[1]=d_amplitudes[Label-1]->get_single_pole();
	result[3]=1.0;//d_amplitudes[Label-1]->get_born();
}

void EvalSubprocess(int Label,double* Momenta,double mu2,double alpha_s,double alpha_ew,double* result){

	vector<vector<double> > moms;
	for (int i=0;i<LH_interface::Instance()->nbr_particles(Label);i++){
		moms.push_back(vector<double>());
		moms.back().push_back(Momenta[i*5+0]);
		moms.back().push_back(Momenta[i*5+1]);
		moms.back().push_back(Momenta[i*5+2]);
		moms.back().push_back(Momenta[i*5+3]);
		// don't need the mass
		//mons.back().push_back(Momenta[i*5+0]);
	}
	BHinput bhi(moms,sqrt(mu2));
	double result_array[4];
	double fac=1; // another valid choice would be fac(2.0*M_PI);
	LH_interface::Instance()->eval(Label,bhi,result_array);
	result[0]=result_array[0]/fac;
	result[1]=result_array[1]/fac;
	result[2]=result_array[2]/fac;
	result[3]=result_array[3];
};

std::vector<double> EvalSubprocessForPython(int Label,std::vector<double>& Momenta,double mu,double alpha_s,double alpha_ew){
	std::vector<double> res(4);
	EvalSubprocess(Label,&Momenta[0],mu,alpha_s,alpha_ew,&res[0]);

	return res;
};


int Init(const char* filename){
	LH_interface::Instance()->read_LH_order(filename,filename,false);

	// 1 means success

	return 1;
}

extern "C" {

void olp_evalsubprocess_(int* Label,double* Momenta,double *mu,double *param,double *result){
	return EvalSubprocess(*Label,Momenta,(*mu)*(*mu),param[0],param[1],result);
}

void OLP_EvalSubProcess(int* Label,double* Momenta,double *mu,double *param,double *result){
  return EvalSubprocess(*Label,Momenta,(*mu)*(*mu),param[0],param[1],result);
}

} /* extern C */

extern "C" {

void olp_start_(const char* filename,int *status){
	*status = Init(filename);
}

void OLP_Start(const char* filename,int *status){
	*status = Init(filename);
}

} /* extern C */




void SignContract(char* filename1,char* filename2){
	LH_interface::Instance()->read_LH_order(filename1,filename2,true);
}
void PrintHelp(std::ostream& os){
	LH_interface::Instance()->printHelp(os);
}

}

}
