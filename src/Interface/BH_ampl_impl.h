/*
 * BH_ampl_impl.h
 *
 *  Created on: 28 Jul 2009
 *      Author: daniel
 */

#ifndef BH_AMPL_IMPL_H_
#define BH_AMPL_IMPL_H_

#include "BH_Ampl.h"
#include "BH_typedefs.h"
#include "assembly.h"

namespace BH {

class BH_interface_impl;

class BH_Ampl_impl : public BH_Ampl {
	int d_NbrExtParticles;
	int d_PowersOfAlphaS;
	int d_PowersOfAlphaQED;
	int d_GeVdim;
	double d_catch_large_ME2;
	double d_catch_large_born_ME2;
	double d_scheme_shift;
	double d_renormalization_shift;
	int d_eval_id;
	int d_mc_ID_born;
	int d_mc_ID_virt;
	std::vector<int> d_symmetric_finalstate_partons;
	Virtual_SME* d_VSM_p;
protected:
	C d_last_result_2;
	C d_last_result_1;
	C d_last_result_0;
	C d_last_result_tree;
	C d_last_result_HP_2;
	C d_last_result_HP_1;
	C d_last_result_HP_0;
	C d_last_result_VHP_2;
	C d_last_result_VHP_1;
	C d_last_result_VHP_0;
	CHP d_last_result_HP_tree;
	CVHP d_last_result_VHP_tree;

	std::complex<double> d_factor;
	BH_interface_impl* d_parent_p;
public:
	BH_Ampl_impl(Virtual_SME* VSM,BH_interface_impl* bhi,int NbrExtParticles, int PowersOfAlphaS, int PowersOfAlphaQED,int GevDim,const std::complex<double>& factor,double scheme_shift=0,double renormalization_shift=0,const std::vector<int>& symmetric_finalstate_partons = std::vector<int>())  :
		d_VSM_p(VSM) ,d_parent_p(bhi),
		d_NbrExtParticles(NbrExtParticles),
		d_PowersOfAlphaS(PowersOfAlphaS), d_PowersOfAlphaQED(PowersOfAlphaQED),
		d_GeVdim(GevDim),
		d_factor(factor),
		d_catch_large_ME2(BH::settings::BH_interface_settings::s_catch_large_ME2),
		d_catch_large_born_ME2(0.),
		d_eval_id(0),
		d_scheme_shift(scheme_shift),
		d_renormalization_shift(renormalization_shift),
		d_mc_ID_born(0),
		d_symmetric_finalstate_partons(symmetric_finalstate_partons),
		d_mc_ID_virt(0)
		{
		//approximate equation for size of born ME2
	//	d_catch_large_born_ME2=std::pow((double(NbrExtParticles)-2.),2*(NbrExtParticles-3))*(BH::settings::BH_interface_settings::s_catch_large_ME2);
        if(d_VSM_p!=0) d_VSM_p->set_order(PowersOfAlphaS, PowersOfAlphaQED);

        switch(BH::settings::BH_interface_settings::s_BH_interface_eval_mode){
            case BH::settings::BH_interface_settings::polarization_eval:{

	            std::cout<<"*****************************************"<<std::endl;
            	std::cout<<" WARNING:"<<std::endl;
            	std::cout<<" Using evaluation in terms polarization"<<std::endl;
            	std::cout<<" fractions for virtual amplitudes."<<std::endl;
            	std::cout<<" ASSUMING: lepton in 3rd position and "<<std::endl;
              	std::cout<<" anti-lepton in 4th position in process specification."<<std::endl;
            	std::cout<<" Lepton angle with repect to W-boson parametrized "<<std::endl;
            	std::cout<<" by spherical cooridnates theta and phi."<<std::endl<<std::endl;
                std::cout<<" values of A_n: "<<std::endl;
                std::cout<<" -1:    Sum over all polarization fractions;"<<std::endl;
                std::cout<<"        useful as consistency check."<<std::endl;
                std::cout<<" 0-8:   Projection to polarization fraction "<<std::endl;
                std::cout<<"        multiplied by its own angular function."<<std::endl;
                std::cout<<" 10-18: Projection to polarization fraction " <<std::endl;
                std::cout<<"        multiplied by the A_0 angluar dependence: (1+cos(theta)^2) "<<std::endl<<std::endl;
            	std::cout<<" polarization fraction: A_n with n = ";
                std::cout<<(BH::settings::BH_interface_settings::s_W_polarization_A)<<std::endl;
                std::cout<<std::endl;
                std::cout<<"*****************************************"<<std::endl;
            } break;
            default: break;
        }
        
		};
	virtual double get_born();
	virtual double get_single_pole();
	virtual double get_double_pole();
	virtual double get_finite();
	virtual double get_born_HP();
	virtual double get_single_pole_HP();
	virtual double get_double_pole_HP();
	virtual double get_finite_HP();
	virtual double get_born_VHP();
	virtual double get_single_pole_VHP();
	virtual double get_double_pole_VHP();
	virtual double get_finite_VHP();
    void set_partial_born();
    virtual void get_map(
	std::vector<std::vector<int> >& permutation,
        std::vector<std::vector<int> >& helicity);
    virtual void get_vals(
	std::vector<double* >& re_ampl,
        std::vector<double* >& im_ampl);
    virtual int get_order_qcd(){return d_VSM_p->get_order_qcd();};
    virtual int get_order_qed(){return d_VSM_p->get_order_qed();};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_impl(){delete d_VSM_p; }
	virtual void dry_run(){d_VSM_p->dry_run();};
	virtual void print_process(){};
	virtual void print_virt_info();
	virtual void print_born_info();
	bool vanishing_VSM(){return (d_VSM_p==0);};
};

}


#endif /* BH_AMPL_IMPL_H_ */
