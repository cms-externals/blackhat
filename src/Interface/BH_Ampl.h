#ifndef SHERPA_H_
#define SHERPA_H_

#include <vector>
#include <fstream>
#include <iostream>

namespace BH {

struct BHinput {
	std::vector<std::vector<double> > m_momenta;
	double m_mu;
	BHinput(const std::vector<std::vector<double> >& momenta, double mu): m_momenta(momenta),m_mu(mu) {} ;
};


class BH_Ampl {
public:
		virtual double operator()(BHinput& in);
		virtual double get_single_pole();
		virtual double get_double_pole();
		virtual double get_finite(){return 0.;};
		virtual double get_born(){return 0.;};
		virtual double get_single_pole_HP(){return 0.;};
		virtual double get_double_pole_HP(){return 0.;};
		virtual double get_finite_HP(){return 0.;};
		virtual double get_born_HP(){return 0.;};
		virtual double get_single_pole_VHP(){return 0.;};
		virtual double get_double_pole_VHP(){return 0.;};
		virtual double get_finite_VHP(){return 0.;};
		virtual double get_born_VHP(){return 0.;};
        // for real subtraction term
        virtual void set_partial_born(){return;};
        virtual void get_map(
            std::vector<std::vector<int> >& permutation,
            std::vector<std::vector<int> >& helicity){return;};
        virtual void get_vals(
            std::vector<double* >& re_ampl,
            std::vector<double* >& im_ampl){return;};
        virtual int get_order_qcd(){return 0;};
        virtual int get_order_qed(){return 0;};

		//! returns the coefficient of the log(mu^2/mu_0^2) parametrisation of the scale variation
		/**
		 * the parrametrisation is (note the 1/2)
		 *
		 * dsigma^loop(mu) = c_0 + c_1*log(mu^2/mu_0^2+1/2*log^2(mu^2/mu_0^2))
		 *
		 *	logMuOrder =0 returns c_0,  logMuOrder=1 returns c_1 and logMuOrder=2 returns c_2
		 *
		 * */
		virtual double getScaleVariationCoefficient(int logMuOrder){return 0;};

		virtual ~BH_Ampl(){}
		virtual void dry_run(){};
};

class settings_table;
template <class T> class momentum_configuration;

#if 0

class BH_factory {
	settings_table* d_settings_p;
public:
	BH_factory();
	~BH_factory();

	//! Return a pointer to a new BH_Ampl object.
	BH_Ampl* new_ampl(const std::vector<int>&);
	//! Sets the named setting to the value provided.
	/** The type is deduced from the argument. It is the responsibility of the caller to ensure that the type matches the type of the setting */
	template <class T> void set(const std::string& name,T value);
	//! Prints a table of the settings and their current values.
	void print_settings();

};

#endif


}

#endif /*SHERPA_H_*/
