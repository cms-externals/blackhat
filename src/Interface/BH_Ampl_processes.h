/*
 * BH_Ampl_processes.h
 *
 *  Created on: Aug 26, 2008
 *      Author: daniel
 */

#ifndef BH2SHERPA_PROCESSES_H_
#define BH2SHERPA_PROCESSES_H_

#include "BH_typedefs.h"
#include "BH_Ampl.h"
#include "BH_interface.h"
#include "assembly.h"
#include "settings.h"
#include "Interface/BH_ampl_impl.h"
#ifndef BH_PUBLIC
#include "ME2.h"
#endif
#include "ME2_from_file.h"

namespace BH {

#ifndef BH_PUBLIC
class BH_Ampl_automated : public BH_Ampl_impl {
public:
	BH_Ampl_automated(
		std::vector<std::pair<int,int> > particle_labels,
		ME2_factory* me2_factory,
		approx born_or_virt,
		BH_interface_impl* bhi,
		int Nbr_Ext_Particles,
		int NbrPowersOfAlphaS,
		int NbrPowersOfAlphaQED,
		int GeVdim,
		double scheme_shift,
		double renormalization_shift);
	std::vector<int> momenta_assignment;
	std::vector<int> particles;
	void print_process(){std::cout<<particles<<std::endl; std::cout<<momenta_assignment<<std::endl;}
	virtual ~BH_Ampl_automated(){};
};

#endif

class BH_Ampl_data : public BH_Ampl_impl {
	int d_color;
	double renormalization_shift;
    double d_lc;
    double d_fc;
public:
	BH_Ampl_data(
		std::vector<std::pair<int,int> > particle_labels,
		QCDorder lo_or_nlo,
		BH_interface_impl* bhi,
		int Nbr_Ext_Particles,
		int NbrPowersOfAlphaS,
		int NbrPowersOfAlphaQED,
		int GeVdim,
		double scheme_shift,
		double renormalization_shift,
        double d_lc,       
        double d_fc        
        );
	std::vector<int> momenta_assignment;
	std::vector<int> particles;
	void print_process(){std::cout<<particles<<std::endl; std::cout<<momenta_assignment<<std::endl;}
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_data(){};
};



class BH_Ampl_2q2l : public BH_Ampl_impl {
public:
	BH_Ampl_2q2l(bool up_down_quark, int photonZW,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q qb e eb
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q qb l lb; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2l(){};
};

class BH_Ampl_2q1g1y : public BH_Ampl_impl {
public:
	BH_Ampl_2q1g1y(bool up_down_quark,int color, const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g qb photon
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q g qb y; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q1g1y(){};
};

class BH_Ampl_2q2g1y : public BH_Ampl_impl {
public:
	BH_Ampl_2q2g1y(bool up_down_quark,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g g qb photon
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q g g qb y; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2g1y(){};
};

class BH_Ampl_2q2Q1y : public BH_Ampl_impl {
public:
	BH_Ampl_2q2Q1y(bool up_down_quark,int case4q,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb photon
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q Qb Q qb y; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2Q1y(){};
};


class BH_Ampl_2q3g1y : public BH_Ampl_impl {
	int d_color;
public:
	BH_Ampl_2q3g1y();
	BH_Ampl_2q3g1y(bool up_down_quark, int photonZW,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g g g qb y
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q g g g qb y; "<<momenta_assignment<<std::endl;};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_2q3g1y(){};
};


class BH_Ampl_2q2Q1g1y : public BH_Ampl_impl {
	int d_color;
public:
	BH_Ampl_2q2Q1g1y(int photonZW,int case4q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb g y
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q Qb Q qb g y; "<<momenta_assignment<<std::endl;};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_2q2Q1g1y(){};
};



class BH_Ampl_2q1g2l : public BH_Ampl_impl {
public:
	BH_Ampl_2q1g2l(bool up_down_quark, int photonZW,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g qb l lb
	void print_process(){std::cout<<"q g qb l lb; "<<momenta_assignment<<std::endl;};
	std::vector<int> momenta_assignment;
	virtual ~BH_Ampl_2q1g2l(){};
};

class BH_Ampl_2q2g2l : public BH_Ampl_impl {
public:
	BH_Ampl_2q2g2l(bool up_down_quark, int photonZW,int color, int tree_color, const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g g qb l lb
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q g g qb l lb; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2g2l(){};
};

class BH_Ampl_2q2Q2l : public BH_Ampl_impl {
public:
	BH_Ampl_2q2Q2l(int photonZW,int color, int tree_color,int case4q,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb photon
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q Qb Q qb l lb; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2Q2l(){};
};

class BH_Ampl_2q3g2l : public BH_Ampl_impl {
	int d_color;
public:
	BH_Ampl_2q3g2l();
	BH_Ampl_2q3g2l(bool up_down_quark, int photonZW,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g g g qb l lb
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q g g g qb l lb; "<<momenta_assignment<<std::endl;};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_2q3g2l(){};
};




class BH_Ampl_2q2Q1g2l : public BH_Ampl_impl {
	int d_color;
protected:
public:
	BH_Ampl_2q2Q1g2l(int photonZW,int case4q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb g l lb
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q Qb Q qb g l lb; "<<momenta_assignment<<std::endl;};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_2q2Q1g2l(){};
};


class BH_Ampl_2q4g2l : public BH_Ampl_impl {
	int d_color;
public:
	BH_Ampl_2q4g2l(bool up_down_quark, int photonZW,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g g g g qb l lb 
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q g g g g qb l lb; "<<momenta_assignment<<std::endl;};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_2q4g2l(){};
};


class BH_Ampl_2q2Q2g2l : public BH_Ampl_impl {
	int d_color;
protected:
public:
	BH_Ampl_2q2Q2g2l(int photonZW,int case4q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb g g l lb
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q Qb Q g g qb l lb; "<<momenta_assignment<<std::endl;};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_2q2Q2g2l(){};
};


class BH_Ampl_6q2l : public BH_Ampl_impl {
	int d_color;
protected:
public:
	BH_Ampl_6q2l(int photonZW,int case6q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q q1b q1 q2b q2 qb l lb
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q q1b q1 q2b q2 qb l lb; "<<momenta_assignment<<std::endl;};
	virtual double getScaleVariationCoefficient(int logMuOrder);
	virtual ~BH_Ampl_6q2l(){};
};


class BH_Ampl_2q5g2l : public BH_Ampl_impl {
public:
	BH_Ampl_2q5g2l(bool up_down_quark, int photonZW,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q g g g g g qb l lb 
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q g g g g g qb l lb; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q5g2l(){};
};


class BH_Ampl_2q2Q3g2l : public BH_Ampl_impl {
protected:
public:
	BH_Ampl_2q2Q3g2l(int photonZW,int case4q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb g g g l lb
	void print_process(){std::cout<<"q Qb Q qb g g g l lb; "<<momenta_assignment<<std::endl;};
	std::vector<int> momenta_assignment;

	virtual ~BH_Ampl_2q2Q3g2l(){};
};


class BH_Ampl_6q1g2l : public BH_Ampl_impl {
protected:
public:
	BH_Ampl_6q1g2l(int photonZW,int case6q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q q1b q1 q2b q2 qb photon
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q q1b q1 q2b q2 qb g l lb; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_6q1g2l(){};
};


class BH_Ampl_2q2Q1g : public BH_Ampl_impl {
public:
	BH_Ampl_2q2Q1g(int case4q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb g
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q Qb Q qb g; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2Q1g(){};
};



class BH_Ampl_2q2Q2g : public BH_Ampl_impl {
public:
	BH_Ampl_2q2Q2g(int case4q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q Qb Q qb g g
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q Qb Q qb g g; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2Q2g(){};
};



class BH_Ampl_2q2Q2P : public BH_Ampl_impl {
public:
	BH_Ampl_2q2Q2P(int case6q,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi);
	// storing the momenta in the BH canonical order q q1b q1 q2b q2 qb 
	std::vector<int> momenta_assignment;
	void print_process(){std::cout<<"q qb1 q1 qb2 q2 qb; "<<momenta_assignment<<std::endl;};
	virtual ~BH_Ampl_2q2Q2P(){};
};


}

#endif /* BH2SHERPA_PROCESSES_H_ */
