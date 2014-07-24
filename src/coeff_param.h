/*
 * coeff_param.h
 *
 *  Created on: 24 Jul 2009
 *      Author: daniel
 */

#ifndef COEFF_PARAM_H_
#define COEFF_PARAM_H_

#include <complex>
#include "spinor.h"


namespace BH {

#define RATPOINTS_STD 3
#define RATPOINTS_SUB_STD 2

#define BUBPOINTS_EXT_STD 3
#define TRIPOINTS_EXT_STD 7


template <class T> class coeffparam_old
{
public:
	size_t K1t, K2t, K3t, K4t;
	size_t K1flat, K2flat, K1flatb;
	size_t chi;
	size_t vec1, vec2;
	std::complex<T> gamma, gammab, S1, S2, S3, S4;
	std::complex<T> f1, f2, alp1, alp2;
	std::complex<T> triC0res[7];
	std::complex<T> Nboxpoles[2], boxcoeffs[2], denfac1, denfac2;
	int masslessleg, tri_corner, box_corner, box_type, ltype;
	int nsols, polepos;
	std::complex<T> nboxpart1;
	std::complex<T> alpha1,alpha2;
	std::complex<T> vec_alp[2];
	std::complex<T> alp_11, alp_12, alp_21, alp_22;

	std::vector<int> tri_legs, box_legs;

	long int mcID;
};

template <class T> class coeffparam_old_mass : public coeffparam_old<T>
{
public:
	std::complex<T> Nboxpoles_mass[2*RATPOINTS_STD], boxcoeffs_mass[2*RATPOINTS_STD];
	std::complex<T> coeffsret_mass[(TRIPOINTS_EXT_STD+1)*RATPOINTS_STD];
	std::complex<T> triC0res_mass[7*RATPOINTS_SUB_STD];
	std::complex<T> pentcoeff_mass;
	std::complex<T> pentpole_mass;
	std::complex<T> coeffsret_mu;
	std::complex<T> bub_acc_pieces[3];
	double accuracy;
	size_t K5t;
	std::complex<T> S5;
	std::complex<T> m0sol;
	std::complex<T> gammap;
	size_t pent_corner, pent_type;
	size_t vertex_ref;

	std::vector<int> pent_legs;
};


template <class T, int CPOINTS=7> class coeffparam
{
public:
	std::complex<T> gamma, gammab, gammap, S1, S2, S4;
	std::complex<T> f1, f2, alp1, alp2;
	std::complex<T> triC0res[CPOINTS];
	std::complex<T> Nboxpoles[2], boxcoeffs[2], denfac1, denfac2;
	int masslessleg, tri_corner, box_corner, box_type, ltype;
	int nsols, polepos;
	std::complex<T> nboxpart1;
	std::complex<T> vec_alp[6];

	momentum<std::complex<T> > K1, K2, K4;
	Cmom<T> vec1c, vec2c;
	Cmom<T> K1flatc, K2flatc, chic, K1flatbc;

	std::complex<T> amp_err;

	std::vector<int> tri_legs, box_legs;

	size_t legs[4];

	int reverse;

	long int mcID;

	RVHP tri_sub_acc;
};

template <class T, int EXTMUPOINTS=RATPOINTS_STD, int EXTMUPOINTS_SUB=RATPOINTS_SUB_STD, int CPOINTS=7, int TRIPOINTS_EXT=TRIPOINTS_EXT_STD> class coeffparam_mass : public coeffparam<T,CPOINTS>
{
public:
	std::complex<T> Nboxpoles_mass[2*EXTMUPOINTS], boxcoeffs_mass[2*EXTMUPOINTS];
	std::complex<T> coeffsret_mass[(TRIPOINTS_EXT+1)*EXTMUPOINTS];
	std::complex<T> triC0res_mass[CPOINTS*EXTMUPOINTS_SUB];
	std::complex<T> pentcoeff_mass;
	std::complex<T> pentpole_mass;
	std::complex<T> coeffsret_mu;
	std::complex<T> m0sol;

	std::complex<T> bub_acc_pieces[3];
	double accuracy;

	momentum<std::complex<T> > vec1c_mass, vec2c_mass;
	std::complex<T> S5;
	momentum<std::complex<T> > K5;

	size_t pent_corner, pent_type;
	size_t vertex_ref;

	std::vector<int> pent_legs;
};



}

#endif /* COEFF_PARAM_H_ */
