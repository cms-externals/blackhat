/*!\file tree_amp.h
\brief Header for the rational vertices of the recursion for the rational terms
*/
#ifndef TREE_AMP_H_
#define TREE_AMP_H_

#include <complex>
#include <vector>
#include "mom_conf.h"
#include "eval_param.h"
#include "qd_suppl.h"

namespace BH {

// Enum taken from tree.cpp
enum zero_amplitudes { odd_nbr_q = -1, odd_nbr_q2 = -2, odd_nbr_l = -3, odd_nbr_l1 = -4, odd_nbr_l2 = -5 };

template <class T> std::complex<T> ZeroF(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	//_MESSAGE("Zero");
	return std::complex<T>(0.,0.);
}

template <int i1, int i2, class T> std::complex<T> A03gm(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	//_MESSAGE3("3-pt minus ",i1,i2);
	return (std::complex<T>(0.,1.)*pow(mc.spa(ind.at(i1),ind.at(i2)),4))/(mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(2),ind.at(0)));
}

template <int i1, int i2, class T> std::complex<T> A03gp(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	//_MESSAGE3("3-pt plus ",i1,i2);
	return (std::complex<T>(0.,1.)*pow(mc.spb(ind.at(i1),ind.at(i2)),4))/(mc.spb(ind.at(0),ind.at(2))*mc.spb(ind.at(2),ind.at(1))*mc.spb(ind.at(1),ind.at(0)));
}

template <int i1, int i2, class T> std::complex<T> A04g(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	return (std::complex<T>(0.,1.)*pow(mc.spa(ind.at(i1),ind.at(i2)),4))/(mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(2),ind.at(3))*mc.spa(ind.at(3),ind.at(0)));
}

template <int i1, int i2, class T> std::complex<T> A05g2m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	//_MESSAGE3("5 pt minus",i1,i2);
	return (std::complex<T>(0.,1.)*pow(mc.spa(ind.at(i1),ind.at(i2)),4))/(mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(2),ind.at(3))*mc.spa(ind.at(3),ind.at(4))*mc.spa(ind.at(4),ind.at(0)));
}

template <int i1, int i2, class T> std::complex<T> A05g2p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	//_MESSAGE3("5 pt plus",i1,i2);
	return (std::complex<T>(0.,-1.)*pow(mc.spb(ind.at(i1),ind.at(i2)),4))/(mc.spb(ind.at(0),ind.at(1))*mc.spb(ind.at(1),ind.at(2))*mc.spb(ind.at(2),ind.at(3))*mc.spb(ind.at(3),ind.at(4))*mc.spb(ind.at(4),ind.at(0)));
}

template <int i1, int i2, class T> std::complex<T> A06g2m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	return (std::complex<T>(0.,1.)*pow(mc.spa(ind.at(i1),ind.at(i2)),4))/(mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(2),ind.at(3))*mc.spa(ind.at(3),ind.at(4))*mc.spa(ind.at(4),ind.at(5))*mc.spa(ind.at(5),ind.at(0)));
}

template <int i1, int i2, class T> std::complex<T> A06g2p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	return (std::complex<T>(0.,1.)*pow(mc.spb(ind.at(i1),ind.at(i2)),4))/(mc.spb(ind.at(0),ind.at(1))*mc.spb(ind.at(1),ind.at(2))*mc.spb(ind.at(2),ind.at(3))*mc.spb(ind.at(3),ind.at(4))*mc.spb(ind.at(4),ind.at(5))*mc.spb(ind.at(5),ind.at(0)));
}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> std::complex<T> A06g2mt(momentum_configuration<T>& mc, std::vector<int>& ind)
{
	return A06g2m<i1,i2>(mc,ind);
}


template <class T> std::complex<T> ZeroF_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	//_MESSAGE("Zero");
	return std::complex<T>(0.,0.);
}

template <int i1, int i2, class T> std::complex<T> A03gm_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	//_MESSAGE3("3-pt minus ",i1,i2);
	return (std::complex<T>(0.,1.)*pow(ep.spa(i1,i2),4))/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,0));
}

template <int i1, int i2, class T> std::complex<T> A03gp_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	//_MESSAGE3("3-pt plus ",i1,i2);
	return (std::complex<T>(0.,1.)*pow(ep.spb(i1,i2),4))/(ep.spb(0,2)*ep.spb(2,1)*ep.spb(1,0));
}

template <int i1, int i2, class T> std::complex<T> A04g_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	return (std::complex<T>(0.,1.)*pow(ep.spa(i1,i2),4))/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,0));
}

template <int i1, int i2, class T> std::complex<T> A05g2m_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	//_MESSAGE3("5 pt minus",i1,i2);
	return (std::complex<T>(0.,1.)*pow(ep.spa(i1,i2),4))/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,0));
}

template <int i1, int i2, class T> std::complex<T> A05g2p_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	//_MESSAGE3("5 pt plus",i1,i2);
	return (std::complex<T>(0.,-1.)*pow(ep.spb(i1,i2),4))/(ep.spb(0,1)*ep.spb(1,2)*ep.spb(2,3)*ep.spb(3,4)*ep.spb(4,0));
}

template <int i1, int i2, class T> std::complex<T> A06g2m_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	return (std::complex<T>(0.,1.)*pow(ep.spa(i1,i2),4))/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)*ep.spa(5,0));
}

template <int i1, int i2, class T> std::complex<T> A06g2p_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	return (std::complex<T>(0.,1.)*pow(ep.spb(i1,i2),4))/(ep.spb(0,1)*ep.spb(1,2)*ep.spb(2,3)*ep.spb(3,4)*ep.spb(4,5)*ep.spb(5,0));
}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> std::complex<T> A06g2mt_eval(const eval_param<T>& ep, const mass_param_coll& mpc)
{
	return A06g2m_eval<i1,i2>(ep,mpc);
}

}

#endif /*TREE_AMP_H_*/
