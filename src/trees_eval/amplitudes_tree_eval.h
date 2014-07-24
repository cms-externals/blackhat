#ifndef AMPLITUDES_TREE_EVAL_H_
#define AMPLITUDES_TREE_EVAL_H_



#define _FAST_COMPILE 1  /* if 1 the 7 and 8 point gluons are not compiled (they return 0) */

#include <complex>
#include <vector>

namespace BH {

class process;
class particle;
class mass_param;
class mass_param_coll;
template <class T> class eval_param;


template <class T> std::complex<T> (*A_Tree_Ptr_eval(const process& pro))(const eval_param<T>&, const mass_param_coll&);
class particle_ID;


template <class T> std::complex<T> (*A2q1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q3g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q4g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q2Q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q1g2Q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q2g2Q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q3g2Q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A4q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A4q1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A4q2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A4q3g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A3g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A4g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A5g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A6g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A7g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A8g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2s1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2s2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2s3g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2s4g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2s2q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2s2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2QM2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2q2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q1g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q2g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q3g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q4g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q5g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2QM1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2QM2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2QM2q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A1QM1qs_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A1QM1q1gs_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T>  (*A2q2l2G_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) ;

template <class T> std::complex<T>  (*A2q2Q2l_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) ;
template <class T> std::complex<T>  (*A2q2Q1g2l_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) ;
template <class T> std::complex<T>  (*A2q2Q2g2l_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) ;
template <class T> std::complex<T>  (*A2q2Q3g2l_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) ;

template <class T> std::complex<T>  (*A2q1g1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2s2q2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2QM1g1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2QM1g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);


template <class T> std::complex<T>  (*A2q2g1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A2q2Q1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T>  (*A2q2G1y_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2LM1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2LM2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2LM2G_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A1LM1Gs_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A1LM1G1gs_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A1QM1q1G1L_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2LM2q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2LM2q1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2QM2q1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2s2G_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T>  (*A1ph2g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph3g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph1g2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2g2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2q2Q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2q2QM_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2sc_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2sc1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2QM1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph1QM1sc1q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2sc2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);

template <class T> std::complex<T> (*A2sc1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2sc2g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2sc3g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2scm1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2SM_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T>  (*A1ph2SM1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A1ph1QM1SM1q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A1ph2SM2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
    
template <class T> std::complex<T>  (*A1ph2Gsc_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2Gsc1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A1Gsc1gM1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);

/* special purpose trees */
template <class T> std::complex<T> (*A4g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A5g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q2g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q3g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);
template <class T> std::complex<T> (*A2q1g2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&);



int helcode_g(const process& pro);
int helcode_2q(const process& pro);
int helcode_4q(const process& pro);
int helcode_2q2l(const process& pro);
int helcode_2q2Q(const process& pro);

int helcode_2q1y(const process& pro);

int helcode_2q1_2q2_2q3(const process& pro);
int helcode_4q1_2q2(const process& pro);

int helcode_2s(const process& pro);
int helcode_2qs_massive(const process& pro);
int helcode_2Q2qs_massive(const process&);
int helcode_2Q1g1y_massive(const process& pro);
int helcode_2Q2qs_lepton_massive(const process& pro);
long helcode_2q2l2Q(const process& pro);
long helcode_2q2l2G(const process& pro);
long helcode_2q2G1y(const process& pro);

int helcode_2Ls_massive(const process& pro);
int helcode_2L2qs_massive(const process& pro);
int helcode_2L2Gs_massive(const process& pro);
int helcode_2L2Gs_lepton_massive(const process& pro);
int helcode_1Q1q1G1L_massive(const process& pro);
int helcode_2L2qs1y_massive(const process& pro);
int helcode_2Q2qs1y_massive(const process& pro);

long helcode_phi_1q(const process& pro);
long helcode_phi_SM_1q(const process& pro);
long helcode_phi_2q(const process& pro);
	
long helcode_2qs_massless(const process& pro);

size_t nbr_of_flavors(const process& pro,const particle& type);
process fix_flavors(const process& pro);

}
#endif /*AMPLITUDES_TREE_EVAL_H_*/
