/*
 * helcode.h
 *
 *  Created on: 24-Nov-2008
 *      Author: daniel
 */

#ifndef HELCODE_H_
#define HELCODE_H_



namespace BH {

class process;

int helcode_g(const process& pro);
int helcode_2q(const process& pro);
int helcode_4q(const process& pro);
int helcode_2q2l(const process& pro);
int helcode_2q2Q(const process& pro);

int helcode_2q1y(const process& pro);

int helcode_2q1_2q2_2q3(const process& pro);
int helcode_4q1_2q2(const process& pro);

int helcode_2s(const process& pro);
long helcode_2Gsc(const process& pro);
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

int helcode_Ng1ph(const process& pro);
long helcode_phi_1q(const process& pro);
long helcode_phi_2q(const process& pro);
int helcode_phi_2q2Q(const process& pro);

}



#endif /* HELCODE_H_ */
