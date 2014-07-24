/*
 * triangle_subtraction.h
 *
 *  Created on: 13 Nov 2009
 *      Author: daniel
 */

#ifndef TRIANGLE_SUBTRACTION_H_
#define TRIANGLE_SUBTRACTION_H_

namespace BH {
namespace cut {
namespace Darren {

template <class T> struct tri_sub_info;

template <class T> struct tri_sub_info {
	vector<complex<T> > *triboxcoeff, *triboxpoles, *denfac;
	Cmom<T> v1old, v2old;
    complex<T> gamma_old, *orig_coeffs;
    const complex<T>* eval_pts;
	const complex<T>* eval_ypts;
	const complex<T>* ymatrixpoint;
    complex<T> *Nysolp;
    complex<T> *Nysolm;
    complex<T> alp1, alp2;
    int reverse;
};


template <class cutDType,int CPOINTS,int TPOINTSBUB> class Normal_Triangle_Subtraction {
public:
	template <class T> complex<T> get_sub_terms_work_3m(const cutDType& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
	template <class T> complex<T> get_sub_terms_work_pm(const cutDType& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
	template <class T> complex<T> get_sub_terms_work_3m(const cutDType& cd,const eval_param<T>&, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
	template <class T> complex<T> get_sub_terms_work_pm(const cutDType& cd,const eval_param<T>&, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
};

template <class cutDType,int CPOINTS,int TPOINTSBUB,int YPOINTS> class General_Triangle_Subtraction {
public:
	template <class T> complex<T> get_sub_terms_work_3m(const cutDType& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
	template <class T> complex<T> get_sub_terms_work_pm(const cutDType& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
	template <class T> complex<T> get_sub_terms_work_3m(const cutDType& cd,const eval_param<T>&, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
	template <class T> complex<T> get_sub_terms_work_pm(const cutDType& cd,const eval_param<T>&, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi);
};



} /* Darren */

} /* cut */

} /* BH */


#endif /* TRIANGLE_SUBTRACTION_H_ */
