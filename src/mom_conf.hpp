/*
 * mom_conf.hpp
 *
 *  Created on: 8 Mar 2010
 *      Author: daniel
 */

#ifndef MOM_CONF_HPP_
#define MOM_CONF_HPP_

#include "mom_conf.h"
#include "BH_error.h"
#include <iostream>

using std::complex;

namespace BH {

template <class key_type, class T> void find_and_stat(const std::string& tag, key_type& key,momentum_configuration<T>* mc ,const  momentum<complex<T> >& mom, size_t& res ) {


	typename std::map<std::string,typename momentum_configuration<T>::hash_stat>::iterator pos_map;
	pos_map=momentum_configuration<T>::hash_stat_table.find(tag);

	if (pos_map == momentum_configuration<T>::hash_stat_table.end()){
		momentum_configuration<T>::hash_stat_table.insert(std::pair<std::string,typename momentum_configuration<T>::hash_stat>(tag,typename momentum_configuration<T>::hash_stat()));
		pos_map=momentum_configuration<T>::hash_stat_table.find(tag);
	}

	if ( ! mc->get_label(key,res) ){
		pos_map->second.nbr_no_hit++;

		res=mc->insert(mom);
		mc->put_label(key,res);
	}
	pos_map->second.nbr_try++;
}

//-------------------------------------------------------------------------------------------------------


template <class T> momentum_configuration<T>::momentum_configuration(size_t i): nbr(0), _offset(0), _parent(0) {
	ps.reserve(i);
	ms.reserve(i);
}

template <class T> momentum_configuration<T>::momentum_configuration(const std::vector<Cmom<T> >& pl ): _offset(0), _parent(0) {
	for (size_t j=0;j<pl.size();j++){
		ps.push_back(pl[j]);ms.push_back(pl[j]*pl[j]);
	}
	nbr=pl.size();
}

template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2): _offset(0), _parent(0)  {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3): _offset(0), _parent(0)  {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5): _offset(0), _parent(0)  {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5,const Cmom<T>& p6): _offset(0), _parent(0)  {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	ps.push_back(p6);ms.push_back(p6*p6);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5,const Cmom<T>& p6,const Cmom<T>& p7): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	ps.push_back(p6);ms.push_back(p6*p6);
	ps.push_back(p7);ms.push_back(p7*p7);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5,const Cmom<T>& p6,const Cmom<T>& p7,const Cmom<T>& p8): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	ps.push_back(p6);ms.push_back(p6*p6);
	ps.push_back(p7);ms.push_back(p7*p7);
	ps.push_back(p8);ms.push_back(p8*p8);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5,const Cmom<T>& p6,const Cmom<T>& p7,const Cmom<T>& p8,const Cmom<T>& p9): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	ps.push_back(p6);ms.push_back(p6*p6);
	ps.push_back(p7);ms.push_back(p7*p7);
	ps.push_back(p8);ms.push_back(p8*p8);
	ps.push_back(p9);ms.push_back(p9*p9);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5,const Cmom<T>& p6,const Cmom<T>& p7,const Cmom<T>& p8,const Cmom<T>& p9,const Cmom<T>& p10): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	ps.push_back(p6);ms.push_back(p6*p6);
	ps.push_back(p7);ms.push_back(p7*p7);
	ps.push_back(p8);ms.push_back(p8*p8);
	ps.push_back(p9);ms.push_back(p9*p9);
	ps.push_back(p10);ms.push_back(p10*p10);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5,const Cmom<T>& p6,const Cmom<T>& p7,const Cmom<T>& p8,const Cmom<T>& p9,const Cmom<T>& p10,const Cmom<T>& p11): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	ps.push_back(p6);ms.push_back(p6*p6);
	ps.push_back(p7);ms.push_back(p7*p7);
	ps.push_back(p8);ms.push_back(p8*p8);
	ps.push_back(p9);ms.push_back(p9*p9);
	ps.push_back(p10);ms.push_back(p10*p10);
	ps.push_back(p11);ms.push_back(p11*p11);
	nbr=ps.size();
}
template <class T> momentum_configuration<T>::momentum_configuration(const Cmom<T>& p1,const Cmom<T>& p2,const Cmom<T>& p3,const Cmom<T>& p4,const Cmom<T>& p5,const Cmom<T>& p6,const Cmom<T>& p7,const Cmom<T>& p8,const Cmom<T>& p9,const Cmom<T>& p10,const Cmom<T>& p11,const Cmom<T>& p12): _offset(0), _parent(0) {
	ps.push_back(p1);ms.push_back(p1*p1);
	ps.push_back(p2);ms.push_back(p2*p2);
	ps.push_back(p3);ms.push_back(p3*p3);
	ps.push_back(p4);ms.push_back(p4*p4);
	ps.push_back(p5);ms.push_back(p5*p5);
	ps.push_back(p6);ms.push_back(p6*p6);
	ps.push_back(p7);ms.push_back(p7*p7);
	ps.push_back(p8);ms.push_back(p8*p8);
	ps.push_back(p9);ms.push_back(p9*p9);
	ps.push_back(p10);ms.push_back(p10*p10);
	ps.push_back(p11);ms.push_back(p11*p11);
	ps.push_back(p12);ms.push_back(p12*p12);
	nbr=ps.size();
}

template <class T> void momentum_configuration<T>::lance_print() {
	std::cout << "setkin := proc()\nglobal K;\n";
	for (size_t i=1;i<=nbr;i++){
		std::cout <<"K[" << i << "]:=vector(["
		<< this->p(i).E() << ","
		<< this->p(i).Z() << ","
		<< this->p(i).X() << ","
		<< this->p(i).Y() << "]);"
		<< std::endl;
	}
	std::cout << "RETURN();\n end:" << std::endl ;

}


template <class T> std::ostream& operator<<(std::ostream& s, const momentum_configuration<T>& mc){
	s<<"(";
	for (size_t j=0;j< mc.nbr-1;j++){
		s << mc.ps[j].P() << ",";
	}
	return s << mc.ps[mc.nbr-1].P() <<")";
}


#if _WITH_CACHING

//! Prints the cache of T values
template <class T> void momentum_configuration<T>::DisplayCache() {
	_MESSAGE3("Cache contains ",cache.size()," elements.");
#if _USE_GCC
	typename std::hash_map<std::string, std::complex<T>, std::hash<std::string> >::iterator pos;
#endif
#if _USE_SUN_CC || _USE_PGCC
	typename std::hash_map<string, complex<T> >::iterator pos;
#endif
	for ( pos = cache.begin();  pos != cache.end();  pos++){
		std::cout << pos->first << ": " << pos->second << std::endl;
	}
}
//! Prints the cache for integer labels
template <class T> void momentum_configuration<T>::DisplayLabelCache() {
	_MESSAGE3("Label Cache contains ",labelscache.size()," elements.");
	LabelCache::iterator pos;
	for ( pos = labelscache.begin();  pos != labelscache.end();  pos++){
		std::cout << pos->first << ": " << pos->second << " " << p(pos->second) << std::endl;
	}
}

#endif /*_WITH_CACHING */

//! sum of momenta
/**
The sum of the momenta is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param v vector of plabel of momenta
\return integer index in the momentum_configuration<T> of the sum.
*/
template <class T> int momentum_configuration<T>::Sum(const std::vector<plabel>& v){
#if _WITH_CACHING == 1
	vector<int> indices;

	for (size_t i=0;i<v.size();i++){indices.push_back(v[i].ind());}
	sort (indices.begin(),indices.end());
	string key=GenKey("Sum",indices);
	size_t res;
#if _CACHE_STATISTICS
	momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
			for (size_t i=0;i< indices.size();i++){
				sum+=mom(indices[i]);
			}
	find_and_stat(string("Sum_plabel"),  key,this ,sum, res );
#else
	if ( !get_label(key,res) ){
		momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
		for (size_t i=0;i< indices.size();i++){
			sum+=mom(indices[i]);
		}

		res= insert(sum);
		put_label(key,res);
	}
#endif
	return res;
#else
	momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
	for (size_t i=0;i< v.size();i++){
		sum+=mom(v[i].ind());
	}
	if (v.size() == 1){
		return insert(sum,_mt_unknown);
	} else {
		return insert(sum,_mt_massive);
	}
	#endif
}
//! sum of momenta
/**
The sum of the momenta in the two vectors is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param v1 vector of plabel of momenta
\param v2 vector of plabel of momenta
\return integer index in the momentum_configuration<T> of the sum.
*/
template <class T> int momentum_configuration<T>::Sum(const std::vector<plabel>& v1,const std::vector<plabel>& v2){
#if _WITH_CACHING
	vector<int> indices;

	for (size_t i=0;i<v1.size();i++){indices.push_back(v1[i].ind());}
	for (size_t i=0;i<v2.size();i++){indices.push_back(v2[i].ind());}
	sort (indices.begin(),indices.end());
	string key=GenKey("Sum",indices);
	size_t res;

#if _CACHE_STATISTICS
	momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
//		Cmom<T> sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
	for (size_t i=0;i< indices.size();i++){
		sum+=mom(indices[i]);
	}

	find_and_stat(string("Sum_2plabel"),  key,this ,sum, res );
#else

	if ( !get_label(key,res) ){
		momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
//		Cmom<T> sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
		for (size_t i=0;i< indices.size();i++){
			sum+=mom(indices[i]);
		}

		res= insert(sum);
		put_label(key,res);
	}
#endif
	return res;
#else
	momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
	for (size_t i=0;i< v1.size();i++){
			sum+=mom(v1[i].ind());
		}
	for (size_t i=0;i< v2.size();i++){
			sum+=mom(v2[i].ind());
		}

	if (v1.size() + v2.size() == 1){
		return insert(sum,_mt_unknown);
	} else {
		return insert(sum,_mt_massive);
	}
	#endif
}

//! sum of momenta
/**
The sum of the momenta in the vector is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param v vector of integer labels of momenta
\return integer index in the momentum_configuration<T> of the sum.
*/
template <class T> int momentum_configuration<T>::Sum(const std::vector<int>& vi){
	momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
#if _WITH_CACHING
	vector<int> v(vi);
	sort (v.begin(),v.end());
	std::string key=GenKey1("Sum",v);
	size_t res;


#if _CACHE_STATISTICS

	for (size_t i=0;i< v.size();i++){
				sum+=mom(v[i]);
			}
	find_and_stat(string("Sumvect"),  key,this ,sum, res );
#else

	if ( !get_label(key,res) ){

		for (size_t i=0;i< v.size();i++){
			sum+=mom(v[i]);
		}
		res=insert(sum);
		put_label(key,res);
	}

#endif
	return res;
#else
		for (size_t i=0;i< vi.size();i++){
			sum+=mom(vi[i]);
		}
		if (vi.size() == 1){
			return insert(sum,_mt_unknown);
		} else {
			return insert(sum,_mt_massive);
		}
#endif
}
//! sum of momenta
/**
The sum of the momenta in the two vectors is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param v1 vector of integer label of momenta
\param v2 vector of integer label of momenta
\return integer index in the momentum_configuration<T> of the sum.
*/
template <class T> int momentum_configuration<T>::Sum(const std::vector<int>& v1,const std::vector<int>& v2){
	momentum<complex<T> > sum(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
#if _WITH_CACHING
	vector<int> v=v1;
	for (size_t i=0;i<v2.size();i++){v.push_back(v2[i]);}
	sort (v.begin(),v.end());
	string key=GenKey("Sum",v);
	size_t res;

#if _CACHE_STATISTICS

	for (size_t i=0;i< v.size();i++){
				sum+=mom(v[i]);
			}
	find_and_stat(string("Sum2vect"),  key,this ,sum, res );
#else

	if ( !get_label(key,res) ){
		for (size_t i=0;i< v.size();i++){
			sum+=mom(v[i]);
		}
		res=insert(sum);
		put_label(key,res);
	}
#endif
	return res;

#else
	for (size_t i=0;i< v1.size();i++){
		sum+=mom(v1[i]);
	}
	for (size_t i=0;i< v2.size();i++){
		sum+=mom(v2[i]);
	}
	if (v1.size() +v2.size() == 1){
		return insert(sum,_mt_unknown);
	} else {
		return insert(sum,_mt_massive);
	}
#endif
}
//! sum of momenta
/**
The sum of the momenta i1 and i2 is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param i1 label of first momentum
\param i2 label of second momenta
\return integer index in the momentum_configuration<T> of the sum.
*/

template <class T> int momentum_configuration<T>::Sum(int i1,int i2){
#if _WITH_CACHING==1
	string key;
	if (i1<i2){
		key=GenKey1("Sum",i1,i2);
	}
	else{
		key=GenKey1("Sum",i2,i1);
	}
	size_t res;
#if _CACHE_STATISTICS
	find_and_stat(string("Sum2"),  key,this ,mom(i1)+mom(i2), res );
#else
	if ( !get_label(key,res) ){
		res=insert(mom(i1)+mom(i2));
		put_label(key,res);
	}
#endif
	return res;
#else
	return insert(mom(i1)+mom(i2),_mt_massive);
	#endif
}
//! sum of momenta
/**
The sum of the momenta with the specified integer labels is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param i1 label of first momentum
\param i2 label of second momenta
\param i3 label of third momenta
\return integer index in the momentum_configuration<T> of the sum.
*/

template <class T> int momentum_configuration<T>::Sum(int i1,int i2, int i3){
#if _WITH_CACHING ==1
	vector<int> indices(3);

	indices[0]=i1;
	indices[1]=i2;
	indices[2]=i3;

	sort (indices.begin(),indices.end());
	string key=GenKey1("Sum",i1,i2,i3);
	size_t res;

#if _CACHE_STATISTICS
	find_and_stat(string("Sum3"),  key,this ,mom(i1)+mom(i2)+mom(i3), res );
#else

	if ( !get_label(key,res) ){
		res=insert(mom(i1)+mom(i2)+mom(i3));
		put_label(key,res);

	}
#endif
	return res;
#else
	return insert(mom(i1)+mom(i2)+mom(i3),_mt_massive);
#endif
}
//! sum of momenta
/**
The sum of the momenta with the specified integer labels is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param i1 label of first momentum
\param i2 label of second momenta
\param i3 label of third momenta
\param i4 label of fourth momenta
\return integer index in the momentum_configuration<T> of the sum.
*/
template <class T> int momentum_configuration<T>::Sum(int i1,int i2, int i3, int i4){
#if _WITH_CACHING
	vector<int> indices(4);

	indices[0]=i1;
	indices[1]=i2;
	indices[2]=i3;
	indices[3]=i4;

	sort (indices.begin(),indices.end());
	string key=GenKey1("Sum",i1,i2,i3,i4);
	size_t res;


#if _CACHE_STATISTICS
	find_and_stat(string("Sum4"),  key,this ,mom(i1)+mom(i2)+mom(i3), res );
#else

	if ( !get_label(key,res) ){
		res=insert(mom(i1)+mom(i2)+mom(i3)+mom(i4));
		put_label(key,res);
	}
#endif
	return res;
#else
	return insert(mom(i1)+mom(i2)+mom(i3)+mom(i4),_mt_massive);
#endif
	}
//! sum of momenta
/**
The sum of the momenta with the specified integer labels is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param i1 label of first momentum
\param i2 label of second momenta
\param i3 label of third momenta
\param i4 label of fourth momenta
\param i5 label of fourth momenta
\return integer index in the momentum_configuration<T> of the sum.
*/

template <class T> int momentum_configuration<T>::Sum(int i1,int i2, int i3, int i4, int i5){
#if _WITH_CACHING
	vector<int> indices(5);

	indices[0]=i1;
	indices[1]=i2;
	indices[2]=i3;
	indices[3]=i4;
	indices[3]=i5;

	sort (indices.begin(),indices.end());
	string key=GenKey("Sum",i1,i2,i3,i4,i5);
	size_t res;
#if _CACHE_STATISTICS
	find_and_stat(string("Sum5"),  key,this ,mom(i1)+mom(i2)+mom(i3), res );
#else

	if ( !get_label(key,res) ){
		res=insert(mom(i1)+mom(i2)+mom(i3)+mom(i4)+mom(i5));
		put_label(key,res);
	}
#endif
	return res;
#else
	return insert(mom(i1)+mom(i2)+mom(i3)+mom(i4)+mom(i5),_mt_massive);
#endif
}

//! sum of momenta
/**
The sum of the momenta with the specified integer labels is build and saved in the momentum configuration. the function returns the index of this sum of vectors. The sum is only computed once and cached for later use.
\param i1 label of first momentum
\param i2 label of second momenta
\param i3 label of third momenta
\param i4 label of fourth momenta
\param i5 label of fourth momenta
\param i6 label of fourth momenta
\return integer index in the momentum_configuration<T> of the sum.
*/
template <class T> int momentum_configuration<T>::Sum(int i1,int i2, int i3, int i4, int i5,int i6){
#if _WITH_CACHING
	vector<int> indices(6);

	indices[0]=i1;
	indices[1]=i2;
	indices[2]=i3;
	indices[3]=i4;
	indices[3]=i5;
	indices[3]=i6;

	sort (indices.begin(),indices.end());
	string key=GenKey("Sum",i1,i2,i3,i4,i5,i6);
	size_t res;
#if _CACHE_STATISTICS
	find_and_stat(string("Sum6"),  key,this ,mom(i1)+mom(i2)+mom(i3), res );
#else

	LabelCache::iterator pos=labelscache.find(key);
	if (pos != labelscache.end()) {
		res=pos->second;
	}
	else {
		res=insert(mom(i1)+mom(i2)+mom(i3)+mom(i4)+mom(i5)+mom(i6));
		labelscache[key]=res;
	}
#endif
	return res;
#else
	return insert(mom(i1)+mom(i2)+mom(i3)+mom(i4)+mom(i5)+mom(i6),_mt_massive);
#endif
	}

template <class T> mom_conf_reader<T>::mom_conf_reader(const char* path,size_t nbr_p):  BH::momentum_configuration<T>(0), nth(0), nbr_particles(nbr_p)  {
	input.open(path);
	if (! input) {
		std::string errorStr="No file ";
		errorStr+=path;
		errorStr+=" for the constructor mc_reader::mc_reader.";
		throw BH::BHerror(errorStr.c_str());
	}

}

template <class T> bool mom_conf_reader<T>::go_to_pos(std::ios::pos_type pos,size_t n) {
	if (nth==n) return true;
	input.seekg(pos);
	nth=n-1;
	return next();
}


template <class T> bool mom_conf_reader<T>::go_to(size_t n) {
	T  dummy;
	if (nth==n) return true;
	if (n>nth) {
		for (size_t i=1;i<n-nth;i++){
			for (size_t j=1;j<=nbr_particles;j++){
				if (!(input >> dummy && input >> dummy && input >> dummy && input >> dummy)) return false;
			}
		}
		nth=n-1;
	return next();
	}
	else if (n<nth) {
		input.seekg(0,std::ios::beg);
		for (size_t i=1;i<n;i++){
			if (!(input >> dummy && input >> dummy && input >> dummy && input >> dummy)) return false;
		}
		nth=n-1;

	}
	return next();
}

template <class T> bool mom_conf_reader<T>::next() {T  e,x,y,z;
	start_pos=input.tellg();
	momentum_configuration<T>::clear();
	momentum_configuration_base::renew_ID();
	for (size_t j=1;j<=nbr_particles;j++){
		if (!(input >> e && input >> x && input >> y && input >> z)) return false;
		this->insert(BH::Cmom<T>(complex<T>(e,0.),complex<T>(x,0.),complex<T>(y,0.),complex<T>(z,0.)),_mt_massless);
	}
	nth++;
	return true;
}

template <class T> void momentum_configuration<T>::clear(){
	ps.clear();
	ms.clear();
#if _WITH_CACHING
	cache.clear();
	labelscache.clear();
	Mspa.clear();
	Mspb.clear();
#endif /* _WITH_CACHING */
	nbr=0;
	renew_ID();
}

//! gets a value from the cache
/** \param key string representing the key \param res complex variable passed by reference set to the value cached if it exists in the cache. \return true if the value correspoonding to the key is in the cache, false if not.  The value of res is set to the cached value if true is returned. The typical call would be
 *  \n T res; if ( !get_value(the_key,res) ) {... compute result ....; res= the_result; mc.put_value(the_key,res)} ... continue using res */
template <class T> bool momentum_configuration<T>::get_value(const std::string& key,complex<T> &res){
#if _USE_GCC
	typename __gnu_cxx::hash_map<std::string, std::complex<T>, __gnu_cxx::hash<std::string> >::iterator pos=cache.find(key);
#endif
#if _USE_SUN_CC || _USE_PGCC
	typename std::hash_map<string, complex<T> >::iterator pos=cache.find(key);
#endif
	if (pos != cache.end()) {
		res=pos->second;
		return true;
	}
	else return false;
}
//! gets an integer label from the label cache
/** \param key string representing the key \param res integer variable passed by reference set to the label cached if it exists in the cache. \return true if the integer label corresponding to the key is in the cache, false if not.  The value of res is set to the cached value if true is returned. The typical call would be
 *  \n T res; if ( !get_label(the_key,res) ) { compute result .... res= the_result; mc.put_label(the_key,res)} ... continue using res */
template <class T> bool momentum_configuration<T>::get_label(const std::string& key,size_t &res){
	LabelCache::iterator pos=labelscache.find(key);
	if (pos != labelscache.end()) {
		res=pos->second;
		return true;
	}
	else return false;
}

template <class T> std::vector<int> momentum_configuration<T>::shiftBA(const std::vector<int>& ind, size_t b, size_t a,const complex<T>& z){
	size_t a_pos,b_pos;
	for (size_t i=0;i<ind.size();i++){
		if (ind[i]==a) a_pos=i;
		if (ind[i]==b) b_pos=i;
	}
	std::vector<int> res=ind;
	res[b_pos]=insert(L(b),Lt(b)-z*Lt(a));
	res[a_pos]=insert(L(a)+z*L(b),Lt(a));
	return res;
}

//template <class T> const momentum<complex<T> >& sub_momentum_configuration<T>::p(size_t n) const {
//	if (n<=momentum_configuration<T>::nbr){
//		if ( n>_offset ) {
//
//			return momentum_configuration<T>::ps[n-_offset-1].P();
//		}
//		else
//		{
//			return _parent->p(n);
//		}
//	}
//	else{_WARNING5("Too large momentum index in sub_momentum_configuration::p: ",n," (max=",momentum_configuration<T>::nbr,")"  );throw; }
//}


template <class T> sub_momentum_configuration<T>::sub_momentum_configuration(momentum_configuration<T>& mc) : momentum_configuration<T>(1000) {
	momentum_configuration<T>::nbr=mc.n();
	momentum_configuration<T>::_offset=(mc.n());
	momentum_configuration<T>::_parent=&mc;
}


template <class T> bool sub_momentum_configuration<T>::get_value(const std::string& key,std::complex<T> &res){
#if _USE_GCC
	typename __gnu_cxx::hash_map<std::string, std::complex<T>, __gnu_cxx::hash<std::string> >::iterator pos=momentum_configuration<T>::cache.find(key);
#endif
#if _USE_SUN_CC || _USE_PGCC
	typename std::hash_map<string, complex<T> >::iterator pos=momentum_configuration<T>::cache.find(key);
#endif
	if (pos != momentum_configuration<T>::cache.end()) {
			res=pos->second;
			return true;
		}
	if (momentum_configuration<T>::_parent->get_value(key,res) ) {
		return true;
		}

	else return false;
}

template <class T> bool sub_momentum_configuration<T>::get_label(const std::string& key,size_t &res){
	LabelCache::iterator pos=momentum_configuration<T>::labelscache.find(key);
	if (pos != momentum_configuration<T>::labelscache.end()) {
		res=pos->second;
		return true;
	}
	if (momentum_configuration<T>::_parent->get_label(key,res) ) {
		if (res>momentum_configuration<T>::_offset) {
//			_MESSAGE("vector available in cache of parent, but its definition came after branching to current sub_momentum_configuration. Duplicating vector in current sub_momentum_configuration (no reason to panic, you just spoiled cpu time...).");
					return false;
		}
//		_PRINT("using stored value in parent mc");
		return true;
		}
	else return false;
}

template<class Te,class T0> momentum_configuration<Te>
  extend(momentum_configuration<T0>& k,
         const std::vector<int>& indices /* of the momenta within "k" */,
         const Te& dummy /* To signal target type */);




template <class T> int momentum_configuration<T>::insert(const Cmom<T>& m){
	#if BH_USE_OMP
	get_lock();
	#endif
	ps.push_back(m);if (m.type() == _mt_massless) { ms.push_back(complex<T>(T(0),T(0)));} else { ms.push_back(m*m);}
	int new_position = ++nbr;
	#if BH_USE_OMP
	#pragma omp flush(nbr)
	release_lock();
	#endif
	return new_position;
}
template <class T> int momentum_configuration<T>::insert(const Cmom<T>& m,momentum_type ty){
	#if BH_USE_OMP
	get_lock();
	#endif
	ps.push_back(m);if (ty == _mt_massless) { ms.push_back(complex<T>(T(0),T(0)));} else { ms.push_back(m*m);}
	int new_position = ++nbr;
	#if BH_USE_OMP
	#pragma omp flush(nbr)
	release_lock();
	#endif
	return new_position;
}


#if _CACHE_STATISTICS

template <class T> void momentum_configuration<T>::cache_statistics(){
	_MESSAGE("------------------------------------------------------------");
	_MESSAGE("CACHING STATISTICS");
	_MESSAGE("------------------------------------------------------------");
	for (typename std::map<std::string,hash_stat>::iterator pos=momentum_configuration<T>::hash_stat_table.begin(); pos != momentum_configuration<T>::hash_stat_table.end();pos++){
		_MESSAGE8(pos->first," ] hit rate: ",pos->second.hit_rate(), " reuse rate: ",pos->second.reuse_rate(), " ( total: ",pos->second.nbr_try, ")");


}
	_MESSAGE("------------------------------------------------------------");
}

template <class key_type, class T> void find_and_stat(const std::string& tag, key_type& key,std::map<key_type,std::complex<T> >& value_map, complex<T> insert_value,complex<T> return_value, complex<T>& res ) {

	typename std::map<int,std::complex<T> >::iterator pos=value_map.find(key);
	typename std::map<std::string,typename momentum_configuration<T>::hash_stat>::iterator pos_map;

	if (pos_map == momentum_configuration<T>::hash_stat_table.end()){
		momentum_configuration<T>::hash_stat_table.insert(std::pair<std::string,typename momentum_configuration<T>::hash_stat>(tag,typename momentum_configuration<T>::hash_stat()));
		pos_map=momentum_configuration<T>::hash_stat_table.find(tag);
	}

	pos_map=momentum_configuration<T>::hash_stat_table.find(tag);
	if (pos != value_map.end()){
		res =  return_value ;
	}
	else {
//		_MESSAGE5(key,": return_value ", return_value,"  insert_value ", insert_value);
 		res=return_value;
		value_map.insert(pair<key_type,complex<T> >(key,insert_value));
		pos_map->second.nbr_no_hit++;
	}
	pos_map->second.nbr_try++;
}

template <class key_type, class T> void find_and_stat(const std::string& tag, key_type& key,hash_map<std::string, complex<T>, hash<std::string> >& value_map, complex<T> insert_value,complex<T> return_value, complex<T>& res ) {


	typename std::map<std::string,typename momentum_configuration<T>::hash_stat>::iterator pos_map;
	pos_map=momentum_configuration<T>::hash_stat_table.find(tag);

	if (pos_map == momentum_configuration<T>::hash_stat_table.end()){
		momentum_configuration<T>::hash_stat_table.insert(std::pair<std::string,typename momentum_configuration<T>::hash_stat>(tag,typename momentum_configuration<T>::hash_stat()));
		pos_map=momentum_configuration<T>::hash_stat_table.find(tag);
	}

	typename hash_map<std::string, complex<T>, hash<std::string> >::iterator pos=value_map.find(key);
	if (pos != value_map.end()){

		res =  return_value;
	}
	else {
 		res=return_value;
		value_map.insert(pair<key_type,complex<T> >(key,insert_value));
		pos_map->second.nbr_no_hit++;

	}
	pos_map->second.nbr_try++;
}

template <class key_type, class T> void find_and_stat(const std::string& tag, key_type& key,momentum_configuration<T>* mc , complex<T> insert_value,complex<T> return_value, complex<T>& res ) {


	typename std::map<std::string,typename momentum_configuration<T>::hash_stat>::iterator pos_map;
	pos_map=momentum_configuration<T>::hash_stat_table.find(tag);

	if (pos_map == momentum_configuration<T>::hash_stat_table.end()){
		momentum_configuration<T>::hash_stat_table.insert(std::pair<std::string,typename momentum_configuration<T>::hash_stat>(tag,typename momentum_configuration<T>::hash_stat()));
		pos_map=momentum_configuration<T>::hash_stat_table.find(tag);
	}

	if ( ! mc->get_value(key,res) ){
		mc->put_value(key,insert_value);
		pos_map->second.nbr_no_hit++;
	}
	pos_map->second.nbr_try++;
	res=return_value;
}



#endif

}

#endif /* MOM_CONF_HPP_ */
