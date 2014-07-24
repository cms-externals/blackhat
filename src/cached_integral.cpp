/*
 * cached_integral.cpp
 *
 *  Created on: 28-Nov-2008
 *      Author: daniel
 */

#include "cached_integral.h"
#include <algorithm>
#include <set>
#include <functional>
#include "BH_typedefs.h"
#include "helcode.h"
#include "integrals.h"
#include "index_vector.h"

using namespace std;

namespace BH {


namespace CachedIntegral {




Cached_Bubble_Integral::Cached_Bubble_Integral(const std::vector<int>& v1,const std::vector<int>& v2) : d_index_vector_1(v1.begin(),v1.end()), d_index_vector_2(v2.begin(),v2.end())  {
	size_t n=v1.size()+v2.size();
	d_scode=get_invariant_code(v1,n);
	d_nbr_indices=n;
}


const SeriesC<R>* Cached_Integral_impl::get() const {
	return &d_value;
}
const SeriesC<RHP>* Cached_Integral_impl::get_HP() const {
	return &d_value_HP;
}
const SeriesC<RVHP>* Cached_Integral_impl::get_VHP() const {
	return &d_value_VHP;
}
#if BH_USE_GMP
const SeriesC<RGMP>* Cached_Integral_impl::get_GMP() const {
	return &d_value_GMP;
}
#endif




void Cached_Bubble_Integral::prepare(mom_conf& mc,int mu){
	set_mcID( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2));
}
void Cached_Bubble_Integral::prepare(mom_conf_HP& mc,int mu){
	set_mcID_HP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2));
}
void Cached_Bubble_Integral::prepare(mom_conf_VHP& mc,int mu){
	set_mcID_VHP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2));
}

#if BH_USE_GMP
void Cached_Bubble_Integral::prepare(momentum_configuration<RGMP>& mc,int mu){
	set_mcID_GMP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2));
}
#endif


std::ostream& operator<<(std::ostream& s,const Cached_Bubble_Integral& cb){
	s << "I2(" << cb.d_scode  << ") used " << cb.get_load() << " times. Last value:" <<  *cb.get(); return s << endl;
}

bool operator<(const Cached_Bubble_Integral& c1,const Cached_Bubble_Integral& c2){
	return c1.d_scode < c2.d_scode;
}

//const SeriesC<R>* Cached_Triangle_Integral::get() const {
//	return &d_value;
//}
//const SeriesC<RHP>* Cached_Triangle_Integral::get_HP() const {
//	return &d_value_HP;
//}
//const SeriesC<RVHP>* Cached_Triangle_Integral::get_VHP() const {
//	return &d_value_VHP;
//}

void Cached_Triangle_Integral::prepare(mom_conf& mc,int mu){
	set_mcID( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3));
}
void Cached_Triangle_Integral::prepare(mom_conf_HP& mc,int mu){
	set_mcID_HP( mc.get_ID());
	set_value( Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3));
}
void Cached_Triangle_Integral::prepare(mom_conf_VHP& mc,int mu){
	set_mcID_VHP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3));
}

#if BH_USE_GMP
void Cached_Triangle_Integral::prepare(momentum_configuration<RGMP>& mc,int mu){
	set_mcID_GMP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3));
}
#endif


Cached_Triangle_Integral::Cached_Triangle_Integral(const std::vector<int>& v1,const std::vector<int>& v2,const std::vector<int>& v3) : d_index_vector_1(v1.begin(),v1.end()), d_index_vector_2(v2.begin(),v2.end()) , d_index_vector_3(v3.begin(),v3.end())  {
	size_t n=v1.size()+v2.size()+v3.size();
	d_scode1=get_invariant_code(v1,n);
	d_scode2=get_invariant_code(v2,n);
	d_scode3=get_invariant_code(v3,n);
	d_nbr_indices=n;
}

std::ostream& operator<<(std::ostream& s,const Cached_Triangle_Integral& ct){
	s << "I3(" << ct.d_scode1 << "," <<ct.d_scode2 << "," << ct.d_scode3 <<  ") used " << ct.get_load() << " times" ; return s << endl;
}

bool operator<(const Cached_Triangle_Integral& c1,const Cached_Triangle_Integral& c2){
	if (c1.d_scode1 < c2.d_scode1) return true;
	if (c1.d_scode2 < c2.d_scode2) return true;
	if (c1.d_scode3 < c2.d_scode3) return true;
	return false;
}
//
//const SeriesC<R>* Cached_Box_Integral::get() const {
//	return &d_value;
//}
//const SeriesC<RHP>* Cached_Box_Integral::get_HP() const {
//	return &d_value_HP;
//}
//const SeriesC<RVHP>* Cached_Box_Integral::get_VHP() const {
//	return &d_value_VHP;
//}

void Cached_Box_Integral::prepare(mom_conf& mc,int mu){
	set_mcID( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3,d_index_vector_4));
}
void Cached_Box_Integral::prepare(mom_conf_HP& mc,int mu){
	set_mcID_HP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3,d_index_vector_4));
}
void Cached_Box_Integral::prepare(mom_conf_VHP& mc,int mu){
	set_mcID_VHP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3,d_index_vector_4));
}

#if BH_USE_GMP
void Cached_Box_Integral::prepare(momentum_configuration<RGMP>& mc,int mu){
	set_mcID_GMP( mc.get_ID());
	set_value(Int(mc, mu,d_index_vector_1,d_index_vector_2,d_index_vector_3,d_index_vector_4));
}
#endif

Cached_Box_Integral::Cached_Box_Integral(const std::vector<int>& v1,const std::vector<int>& v2,const std::vector<int>& v3,const std::vector<int>& v4) :  d_index_vector_1(v1.begin(),v1.end()), d_index_vector_2(v2.begin(),v2.end()) , d_index_vector_3(v3.begin(),v3.end()), d_index_vector_4(v4.begin(),v4.end())  {
	size_t n=v1.size()+v2.size()+v3.size()+v4.size();
	d_scode1=get_invariant_code(v1,n);
	d_scode2=get_invariant_code(v2,n);
	d_scode3=get_invariant_code(v3,n);
	d_scode4=get_invariant_code(v4,n);
	d_nbr_indices=n;
}

std::ostream& operator<<(std::ostream& s,const Cached_Box_Integral& ct){
	s << "I4(" << ct.d_scode1 << "," <<ct.d_scode2 << "," << ct.d_scode3 << "," << ct.d_scode4 <<  ") used " << ct.get_load() << " times" ; return s << endl;
}

bool operator<(const Cached_Box_Integral& c1,const Cached_Box_Integral& c2){
	if (c1.d_scode1 < c2.d_scode1) return true;
	if (c1.d_scode2 < c2.d_scode2) return true;
	if (c1.d_scode3 < c2.d_scode3) return true;
	if (c1.d_scode4 < c2.d_scode4) return true;
	return false;
}




void Cached_Integral_Factory::prepare(mom_conf& mc,int mu){
	for (std::vector<Cached_Bubble_Integral*>::iterator it=d_bubble_integrals.begin();it != d_bubble_integrals.end();++it){
		(*it)->prepare(mc,mu);
	}
	for (std::vector<Cached_Triangle_Integral*>::iterator it=d_triangle_integrals.begin();it != d_triangle_integrals.end();++it){
		(*it)->prepare(mc,mu);
	}
	for (std::vector<Cached_Box_Integral*>::iterator it=d_box_integrals.begin();it != d_box_integrals.end();++it){
		(*it)->prepare(mc,mu);
	}
}

Cached_Integral_Factory Cached_Integral_Factory::s_default_CIF;


Cached_Integral_Factory::~Cached_Integral_Factory(){
	for_each(d_bubble_integrals.begin(),d_bubble_integrals.end(),do_delete<Cached_Bubble_Integral>());
	for_each(d_triangle_integrals.begin(),d_triangle_integrals.end(),do_delete<Cached_Triangle_Integral>());
	for_each(d_box_integrals.begin(),d_box_integrals.end(),do_delete<Cached_Box_Integral>());
}

void Cached_Integral_Factory::print_state(){
	std::cout << "=-=-=-=-=-=-=-=-=-=-= integral_factory =-=-=-=-=-=-=-=-=-=-= " << std::endl;
	size_t nbr_bubbles=d_bubble_integrals.size();
	cout << nbr_bubbles << " bubbles: " << endl;
	for (int j=0;j< nbr_bubbles;j++){

		std::cout << *d_bubble_integrals[j] ;
	}
	std::cout << std::endl;
	size_t nbr_triangles=d_triangle_integrals.size();
	cout << nbr_triangles << " triangles: " << endl;
	for (int j=0;j< nbr_triangles;j++){

		std::cout << *d_triangle_integrals[j] ;
	}
	std::cout << std::endl;
	size_t nbr_boxes=d_box_integrals.size();
	cout << nbr_boxes << " boxes: " << endl;
	for (int j=0;j< nbr_boxes;j++){

		std::cout << *d_box_integrals[j] ;
	}
	std::cout << std::endl;
	std::cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= " << std::endl;

}

struct matching_bubble_integral : public std::unary_function<Cached_Bubble_Integral*,bool> {
	long d_scode;
	size_t d_nbr;
	matching_bubble_integral(long code,size_t n): d_scode(code),d_nbr(n) {};
	bool operator()(Cached_Bubble_Integral* cbi){ return ( (cbi->get_scode() == d_scode) && ( cbi->get_nbr() == d_nbr ) ); };
};

struct matching_triangle_integral : public std::unary_function<Cached_Triangle_Integral*,bool> {
	long d_scode1;
	long d_scode2;
	long d_scode3;
	size_t d_nbr;
	matching_triangle_integral(long code1,long code2,long code3,size_t n): d_scode1(code1),d_scode2(code2),d_scode3(code3),d_nbr(n) {};
	bool operator()(Cached_Triangle_Integral* cbi){ return ( (cbi->get_scode(1) == d_scode1) &&(cbi->get_scode(2) == d_scode2) &&(cbi->get_scode(3) == d_scode3) && ( cbi->get_nbr() == d_nbr ) ); };
};

struct matching_box_integral : public std::unary_function<Cached_Box_Integral*,bool> {
	long d_scode1;
	long d_scode2;
	long d_scode3;
	long d_scode4;
	size_t d_nbr;
	matching_box_integral(long code1,long code2,long code3,long code4,size_t n): d_scode1(code1),d_scode2(code2),d_scode3(code3),d_scode4(code4),d_nbr(n) {};
	bool operator()(Cached_Box_Integral* cbi){ return ( (cbi->get_scode(1) == d_scode1) &&(cbi->get_scode(2) == d_scode2) &&(cbi->get_scode(3) == d_scode3) && (cbi->get_scode(4) == d_scode4) && ( cbi->get_nbr() == d_nbr ) ); };
};



Cached_Bubble_Integral* Cached_Integral_Factory::new_integral(const vector<int>& v1,const vector<int>& v2){
	vector<Cached_Bubble_Integral*>::iterator it;
	size_t n=v1.size()+v2.size();
	it=find_if(d_bubble_integrals.begin(),d_bubble_integrals.end(),matching_bubble_integral(get_invariant_code(v1,n),n));
	if ( it != d_bubble_integrals.end() ){
		(*it)->increase_load();
		return *it;
	} else {
		d_bubble_integrals.push_back(new Cached_Bubble_Integral(v1,v2));
		d_bubble_integrals.back()->increase_load();
		return d_bubble_integrals.back();
	}
}
/** The triangles are ordered using their invariants. Smaller invariants come first. */
Cached_Triangle_Integral* Cached_Integral_Factory::new_integral(const vector<int>& v1,const vector<int>& v2,const vector<int>& v3){
	vector<Cached_Triangle_Integral*>::iterator it;
	size_t n=v1.size()+v2.size()+v3.size();
	long s1=get_invariant_code(v1,n);
	long s2=get_invariant_code(v2,n);
	long s3=get_invariant_code(v3,n);
	long s[]={s1,s2,s3};
	sort(s,s+3);
	it=find_if(d_triangle_integrals.begin(),d_triangle_integrals.end(),matching_triangle_integral(s[0],s[1],s[2],n));
	if ( it != d_triangle_integrals.end() ){
		(*it)->increase_load();
		return *it;
	} else {
		size_t pos1=find(s,s+3,s1)-s;
		size_t pos2=find(s,s+3,s2)-s;
		size_t pos3=find(s,s+3,s3)-s;
		const vector<int>* vs[3];
		vs[pos1]=&v1;
		vs[pos2]=&v2;
		vs[pos3]=&v3;
		d_triangle_integrals.push_back(new Cached_Triangle_Integral(*vs[0],*vs[1],*vs[2]));
		return d_triangle_integrals.back();
	}
}

/** The boxes are ordered using their invariants. The procedure is the following. The first invariant S[0] in the list is the smallest one. Then we put all the other invariants
 * starting from there, in the direction determined by the following procedure. If the element left (cyclically) of S[0] is smaller than the invariant right of S[0], then we order the invariants from the right to the left
 *  starting with S[0], otherwise (also when both are equal) we order them from the left to the right, starting with s1.*/
Cached_Box_Integral* Cached_Integral_Factory::new_integral(const vector<int>& v1,const vector<int>& v2,const vector<int>& v3,const vector<int>& v4){
	size_t n=v1.size()+v2.size()+v3.size()+v4.size();
	long s1=get_invariant_code(v1,n);
	long s2=get_invariant_code(v2,n);
	long s3=get_invariant_code(v3,n);
	long s4=get_invariant_code(v4,n);
	long S[4];
	long s[]={s1,s2,s3,s4};
	size_t min_pos=min_element(s,s+4)-s;
	S[0]=s[min_pos];
	if (s[(min_pos+3)%4]<s[(min_pos+1)%4]) {
		S[1]=s[(min_pos+3)%4];
		S[2]=s[(min_pos+2)%4];
		S[3]=s[(min_pos+1)%4];
	} else {
		S[1]=s[(min_pos+1)%4];
		S[2]=s[(min_pos+2)%4];
		S[3]=s[(min_pos+3)%4];
	}

	vector<Cached_Box_Integral*>::iterator it;
	it=find_if(d_box_integrals.begin(),d_box_integrals.end(),matching_box_integral(S[0],S[1],S[2],S[3],n));
	if ( it != d_box_integrals.end() ){
		(*it)->increase_load();
		return *it;
	} else {
		size_t pos1=find(S,S+4,s1)-S;
		size_t pos2=find(S,S+4,s2)-S;
		size_t pos3=find(S,S+4,s3)-S;
		size_t pos4=find(S,S+4,s4)-S;
		const vector<int>* vs[4];
		vs[pos1]=&v1;
		vs[pos2]=&v2;
		vs[pos3]=&v3;
		vs[pos4]=&v4;
		d_box_integrals.push_back(new Cached_Box_Integral(*vs[0],*vs[1],*vs[2],*vs[3]));
		return d_box_integrals.back();
	}
}


void Cached_Bubble_Integral_User::print_info(){
	_MESSAGE("-------------------------------------------------------------");
	_MESSAGE2("Invariant: s",get_invariant_code(d_corner1,d_corner1.size()+d_corner2.size()));
	_MESSAGE2("Nbr of cached permutations: ",d_values.size());
	for (iter it=d_values.begin();it!=d_values.end();++it){
		cout << it->first << ": " <<  *it->second << endl;
	}
	_MESSAGE("-------------------------------------------------------------");
}

Cached_Bubble_Integral_User::Cached_Bubble_Integral_User(const vector<int>& v1,const vector<int>& v2){
	d_corner1=v1;
	d_corner2=v2;
}

template <class T> const SeriesC<T>* Cached_Bubble_Integral_User::get_value_fn(momentum_configuration<T>& mc,const Index_Vector& IV,int mu){
	iter it=d_values.find(IV.get_permutation_code());
	if (it != d_values.end()){
		if ( (*it).second->get_mc_ID<T>() == mc.get_ID() ){
			return ((*it).second->get_value<T>());
		} else {
			(*it).second->prepare(mc,mu);
			return ((*it).second->get_value<T>());
		}
	} else {
		vector<int> v1,v2;
		for (int i=0;i<d_corner1.size();i++){
			v1.push_back(IV[d_corner1[i]-1]);
		}
		for (int i=0;i<d_corner2.size();i++){
			v2.push_back(IV[d_corner2[i]-1]);
		}
		Cached_Bubble_Integral* CBI=Cached_Integral_Factory::s_default_CIF.new_integral(v1,v2);
		CBI->prepare(mc,mu);
		it= (d_values.insert(std::make_pair(IV.get_permutation_code(),CBI))).first ;
		return ((*it).second->get_value<T>());
	}
}


template const SeriesC<R>* Cached_Bubble_Integral_User::get_value_fn(mom_conf& mc,const Index_Vector& IV,int mu);
template const SeriesC<RHP>* Cached_Bubble_Integral_User::get_value_fn(mom_conf_HP& mc,const Index_Vector& IV,int mu);
template const SeriesC<RVHP>* Cached_Bubble_Integral_User::get_value_fn(mom_conf_VHP& mc,const Index_Vector& IV,int mu);

#if BH_USE_GMP
template const SeriesC<RGMP>* Cached_Bubble_Integral_User::get_value_fn(momentum_configuration<RGMP>& mc,const Index_Vector& IV,int mu);
#endif

void Cached_Triangle_Integral_User::print_info(){
	size_t n=d_corner1.size()+d_corner2.size()+d_corner3.size();
	_MESSAGE("-------------------------------------------------------------");
	_MESSAGE6("Invariants: s",get_invariant_code(d_corner1,n),", s",get_invariant_code(d_corner2,n),", s",get_invariant_code(d_corner3,n));
	_MESSAGE2("Nbr of cached permutations: ",d_values.size());
	for (iter it=d_values.begin();it!=d_values.end();++it){
		cout << it->first << ": " <<  *it->second << endl;
	}
	_MESSAGE("-------------------------------------------------------------");
}

Cached_Triangle_Integral_User::Cached_Triangle_Integral_User(const vector<int>& v1,const vector<int>& v2,const vector<int>& v3) : d_corner1(v1.begin(),v1.end()) , d_corner2(v2.begin(),v2.end()), d_corner3(v3.begin(),v3.end()){
}




template <class T> const SeriesC<T>* Cached_Triangle_Integral_User::get_value_fn(momentum_configuration<T>& mc,const Index_Vector& IV,int mu){
	iter it=d_values.find(IV.get_permutation_code());
	if (it != d_values.end()){
		if ( (*it).second->get_mc_ID<T>() == mc.get_ID() ){
			return ((*it).second->get_value<T>());
		} else {
			(*it).second->prepare(mc,mu);
			return ((*it).second->get_value<T>());
		}
	} else {
		vector<int> v1,v2,v3;
		for (int i=0;i<d_corner1.size();i++){
			v1.push_back(IV[d_corner1[i]-1]);
		}
		for (int i=0;i<d_corner2.size();i++){
			v2.push_back(IV[d_corner2[i]-1]);
		}
		for (int i=0;i<d_corner3.size();i++){
			v3.push_back(IV[d_corner3[i]-1]);
		}
		Cached_Triangle_Integral* CBI=Cached_Integral_Factory::s_default_CIF.new_integral(v1,v2,v3);
		CBI->prepare(mc,mu);
		it= (d_values.insert(std::make_pair(IV.get_permutation_code(),CBI))).first ;
		return ((*it).second->get_value<T>());
	}
}

template const SeriesC<R>* Cached_Triangle_Integral_User::get_value_fn(mom_conf& mc,const Index_Vector& IV,int mu);
template const SeriesC<RHP>* Cached_Triangle_Integral_User::get_value_fn(mom_conf_HP& mc,const Index_Vector& IV,int mu);
template const SeriesC<RVHP>* Cached_Triangle_Integral_User::get_value_fn(mom_conf_VHP& mc,const Index_Vector& IV,int mu);
#if BH_USE_GMP
template const SeriesC<RGMP>* Cached_Triangle_Integral_User::get_value_fn(momentum_configuration<RGMP>& mc,const Index_Vector& IV,int mu);
#endif

void Cached_Box_Integral_User::print_info(){
	size_t n=d_corner1.size()+d_corner2.size()+d_corner3.size()+d_corner4.size();
	_MESSAGE("-------------------------------------------------------------");
	_MESSAGE8("Invariants: s",get_invariant_code(d_corner1,n),", s",get_invariant_code(d_corner2,n),", s",get_invariant_code(d_corner3,n),", s",get_invariant_code(d_corner3,n));
	_MESSAGE2("Nbr of cached permutations: ",d_values.size());
	for (iter it=d_values.begin();it!=d_values.end();++it){
		cout << it->first << ": " <<  *it->second << endl;
	}
	_MESSAGE("-------------------------------------------------------------");
}

Cached_Box_Integral_User::Cached_Box_Integral_User(const vector<int>& v1,const vector<int>& v2,const vector<int>& v3,const vector<int>& v4) : d_corner1(v1.begin(),v1.end()) , d_corner2(v2.begin(),v2.end()), d_corner3(v3.begin(),v3.end()), d_corner4(v4.begin(),v4.end()){
#if _VERBOSE
	cout << "Cached_Box_Integral_User constructed with:" << endl;
	cout << "corner 1:" ; copy(d_corner1.begin(),d_corner1.end(),ostream_iterator<int>(cout ," ")); cout << endl;
	cout << "corner 2:" ; copy(d_corner2.begin(),d_corner2.end(),ostream_iterator<int>(cout ," ")); cout << endl;
	cout << "corner 3:" ; copy(d_corner3.begin(),d_corner3.end(),ostream_iterator<int>(cout ," ")); cout << endl;
	cout << "corner 4:" ; copy(d_corner4.begin(),d_corner4.end(),ostream_iterator<int>(cout ," ")); cout << endl;
#endif
}

template <class T> const SeriesC<T>* Cached_Box_Integral_User::get_value_fn(momentum_configuration<T>& mc,const Index_Vector& IV,int mu){

	iter it=d_values.find(IV.get_permutation_code());
	if (it != d_values.end()){
		if ( (*it).second->get_mc_ID<T>() == mc.get_ID() ){
			return ((*it).second->get_value<T>());
		} else {
			(*it).second->prepare(mc,mu);
			return ((*it).second->get_value<T>());
		}
	} else {
		vector<int> v1,v2,v3,v4;
		for (int i=0;i<d_corner1.size();i++){
			v1.push_back(IV[d_corner1[i]-1]);
		}
		for (int i=0;i<d_corner2.size();i++){
			v2.push_back(IV[d_corner2[i]-1]);
		}
		for (int i=0;i<d_corner3.size();i++){
			v3.push_back(IV[d_corner3[i]-1]);
		}
		for (int i=0;i<d_corner4.size();i++){
			v4.push_back(IV[d_corner4[i]-1]);
		}
		Cached_Box_Integral* CBI=Cached_Integral_Factory::s_default_CIF.new_integral(v1,v2,v3,v4);
		CBI->prepare(mc,mu);
		it= (d_values.insert(std::make_pair(IV.get_permutation_code(),CBI))).first ;
		return ((*it).second->get_value<T>());
	}


}

template const SeriesC<R>* Cached_Box_Integral_User::get_value_fn(momentum_configuration<R>& mc,const Index_Vector& IV,int mu);
template const SeriesC<RHP>* Cached_Box_Integral_User::get_value_fn(momentum_configuration<RHP>& mc,const Index_Vector& IV,int mu);
template const SeriesC<RVHP>* Cached_Box_Integral_User::get_value_fn(momentum_configuration<RVHP>& mc,const Index_Vector& IV,int mu);
#if BH_USE_GMP
template const SeriesC<RGMP>* Cached_Box_Integral_User::get_value_fn(momentum_configuration<RGMP>& mc,const Index_Vector& IV,int mu);
#endif






Known_Cut_wCI::Known_Cut_wCI(const process& PRO,color_structure cs) : Cut_Part_base(PRO) {
	vector<int> ind; for (int i=1;i<=PRO.n();i++){ind.push_back(i);};
	if ( ! (d_CPwCI_p=CwCI_Ptr(PRO,cs,ind)) ) {_WARNING("Known_Cut_wCI construction failed" );}
}

Known_Cut_wCI::Known_Cut_wCI(const process& PRO,Cut_Part_wCI* CPp) : Cut_Part_base(PRO), d_CPwCI_p(CPp) {
}

Known_Cut_wCI::~Known_Cut_wCI(){
	delete d_CPwCI_p;
}

SeriesC<R> Known_Cut_wCI::eval(mom_conf& mc,const std::vector<int>& ind){
	if (d_mu_index==0){d_mu_index = DefineMu<R>(mc,m_MU);};
	return d_CPwCI_p->eval(mc,static_cast<const Index_Vector&>(ind),d_mu_index);
}
SeriesC<RHP> Known_Cut_wCI::eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind){
	if (d_mu_index_HP==0){d_mu_index_HP = DefineMu<RHP>(mc,m_MU);};
	return d_CPwCI_p->eval(mc,static_cast<const Index_Vector&>(ind),d_mu_index_HP);
//	_WARNING("RHP for Known_Cut_wCI::eval Not yet implemented");
}
SeriesC<RVHP> Known_Cut_wCI::eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){
	if (d_mu_index_VHP==0){d_mu_index_VHP = DefineMu<RVHP>(mc,m_MU);};
	return d_CPwCI_p->eval(mc,static_cast<const Index_Vector&>(ind),d_mu_index_VHP);
//	_WARNING("RVHP for Known_Cut_wCI::eval Not yet implemented");
}
#if BH_USE_GMP
SeriesC<RGMP> Known_Cut_wCI::eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){
	if (d_mu_index_GMP==0){d_mu_index_GMP = DefineMu<RGMP>(mc,m_MU);};
	return d_CPwCI_p->eval(mc,static_cast<const Index_Vector&>(ind),d_mu_index_GMP);
//	_WARNING("RVHP for Known_Cut_wCI::eval Not yet implemented");
}
#endif


SeriesC<R> Known_Cut_wCI::eval(const eval_param<R>& ep){
//	if (d_mu_index==0){d_mu_index = DefineMu<R>(mc,m_MU);};
//	return d_CPwCI_p->eval(mc,static_cast<const Index_Vector&>(ind),d_mu_index);
	_WARNING("R for Known_Cut_wCI::eval Not yet implemented");
}
SeriesC<RHP> Known_Cut_wCI::eval(const eval_param<RHP>& ep){
	_WARNING("RHP for Known_Cut_wCI::eval Not yet implemented");
}
SeriesC<RVHP> Known_Cut_wCI::eval(const eval_param<RVHP>& ep){
	_WARNING("RVHP for Known_Cut_wCI::eval Not yet implemented");
}


Known_Cut_wCI_offset::Known_Cut_wCI_offset(const process& PRO,color_structure cs,size_t offset) : Cut_Part_base(PRO) , d_offset(offset) {
	vector<particle_ID> pp;
	for (int i=1;i<=d_process.n();i++){ pp.push_back(d_process.p(( d_offset+(i-2))%d_process.n()+1)) ; }
	process new_pro(pp);
	vector<int> ind; for (int i=1;i<=d_process.n();i++){ ind.push_back( i ); }
//	vector<int> ind; for (int i=1;i<=pro.n();i++){ ind.push_back( (pro.n()-offset+i)%pro.n()+1 ); }
//	vector<int> ind; for (int i=1;i<=pro.n();i++){ ind.push_back( (d_offset+(i-2))%pro.n()+1) ; }
	if ( ! (d_CPwCI_p=CwCI_Ptr(new_pro,cs,ind)) ) {_WARNING("Known_Cut_wCI construction failed" );}
}
Known_Cut_wCI_offset::Known_Cut_wCI_offset(const process& pro,CachedIntegral::Cut_Part_wCI* CPp,size_t offset) : Cut_Part_base(pro) ,d_CPwCI_p(CPp), d_offset(offset) {
}

Known_Cut_wCI_offset::~Known_Cut_wCI_offset(){
	delete d_CPwCI_p;
}

SeriesC<R> Known_Cut_wCI_offset::eval(mom_conf& mc,const std::vector<int>& ind){
	if (d_mu_index==0){d_mu_index = DefineMu<R>(mc,m_MU);};
	vector<int> new_ind(ind.begin(),ind.end());
	rotate(new_ind.begin(),new_ind.begin()+(d_offset-1),new_ind.end());
	Index_Vector IV(new_ind);
	return d_CPwCI_p->eval(mc,IV,d_mu_index);
}
SeriesC<RHP> Known_Cut_wCI_offset::eval(momentum_configuration<RHP>& mc,const std::vector<int>& ind){
	_WARNING("RHP for Known_Cut_wCI::eval Not yet implemented");
}
SeriesC<RVHP> Known_Cut_wCI_offset::eval(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){
	_WARNING("RVHP for Known_Cut_wCI::eval Not yet implemented");
}
#if BH_USE_GMP
SeriesC<RGMP> Known_Cut_wCI_offset::eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){
	_WARNING("RGMP for Known_Cut_wCI::eval Not yet implemented");
}
#endif

SeriesC<R> Known_Cut_wCI_offset::eval(const eval_param<R>& ep){
//	if (d_mu_index==0){d_mu_index = DefineMu<R>(mc,m_MU);};
//	vector<int> new_ind(ind.begin(),ind.end());
//	rotate(new_ind.begin(),new_ind.begin()+(d_offset-1),new_ind.end());
//	Index_Vector IV(new_ind);
//	return d_CPwCI_p->eval(mc,IV,d_mu_index);
	_WARNING("R for Known_Cut_wCI::eval Not yet implemented");
}
SeriesC<RHP> Known_Cut_wCI_offset::eval(const eval_param<RHP>& ep){
	_WARNING("RHP for Known_Cut_wCI::eval Not yet implemented");
}
SeriesC<RVHP> Known_Cut_wCI_offset::eval(const eval_param<RVHP>& ep){
	_WARNING("RVHP for Known_Cut_wCI::eval Not yet implemented");
}

Cut_Part_wCI::~Cut_Part_wCI(){
	for_each(CI_users.begin(),CI_users.end(),do_delete<Cached_Integral_User>());
}
#ifdef USE_OLD_CUT_PART

Unknown_Cut_Part_wCI::~Unknown_Cut_Part_wCI(){
	for_each(d_box_integrals.begin(),d_box_integrals.end(),do_delete<CachedIntegral::Cached_Box_Integral_User>());
	for_each(d_triangle_integrals.begin(),d_triangle_integrals.end(),do_delete<CachedIntegral::Cached_Triangle_Integral_User>());
	for_each(d_bubble_integrals.begin(),d_bubble_integrals.end(),do_delete<CachedIntegral::Cached_Bubble_Integral_User>());
}

void Unknown_Cut_Part_wCI::construct_integral_table(){
	for (vector<bubbleD*>::iterator it=bubbles.begin();it!=bubbles.end();++it){
//		_PRINT(*(*it));
		vector<int> cor1,cor2;
		transform((*it)->c(1).begin(),(*it)->c(1).end(),back_inserter(cor1),std::mem_fun_ref(&plabel::ind));
		transform((*it)->c(2).begin(),(*it)->c(2).end(),back_inserter(cor2),mem_fun_ref(&plabel::ind));
		d_bubble_integrals.push_back(new Cached_Bubble_Integral_User(cor1,cor2));
	}
	for (vector<triangleD*>::iterator it=triangles.begin();it!=triangles.end();++it){
//		_PRINT(*(*it));
		vector<int> cor1,cor2,cor3;
		transform((*it)->c(1).begin(),(*it)->c(1).end(),back_inserter(cor1),mem_fun_ref(&plabel::ind));
		transform((*it)->c(2).begin(),(*it)->c(2).end(),back_inserter(cor2),mem_fun_ref(&plabel::ind));
		transform((*it)->c(3).begin(),(*it)->c(3).end(),back_inserter(cor3),mem_fun_ref(&plabel::ind));
		d_triangle_integrals.push_back(new Cached_Triangle_Integral_User(cor1,cor2,cor3));
	}
	for (vector<boxD*>::iterator it=boxes.begin();it!=boxes.end();++it){
//		_PRINT(*(*it));
		vector<int> cor1,cor2,cor3,cor4;
		transform((*it)->c(1).begin(),(*it)->c(1).end(),back_inserter(cor1),mem_fun_ref(&plabel::ind));
		transform((*it)->c(2).begin(),(*it)->c(2).end(),back_inserter(cor2),mem_fun_ref(&plabel::ind));
		transform((*it)->c(3).begin(),(*it)->c(3).end(),back_inserter(cor3),mem_fun_ref(&plabel::ind));
		transform((*it)->c(4).begin(),(*it)->c(4).end(),back_inserter(cor4),mem_fun_ref(&plabel::ind));
		d_box_integrals.push_back(new Cached_Box_Integral_User(cor1,cor2,cor3,cor4));
	}
}

template <class T> SeriesC<T> Unknown_Cut_Part_wCI::Unknown_Cut_Part_wCI_eval(momentum_configuration<T>& mc,const vector<int>& ind){
	complex<T> fudgetri(1.,0.);    // fudge in triangle coefficients
	complex<T> fudgebox(1.,0.);        // fudge in box coefficients
	complex<T> fudgebubble(1.,0.);        // fudge in box coefficients

	SeriesC<T> intval(-2,0);


#if _VERBOSE
	complex<T> Atree;
//needed since some amplitudes throw exceptions for zero amplitudes
	if (TreeHelAmpl(pro).is_zero()) {Atree=complex<T>(0,0);} else Atree=((*A_Tree_Ptr<T>(pro))(mc,ind));
	_PRINT(Atree);
#endif
	// Add up the bubbles
	SeriesC<T> totalbubble(-2,0);
	complex<T> coeffbubble(0.,0.);

	if (get_mu_index<T>() == 0){
		set_mu_index<T>( DefineMu<T>(mc,m_MU));
	}


	 for (size_t i=1;i<=nbr_bubbles();i++){
		 	bubbleD* thebubble=bubble(i);

	       intval =(* d_bubble_integrals[i-1]->get_value(mc,static_cast<const Index_Vector&>(ind),get_mu_index<T>()));

	       coeffbubble = bubbles[i-1]->get_value(mc,ind);
//	       _PRINT(*bubble(i));
//	       _PRINT(coeffbubble);
	       totalbubble += coeffbubble*intval;
	  }
#if _VERBOSE
	 _PRINT(totalbubble/Atree);
#endif
	 SeriesC<T> totaltri(-2,0);
	 complex<T> coefftri(0.,0.);

	 for (size_t i=1;i<=nbr_triangles();i++){


	       intval =(* d_triangle_integrals[i-1]->get_value(mc,static_cast<const Index_Vector&>(ind),get_mu_index<T>()));

	       // The contribution to the coefficient divided by tree
	       coefftri = triangles[i-1]->get_value(mc,ind);
//	      _PRINT(intval);
	       totaltri += coefftri*intval;
//	       if( verbose ) {
	  }

#if _VERBOSE
	 _PRINT(totaltri/Atree);
#endif
	// Now add up the boxes

	 SeriesC<T> totalbox(-2,0);
	 complex<T> coeffbox;

	 for (unsigned int i=1;i<=nbr_boxes();i++){

	       intval =(* d_box_integrals[i-1]->get_value(mc,static_cast<const Index_Vector&>(ind),get_mu_index<T>()));

	       coeffbox = boxes[i-1]->get_value(mc,ind);
	       totalbox += coeffbox*intval;
	  }
#if _VERBOSE
	 _PRINT(totalbox/Atree);
#endif
	 // Sum together the triangle and box results


	 SeriesC<T> result;
	 result = fudgebubble*totalbubble + fudgetri*totaltri + fudgebox*totalbox ;    // still need to add
	                                                    // bubble and triangle
	  return(result);

}


SeriesC<R> Unknown_Cut_Part_wCI::eval(mom_conf& mc,const vector<int>& ind){
	  return Unknown_Cut_Part_wCI_eval(mc,ind);
}
SeriesC<RHP> Unknown_Cut_Part_wCI::eval(mom_conf_HP& mc,const vector<int>& ind){
	return Unknown_Cut_Part_wCI_eval(mc,ind);
}
SeriesC<RVHP> Unknown_Cut_Part_wCI::eval(mom_conf_VHP& mc,const vector<int>& ind){
	return Unknown_Cut_Part_wCI_eval(mc,ind);
}
#if BH_USE_GMP
SeriesC<RGMP> Unknown_Cut_Part_wCI::eval(momentum_configuration<RGMP>& mc,const vector<int>& ind){
	return Unknown_Cut_Part_wCI_eval(mc,ind);
}
#endif

#endif

}

}
