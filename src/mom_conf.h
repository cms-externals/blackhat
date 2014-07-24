/*!\file mom_conf.h
\brief Header for the momentum configurations
*/

#ifndef MOM_CONF_H_
#define MOM_CONF_H_


#define _WITH_CACHING 0
// with _CACHE_STATISTICS set to 1 (and _WITH_CACHING == 1), information about the caching efficiency is collected. The assumption in programming this was that the speed was not the goal when collecting this information. Please do not use it if you are interested in speed (or don't complain... )
#define _CACHE_STATISTICS 0
#include<iostream>
#include<fstream>
#include<iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <complex>
#include <string>
#include "spinor.h"
#include "particles.h"
#include "BH_typedefs.h"
#include "BH_error.h"
#include "BH_utilities.h"
#include "mc_storage.h"


#if BH_USE_OMP
#include <omp.h>
#include "BH_omp.h"
#endif

#if _USE_GCC
#include <ext/hash_map>

using namespace __gnu_cxx;

namespace __gnu_cxx {
template<> struct hash<std::string>{
  size_t operator()(const std::string& s) const
     {return hash<char const*>()(s.c_str());}
 };
}

#endif
#if SWIG
using namespace std;
#endif

#if _USE_SUN_CC || _USE_PGCC
#include <hash_map>
#endif


namespace BH {

#if _USE_GCC
typedef hash_map<std::string, C, hash<std::string> > ValueCache;
typedef hash_map<std::string, size_t, hash<std::string> > LabelCache;
typedef hash_map<std::string, CHP, hash<std::string> > ValueCacheHP;
#endif

#if _USE_SUN_CC || _USE_PGCC
typedef std::hash_map<std::string, C> ValueCache;
typedef std::hash_map<std::string, size_t> LabelCache;
typedef std::hash_map<std::string, CHP > ValueCacheHP;

#endif


template <class T> class momentum_configuration;
template <class T_low,class T_high> momentum_configuration<T_high> extend_daniel(const momentum_configuration<T_low>& mc,const std::vector<int>& ind);
template <class T_low,class T_high> momentum_configuration<T_high> extend_real(const momentum_configuration<T_low>& mc,const std::vector<int>& ind);

template<class Te,class T0> momentum_configuration<Te>
  extend(const momentum_configuration<T0>& k,
         const std::vector<int>& indices /* of the momenta within "k" */,
         const Te& dummy /* To signal target type */);

template <class T> std::ostream& operator<<(std::ostream&, const momentum_configuration<T>&);

//! template base class for momentum configurations
/** The momentum_configuration_base carries an ID number. This number changes for each new momentum configuration so that it can be used to check that the mom_conf is still the same as the one cached previously.   */
class momentum_configuration_base {

	long int _mc_ID;
	static long int mom_conf_next_ID;

public:
	long int get_ID() const {return _mc_ID;}
	void renew_ID(){_mc_ID=mom_conf_next_ID;mom_conf_next_ID++;};
	momentum_configuration_base() : _mc_ID(mom_conf_next_ID) {mom_conf_next_ID++;};
	virtual ~momentum_configuration_base(){};
};
//! mom_conf contains a collection of (complex) momenta
/** A mom_conf object contains a collection of momenta. In addition, mom_conf objects contain a cache
  for values of expression that are related to the momenta in the momentum configuration. The spinor
  products, invariants and scalar products are accessed through the momentum configuration. The i-th momenta of the
  mom_conf mc can be accessed using mc[i] or mc(i).

   */
template <class T> class momentum_configuration : public momentum_configuration_base 
#if BH_USE_OMP
, private Tools::HasOMPLock
#endif
{
#ifndef SWIG
	friend std::ostream& operator<<<>(std::ostream&, const momentum_configuration<T>&);
#endif
#if _WITH_CACHING
	std::map<int,std::complex<T> > Mspa;
	std::map<int,std::complex<T> > Mspb;

#if _WITH_CACHING && _CACHE_STATISTICS
public:
	struct hash_stat {
		long nbr_try;
		long nbr_no_hit;
		hash_stat() : nbr_try(0) , nbr_no_hit(0) {};
		double hit_rate() {return double(nbr_try-nbr_no_hit)/double (nbr_try);};
		//! says how many times a value is reused in average
		double reuse_rate() {return double(nbr_try-nbr_no_hit)/double (nbr_no_hit);};
	};
	static std::map<std::string,hash_stat> hash_stat_table;
	void cache_statistics();


#endif
#endif
protected:
	size_t nbr;
#if BH_USE_OMP
	Tools::FSArray<Cmom<T>,1000> ps;
	Tools::FSArray<std::complex<T>,1000> ms;
#else
	std::vector<Cmom<T> > ps;
	std::vector<std::complex<T> > ms;
#endif

	size_t _offset;
	momentum_configuration<T>* _parent;



public:

#if _USE_GCC
	hash_map<std::string, std::complex<T>, hash<std::string> > cache;
#endif

#if _USE_SUN_CC || _USE_PGCC
	std::hash_map<std::string, std::complex<T> > cache;
#endif

#ifndef SWIG
	LabelCache labelscache;
#endif
	/*! \return number of momenta  */
	size_t n(){return nbr;};
	const Cmom<T>& p(size_t n) const ;
	const momentum<std::complex<T> >& mom(size_t n) const {
		return p(n).P();
	};
	const Cmom<T>& operator[](const size_t n){
			return p(n);
	};
	const Cmom<T>& operator()(const size_t n){
			return this->p(n);
	};
	momentum_configuration(size_t i=1000);
#ifndef SWIG
	momentum_configuration(const std::vector<plabel>&);
#endif
	momentum_configuration(const std::vector<Cmom<T> >&);
	momentum_configuration(const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	momentum_configuration(const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&,const Cmom<T>&);
	int insert(const Cmom<T>&);
	int insert(const Cmom<T>&,momentum_type);
	int insert(const momentum<std::complex<T> >&);
	int insert(const momentum<std::complex<T> >&,momentum_type type);
	int insert(const lambda<T>&,const lambdat<T>&);
	int insert(const lambdat<T>&,const lambda<T>&);
	void lance_print() ;
	std::complex<T> spa(int, int);
	std::complex<T> spb(int, int);
	std::complex<T> spab(int, int,int);
	std::complex<T> spba(int, int,int);
	std::complex<T> spab(int i, const std::vector<int>& ,int k);
	std::complex<T> spba(int i, const std::vector<int>& , int k);
	std::complex<T> spab(int i, const std::vector<plabel>& ,int k);
	std::complex<T> spba(int i, const std::vector<plabel>& , int k);
	std::complex<T> spaa(int, int,int,int);
	std::complex<T> spbb(int, int,int,int);
	std::complex<T> spab(int, int,int,int,int);
	std::complex<T> spba(int, int,int,int,int);
	std::complex<T> spaa(int, int,int,int,int,int);
	std::complex<T> spbb(int, int,int,int,int,int);
	std::complex<T> s(const std::vector<int>&);
	std::complex<T> s(const std::vector<plabel>&);
	std::complex<T> s(const std::vector<int>&,const std::vector<int>&);
	std::complex<T> s(const std::vector<plabel>&,const std::vector<plabel>&);
	std::complex<T> s(int, int);
	std::complex<T> s(int, int,int);
	std::complex<T> s(int, int,int,int);
	std::complex<T> s(int, int,int,int,int);
	std::complex<T> sp(int, int);
	//! Squared mass of the momentum
	/** \param i index of the momentum \return p_i^2  */
	std::complex<T> m2(size_t i) const;
	//! Lambda spinor of the momentum
	/** \param i index of the momentum \return Lambda spinor corresponding to i (shortcut for p(i).L()  */
	const lambda<T>& L(int i) const  {return this->p(i).L();};
	//! Lambda tilde spinor of the momentum
	/** \param i index of the momentum \return Lambda tilde spinor corresponding to i (shortcut for p(i).Lt()  */
	const lambdat<T>& Lt(int i) const {return this->p(i).Lt();};
	//! Slashed matrix of the momentum
	/** \param i index of the momentum \return slashed matrix corresponding to i (shortcut for p(i).Sm()  */
	smatrix<T> Sm(int i) {return this->p(i).Sm();};
	int Sum(const std::vector<plabel>&);
	int Sum(const std::vector<int>&);
	int Sum(const std::vector<plabel>&,const std::vector<plabel>&);
	int Sum(const std::vector<int>&,const std::vector<int>&);
	int Sum(int,int);
	int Sum(int,int,int);
	int Sum(int,int,int,int);
	int Sum(int,int,int,int,int);
	int Sum(int,int,int,int,int,int);
//	static momentum_configuration* default_mom_conf;
	void clear();
#ifndef SWIG
	void DisplayCache() ;
	void DisplayLabelCache() ;
#endif
	void put_value(const std::string&,std::complex<T> &);
	virtual bool get_value(const std::string&,std::complex<T> &);
	void put_label(const std::string&,size_t&);
	virtual bool get_label(const std::string&,size_t&);
	//! Shifted momenta
	std::vector<int> shiftBA(const std::vector<int>& ind, size_t b, size_t a,const std::complex<T>& z);
	template <class new_T> momentum_configuration<new_T> extend(const std::vector<int>& indices){return BH::extend_real<T,new_T>(*this,indices);}


};

typedef momentum_configuration<R> mom_conf;
typedef momentum_configuration<RHP> mom_conf_HP;
typedef momentum_configuration<RVHP> mom_conf_VHP;

/*! Generates Key to be used in the mom_conf cache*/
std::string GenKey(const char* tag,int);
/*! Generates Key to be used in the mom_conf cache*/
std::string GenKey(const char* tag,int,int);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,int);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,int,int);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,int,int,int);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,int,int,int,int);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,const std::vector<int>&);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,const std::vector<int>&,const std::vector<int>&);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int, const std::vector<int>&);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,int, const std::vector<int>&);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,const std::vector<int>&,const std::vector<int>&);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
/*! Generates Key to be used in the mom_conf cache */
std::string GenKey(const char* tag,int,int,const std::vector<int>&,const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);





template <class T> class sub_momentum_configuration : public momentum_configuration<T> {
public:
	sub_momentum_configuration(momentum_configuration<T>& mc);
	virtual bool get_value(const std::string& key,std::complex<T> &res);
	virtual bool get_label(const std::string& key,size_t &res);
	virtual ~sub_momentum_configuration(){};
};




//! Abstract class for reading momentum configurations from a file. This allows to treat stepping through a file independently from the type (precision) of the points read.
class mom_conf_reader_base  {
public :
	virtual bool next()=0;
	//! reads the n th configuration from the file
	/** \param n One based index, for  n=1 go_to reads the first entry \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read) This member function reads through the file until it arrives at the required position and is therefore slow. If the position where the momentum configuration to be read is known, one should use go_to_pos instead. \sa go_to_pos*/
	virtual bool go_to(size_t n)=0;
	//! reads the n th configuration from the file knowing the position from which to start reading
	/** This function is more efficient than go_to \sa go_to \param pos is the position in the file of the momentum to be read. \param n One based index, for  n=1 go_to reads the first entry. This is necessary so that the mc reader knows where it is in the file, in case the go_to member is called. \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read)  */
	virtual bool go_to_pos(std::ios::pos_type pos,size_t n)=0;
	virtual ~mom_conf_reader_base(){};
};



//! Class for reading momentum configurations with massless momenta from a file
/** A mc_reader object is constructed by giving it the name of the file to read and the number of momenta to read in for each new momentum configuration. When constructed, a new momentum configuration is read from the file by using the next() member function. */
template <class T> class mom_conf_reader : public mom_conf_reader_base, public momentum_configuration<T> {
protected:
	std::ifstream input;
	size_t nth;
	size_t nbr_particles;
	std::ios::pos_type start_pos;
public :
	//! constructor
	mom_conf_reader(const char* path,size_t nbr_p);
	//! reads the next configuration from the file
	/** \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read) */
	virtual bool next();
	//! reads the n th configuration from the file
	/** \param n One based index, for  n=1 go_to reads the first entry \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read) This member function reads through the file until it arrives at the required position and is therefore slow. If the position where the momentum configuration to be read is known, one should use go_to_pos instead. \sa go_to_pos*/
	virtual bool go_to(size_t n);
	//! reads the n th configuration from the file knowing the position from which to start reading
	/** This function is more efficient than go_to \sa go_to \param pos is the position in the file of the momentum to be read. \param n One based index, for  n=1 go_to reads the first entry. This is necessary so that the mc reader knows where it is in the file, in case the go_to member is called. \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read)  */
	virtual bool go_to_pos(std::ios::pos_type pos,size_t n);
};


//! class for multiple precision momentum configuration readers
//** multi_precision_reader_HP objects are derived from mom_conf_reader<RHP> so they can be used as such, but they contain a VHP mom_conf that allow to "jump" to a higher precision when required. The VHP mom_conf only uses next() for the requested momentum configuration. */
#ifndef SWIG
class multi_precision_reader_HP : public mom_conf_reader<RHP> {
	mom_conf_reader<RVHP> _mc_VHP;
public :
	//! constructor
	multi_precision_reader_HP(const char* path,size_t nbr_p): mom_conf_reader<RHP>(path,nbr_p), _mc_VHP(path,nbr_p){};
	//! reference to the corresponding VHP momentum configuration
	mom_conf_VHP& mc_VHP(){_mc_VHP.go_to_pos(start_pos,nth);return _mc_VHP;};
};


//! class for multiple precision momentum configuration readers
//** multi_precision_reader objects are derived from mom_conf_reader<R> so they can be used as such, but they contain a HP and VHP mom_conf that allow to "jump" to a higher precision when required. The HP and VHP mom_conf only use next() for the requested momentum configuration. */

class multi_precision_reader : public mom_conf_reader<R> {
	multi_precision_reader_HP _mc_HP;
public :
	//! constructor
	multi_precision_reader(const char* path,size_t nbr_p): mom_conf_reader<R>(path,nbr_p), _mc_HP(path,nbr_p){};
	//! reference to the corresponding HP momentum configuration
	mom_conf_HP& mc_HP(){_mc_HP.go_to_pos(start_pos,nth);return _mc_HP;};
};


typedef mom_conf_reader<R> mc_reader;
typedef mom_conf_reader<RHP> mc_reader_HP;
typedef mom_conf_reader<RVHP> mc_reader_VHP;


typedef sub_momentum_configuration<R> sub_mom_conf;
typedef sub_momentum_configuration<RHP> sub_mom_conf_HP;
typedef sub_momentum_configuration<RVHP> sub_mom_conf_VHP;


#if _CACHE_STATISTICS

template <class key_type, class T> void find_and_stat(const std::string& tag, key_type& key,std::map<key_type,std::complex<T> >& value_map, std::complex<T> insert_value,std::complex<T> return_value, std::complex<T>& res ) ;
template <class key_type, class T> void find_and_stat(const std::string& tag, key_type& key,hash_map<std::string, std::complex<T>, hash<std::string> >& value_map, std::complex<T> insert_value,std::complex<T> return_value, std::complex<T>& res ) ;
template <class key_type, class T> void find_and_stat(const std::string& tag, key_type& key,momentum_configuration<T>* mc , std::complex<T> insert_value,std::complex<T> return_value, std::complex<T>& res ) ;


#endif

void mathprint(momentum_configuration<RVHP>& momc);

#include "mom_conf_inline.h"
#endif



}
#endif /*MOM_CONF_H_*/
