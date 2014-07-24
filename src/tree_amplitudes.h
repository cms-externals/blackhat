/*
 * tree_amplitudes.h
 *
 *  Created on: 16-Apr-2009
 *      Author: daniel
 */

#ifndef TREE_AMPLITUDES_H_
#define TREE_AMPLITUDES_H_

#include <complex>
#include <vector>
#include <iosfwd>
#include "BH_typedefs.h"
#include "process.h"
#include "mode_dependent_typedefs.h"


#if BH_USE_GMP
#include "gmp_r.h"
#endif

namespace BH {

class process;
template <class T> class momentum_configuration;
template <class T> class eval_param;


//! Class for helicity amplitudes
class HelAmpl {
protected:
	process d_process;
public:
//	virtual C eval(mom_conf);
	virtual void print();
	virtual bool is_zero() const;
	virtual bool is_split_helicity();
	HelAmpl()  {};
	HelAmpl(const process& p) : d_process(p) {};
	virtual ~HelAmpl(){};
	friend std::ostream& operator<<(std::ostream& s, HelAmpl& p);
	const process& get_process() const {return d_process;};
};


//! Class for tree helicity amplitudes
class TreeHelAmpl : public HelAmpl {
	TREE_TYPE* d_tree_ptr;
public :
	virtual std::complex<R> eval(momentum_configuration<R>&,const std::vector<int>& ind);
	virtual std::complex<RHP> eval(momentum_configuration<RHP>&,const std::vector<int>& ind);
	virtual std::complex<RVHP> eval(momentum_configuration<RVHP>&,const std::vector<int>& ind);
#if BH_USE_GMP
	virtual std::complex<RGMP> eval(momentum_configuration<RGMP>&,const std::vector<int>& ind);
#endif
	virtual std::complex<R> eval(const eval_param<R>&);
	virtual std::complex<RHP> eval(const eval_param<RHP>&);
	virtual std::complex<RVHP> eval(const eval_param<RVHP>&);

	virtual bool is_zero () const ;
	TreeHelAmpl(const process& p);
	TreeHelAmpl(const TreeHelAmpl& T);
	TreeHelAmpl operator=(const TreeHelAmpl& T);
	virtual ~TreeHelAmpl();
	TREE_TYPE* pointee(){return d_tree_ptr;};
};


}

#endif /* TREE_AMPLITUDES_H_ */
