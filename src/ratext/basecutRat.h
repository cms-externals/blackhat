/*
 * basecutRat.h
 *
 *  Created on: 7 Jul 2009
 *      Author: daniel
 */

#ifndef BASECUTRAT_H_
#define BASECUTRAT_H_

#ifndef BH_PUBLIC
#include "cutRat_comp.h"
#endif
#include "rec_tree.h"
#if __USE_NEW_REC_TREE_REC_RAT==1
#include "rec_tree_new.h"
#endif

//#define __RAT_TREE_TYPE BG_TREE_FACTORY_V
#define __RAT_TREE_TYPE TREE_FACTORY_TYPE


namespace BH {

#ifndef BH_PUBLIC


/*
 *
 *
 * The basecutD objects from which all other base objects should be derived
 *
 *
 */

class baseCutD : public part{
private:
	symmetry_factor _symmetry_factor;

protected:
#if __USE_NEW_REC_TREE_REC_RAT==1
    std::vector<rec_tree_eval*> d_trees; // Stores the trees for each corner of each contributing diagram
#else
    std::vector<TREE_TYPE*> d_trees; // Stores the trees for each corner of each contributing diagram
#endif
    size_t _width; // Stores the number of legs for the CutD so that we can extract the trees from d_trees correctly.
	size_t _iclass; // Stores the class of this cut object we can use this to store information that distinguishes two otherwise mathcing cuts, e.g. the sign of the momentum param of the triangles

public:
	baseCutD() : part() {_symmetry_factor=symmetry_factor(1,1);};
	baseCutD(const cutD& cd) : part(cd) {_symmetry_factor=symmetry_factor(1,1);};
	baseCutD(const cutD* cd) : part(*cd) {_symmetry_factor=symmetry_factor(1,1);};
	virtual ~baseCutD();

	// Checks to see if the cutD being passed to this has the same corners irrespective of cut propagators
	bool is_match(cutD* cmpr);
	// Checks to see if the cutD being passed to this has the same corners irrespective of cut propagators and was opened from the same object
	bool is_match(cutD* cmpr, size_t iclass);

	// Returns the symmetry factor
	template <class T> T get_symmetry_factor(){return _symmetry_factor.eval<T>();}
	symmetry_factor get_symmetry_factor_no_eval() const {return _symmetry_factor;}
	void set_symmetry_factor(const symmetry_factor& sym_fac){_symmetry_factor=sym_fac;};

//	virtual const cutD* get_decendant(int i) const {_WARNING("get_decendant not defined");return 0;};
//	virtual size_t decendant_nbr() const {_WARNING("decendant_nbr not defined");return 0;};
	virtual const cutD* get_decendant(int i) const =0;
	virtual size_t decendant_nbr() const =0;

	// Used for computing the trees at each corner
    std::complex<R> eval_tree(int descendant_n,int tree_n,momentum_configuration<R>& mc,const std::vector<int>& ind);
    std::complex<RHP> eval_tree(int descendant_n,int tree_n,momentum_configuration<RHP>& mc,const std::vector<int>& ind);
    std::complex<RVHP> eval_tree(int descendant_n,int tree_n,momentum_configuration<RVHP>& mc,const std::vector<int>& ind);
    std::complex<R> eval_tree(int descendant_n,int tree_n,const eval_param<R>& ep);
    std::complex<RHP> eval_tree(int descendant_n,int tree_n,const eval_param<RHP>& ep);
    std::complex<RVHP> eval_tree(int descendant_n,int tree_n,const eval_param<RVHP>& ep);
#if BH_USE_GMP
    std::complex<RGMP> eval_tree(int descendant_n,int tree_n,momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
    std::complex<RGMP> eval_tree(int descendant_n,int tree_n,const eval_param<RGMP>& ep);
#endif
    int corner_size(int cor) const {return c(cor).size();}
    int corner_ind(int cor,int ind) const {return c(cor)[ind-1].ind();}
	
private:
    baseCutD(const baseCutD& cd) : part(cd) {};
    baseCutD& operator=(const baseCutD& cd) {};
};



/*
 *
 *
 * The base bubble, triangle, box and pentagon objects these should be used for derivation only
 *
 *
 */

class basetriangleRat;
class baseboxRat;
class basepentagonRat;

class basebubbleRat : public baseCutD, public computable<std::complex>  {
	std::vector<basetriangleRat* > _daughters;
	std::vector<size_t> _opened_corner;

	double _accuracy;

protected:
	std::vector<const bubbleRat_comp*> _decendantbubbleRat;

	// We store a number that gives the accuracy of the computation
	void set_accuracy(double acc){_accuracy=acc;};

public:
	basebubbleRat(const bubbleRat_comp* pcd);
	virtual ~basebubbleRat(){};

	//! adds a new daughter
	/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_daughter(basetriangleRat* p){_daughters.push_back(p);};
	void add_daughter(basetriangleRat* p,size_t opened){_daughters.push_back(p);_opened_corner.push_back(opened);};

	/** \param i (1-based) parent label \return pointer to the ith parent */
	basetriangleRat* get_daughter(size_t i){return _daughters[i-1];}
	//! number of parents
	size_t daughters_nbr(){return _daughters.size();};
	size_t get_opened_corner(size_t i){return _opened_corner[i-1];}
	void set_opened_corner(size_t i,size_t oc){_opened_corner[i-1]=oc;}

	// We store a number that gives the accuracy of the computation
	double get_accuracy(){return _accuracy;};

	// Add a cutD to the list
	void add(const bubbleRat_comp* pcd);
	const cutD* get_decendant(int i) const {return _decendantbubbleRat[i-1];};
	size_t decendant_nbr() const {return _decendantbubbleRat.size();};
};

class basetriangleRat : public baseCutD, public computable<std::complex>  {
	std::vector<basebubbleRat*> _parents;
	std::vector<baseboxRat* > _daughters;
	std::vector<size_t> _opened_corner;
	std::vector<size_t> _closed_corner;

protected:
	std::vector<const triangleRat_comp*> _decendanttriangleRat;

public:
	basetriangleRat(const triangleRat_comp* pcd);
	virtual ~basetriangleRat(){};

	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(basebubbleRat* p){_parents.push_back(p);};
	void add_parent(basebubbleRat* p,size_t closed){_parents.push_back(p);_closed_corner.push_back(closed);};
	//! adds a new daughter
	/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_daughter(baseboxRat* p){_daughters.push_back(p);};
	void add_daughter(baseboxRat* p,size_t opened){_daughters.push_back(p);_opened_corner.push_back(opened);};
	//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	basebubbleRat* get_parent(size_t i){return _parents[i-1];}
	//! pointer to a daughter
	/** \param i (1-based) daughter label \return pointer to the ith daughter */
	baseboxRat* get_daughter(size_t i){return _daughters[i-1];}
	//! number of parents
	size_t parents_nbr(){return _parents.size();};
	//! number of daughters
	size_t daughters_nbr(){return _daughters.size();};

	size_t get_opened_corner(size_t i){return _opened_corner[i-1];}
	size_t get_closed_corner(size_t i){return _closed_corner[i-1];}
	void set_opened_corner(size_t i,size_t oc){_opened_corner[i-1]=oc;}
	void set_closed_corner(size_t i,size_t cc){_closed_corner[i-1]=cc;}

	// Add a cutD to the list
	void add(const triangleRat_comp* pcd);
	const cutD* get_decendant(int i) const {return _decendanttriangleRat[i-1];};
	size_t decendant_nbr() const {return _decendanttriangleRat.size();};

};


class baseboxRat : public baseCutD, public computable<std::complex>  {
	std::vector<basetriangleRat*> _parents;
	std::vector<basepentagonRat* > _daughters;
	std::vector<size_t> _opened_corner;
	std::vector<size_t> _closed_corner;

protected:
	std::vector<const boxRat_comp*> _decendantboxRat;

public:
	baseboxRat(const boxRat_comp* pcd);
	virtual ~baseboxRat(){};

	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(basetriangleRat* p){_parents.push_back(p);};
	void add_parent(basetriangleRat* p,size_t closed){_parents.push_back(p);_closed_corner.push_back(closed);};
	//! adds a new daughter
	/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_daughter(basepentagonRat* p){_daughters.push_back(p);};
	void add_daughter(basepentagonRat* p,size_t opened){_daughters.push_back(p);_opened_corner.push_back(opened);};
	//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	basetriangleRat* get_parent(size_t i){return _parents[i-1];}
	//! pointer to a daughter
	/** \param i (1-based) daughter label \return pointer to the ith daughter */
	basepentagonRat* get_daughter(size_t i){return _daughters[i-1];}
	//! number of parents
	size_t parents_nbr(){return _parents.size();};
	//! number of daughters
	size_t daughters_nbr(){return _daughters.size();};

	size_t get_opened_corner(size_t i){return _opened_corner[i-1];}
	size_t get_closed_corner(size_t i){return _closed_corner[i-1];}
	void set_opened_corner(size_t i,size_t oc){_opened_corner[i-1]=oc;}
	void set_closed_corner(size_t i,size_t cc){_closed_corner[i-1]=cc;}

	// Add a cutD to the list
	void add(const boxRat_comp* pcd);
	const cutD* get_decendant(int i) const {return _decendantboxRat[i-1];};
	size_t decendant_nbr() const {return _decendantboxRat.size();};
};


class basepentagonRat : public baseCutD, public computable<std::complex>  {
	std::vector<baseboxRat*> _parents;
	std::vector<size_t> _closed_corner;

protected:
	std::vector<const pentagonRat_comp*> _decendantpentagonRat;

public:
	basepentagonRat(const pentagonRat_comp* pcd);
	virtual ~basepentagonRat(){};

	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(baseboxRat* p){_parents.push_back(p);};
	void add_parent(baseboxRat* p,size_t closed){_parents.push_back(p);_closed_corner.push_back(closed);};
	/** \param i (1-based) parent label \return pointer to the ith parent */
	baseboxRat* get_parent(size_t i){return _parents[i-1];}
	//! number of parents
	size_t parents_nbr(){return _parents.size();};

	size_t get_closed_corner(size_t i){return _closed_corner[i-1];}
	void set_closed_corner(size_t i,size_t cc){_closed_corner[i-1]=cc;}

	// Add a cutD to the list
	virtual void add(const pentagonRat_comp* pcd);
	virtual const cutD* get_decendant(int i) const {return _decendantpentagonRat[i-1];};
	virtual size_t decendant_nbr() const {return _decendantpentagonRat.size();};
};


#endif

}


#endif /* BASECUTRAT_H_ */
