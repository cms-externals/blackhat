/*
 * rat_worker.h
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef RAT_WORKER_H_
#define RAT_WORKER_H_

#include <vector>
#include "partitions.h" // for symmetry_factor
#include "mode_dependent_typedefs.h"

namespace BH {
class Rec_Tree;
class Tree_factory;
}

namespace BH {

namespace ratext {



class rat_worker {
    std::vector<std::vector<int> > d_corners;
    std::vector<rat_worker*> _parents;
	std::vector<rat_worker* > _daughters;
	std::vector<size_t> _opened_corner;
	std::vector<size_t> _closed_corner;
	symmetry_factor d_symmetry_factor;
protected:
	std::vector<std::vector<TREE_TYPE*> > d_trees;
public:
	rat_worker(std::istream& is);
	void add_parent(rat_worker* p){_parents.push_back(p);};
	void add_parent(rat_worker* p,size_t closed){_parents.push_back(p);_closed_corner.push_back(closed);};
	void add_daughter(rat_worker* p){_daughters.push_back(p);};
	void add_daughter(rat_worker* p,size_t opened){_daughters.push_back(p);_opened_corner.push_back(opened);};
	rat_worker* get_parent(size_t i){return _parents[i-1];}
	rat_worker* get_daughter(size_t i){return _daughters[i-1];}
	//! number of parents
	size_t parents_nbr(){return _parents.size();};
	//! number of daughters
	size_t daughters_nbr(){return _daughters.size();};
	void add_opened_corner(int i){_opened_corner.push_back(i);}
	size_t get_opened_corner(size_t i){return _opened_corner[i-1];}
	size_t get_closed_corner(size_t i){return _closed_corner[i-1];}
    void set_symmetry_factor(int num,int den){d_symmetry_factor=symmetry_factor(num,den);}
    template <class T> T get_symmetry_factor(){return d_symmetry_factor.eval<T>();}
    int corner_nbr() const {return d_corners.size();}
    int corner_size(int cor) const {return d_corners[cor-1].size();}
    int corner_ind(int cor,int ind) const {return d_corners[cor-1][ind-1];}
	size_t decendant_nbr() const {return d_trees.size();};

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
	virtual ~rat_worker();

};

std::ostream&  operator<<(std::ostream& os,const rat_worker&);


}
}
#endif /* RAT_WORKER_H_ */
