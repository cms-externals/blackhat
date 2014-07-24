#ifndef CUT_WORKER_H
#define CUT_WORKER_H

#include "partitions.h"
#include <vector>
#include <map>
#include <complex>
#include <iosfwd>
#include "cut_part.h"
#include "mode_dependent_typedefs.h"

namespace BH {
    class Rec_Tree;
}

using namespace std;

namespace BH {

namespace cut {

namespace worker {


// TODO: split the functionalities into FileReadable, HasDaughters, HasParents
class worker_cutD {
protected:
	std::vector<std::vector<int> > d_corners;
        std::vector<TREE_TYPE*> d_trees;
        std::vector<int> d_opened_corner;
        std::vector<worker_cutD*> d_parents;
        std::vector<worker_cutD*> d_daughters;
        symmetry_factor d_symmetry_factor;
        std::vector<corner_type> d_corner_types;
    	double _accuracy; // A variable for storing the numerical accuracy of the computation
public:
        worker_cutD(std::istream& os);
        void print(std::ostream& os) const ;
        void add_parent(worker_cutD* wcp){ d_parents.push_back(wcp);}
        void add_daughter(worker_cutD* wcp){ d_daughters.push_back(wcp);}
        worker_cutD* get_parent(int n){ return d_parents[n-1];}
        worker_cutD* get_daughter(int n){ return d_daughters[n-1];}
        std::complex<R> eval_tree(int n,momentum_configuration<R>& mc,const std::vector<int>& ind);
        std::complex<RHP> eval_tree(int n,momentum_configuration<RHP>& mc,const std::vector<int>& ind);
        std::complex<RVHP> eval_tree(int n,momentum_configuration<RVHP>& mc,const std::vector<int>& ind);
#if BH_USE_GMP
        std::complex<RGMP> eval_tree(int n,momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
#endif
        std::complex<R> eval_tree(int n,eval_param<R>& ep);
        std::complex<RHP> eval_tree(int n,eval_param<RHP>& ep);
        std::complex<RVHP> eval_tree(int n,eval_param<RVHP>& ep);
#if BH_USE_GMP
        std::complex<RGMP> eval_tree(int n,eval_param<RGMP>& ep);
#endif
        const std::vector<int>& corner(int n) const {return d_corners[n-1];}
        int daughters_nbr() const {return d_daughters.size();}
        int parents_nbr() const {return d_parents.size();}
        void add_opened_corner(int opened_corner){d_opened_corner.push_back(opened_corner);}
        int get_opened_corner(int opened_corner) const { return d_opened_corner[opened_corner-1];}
        template <class T> T get_symmetry_factor() const {return d_symmetry_factor.eval<T>();}
        corner_type c_type(size_t cor) const ;
        int corner_size(int cor) const {return d_corners[cor-1].size();}
        int corner_ind(int cor,int ind) const {return d_corners[cor-1][ind-1];}
    	// We store a number that gives the accuracy of the computation
    	void set_accuracy(double acc) {_accuracy=acc;};
    	double get_accuracy() const {return _accuracy;};

        virtual ~worker_cutD();

};


ostream& operator<<(ostream& os,const worker_cutD& wc);
#ifndef BH_PUBLIC

void write(const part& pa,ostream& os);
void write(const cutD& cut,ostream& os);

#endif
}
}
}

#endif // CUT_WORKER_H


