#include "cut_worker.h"
#include "worker_utils.h"
#include <iostream>
#include <cassert>
#include <memory>
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "cut_Darren.h"
#include "process_utils.h"
#include "cut_part_worker.h"
#include "mode_dependent_typedefs.h"

#ifdef BH_PUBLIC
#include "worker_tree.h"
#endif


using namespace std;

using BH::worker::write;
using BH::worker::read_process_from_stream;

namespace BH {

namespace cut {

namespace worker {

#ifndef BH_PUBLIC

void write_to_file(const process& PRO,const string& filename){

        ofstream file;
        file.open(filename.c_str());

        BH::worker::write(PRO,file);

        file.close();
}

void write(const part& pa,ostream& os){

//	 os << pa.get_code() << " ";

         int nc=pa.nc();
         os << nc << " ";
         for (int i=1;i<=nc;i++){
                 const std::vector<plabel>& corner=pa.c(i);
                 os << corner.size() << " " ;
                 for (int j=0;j<corner.size();j++){
                         os << corner[j].ind() << " ";
                 }
         }
}

void write(const cutD& cut,ostream& os){

                write(static_cast<const part&>(cut),os);
         for (int i=1;i<=cut.nc();i++){
                 BH::worker::write(cut.get_process(i),os);
         }
         os << "SF ";
         os << cut.get_symmetry_factor_no_eval().get_num() << " " << cut.get_symmetry_factor_no_eval().get_den() << " ";
         os << "CT ";
         for (int i=1;i<=cut.nc();i++){
                switch (cut.c_type(i)){
                case a_type: os << "a ";break;
                case b_type: os << "b ";break;
                case m_type: os << "m ";break;
                case zero: os << "z ";break;
                case massive: os << "M ";break;
                default: os << "u ";break;
                };
         }

}

#endif

bool read_process_from_file(process& ret,const string& filename){
        ifstream ifile(filename.c_str());
        if (!ifile) return false;
        read_process_from_stream(ret,ifile);
        ifile.close();
        return true;
}


worker_cutD::~worker_cutD(){
        for (int j=0;j<d_trees.size();j++){
                delete d_trees[j];
        }
}

corner_type worker_cutD::c_type(size_t cor) const {
        return d_corner_types[cor-1];
}

std::complex<R> worker_cutD::eval_tree(int n,momentum_configuration<R>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}
std::complex<RHP> worker_cutD::eval_tree(int n,momentum_configuration<RHP>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}
std::complex<RVHP> worker_cutD::eval_tree(int n,momentum_configuration<RVHP>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}
#if BH_USE_GMP
std::complex<RGMP> worker_cutD::eval_tree(int n,momentum_configuration<RGMP>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}
#endif
std::complex<R> worker_cutD::eval_tree(int n,eval_param<R>& ep){
        return d_trees[n-1]->eval(ep);
}
std::complex<RHP> worker_cutD::eval_tree(int n,eval_param<RHP>& ep){
        return d_trees[n-1]->eval(ep);
}
std::complex<RVHP> worker_cutD::eval_tree(int n,eval_param<RVHP>& ep){
        return d_trees[n-1]->eval(ep);
}
#if BH_USE_GMP
std::complex<RGMP> worker_cutD::eval_tree(int n,eval_param<RGMP>& ep){
        return d_trees[n-1]->eval(ep);
}
#endif


ostream& operator<<(ostream& os,const worker_cutD& wc){
#ifndef BH_PUBLIC
        wc.print(os);
#endif
        return os;
}

#ifndef BH_PUBLIC


void worker_cutD::print(std::ostream& os) const {
        for (int i=0;i<d_corners.size();i++){
                os << "{ ";
                for (int j=0;j<d_corners[i].size();j++){
                        os << d_corners[i][j] << " ";
                }
                os << "}";
        }
        for (int i=0;i<d_trees.size();i++){
                os << d_trees[i]->get_process() << " ";
        }
}

#endif


worker_cutD::worker_cutD(std::istream& is){
        int nbr_corner;
        is >> nbr_corner;
        assert(nbr_corner > 1);
        assert(nbr_corner < 6);
        d_corners.reserve(nbr_corner);
        for (int i=0;i<nbr_corner;i++){
                d_corners.push_back(vector<int>());
                int nbr_entries;
                is >> nbr_entries;
                for (int j=0;j<nbr_entries;j++){
                        int entry;
                        is >> entry ;
                        assert(entry);
                        d_corners[i].push_back(entry);
                }
        }
        TREE_FACTORY_TYPE TF;
        for (int i=0;i<nbr_corner;i++){
                process PRO;
                read_process_from_stream(PRO,is);
                d_trees.push_back(TF.new_tree(fix_flavors(PRO)));
        }
        string title;
        is >> title;
        assert(title=="SF");
        int num,den;
        is >> num; is >> den;
        d_symmetry_factor=symmetry_factor(num,den);
        is >> title;
        assert(title=="CT");
         for (int i=1;i<=nbr_corner;i++){
                char ctype;
                is >> ctype;
                switch (ctype){
                case 'a': d_corner_types.push_back(a_type);break;
                case 'b': d_corner_types.push_back(b_type);break;
                case 'z': d_corner_types.push_back(zero);break;
                case 'M': d_corner_types.push_back(massive);break;
                case 'm': d_corner_types.push_back(m_type);break;
                };
         }

}





}



}
}
