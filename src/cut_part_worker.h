/*
 * cut_part_worker.h
 *
 *  Created on: 30-Apr-2009
 *      Author: daniel
 */

#ifndef CUT_PART_WORKER_H_
#define CUT_PART_WORKER_H_

#include <vector>
#include <map>
#include <complex>
#include <iosfwd>
#include "cut_part.h"
#include "cut_Darren.h"
#include "standard_cut_part.h"
#include "cut_part_factory.h"

#include "bubble_specification.h"


using namespace std;

namespace BH {

namespace cut {

template <class CutDFactory> class normal_cut_part;


namespace worker {

typedef Normal_Triangle_Specification<BH::cut::worker::worker_cutD> NTS;
typedef Normal_Bubble_Specification<BH::cut::worker::worker_cutD> NBS;
//typedef General_Bubble_Specification<BH::cut::worker::worker_cutD,CTRIPOINTS_STD,BUBPOINTS_STD,TRIPOINTS_STD,BUBYPOINTS_STD> NBS;
//typedef General_Triangle_Specification<BH::cut::worker::worker_cutD::worker_cutD,CTRIPOINTS_STD,BUBPOINTS_STD,TRIPOINTS_STD,BUBYPOINTS_STD> NTS;
	
	
typedef BH::cut::Darren::box_Darren<BH::cut::worker::worker_cutD,CTRIPOINTS_STD,TRIPOINTS_STD> worker_boxDarren;
typedef BH::cut::Darren::triangle_Darren<BH::cut::worker::worker_cutD,NTS> worker_triangleDarren;
typedef BH::cut::Darren::bubble_Darren<BH::cut::worker::worker_cutD,NBS> worker_bubbleDarren;

#ifndef BH_PUBLIC

//void write(Cut_Part& CP,ostream& os);
template <class CutDFactory> void write(normal_cut_part<CutDFactory>& CP,ostream& os);

#endif

class worker_cut_part: public standard_cut_part<worker_boxDarren,worker_triangleDarren,worker_bubbleDarren> {
public:
	worker_cut_part(const process& PRO, std::istream& is);
	virtual SeriesC<R> eval(mom_conf&,const std::vector<int>&);
    virtual SeriesC<RHP> eval(mom_conf_HP&,const std::vector<int>&);
    virtual SeriesC<RVHP> eval(mom_conf_VHP&,const std::vector<int>&);

	virtual SeriesC<R> eval(const eval_param<R>&);
    virtual SeriesC<RHP> eval(const eval_param<RHP>&);
    virtual SeriesC<RVHP> eval(const eval_param<RVHP>&);

private:
        //no copy or operator assignement
        worker_cut_part(const worker_cut_part&);
        worker_cut_part& operator=(const worker_cut_part&);
};


class worker_cut_part_factory : public cut_part_factory<Cut_Part_base> {
public:
	virtual Cut_Part_base* new_cut_part(const process&,color_structure);
	static cut_part_factory<Cut_Part_base>* s_default_WCPF;
};


}

}

typedef cut::worker::worker_cut_part_factory cut_part_factory_worker  ;
typedef cut::worker::worker_cut_part cut_part_worker  ;

}



#endif /* CUT_PART_WORKER_H_ */
