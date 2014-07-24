/*
 * cut_part_worker.cpp
 *
 *  Created on: 30-Apr-2009
 *      Author: daniel
 */

#include "cut_part_worker.h"
#ifndef BH_PUBLIC
#include "cut_part_normal.h"
#endif
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <unistd.h>
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "cut_Darren.h"
#ifndef BH_PUBLIC
#include "partial_order.h"
#endif
#include "settings.h"
#include "worker_utils.h"
#include "BH_error.h"
#include "from_file.h"


namespace BH {
	string string_name(const process& pro);
}

using namespace std;

namespace BH {

namespace cut {

namespace worker {

template <class T> struct do_delete : public std::unary_function<T*,void> {
	void operator()(T* ptr) { delete ptr;}
};





worker_cut_part::worker_cut_part(const process& PRO,std::istream& is): standard_cut_part<worker_boxDarren,worker_triangleDarren,worker_bubbleDarren>(PRO){
        string title;
        is >> title;
        assert(title=="boxes");
        int nbr_boxes;
        is >> nbr_boxes;
        for (int i=0;i<nbr_boxes;i++){
//		worker_cutD* new_wcd=new worker_cutD(is);
                d_boxes.push_back(new worker_boxDarren(is));
        }
        is >> title;
        assert(title=="triangles");
        int nbr_triangles;

	// TODO: FIX THIS SO THAT WE DO NOT NEED Normal_Triangle_Specification<worker_cutD> and can use any specalisation we want.
        BH::cut::Darren::triangle_Darren_factory<worker_cutD,Normal_Triangle_Specification<worker_cutD> > wf;

        is >> nbr_triangles;
        for (int i=0;i<nbr_triangles;i++){
                d_triangles.push_back(wf.new_triangle(is));
        }
        is >> title;
        assert(title=="bubbles");
        int nbr_bubbles;
        is >> nbr_bubbles;
        for (int i=0;i<nbr_bubbles;i++){
                d_bubbles.push_back(new worker_bubbleDarren(is));
        }

        is >> title;
        assert(title=="box_triangle_links");
        for ( int i=1;i<=nbr_triangles;i++){
                int daughters_nbr;
                is >> daughters_nbr;
                worker_cutD* tri=d_triangles[i-1];
                for (int j=1;j<=daughters_nbr;j++){
                        int daughter_index;
                        is >> daughter_index;
                        assert(daughter_index>0);
                        assert(daughter_index<=nbr_boxes);
                        worker_cutD* daughter=d_boxes[daughter_index-1];
                        tri->add_daughter(daughter);
                        int opened_corner;
                        is >> opened_corner;
                        tri->add_opened_corner(opened_corner);
                        daughter->add_parent(tri);
                }
        }

        is >> title;
        assert(title=="triangle_bubble_links");
        for ( int i=1;i<=nbr_bubbles;i++){
                int daughters_nbr;
                is >> daughters_nbr;
                worker_cutD* bub=d_bubbles[i-1];
                for (int j=1;j<=daughters_nbr;j++){
                        int daughter_index;
                        is >> daughter_index;
                        assert(daughter_index>0);
                        assert(daughter_index<=nbr_triangles);
                        worker_cutD* daughter=d_triangles[daughter_index-1];
                        bub->add_daughter(daughter);
                        int opened_corner;
                        is >> opened_corner;
                        bub->add_opened_corner(opened_corner);
                        daughter->add_parent(bub);
                }
        }


}




SeriesC<R> worker_cut_part::eval(mom_conf& mc,const vector<int>& ind){
          return eval_fn(mc,ind);
}
SeriesC<RHP> worker_cut_part::eval(mom_conf_HP& mc,const vector<int>& ind){
          return eval_fn(mc,ind);
}
SeriesC<RVHP> worker_cut_part::eval(mom_conf_VHP& mc,const vector<int>& ind){
          return eval_fn(mc,ind);
}

#ifndef BH_PUBLIC
int ttype(const triangleD& cd)
{
//	_MESSAGE8("TRI : ",cd," is type ",cd.c_type(3),", ",cd.c_type(2),", ",cd.c_type(1));
        //Find out what type of triangle this is
        int type[3]={cd.c_type(1),cd.c_type(2),cd.c_type(3)};
        switch(type[2]){
        case a_type://if leg 3 is massless then it means that K2 will be massless
        case b_type:
                return 1;
                break;
        case massive:
                switch(type[1]){
                case a_type://if leg 2 is massless then it means that K3 will be massless
                case b_type:
                        return 1;
                        break;
                case massive:
                        switch(type[0]){
                        case massive:
                                return 3;
                                break;
                        case a_type://If leg 1 is massless then we do not know until we open a bubble leg which type it is as it can be either K1 or K3
                        case b_type:
                                return 1;
                                break;
                        case zero:
                                return 0;
                                break;
                        }
                case zero:
                        return 0;
                        break;
                }
        case zero:
                return 0;
                break;
        }

        //We will never get here but to avoid compiler warnings
        return 0;
}


template <class cutDbase, int CPOINTS, int TPOINTSTRI>  void write(const cut::Darren::box_Darren<cutDbase,CPOINTS,TPOINTSTRI>& cut,ostream& os){
    write(static_cast<const cutDbase&>(cut),os);
    os << "CDspecific " ;
    os << cut.get_kleg(1)  << " ";
    os << cut.get_kleg(2)  << " ";
    os << cut.get_kleg(3)  << " ";
    os << cut.get_kleg(4)  << " ";
    os << cut.get_masslessleg_type()  << " ";
    os << cut.get_massless_K1() << " ";
}

template <class cutDbase, class TriangleSpecs> void write(const cut::Darren::triangle_Darren<cutDbase,TriangleSpecs>& cut,ostream& os){
    write(static_cast<const cutDbase&>(cut),os);
    os << "CDspecific " ;
    os << cut.get_kleg(1)  << " ";
    os << cut.get_kleg(2)  << " ";
    os << cut.get_kleg(3)  << " ";
    os << cut.get_masslessleg_type()  << " ";
}


template <class CutDFactory> void write(normal_cut_part<CutDFactory>& NCP,ostream& os){


        os << "\nboxes\n";
        os << NCP.nbr_boxes() << " ";
        for ( int i=1;i<=NCP.nbr_boxes();i++){
                write(*(NCP.box(i)),os);
        }
        os << "\ntriangles\n";
        os << NCP.nbr_triangles() << " ";
        for ( int i=1;i<=NCP.nbr_triangles();i++){
                os << "TType ";
                os << ttype(*NCP.triangle(i)) << " ";
                write(*(NCP.triangle(i)),os);
        }
        os << "\nbubbles\n";
        os << NCP.nbr_bubbles() << " ";
        for ( int i=1;i<=NCP.nbr_bubbles();i++){
        	write(*(NCP.bubble(i)),os);
        }


        os << "\nbox_triangle_links\n";
        for ( int i=1;i<=NCP.nbr_triangles();i++){
                os << NCP.triangle(i)->daughters_nbr() << " ";
                for (int j=1;j<=NCP.triangle(i)->daughters_nbr();j++){
                        boxD* daughter=NCP.triangle(i)->get_daughter(j);
                        long boxID=daughter->get_ID();
                        int offset=-1;
                        for (int k=1;k<=NCP.nbr_boxes();k++){
                                if (NCP.box(k)->get_ID() == boxID ){
                                        offset=k;
                                        break;
                                }
                        }
                        assert(offset>=0);
                        os << offset << " ";
                        os << NCP.triangle(i)->get_opened_corner(j) <<" ";
                }
        }
        os << "\ntriangle_bubble_links\n";
        for ( int i=1;i<=NCP.nbr_bubbles();i++){
                os << NCP.bubble(i)->daughters_nbr() << " ";
                for (int j=1;j<=NCP.bubble(i)->daughters_nbr();j++){
                        triangleD* daughter=NCP.bubble(i)->get_daughter(j);
                        long triangleID=daughter->get_ID();
                        int offset=-1;
                        for (int k=1;k<=NCP.nbr_triangles();k++){
                                if (NCP.triangle(k)->get_ID() == triangleID ){
                                        offset=k;
                                        break;
                                }
                        }
                        assert(offset>=0);
                        os << offset << " ";
                        os << NCP.bubble(i)->get_opened_corner(j) << " ";
                }
        }
}


template void write(normal_cut_part<Darren_CutD_Factory>& NCP,ostream& os);
#endif

string cut_filename(const process& pro,color_structure cs){
	stringstream ss;
	ss << get_worker_dir("cut");
	ss << "/" << pro.n() ;
	if ( access( ss.str().c_str(), 0 ) != 0 ){
			_WARNING3("Data path ",ss.str(),"not present. Please create it. ");
			throw BHerror("Missing path");
		}	ss << "/cut_";
	ss << pro;
	ss << "_";
	ss << cs;
	ss << ".dat";
	return ss.str();
}
template <class CutDFactory> void add_cut_to_lib(const process& PRO,normal_cut_part<CutDFactory>& NCP,color_structure cs){
		string filename=cut_filename(PRO,cs);
		ifstream is;
		is.open(filename.c_str(),ios::in);
		if( is.is_open() ) {
			_MESSAGE3("File ", filename ," already present. Did nothing.");
		} else {
			ofstream file;
			file.open(filename.c_str());
			_MESSAGE4("Creating cut file for ",NCP.get_process()," ",cs );
			write(NCP,file);
		}
}




Cut_Part_base* worker_cut_part_factory::new_cut_part(const process& PRO,color_structure cs){
	// first try to see if the cut part is known
	if (settings::general::s_use_known_formulae){
		Cut_Part_base* CP_known=Known_cut_part_factory::s_default_KCPF->new_cut_part(PRO,cs);
		if (CP_known){
			return CP_known;
		}
	}

	// if not go on

	string filename=cut_filename(PRO,cs);
		ifstream ifile;
		ifile.open(filename.c_str(),ios::in);
		if( ifile.is_open() ) {
			return new worker_cut_part(PRO,ifile);
		} else {
#ifndef BH_PUBLIC
			cut::normal_cut_part_factory<cut::Darren_CutD_Factory> NCPF;
			_MESSAGE5("Process ",PRO," not in the database for color structure ",cs," creating it now...");
			cut::normal_cut_part<cut::Darren_CutD_Factory>* NCP=dynamic_cast<cut::normal_cut_part<cut::Darren_CutD_Factory>*>(NCPF.new_cut_part(PRO,cs));
			if ( NCP ){
				add_cut_to_lib(PRO,*NCP,cs);
				delete NCP;
			} else {
				_WARNING2("Expect the type to be normal_cut_part<cut::Darren_CutD_Factory>* but got ",typeid(*NCP).name());
				throw BHerror("Unexpected type");
			}
			_MESSAGE("Done.");
			ifstream ifile;
			ifile.open(filename.c_str(),ios::in);
			if( ifile.is_open() ) {
				return new worker_cut_part(PRO,ifile);
			} else {
				// perhaps the data directory exists
				stringstream directory;
				directory << settings::general::s_data_path  << "/Darren_cut/";
				string dir_str=directory.str();
				if ( access( dir_str.c_str(), 0 ) != 0 )
				    {
				        cout << "The directory " << dir_str << " doesn't exist, please create it." << endl;
				        throw BHerror("missing data directory");
				    }
				_WARNING("Amplitude not present, even though we just created it.");
				throw BHerror("Amplitude not present");
			}
#else
			throw BHerror("Amplitude not present");
#endif
		}
}

}
}
}


