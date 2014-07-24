/*
 * ratext_part_worker.cpp
 *
 *  Created on: 7 Jul 2009
 *      Author: daniel
 */

#include "ratext_part.h"
#include "ratext_part_worker.h"
#ifndef BH_PUBLIC
#include "ratext_part_normal.h"
#endif
#include "ratext_part.hpp"
#ifndef BH_PUBLIC
#include "partial_order.h"
#endif
#include "cut_worker.h"  // for read_process_from_stream
#include <iostream>
#include <string>
#include <cassert>
#include "OneLoopHelAmpl.h"
#include <fstream>
#include <typeinfo>
#ifndef BH_PUBLIC
#include "ratext/data_files.h"
#endif
#include "ratext/filename.h"
#include "BH_typedefs.h"
#include "known_rational.h"


using namespace std;

namespace BH {

namespace ratext {


template <class RatBubSpecs, class RatTriSpecs, class RatBoxSpecs, class RatPentSpecs> general_worker_ratext<RatBubSpecs,RatTriSpecs,RatBoxSpecs,RatPentSpecs>::general_worker_ratext(std::istream& is): Rational_base(read_process_from_stream(is)) ,RP(d_process) {
    string title;
    is >> title;
    assert(title=="pentagons");
    int nbr_pentagons;
    is >> nbr_pentagons;
    for (int i=0;i<nbr_pentagons;i++){
    	pent_type* newPent=new pent_type(is);
		RP::d_pentagons.push_back(newPent);
    	string title;
    	is >> title;
    	assert(title == "massLabel" );
    	int m1,m2,m3,m4,m5;
    	is >> m1; is >> m2; is >> m3; is >> m4; is >> m5;
    	newPent->add_mass(m1,m2,m3,m4,m5);
    }

    is >> title;
    assert(title=="boxes");
    int nbr_boxes;
    is >> nbr_boxes;
    for (int i=0;i<nbr_boxes;i++){
//		worker_cutD* new_wcd=new worker_cutD(is);
    	box_type* newBox=new box_type(is);
    	RP::d_boxes.push_back(newBox);
    	string title;
    	is >> title;
    	assert(title == "massLabel" );
    	int m1,m2,m3,m4;
    	is >> m1; is >> m2; is >> m3; is >> m4;
    	newBox->add_mass(m1,m2,m3,m4);

    }
    is >> title;
    assert(title=="triangles");
    int nbr_triangles;


    is >> nbr_triangles;
    for (int i=0;i<nbr_triangles;i++){
    	string title;
    	is >> title;
    	assert(title=="TType");
    	int tri_type;
    	is >> tri_type;
    	triangle_Rat<rat_worker,RatTriSpecs>* newTri;
    	switch (tri_type){
    	case 1: case 2: newTri= new triangle_Rat_plusminus<rat_worker,RatTriSpecs>(is); break;
    	case 3: newTri=new triangle_Rat_3mass<rat_worker,RatTriSpecs>(is);break;
    	}
    	RP::d_triangles.push_back (newTri);

    	is >> title;
    	assert(title == "massLabel" );
    	int m1,m2,m3;
    	is >> m1; is >> m2; is >> m3;
    	newTri->add_mass(m1,m2,m3);

    }
    is >> title;
    assert(title=="bubbles");
    int nbr_bubbles;
    is >> nbr_bubbles;
    for (int i=0;i<nbr_bubbles;i++){
    	bub_type* newBubble=new bub_type(is);
    	RP::d_bubbles.push_back(newBubble);
    	string title;
    	is >> title;
    	assert(title == "massLabel" );
    	int m1,m2;
    	is >> m1; is >> m2;
    	newBubble->add_mass(m1,m2);

    }

    is >> title;
    assert(title=="pentagon_box_links");
    for ( int i=1;i<=nbr_boxes;i++){
            int daughters_nbr;
            is >> daughters_nbr;
            rat_worker* box=RP::d_boxes[i-1];
            for (int j=1;j<=daughters_nbr;j++){
                    int daughter_index;
                    is >> daughter_index;
                    assert(daughter_index>0);
                    assert(daughter_index<=nbr_pentagons);
                    rat_worker* daughter=RP::d_pentagons[daughter_index-1];
                    box->add_daughter(daughter);
                    int opened_corner;
                    is >> opened_corner;
                    box->add_opened_corner(opened_corner);
                    daughter->add_parent(box);
            }
    }

    is >> title;
    assert(title=="box_triangle_links");
    for ( int i=1;i<=nbr_triangles;i++){
            int daughters_nbr;
            is >> daughters_nbr;
            rat_worker* tri=RP::d_triangles[i-1];
            for (int j=1;j<=daughters_nbr;j++){
                    int daughter_index;
                    is >> daughter_index;
                    assert(daughter_index>0);
                    assert(daughter_index<=nbr_boxes);
                    rat_worker* daughter=RP::d_boxes[daughter_index-1];
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
            rat_worker* bub=RP::d_bubbles[i-1];
            for (int j=1;j<=daughters_nbr;j++){
                    int daughter_index;
                    is >> daughter_index;
                    assert(daughter_index>0);
                    assert(daughter_index<=nbr_triangles);
                    rat_worker* daughter=RP::d_triangles[daughter_index-1];
                    bub->add_daughter(daughter);
                    int opened_corner;
                    is >> opened_corner;
                    bub->add_opened_corner(opened_corner);
                    daughter->add_parent(bub);
            }
    }

    is >> title;

    assert (title == "CF");
    int num,den;
    is >> num; is >> den;
    RP::set_colour_fac(num,den);

	// We need to set up some storage in the triangles
	for(int fin=0;fin<RP::d_triangles.size();fin++){
		RP::d_triangles[fin]->final_initialisation();
	}


}




Rational_base* worker_rational_factory::new_rational(const process& PRO,color_structure cs){
	if ( settings::rational_settings::s_set_all_zero ) {
		 return new Known_Rec_Rational(PRO,cs);
	};

	// first try to see if the rational term is known
	if (settings::general::s_use_known_formulae){
		Rational_base* Rb_known=Known_Rational_factory::s_default_KRF->new_rational(PRO,cs);
		if (Rb_known){
			return Rb_known;
		}
	}

	// if not go on

	string filename=rat_filename(PRO,cs);
		ifstream ifile;
		ifile.open(filename.c_str(),ios::in);
		if( ifile.is_open() ) {
			return new worker_ratext (ifile);
		} else {
#ifndef BH_PUBLIC
			normal_ratext_factory NRF;
			_MESSAGE5("Process ",PRO," not in the database for color structure ",cs," creating it now...");
			Rational_base* Rb=NRF.new_rational(PRO,cs);
			normal_ratext* NR=dynamic_cast<normal_ratext*>(Rb);
			if (NR){
				add_rational_to_lib(NR,cs);
				delete NR;
			} else {
				// if the rational term is known, we should not get there
				_WARNING2("Expect the type to be normal_ratext* but got ",typeid(*Rb).name()); throw BHerror("Unexpected type");
			}
			_MESSAGE("Done.");
			ifstream ifile;
			ifile.open(filename.c_str(),ios::in);
			if( ifile.is_open() ) {
				return new worker_ratext (ifile);
			} else {
				_WARNING("Amplitude not present, even though we just created it.");
			}
#else
			_WARNING("Amplitude not present.");
			return 0;
#endif
		}
}


template class ratext_part<pentagon_Rat<rat_worker,Normal_RatPent_Specification<rat_worker> >,box_Rat<rat_worker,Normal_RatBox_Specification<rat_worker> >,triangle_Rat<rat_worker,Normal_RatTri_Specification<rat_worker> >,bubble_Rat<rat_worker,Normal_RatBub_Specification<rat_worker> > >;
template class ratext_part<pentagon_Rat<rat_worker,Higgs_RatPent_Specification<rat_worker> >,box_Rat<rat_worker,Higgs_RatBox_Specification<rat_worker> >,triangle_Rat<rat_worker,Higgs_RatTri_Specification<rat_worker> >,bubble_Rat<rat_worker,Higgs_RatBub_Specification<rat_worker> > >;


}

}
