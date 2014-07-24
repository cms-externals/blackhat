/*
 * print_cutD.h
 *
 *  Created on: 09-Oct-2008
 *      Author: daniel
 */

#ifndef PRINT_CUTD_H_
#define PRINT_CUTD_H_

#include <vector>

namespace  BH {

class Cut_Part;
class OneLoopRawAmplitude;
class Cut_Part_base;
class Cut_Part;
class cutD;

namespace cut {
	template <class T> class normal_cut_part;
	class Darren_CutD_Factory;
}

typedef cut::normal_cut_part<cut::Darren_CutD_Factory> CutType;


template <class T> class momentum_configuration;

template <class CT> void print_tree_graph(CT* A,const char* path);

void print_tree_graph(OneLoopRawAmplitude* A,const char* path);

void print_cut_part_graph(Cut_Part_base* A,const char* path);

void print_tree_graph_non_zero(Cut_Part* A,const char* path,momentum_configuration<double>& mc,const std::vector<int>& ind);

void print_cutD(const cutD& cd,const char* filename);
void print_cutD(const cutD& cd,const char* filename, const std::vector<int>& ind);

void display_cut(const cutD& cd);

}


#endif /* PRINT_CUTD_H_ */
