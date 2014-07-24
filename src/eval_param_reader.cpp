/*
 *  eval_param_reader.cpp
 *  BlackHat
 *
 *  Created by Darren Forde on 25/06/2010.
 *  Copyright 2010 BlackHat Collaboration. All rights reserved.
 *
 */

#include "eval_param_reader.h"
#include "BH_error.h"
#include <iostream>

using std::complex;

namespace BH {

template <class T> eval_param_reader<T>::eval_param_reader(const char* path,int nbr_p):  BH::eval_param<T>(nbr_p), nth(0), nbr_particles(nbr_p)  {
	for(int i=0;i<nbr_p;i++){
		_momenta.push_back(new Cmom<T>());
		this->set(i,_momenta[i]);
	}
	
	input.open(path);
	if (! input) {
		std::string errorStr="No file ";
		errorStr+=path;
		errorStr+=" for the constructor eval_param_reader::eval_param_reader.";
		throw BH::BHerror(errorStr.c_str());
	}
	else{
		// Advance to the first element
		next();
	}
	
}

template <class T> eval_param_reader<T>::~eval_param_reader()  {
	for(int i=0;i<nbr_particles;i++){
		delete _momenta[i];
	}
}
	

template <class T> bool eval_param_reader<T>::go_to_pos(std::ios::pos_type pos,size_t n) {
	if (nth==n) return true;
	input.seekg(pos);
	nth=n-1;
	return next();
}


template <class T> bool eval_param_reader<T>::go_to(size_t n) {
	T  dummy;
	if (nth==n) return true;
	if (n>nth) {
		for (size_t i=1;i<n-nth;i++){
			for (size_t j=1;j<=nbr_particles;j++){
				if (!(input >> dummy && input >> dummy && input >> dummy && input >> dummy)) return false;
			}
		}
		nth=n-1;
		return next();
	}
	else if (n<nth) {
		input.seekg(0,std::ios::beg);
		for (size_t i=1;i<n;i++){
			if (!(input >> dummy && input >> dummy && input >> dummy && input >> dummy)) return false;
		}
		nth=n-1;
		
	}
	return next();
}

template <class T> bool eval_param_reader<T>::next() {T  e,x,y,z;
	start_pos=input.tellg();
	this->renew_ID(); // This is a new eval_param so change the ID, this will automatically update this whereever it is used.
	
	// Add the new momenta to the eval_param
	for (size_t j=0;j<nbr_particles;j++){
		if (!(input >> e && input >> x && input >> y && input >> z)) return false;
		*_momenta[j]=BH::Cmom<T>(complex<T>(e,0.),complex<T>(x,0.),complex<T>(y,0.),complex<T>(z,0.));
	}
	nth++;
	return true;
}
	
//
// Explicit Instantiation
//
	
template class eval_param_reader<R>;
template class eval_param_reader<RHP>;
template class eval_param_reader<RVHP>;
    
#if BH_USE_GMP
template class eval_param_reader<RGMP>;
#endif   

}
