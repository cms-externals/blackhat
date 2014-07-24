#include <iostream>
#include <fstream>
#include "histogram.h"
#include "BH_typedefs.h"


template <class T> void histogram<T>::put(T x){
	_total++;
	if (x < _range_min){ _lower_than_range++; return; }
	for (size_t j=1;j<=_nbr_bins;j++){
		if ( x<_limit[j] ) {bin[j-1]++; return ;}
	}
	_higher_than_range++;
}

template <class T> void histogram<T>::print(){
	for (size_t j=0;j<_nbr_bins;j++){
		std::cout << _limit[j] << " " << int(bin[j]) <<std::endl;
	}
}

template <class T> void histogram<T>::print_to_file(const std::string& path){
	std::ofstream output;
	output.open(path.c_str());
	for (size_t j=0;j<_nbr_bins;j++){
		output << _limit[j] << " " << int(bin[j]) <<std::endl;
	}
	output.close();
}

template <class T> void histogram<T>::print_normalized_to_file(const std::string& path){
	std::ofstream output;
	output.open(path.c_str());
	for (size_t j=0;j<_nbr_bins;j++){
		output << _limit[j] << " " << double(bin[j])/double(_total) <<std::endl;
	}
	output.close();
}


template <class T> linear_histogram<T>::linear_histogram(int n,T min,T max) : histogram<T>(n,min,max){
	for (size_t j=0;j<this->nbr_bins();j++){
		this->_limit[j]= min+T(double(j))*(max-min)/T(double(this->nbr_bins()));
	}
}

template class  histogram<BH::R>;
template class  histogram<BH::RHP>;
template class  linear_histogram<BH::R>;
template class  linear_histogram<BH::RHP>;
