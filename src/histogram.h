#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <vector>

/*
  
 A histogram with n bins for values from _range_min to _range_max. It has n+1 values as limits, _limits[0]=_range_min and _limits[n]=_range_max. 
 The i-th bin (bin[i-1]) contains the number of values satisfying
 
 _limits[i-1] <= x < _limits[i]
 
 If values are outside the range, the value of _out_of_range is increased and the value of _higher_than_range or _lower_than_range is incremented, as appropriate.
 
 */



template <class T> class histogram {
long _total;
protected:
	std::vector<int> bin;
	T _range_min;
	T _range_max;
	size_t _nbr_bins;
	std::vector<T> _limit;
	size_t _lower_than_range;
	size_t _higher_than_range;
public:
	explicit histogram(size_t i, T min , T max): _total(0), bin(i) , _range_min(min), _range_max(max), _nbr_bins(i), _limit(i+1), _lower_than_range(0) , _higher_than_range(0) {};
	void put(T);
	size_t nbr_bins(){return _nbr_bins;};
	void print();
	void print_to_file(const std::string& path);
	void print_normalized_to_file(const std::string& path);
	size_t lower_than_range(){return _lower_than_range;};
	size_t higher_than_range(){return _higher_than_range;};
	size_t out_of_range(){return _lower_than_range + _higher_than_range;};
};
 
template <class T> class linear_histogram : public  histogram<T> {
public:
	linear_histogram(int n,T min,T max);
};


#endif /*HISTOGRAM_H_*/
