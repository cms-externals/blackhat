/*
 * counted.h
 *
 *  Created on: 24-Apr-2009
 *      Author: daniel
 */

#ifndef COUNTED_H_
#define COUNTED_H_

#include <string>
#include <vector>
#include <typeinfo>

// not many types should be handeled at the same time
class counter_manager {
	std::vector<std::string> d_types;
	std::vector<long> d_alive;
	std::vector<long> d_existed;
public:
	void add(const std::string& name);
	void remove(const std::string& name);
	void print() const ;
	static counter_manager s_cm;
	~counter_manager();
};


template <class T> class counted {
public:
	counted(){
		counter_manager::s_cm.add(std::string(typeid(T).name()));
	};
	~counted(){
		counter_manager::s_cm.remove(std::string(typeid(T).name()));
	};
	counted(const counted& c){
		counter_manager::s_cm.add(std::string(typeid(T).name()));
	}
};



#endif /* COUNTED_H_ */
