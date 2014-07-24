/*
 * counted.cpp
 *
 *  Created on: 24-Apr-2009
 *      Author: daniel
 */

#include "counted.h"
#include <iostream>
#include <cassert>
#include <algorithm>
#include "BH_error.h"

using namespace std;

counter_manager counter_manager::s_cm;

void counter_manager::add(const std::string& name){
	vector<string>::iterator it=find(d_types.begin(),d_types.end(),name);
	static vector<int> max_alive;
	if (it==d_types.end()){
		d_types.push_back(name);
		d_alive.push_back(0);
		d_existed.push_back(0);
		max_alive.push_back(0);
		it=d_types.end();
		it--;

	}
	int offset=it-d_types.begin();
	d_alive[offset]++;
	d_existed[offset]++;

	if ( (d_alive[offset]%10000)==0){
		if (d_alive[offset] > max_alive[offset]){
			cout << "Count of alive instances for " << *it << " reached " << d_alive[offset] << endl;
			max_alive[offset]=d_alive[offset];
		}
	}
}

void counter_manager::remove(const std::string& name){
	vector<string>::iterator it=find(d_types.begin(),d_types.end(),name);
	if (it==d_types.end()){
		throw BH::BHerror("Ref counting error");
	}
	int offset=it-d_types.begin();
	d_alive[offset]--;
	if ( d_alive[offset] < 0 ){
	  cerr << "Problem with the counting of instances of type " << name << ": Supposedly " << d_alive[offset] << " alive.";
	}
}

void counter_manager::print() const {
	int offset=0;
	for (vector<string>::const_iterator it=d_types.begin();it!=d_types.end();it++){
		cout << *it << ": alive " << d_alive[offset] << " existed " << d_existed[offset] << "\n";
		offset++;
	}
}

counter_manager::~counter_manager(){
	print();
}
