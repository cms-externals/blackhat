/*!\file tree.cpp
\brief Implementation of the tree level amplitudes
*/



#include "process.h"
#include "BH_A0.h"
#include "mom_conf.h"


using namespace std;

namespace BH {

enum zero_amplitudes { odd_nbr_q = -1, odd_nbr_q2 = -2, odd_nbr_l = -3, odd_nbr_l1 = -4, odd_nbr_l2 = -5, l_no_q = -6, no_ferm_conservation = -7 };


size_t nbr_of_flavors(const process& pro,const particle& type){
	size_t count=0;
  	vector<int> flavors;
  	for(int i=1;i<=pro.n();i++)
  	{
  		if( (*pro.p(i).type()) == type ) {
  			vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
  	    	if ( pos == flavors.end() ){
  	    		count++;
  	    		flavors.push_back(pro.p(i).flavor());
  	    	}
		}
  	}
  	return count;
}

size_t nbr_of_flavors(const process& pro,const particle& type1,const particle& type2){
	size_t count=0;
  	vector<int> flavors;
  	for(int i=1;i<=pro.n();i++)
  	{
  		if( ((*pro.p(i).type()) == type1) || ((*pro.p(i).type()) == type2) ) {
  			vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
  	    	if ( pos == flavors.end() ){
  	    		count++;
  	    		flavors.push_back(pro.p(i).flavor());
  	    	}
		}
  	}
  	return count;
}



}
