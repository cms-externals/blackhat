/*
 * index_vector.cpp
 *
 *  Created on: 01-Dec-2008
 *      Author: daniel
 */

/*
 * cached_integral.cpp
 *
 *  Created on: 28-Nov-2008
 *      Author: daniel
 */

#include "index_vector.h"
#include "cached_integral.h"
#include <algorithm>
#include <set>
#include "BH_typedefs.h"
#include "helcode.h"
#include "integrals.h"

using namespace std;

namespace BH {

Index_Vector::Index_Vector(const std::vector<int>& ind) : std::vector<int>(ind){ d_permutation_code=compute_permutation_code(ind);}


long compute_permutation_code(const vector<int>& inds){

	size_t n=inds.size();
	  long result = 0;
	  long base  = 10;
	  long power = 1;
	  for (int i=1;i<=n;i++){
	    result += power*inds[n-i];
	  //   std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
	    power *= base;
	    }
	  return(result);
}

long index_combination_code(const vector<int>& inds){
	vector<int> sorted(inds.begin(),inds.end());
	sort(sorted.begin(),sorted.end());
	return compute_permutation_code(sorted);
}

vector<int> index_complement(const vector<int>& inds, size_t n){
	vector<int> ref(n),res;
	set<int> sorted;copy(inds.begin(),inds.end(),inserter(sorted,sorted.begin()));
	for (int i=1;i<=n;i++){ref[i-1]=i;};
	set_difference(ref.begin(),ref.end(),sorted.begin(),sorted.end(),back_inserter(res));
	return res;
}

long get_invariant_code(const vector<int>& inds, size_t n){
	if ( 2*inds.size() < n ) {
		return index_combination_code(inds);
	} else if ( 2*inds.size() > n ){
		return index_combination_code(index_complement(inds,n));
	} else {
		if (find(inds.begin(),inds.end(),1) != inds.end()){
			return index_combination_code(inds);
		} else {
			return index_combination_code(index_complement(inds,n));
		}
	}
	return 0;
}


//long compute_permutation_code(const vector<int>& inds){
//
//	size_t n=inds.size();
//	  long result = 0;
//	  long base  = 10;
//	  long power = 1;
//	  for (int i=1;i<=n;i++){
//	    result += power*inds[n-i];
//	  //   std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
//	    power *= base;
//	    }
//	  return(result);
//}

}
