/*
 * cached_integral.h
 *
 *  Created on: 28-Nov-2008
 *      Author: daniel
 */

#ifndef INDEX_VECTOR_H_
#define INDEX_VECTOR_H_

#include <stdlib.h>
#include <vector>

namespace BH {


//! Class for index vectors
/** an Index vector is basically a vector<int> but allows to carry more inforamtion, like the permutation code. */
class Index_Vector: public std::vector<int> {
	long d_permutation_code;
public:
	Index_Vector(const std::vector<int>& ind);
	long get_permutation_code() const {return d_permutation_code;};
};

//! permutation code
/** the permutation code is supposed to be used only for permutation of index vectors that are permutations of ind4 ... ind9. the long type has to be at least able to hold 2147483647
 so it will work. If compute_permutation_code is called with an index vector longer than 9 or containing larger indices than 9, the result will be inconsistent.
 */
long compute_permutation_code(const std::vector<int>& inds);


//! code for a combination of indices
/** returns a long integer representing a group of external momenta. The indices are sorted so that the ordering of
 * the momenta in the group is  irrelevant
 * example {3,2,5} -> 235    {2,6,1,7} -> 1267
 *
 * the indices in the vector are assumed to be in the range 1,..,9. If this is not the case, the code will not be unique
 * 21 = {2,1} or {1,11} ?
 * */
long index_combination_code(const std::vector<int>& inds);

//! return the indices in {1,...,n} that are not in the vector of integers inds
std::vector<int> index_complement(const std::vector<int>& inds, size_t n);


//1 returns a code corresponding to the invariant formed by the specified momenta
/** the second argument specifies the total number of momenta so that momentum conservation can be used. It is used to give the invariant that has the smalest number of momenta.
 * If both combinations have the same number of momenta, the one containing particle 1 will be chosen.
 * example ({1,2},5) -> 12    ({1,2,3},5) -> 45   ({2,3,5},6) -> 146
 * */
long get_invariant_code(const std::vector<int>& inds, size_t n);


}

#endif /* INDEX_VECTOR_H_ */
