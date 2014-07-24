/*
 * partitions_fwd.h
 *
 *  Created on: 22-Apr-2009
 *      Author: daniel
 */

#ifndef PARTITIONS_FWD_H_
#define PARTITIONS_FWD_H_

namespace BH {

//! corner type
/** A corner i of type
  \n - "a" if the corner is massless and it is expressed only in terms of < > products
  \n - "b" if the corner is massless and it is expressed only in terms of [ ] products
  \n - "massive" if the corner is massive
*/
enum corner_type {
	a_type, 		/*!< a-type */
	b_type, 	/*!< b-type */
	m_type, 	/*!< massless external, but massive internal particle */

	massive,		/*!< massive */
	zero    /*!< zero */
};

class raw_box;
class raw_triangle;
class raw_bubble;

class boxD;
class triangleD;
class bubbleD;

}

#endif /* PARTITIONS_FWD_H_ */
