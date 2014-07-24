/*
 * mode_dependent_typedefs.h
 *
 *  Created on: 24 Jul 2009
 *      Author: daniel
 */

#ifndef MODE_DEPENDENT_TYPEDEFS_H_
#define MODE_DEPENDENT_TYPEDEFS_H_

namespace BH {

#ifndef BH_PUBLIC
class Rec_Tree;
class Tree_factory;

typedef Tree_factory TREE_FACTORY_TYPE;
typedef Rec_Tree TREE_TYPE;
#else
class worker_tree;
class worker_tree_factory;

typedef worker_tree_factory TREE_FACTORY_TYPE;
typedef worker_tree TREE_TYPE;
#endif

}

#endif /* MODE_DEPENDENT_TYPEDEFS_H_ */
