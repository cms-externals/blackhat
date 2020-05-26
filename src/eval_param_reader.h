/*
 *  eval_param_reader.h
 *  BlackHat
 *
 *  Created by Darren Forde on 25/06/2010.
 *  Copyright 2010 BlackHat Collaboration. All rights reserved.
 *
 */

#ifndef EVAL_PARAM_READER_H_
#define EVAL_PARAM_READER_H_

#include "eval_param.h"

namespace BH {
	
//! Class for reading momentum configurations from a file. This allows to treat stepping through a file independently from the type (precision) of the points read.
/** An eval_param object is constructed by giving it the name of the file to read and the number of momenta to read in for each new momentum configuration. When constructed, a new momentum configuration is read from the file by using the next() member function. */
template <class T> class eval_param_reader : public eval_param<T> {
protected:
	std::ifstream input;
	size_t nth;
	size_t nbr_particles;
	std::ios::pos_type start_pos;
	
	std::vector<Cmom<T>* > _momenta;
public :
	//! constructor
	eval_param_reader(const char* path,int nbr_p);
	~eval_param_reader();
	//! reads the next configuration from the file
	/** \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read) */
	virtual bool next();
	//! reads the n th configuration from the file
	/** \param n One based index, for  n=1 go_to reads the first entry \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read) This member function reads through the file until it arrives at the required position and is therefore slow. If the position where the momentum configuration to be read is known, one should use go_to_pos instead. \sa go_to_pos*/
	virtual bool go_to(size_t n);
	//! reads the n th configuration from the file knowing the position from which to start reading
	/** This function is more efficient than go_to \sa go_to \param pos is the position in the file of the momentum to be read. \param n One based index, for  n=1 go_to reads the first entry. This is necessary so that the eval_param reader knows where it is in the file, in case the go_to member is called. \return true if the reading has been successful, false if not (for example when not all momentum configurations have been read)  */
	virtual bool go_to_pos(std::ios::pos_type pos,size_t n);
};

}

#endif /* EVAL_PARAM_READER_H_ */


