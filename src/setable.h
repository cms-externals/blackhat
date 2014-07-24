/*
 * setable.h
 *
 *  Created on: 30 Mar 2010
 *      Author: daniel
 */

#ifndef SETABLE_H_
#define SETABLE_H_

namespace BH {

class setable {
public:
	virtual bool set(const std::string& name,double)=0;
	virtual bool set(const std::string& name,int)=0;
	virtual bool set(const std::string& name,bool)=0;
	virtual bool set(const std::string& name,std::string)=0;
	virtual ~setable(){}
};
} /* BH */

#endif /* SETABLE_H_ */
