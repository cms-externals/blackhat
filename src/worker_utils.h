/*
 * worker_utils.h
 *
 *  Created on: 31 Jul 2009
 *      Author: daniel
 */

#ifndef WORKER_UTILS_H_
#define WORKER_UTILS_H_

#include <iosfwd>

namespace BH {

class process;

namespace worker {

void write(const process& PRO,std::ostream& os);
bool read_process_from_stream(process& ret,std::istream& is);
process read_process_from_stream(std::istream& is);

} /* worker */

} /* BH */

#endif /* WORKER_UTILS_H_ */
