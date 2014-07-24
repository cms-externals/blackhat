/*
 * timing.h
 *
 *  Created on: 6 Nov 2009
 *      Author: daniel
 */

#ifndef TIMING_H_
#define TIMING_H_


#include <string>
#include <fstream>


#ifdef BH_TIMING_ON
#define BH_START_TIMER(NAME) static BH::timing::counter c_##NAME (#NAME); c_##NAME.start();
#define BH_STOP_TIMER(NAME)  c_##NAME.stop();
#else
#define BH_START_TIMER(NAME)
#define BH_STOP_TIMER(NAME)
#endif

namespace BH {

namespace timing {

class counter {
	clock_t d_start;
	long d_nbr_calls;
	long d_total_time;
	std::string d_name;
	static std::ofstream ofile;
	static bool s_is_open;
	bool is_open(){ if (! s_is_open){ ofile.open("BHstatistics");s_is_open=true;}  return s_is_open;};
public:
	counter(const std::string& name): d_name(name),d_nbr_calls(0),d_total_time(0) {}
	void add(long t){ ++d_nbr_calls; d_total_time+=t; };
	void start(){d_start=clock();};
	void stop(){ add(clock()-d_start);}
	// upon destruction an entry in the statistics file is made.
	~counter();
};

} /* timing */
} /* BH */

#endif /* TIMING_H_ */
