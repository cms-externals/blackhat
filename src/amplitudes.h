/*!\file amplitudes.h
\brief Header file for amplitudes: OneLoopHelAmpl, OneLoopAmplitude
*/
#ifndef AMPLITUDES_H_
#define AMPLITUDES_H_

#include "partitions.h"
#include "particles.h"



#include "Series.h"
#ifndef BH_PUBLIC
#include "options.h"
#include "ordering.h"
#endif
#include "tree_amplitudes.h"

namespace BH {


class Amplitude;
class OneLoopAmplitude;
struct ordering_constraint;


//! class for one-loop helicity amplitudes
/** The objects of this class contain three std::vectors of cutD objects to represent all the triangle, bubble and box diagrams. The constructor for the one-loop amplitudes generate all non-vanishing bubbles, triangles and boxes configurations.*/


std::ostream& operator<<(std::ostream& s, HelAmpl& p);
//! class for amplitudes
/** The objects of this class contain three vectors of raw_part objects to represent all the triangle, bubble and box partitions. The constructor for the one-loop amplitudes generates all helicity amplitudes for the external particles.*/
class Amplitude {
	std::vector<particle*> _particles;
public:
//	virtual C eval(mom_conf);
//	virtual void print();
//	virtual bool is_zero();
	const std::vector<particle*>& particles(){return _particles;};
	Amplitude(std::vector<particle*> p) : _particles(p){};
	Amplitude() : _particles(0) {};
	virtual ~Amplitude(){};
};
//! class for one-loop helicity amplitudes
/** The objects of this class contain three vectors of cutD objects to represent all the triangle, bubble and box diagrams. The constructor for the one-loop amplitudes generate all non-vanishing bubbles, triangles and boxes configurations.*/


std::vector<part> all_bubble_partitions(const process& p);

std::vector<cutD> all_bubble_cutD(part P);


class OneLoopRawAmplitude : public Amplitude {
	std::vector<raw_bubble> r_bubbles;
	std::vector<raw_triangle> r_triangles;
	std::vector<raw_box> r_boxes;
	std::map<int,raw_part*> code_to_rpart;
public :
	//! constructor
	OneLoopRawAmplitude(std::vector<particle*> p);
	//! constructor
	OneLoopRawAmplitude(const std::vector<particle*>& p,const ordering_constraint& oc);
	/** \param i (1-based) index \return raw_box object representing the ith box diagram */
	raw_box* box(size_t i) ;
	//! triangle cut diagram
	/** \param i (1-based) index \return raw_triangle object representing the ith triangle diagram */
	raw_triangle* triangle(size_t i);
	//! bubble cut diagram
	/** \param i (1-based) index \return raw_bubble object representing the ith bubble diagram */
	raw_bubble* bubble(size_t i);
	/** \return number of box raw partitions */
	size_t nbr_boxes() const ;
	/** \return number of triangle raw partitions */
	size_t nbr_triangles() const ;
	/** \return number of bubble raw partitions */
	size_t nbr_bubbles() const ;
	//! raw partition from code
	/** \param code integer code \return pointer to the corresponding raw partition*/
	raw_part* rpart_from_code(int code);
	//! helicity amplitude i
	/** \param i index of the helicity amplitude \return pointer to the helicity amplitude. */
	virtual ~OneLoopRawAmplitude(){};
};


}
#endif /*AMPLITUDES_H_*/
