/*!\file partitions.h
\brief Header file for part, raw_part and cutD
*/
#ifndef PARTITIONS_H_
#define PARTITIONS_H_

#include <vector>
#include <iostream>
#include "BH_typedefs.h"
#include "BH_error.h"
#include "particles.h"
#include "process.h"
#include "mom_conf.h"
#include "computable.h"
#include "partitions_fwd.h"


namespace BH {

class cutD;
//! Class for partitions
/**
 The class part contains an integer label specifying the partition and a vector of vectors of plabel that represent particles on each corner.
 The integer code is constructd as follows.
 - The first corner is always the one that contains the particle with label 1.
 - the code is formed from the labels of the first momenta in the corners, so for example
       \n{12} {3} {456} -> 134
       \n{561} {2} {3} {4} -> 5234 \n
It is the responsability of the procedure constructing the partition to set the code,
   as it can't be reconstructed afterwards.

The first propagator l(1) is the one connecting the last corner with the first. The i-th propagator l(i) is the one connecting corner i-1 with corner i.
A partition P is called the parent of a partition D if D can be obtained by opening a corner in P. In this case D is called a daughter of P.
part objects can be printed to standard streams.
*/

class part {
	std::vector<std::vector<plabel> > corner;
	int code;
#ifndef SWIG
	friend std::ostream& operator<<(std::ostream& ,const  part& );
#endif
	friend class cutD;
	process _process;
public :
	typedef std::vector<std::vector<plabel> >::iterator iterator;
	typedef std::vector<std::vector<plabel> >::const_iterator const_iterator;
	//! inserts a new particle label in a corner
	/** \param cor corner where the plabel should be added (1-based)
	\param p plabel to be inserted */
	void add(size_t cor,plabel p);
	/** \return number of particles in the cut diagram */
	size_t np() const ;
	/** \return number of corners in the cut diagram */
	size_t nc() const ;
	/** \param i corner (1-based)
	 \return vector containing all the plabels of the corner i */
	int corner_size(int i) const {return corner[i-1].size();}
	int corner_ind(int i,int j) const {return corner[i-1][j-1].ind();}
	const std::vector<plabel>& c(int i) const ;
	//! Sets the corner i
	void set_corner(size_t cor,const std::vector<plabel>& new_cor){corner[cor-1]=new_cor;set_process();} ;
	//! Sets the code
	void set_code(int newcode) ;
	//* \return the integer code for the partition
	int get_code() const;
	//* \return the process the partition is a partition of.
	const process& extern_process() const {return _process;};
	//! default constructor
	part() {};
	//!constructor for i corners
	part(size_t i) : corner(i), code(0) {};
//	part(const part&);
//	part& operator=(const part&);
	virtual ~part(){};
	const_iterator begin() const {return corner.begin();};
	const_iterator end() const {return corner.end();};
	size_t size() const {return corner.size();}
#ifdef SWIG
%rename(__getitem__) operator[];
#endif
	const std::vector<plabel>& operator[](size_t i) const {return corner[i];} ;
private:
	void set_process();
};

std::ostream& operator<<(std::ostream& s, const part& p);

bool operator<(const part& c1, const part& c2);
bool operator==(const part& c1, const part& c2);

std::vector<std::vector<size_t> > splits4(size_t n);
std::vector<std::vector<size_t> > splits2(size_t n);
std::vector<std::vector<size_t> > splits1(size_t n);

//!class for partitions without helicity assignment
/** A raw partition is an object containing particle_type organized in corners. There is no information about the helicity of the external particle or the type and helicity of the propagtors.
The labels are integers starting from 1 to the number of external particles. The first corner is the one containing leg labelled 1. A raw_part object also contains a vector of pointers to cutD objects which is intended to contain the reference to cut diagram that share the raw partition partitioning of the corners. A raw partition P is called the parent of a raw partition D if D can be obtained by opening a corner in P. In this case D is called a daughter of P.   \sa part \sa cutD
raw_part objects can be printed to standard streams.
raw_part objects have == and < operations */
class raw_part {
	std::vector<std::vector<particle*> > _corner;
	std::vector<std::vector<size_t> > _indices;
	std::vector<cutD*> _hel_cuts;
	int _code;
#ifndef SWIG
	friend std::ostream& operator<<(std::ostream& , raw_part );
#endif
public :
	//! adds a label to a corner
	/** \param  cor corner to add the particle and label to. \param pt particle_type \param ind index of the new particle to be inserted */
	void add(size_t cor,particle* pt,size_t ind);
	//! number of particles
	/** \return number of particles in the cut diagram */
	size_t np() const ;
	//! number of corners
	/** \return number of corners in the cut diagram */
	size_t nc() const ;
	//! corner
	/** \param i corner (1-based)
	\return a vector containing all the plabels of the corner i */
	const std::vector<particle*>& c(int i) const;
	//! corner
	/** \param i corner (1-based)
	\return a vector containing all the plabels of the corner i */
	const std::vector<particle*>& corner(int i) const;
	//! vector containing the integer label at the corner i.
	const std::vector<size_t>& ind(int i) const;
	//! returns a pointer to the vector containing the integer label at the corner i.
	const std::vector<size_t>& indices(int i) const;
	//! opens a corner after its nth momentum
	/** \param cor corner to be opened
	  \param nth the corner is opened after the nth momenentum, so that we have 1<=nth<cor
	  \param jump is boolean that is passed by reference and that is set true if a jump has happend, false if not.
	*/
	raw_part open(size_t cor, size_t nth,bool& jump);
	//! Sets the code
	void set_code(int newcode);
	//* \return the integer code for the partition
	int get_code() const;
	//! default constructor
	raw_part() {};
	//!constructor for i corners
	raw_part(size_t i) : _corner(i), _indices(i), _code(0) {};
	//! adds a new pointer to an helicity cut
	//* For internal use (for example in the construction of OneLoopHelAmpl(OneLoopAmplitude))*/
	void add_hel_cut(cutD* cd);
	//! number of pointers to cutD objects corresponding to the raw partition
	size_t nbr_hel_cuts();
	//! i th pointer to cutD
	/** \return poiner to the i th cutD object */
	cutD* hel_cut(size_t i);
	virtual ~raw_part(){};
};


//! class for raw box partitions
/** \sa raw_part*/
class raw_box : public raw_part {
	std::vector<raw_triangle* > parents;
	std::vector<size_t> closed_corner;
public:
	//! adds a new parent
	/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(raw_triangle* p){parents.push_back(p);};
	//! adds a new parent and remembers which corner of the box is the left part of the triangle corner that has been opened.
	/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(raw_triangle* p,size_t closed){parents.push_back(p);closed_corner.push_back(closed);};
	//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	raw_triangle* get_parent(size_t i){return parents[i-1];}
	//! corner number of the opened corner
	/** the corner number is that of the corner that is left of the new propagator */
	size_t get_closed_corner(size_t i){return closed_corner[i-1];}
	//! number of parents
	size_t parents_nbr(){return parents.size();};
	//! constructor
	raw_box(raw_part rp): raw_part(rp){};
	~raw_box(){};
};

//! class for raw triangle partitions
/** \sa raw_part*/
class raw_triangle : public raw_part {
	std::vector<raw_box*> daughters;
	std::vector<raw_bubble* > parents;
	std::vector<size_t> opened_corner;
	std::vector<size_t> closed_corner;
public:
	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	//! adds a new daughter
	void add_daughter(raw_box* p){daughters.push_back(p);};
	//! adds a new daughter and remembers which corner has been opened
	/** adds a new parent and remembers which triangle corner has been opened to obtain the bubble.  */
	void add_daughter(raw_box* p,size_t opened){daughters.push_back(p);opened_corner.push_back(opened);};
	//! adds a new parent
	void add_parent(raw_bubble* p){parents.push_back(p);};
	//! adds a new parent and remembers which corner has been openend
	/** adds a new parent and remembers which corner of the triangle is the left part of the bubble corner that has been opened.  */
	void add_parent(raw_bubble* p,size_t closed){parents.push_back(p);closed_corner.push_back(closed);};
	//! pointer to a daughter
	/** \param i (1-based) parent label \return pointer to the ith daughter */
	raw_box* get_daughter(size_t i){return daughters[i-1];}
	//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	raw_bubble* get_parent(size_t i){return parents[i-1];}
	//! corner number of the opened corner for box(i)
	/** the corner number is that of the corner that is left of the new propagator in box(i)   */
	size_t get_opened_corner(size_t i){return opened_corner[i-1];}
	//! corner number of the corner that has to be closed
	/** the corner number is that of the corner left of the propagator that has to be collapsed to obtain the bubble i */
	size_t get_closed_corner(size_t i){return closed_corner[i-1];}
	//! number of daughters
	size_t daughters_nbr(){return daughters.size();};
	//! number of parents
	size_t parents_nbr(){return parents.size();};
	//! constructor
	raw_triangle(raw_part rp): raw_part(rp){};
	~raw_triangle(){};
};

//! class for raw bubble partitions
/** \sa raw_part*/
class raw_bubble : public raw_part {
	std::vector<raw_triangle*> daughters;
	std::vector<size_t> opened_corner;
public:
	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	//! adds a new daughter
	void add_daughter(raw_triangle* p){daughters.push_back(p);};
	//! adds a new daughter and remembers which corner has been opened
	/** adds a new parent and remembers which triangle corner has been opened to obtain the bubble.  */
	void add_daughter(raw_triangle* p,size_t opened){daughters.push_back(p);opened_corner.push_back(opened);};
	//! pointer to a daughter
	/** \param i (1-based) daughter label \return pointer to the ith parent */
	raw_triangle* get_daughter(size_t i){return daughters[i-1];}
	//! corner number of the opened corner for triangle(i)
	/** the corner number is that of the corner that is left of the new propagator in triangle(i)   */
	size_t get_opened_corner(size_t i){return opened_corner[i-1];}
	//! number of daughters
	size_t daughters_nbr(){return daughters.size();};
	//!constructor
	raw_bubble(raw_part rp): raw_part(rp){};
	~raw_bubble(){};
};

std::ostream& operator<<(std::ostream& s, raw_part p);

bool operator<(const raw_part& c1, const raw_part& c2);
bool operator==(const raw_part& c1, const raw_part& c2);


class symmetry_factor {
	int d_num;
	int d_den;
public:
	symmetry_factor(int num=1,int den=1): d_num(num), d_den(den) {};
	int get_num() const {return d_num;}
	int get_den() const {return d_den;}
	template <class T> T eval() const {return T(d_num)/T(d_den);}
};


//! Class for cut diagrams
/** A cutD object contains a partition, and a set of propagators

 */
class cutD : public part {
//	part pa;
	std::vector<particle_ID> lm;
	std::vector<process> process_list;
	std::vector<long> _process_code;
	raw_part* _raw_partition;

	symmetry_factor _symmetry_factor;

	long int _cutD_ID;
	static long int cutD_next_ID;

public:
	//! loop particle
	/** \param i 1-based label of the loop particle. l(2) returns the particle leaving corner 1 and entering corner 2
	 \return ptype of the loop particle
	  */
	particle_ID l(int i) const { return lm[i-1]; };
#ifndef SWIG
	friend std::ostream& operator<<(std::ostream& , const  cutD& );
#endif
	//	cutD open(size_t corner,size_t n_th,ph_type t); //delete that one when finished
//! opening	a corner
	/**
	 \param corner (1-based) label of the corner to open
	 \param n_th particle in the corner after which to open the new propagator
	 \param t particle type of the new propagator (in ph_type).
	 \param jump boolean passed by reference that will be set according to whether or not a jump occured
	\return cutD with the opened propagator
	a "jump" happens when the corner opened is the first corner and it is opened such that leg one is on the right side of the new propagator. Since the corner containing the leg one is always corner 1, the right side of the new propagator is corner one and the left one is the last corner. In all other cases the left corner has the index of the opened corner and the right part has the next label.
  */
	cutD open(size_t corner,size_t n_th,const particle_ID& t,bool &jump);
//! constructor for tadpoles
	cutD(const part& pp,const std::vector<particle_ID>& ls);
	//! constructor for bubbles
	cutD(const part& ,const particle_ID&,const particle_ID&);
	//! constructor for triangles
	cutD(const part& PA,const particle_ID&,const particle_ID&,const particle_ID&);
	//! constructor for boxes
	cutD(const part& PA,const particle_ID&,const particle_ID&,const particle_ID&,const particle_ID&);

	/** \return whether the cut is (obviously) vanishing  */
	bool is_zero();
	//! the corner_type of the corner i
	/** A corner i of type
	  \n - "a" if the corner is massless and it is expressed only in terms of < > products
	  \n - "b" if the corner is massless and it is expressed only in terms of [ ] products
	  \n - "massive" if the corner is massive
	*/
	corner_type c_type(size_t i) const ;
	//! Process associated to the specified corner
	/**
	 The processes associated with the cut diagrams are not constructed by the constructor, but when one is requested with cutD::get_process, the processes for all corners are generated and saved in the cutD object for later use.
	  \param cor corner number (1-based) \return the process associated with the corner cor. The propagator to the right of the corner will be conjugated. */
	const process& get_process (size_t cor) const;
	//! pcode of the process associated to the specified corner
	/**
	 The processes associated with the cut diagrams are not constructed by the constructor, but when one is requested with cutD::get_process, the processes for all corners are generated and saved in the cutD object for later use.
	  \param cor corner number (1-based) \return the process associated with the corner cor. The propagator to the right of the corner will be conjugated. */
	long get_process_code(size_t cor) const;
	//! gets the raw partition of the cutD
	/** It is not the responsability of the constructor to set the link to the raw partition since there is not always a raw partition to link to. \sa set_raw_partition*/
	raw_part* get_raw_partition();
	//! sets the raw partition of the cutD
	void set_raw_partition(raw_part*);
	size_t nbr_props(){return lm.size();};
	virtual ~cutD(){};
	long int get_ID(){return _cutD_ID;}
	template <class T> T get_symmetry_factor(){return _symmetry_factor.eval<T>();}
	const symmetry_factor& get_symmetry_factor_no_eval() const {return _symmetry_factor;};
private:
	void set_symmetry_factor();
	void construct_processes();

};


//! Class for pentagon cut diagrams
class pentagonD : public cutD, public computable<std::complex> {
protected:
	std::vector<boxD*> parents;
	std::vector<size_t> closed_corner;

public:
	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(boxD* p){parents.push_back(p);};
	void add_parent(boxD* p,size_t closed){parents.push_back(p);closed_corner.push_back(closed);};
	//! clears the vector of parents
	void clear_parents(){parents.clear();};
	//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	boxD* get_parent(size_t i){return parents[i-1];}
	size_t get_closed_corner(size_t i){return closed_corner[i-1];}
	void set_closed_corner(size_t i,size_t cc){closed_corner[i-1]=cc;}
	size_t parents_nbr(){return parents.size();};
	pentagonD(cutD cd): cutD(cd){};
	//! Evaluation of pentagon coefficients
	/** \param mc momentum configuration to be used \param ind indices of the momenta to be used \returns pentagon coefficient*/
	virtual C eval(mom_conf& mc,const std::vector<int>&);
	virtual CHP eval(mom_conf_HP&,const std::vector<int>&);
	virtual CVHP eval(mom_conf_VHP&,const std::vector<int>&);
	virtual C eval(const eval_param<R>&);
	virtual CHP eval(const eval_param<RHP>&);
	virtual CVHP eval(const eval_param<RVHP>&);
#ifdef BH_USE_GMP
	virtual CGMP eval(const eval_param<RGMP>&){throw BHerror("Not implemented");};
#endif

	virtual ~pentagonD(){};
};



//! Class for box cut diagrams
class boxD : public cutD, public computable<std::complex> {
protected:
	mutable std::vector<triangleD*> parents;
	std::vector<cutD* > daughters;
	mutable std::vector<size_t> closed_corner;
//	long int _last_ID;
//	long int _last_ID_HP;
//	long int _last_ID_VHP;
//	C _last_value;
//	CHP _last_value_HP;
//	CVHP _last_value_VHP;
//	std::vector<int> _last_ind;
//	std::vector<int> _last_ind_HP;
//	std::vector<int> _last_ind_VHP;
public:
	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(triangleD* p) const {parents.push_back(p);};
	void add_parent(triangleD* p,size_t closed) const {parents.push_back(p);closed_corner.push_back(closed);};
	//! clears the vector of parents
	void clear_parents(){parents.clear();};
	//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	triangleD* get_parent(size_t i) const {return parents[i-1];}
	size_t get_closed_corner(size_t i) const {return closed_corner[i-1];}
	void set_closed_corner(size_t i,size_t cc){closed_corner[i-1]=cc;}
	size_t parents_nbr() const {return parents.size();};
	boxD(const cutD& cd): cutD(cd){};
	//! Evaluation of box coefficients
	/** \param mc momentum configuration to be used \param ind indices of the momenta to be used \returns box coefficient*/
//	C get_coeff(mom_conf& mc,std::vector<int> ind);
//	CHP get_coeff(mom_conf_HP& mc,std::vector<int> ind);
//	CVHP get_coeff(mom_conf_VHP& mc,std::vector<int> ind);
	virtual C eval(mom_conf& mc,const std::vector<int>&);
	virtual CHP eval(mom_conf_HP&,const std::vector<int>&);
	virtual CVHP eval(mom_conf_VHP&,const std::vector<int>&);
#ifdef BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){throw BHerror("Not implemented");};
#endif
	virtual C eval(const eval_param<R>&);
	virtual CHP eval(const eval_param<RHP>&);
	virtual CVHP eval(const eval_param<RVHP>&);
#ifdef BH_USE_GMP
	virtual CGMP eval(const eval_param<RGMP>&){throw  BHerror("Not implemented");};
#endif

	virtual ~boxD(){};
};

//! Class for triangle cut diagrams
class triangleD : public cutD , public computable<std::complex> {
protected:
	std::vector<bubbleD*> parents;
	std::vector<boxD* > daughters;
	std::vector<size_t> opened_corner;
	std::vector<size_t> closed_corner;
//	long int _last_ID;
//	long int _last_ID_HP;
//	long int _last_ID_VHP;
//	C _last_value;
//	CHP _last_value_HP;
//	CVHP _last_value_VHP;
//	std::vector<int> _last_ind;
//	std::vector<int> _last_ind_HP;
//	std::vector<int> _last_ind_VHP;

public:
	//! adds a new parent
/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_parent(bubbleD* p){parents.push_back(p);};
	void add_parent(bubbleD* p,size_t closed){parents.push_back(p);closed_corner.push_back(closed);};
	//! clears the vector of parents
	void clear_parents(){parents.clear();};
	//! adds a new daughter
	/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_daughter(boxD* p){daughters.push_back(p);};
	void add_daughter(boxD* p,size_t opened){daughters.push_back(p);opened_corner.push_back(opened);};
	//! clears the vector of daughters
	void clear_daughters(){daughters.clear();};
	//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	bubbleD* get_parent(size_t i){return parents[i-1];}
	//! pointer to a daughter
	/** \param i (1-based) daughter label \return pointer to the ith daughter */
	boxD* get_daughter(size_t i){return daughters[i-1];}
	//! number of parents
	size_t parents_nbr(){return parents.size();};
	//! number of daughters
	size_t daughters_nbr(){return daughters.size();};
	//! coefficient
//	/** \param mc momentum configuration \param ind index vector \return The coefficient of the corresponding cut. The value is cached for later use. */
//	C get_coeff(mom_conf& mc,std::vector<int> ind);
//	//! coefficient in high precision
//	/** \param mc momentum configuration \param ind index vector \return The coefficient of the corresponding cut. The value is cached for later use. */
//	CHP get_coeff(mom_conf_HP& mc,std::vector<int> ind);
//	//! coefficient in high precision
//	/** \param mc momentum configuration \param ind index vector \return The coefficient of the corresponding cut. The value is cached for later use. */
//	CVHP get_coeff(mom_conf_VHP& mc,std::vector<int> ind);
//	//! coefficient evaluation
//	/** \param mc momentum configuration \param ind index vector \return The coefficient of the corresponding cut.  The value is only evaluated and not cached. \sa get_coeff */
	virtual C eval(mom_conf&,const std::vector<int>&);
	//! coefficient evaluation for high precision
	/** \param mc momentum configuration \param ind index vector \return The coefficient of the corresponding cut. The value is only evaluated and not cached. \sa get_coeff */
	virtual CHP eval(mom_conf_HP&,const std::vector<int>&);
	//! coefficient evaluation for very high precision
	/** \param mc momentum configuration \param ind index vector \return The coefficient of the corresponding cut. The value is only evaluated and not cached. \sa get_coeff */
	virtual CVHP eval(mom_conf_VHP&,const std::vector<int>&);
#ifdef BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){throw  BHerror("Not implemented");};
#endif
	virtual C eval(const eval_param<R>&);
	virtual CHP eval(const eval_param<RHP>&);
	virtual CVHP eval(const eval_param<RVHP>&);
#ifdef BH_USE_GMP
	virtual CGMP eval(const eval_param<RGMP>&){throw BHerror("Not implemented");};
#endif

	size_t get_opened_corner(size_t i){return opened_corner[i-1];}
	size_t get_closed_corner(size_t i){return closed_corner[i-1];}
	void set_opened_corner(size_t i,size_t oc){opened_corner[i-1]=oc;}
	void set_closed_corner(size_t i,size_t cc){closed_corner[i-1]=cc;}
	triangleD(const cutD& cd): cutD(cd){};
	virtual ~triangleD(){};
};

//! Class for bubble cut diagrams
class bubbleD : public cutD, public computable<std::complex> {
protected:
	std::vector<triangleD* > daughters;
	std::vector<size_t> opened_corner;
//	long int _last_ID;
//	long int _last_ID_HP;
//	long int _last_ID_VHP;
//	C _last_value;
//	CHP _last_value_HP;
//	CVHP _last_value_VHP;
//	std::vector<int> _last_ind;
//	std::vector<int> _last_ind_HP;
//	std::vector<int> _last_ind_VHP;

	double _accuracy; // A variable for storing the numerical accuracy of the computation

	// We store a number that gives the accuracy of the computation
	void set_accuracy(double acc){_accuracy=acc;};

public:
	//! adds a new daughter
	/** should be taken care of by the constructor of the amplitude. Since it is a pointer, take has to be taken that the vector is not changed after the pointer is put in.*/
	void add_daughter(triangleD* p){daughters.push_back(p);};
	void add_daughter(triangleD* p,size_t opened){daughters.push_back(p);opened_corner.push_back(opened);};
	//! clears the vector of daughters
	void clear_daughters(){daughters.clear();};
//! pointer to a parent
	/** \param i (1-based) parent label \return pointer to the ith parent */
	triangleD* get_daughter(size_t i){return daughters[i-1];}
	//! number of parents
	size_t daughters_nbr(){return daughters.size();};
	size_t get_opened_corner(size_t i){return opened_corner[i-1];}
	void set_opened_corner(size_t i,size_t oc){opened_corner[i-1]=oc;}
	bubbleD(const cutD& cd): cutD(cd), _accuracy(64) {}; // Set the default accuracy to 64, it will never be worse than this
//	C get_coeff(mom_conf& mc,std::vector<int> ind);
//	CHP get_coeff(mom_conf_HP& mc,std::vector<int> ind);
//	CVHP get_coeff(mom_conf_VHP& mc,std::vector<int> ind);
	virtual C eval(mom_conf&,const std::vector<int>&);
	virtual CHP eval(mom_conf_HP&,const std::vector<int>&);
	virtual CVHP eval(mom_conf_VHP&,const std::vector<int>&);
#ifdef BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){throw BHerror("Not implemented");};
#endif
	virtual C eval(const eval_param<R>&);
	virtual CHP eval(const eval_param<RHP>&);
	virtual CVHP eval(const eval_param<RVHP>&);
#ifdef BH_USE_GMP
	virtual CGMP eval(const eval_param<RGMP>&){throw BHerror("Not implemented");};
#endif

	virtual ~bubbleD(){};

	// We store a number that gives the accuracy of the computation, this retrieves it
	double get_accuracy(){return _accuracy;};
};

// Output routines
std::ostream& operator<<(std::ostream& s, const pentagonD& p);
std::ostream& operator<<(std::ostream& s, const boxD& p);
std::ostream& operator<<(std::ostream& s, const triangleD& p);
std::ostream& operator<<(std::ostream& s, const bubbleD& p);

bool operator<(const cutD& c1, const cutD& c2);
bool operator==(const cutD& c1, const cutD& c2);

//! abstract factory class
class cutD_factory {
public:
	virtual boxD* new_box(const boxD&)=0;
	virtual triangleD* new_triangle(const triangleD&)=0;
	virtual bubbleD* new_bubble(const bubbleD&)=0;
	static cutD_factory* default_CF;
	virtual ~cutD_factory(){};
};

}

#include "partitions_inline.hpp"

#endif /*PARTITIONS_H_*/
