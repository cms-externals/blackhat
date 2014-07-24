/*!\file OneLoopHelAmpl.h
\brief header for the OneLoopHelAmpl objects
*/
#ifndef ONELOOPHELAMPL_H_
#define ONELOOPHELAMPL_H_

#include "amplitudes.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "cut_part.h"
#include "BH_utilities.h"
#include "scheme.h"
#include "cut_part_factory.h"
#include "rational.h"
#include "rational_factory.h"
#include "mode_dependent_typedefs.h"  // for TREE_TYPE and TREE_FACTORY_TYPE
#include "BH_debug.h"
#ifdef BH_PUBLIC
#include "worker_tree.h"
#endif
namespace BH {


// forward declaration
template <class rational_type> class Rational_factory;
template <class cut_part_type> class cut_part_factory;
class Rec_Tree;

class OneLoopAmplitude_base : public HelAmpl , public computable<SeriesC>  {
	color_structure d_cs;
public:
	OneLoopAmplitude_base(const process& p,color_structure cs) : HelAmpl(p) , d_cs(cs){};
	//! returns the cut part
	virtual std::complex<R> get_rational(momentum_configuration<R>& mc, const std::vector<int>& ind)=0;
	virtual std::complex<RHP> get_rational(momentum_configuration<RHP>& mc, const std::vector<int>& ind)=0;
	virtual std::complex<RVHP> get_rational(momentum_configuration<RVHP>& mc, const std::vector<int>& ind)=0;
	//! returns the rational part
	virtual SeriesC<R> get_cut(momentum_configuration<R>& mc, const std::vector<int>& ind)=0;
	virtual SeriesC<RHP> get_cut(momentum_configuration<RHP>& mc, const std::vector<int>& ind)=0;
	virtual SeriesC<RVHP> get_cut(momentum_configuration<RVHP>& mc, const std::vector<int>& ind)=0;
	//! returns the tree
	virtual std::complex<R> get_tree(momentum_configuration<R>& mc, const std::vector<int>& ind)=0;
	virtual std::complex<RHP> get_tree(momentum_configuration<RHP>& mc, const std::vector<int>& ind)=0;
	virtual std::complex<RVHP> get_tree(momentum_configuration<RVHP>& mc, const std::vector<int>& ind)=0;

	virtual double get_accuracy()=0;
    virtual SeriesC<R> get_conjugate_amplitude()=0;
    virtual SeriesC<RHP> get_conjugate_amplitude_HP()=0;
    virtual SeriesC<RVHP> get_conjugate_amplitude_VHP()=0;
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> get_rational(momentum_configuration<RGMP>& mc, const std::vector<int>& ind)=0;
	virtual SeriesC<RGMP> get_cut(momentum_configuration<RGMP>& mc, const std::vector<int>& ind)=0;
	virtual std::complex<RGMP> get_tree(momentum_configuration<RGMP>& mc, const std::vector<int>& ind)=0;
    virtual SeriesC<RGMP> get_conjugate_amplitude_GMP()=0;
#endif
	
    virtual void set_mu(int index)=0;
	virtual void set_mu_HP(int index)=0;
	virtual void set_mu_VHP(int index)=0;
#ifdef BH_USE_GMP
	virtual void set_mu_GMP(int index)=0;
	virtual void set_mu(R val,const RHP& val_HP,const RVHP& val_VHP,const RGMP& val_GMP)=0;
#else
	virtual void set_mu(R val,RHP val_HP,RVHP val_VHP)=0;
#endif
	virtual void set_mu(const multi_precision_constant& multi)=0;
	virtual void set_scheme(scheme sc)=0;

	virtual void dry_run(const std::vector<int>&){};

	virtual ~OneLoopAmplitude_base(){};
	color_structure color_struct(){return d_cs;};
};


//! class for one loop helicity amplitudes
/** \sa Cut_Part OneLoopHelAmpl objects contain a cut part (Cut_Part object) and a rational part (Rec_Rational object) */
class One_Loop_Helicity_Amplitude : public OneLoopAmplitude_base {
private: 
	scheme _scheme;
	
protected:
	TREE_TYPE* d_tree_ptr;
	Rational_base* _rational_part;
	Cut_Part_base* _cut_part;
	
	double _accuracy; //Store the estimated accuracy
    	SeriesC<R>_conjugate_cut_part;
    	SeriesC<RHP>_conjugate_cut_part_HP;
    	SeriesC<RVHP>_conjugate_cut_part_VHP;
    	std::complex<R>_conjugate_rat_part;
    	std::complex<RHP>_conjugate_rat_part_HP;
    	std::complex<RVHP>_conjugate_rat_part_VHP;
#ifdef BH_USE_GMP
    	SeriesC<RGMP>_conjugate_cut_part_GMP;
    	std::complex<RGMP>_conjugate_rat_part_GMP;
#endif
public :
	SeriesC<R> eval(momentum_configuration<double>& mc,const std::vector<int>& ind);
	SeriesC<RHP> eval(momentum_configuration<dd_real>& mc,const std::vector<int>& ind);
	SeriesC<RVHP> eval(momentum_configuration<qd_real>& mc,const std::vector<int>& ind);
#ifdef BH_USE_GMP
	SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
#endif

	SeriesC<R> eval(const eval_param<R>&);
	SeriesC<RHP> eval(const eval_param<RHP>&);
	SeriesC<RVHP> eval(const eval_param<RVHP>&);
#ifdef BH_USE_GMP
	SeriesC<RGMP> eval(const eval_param<RGMP>&);
#endif

	//! pointer to the tree
	TREE_TYPE* tree(){return d_tree_ptr;};
	//! pointer to the cut part
	Cut_Part_base* cut_part(){return _cut_part;};
//! pointer to the rational part
	Rational_base* rational_part(){return _rational_part;};

	//! returns the tree
	virtual std::complex<R> get_tree(momentum_configuration<R>& mc, const std::vector<int>& ind){return d_tree_ptr->get_value(mc,ind);};
	virtual std::complex<RHP> get_tree(momentum_configuration<RHP>& mc, const std::vector<int>& ind){return d_tree_ptr->get_value(mc,ind);};
	virtual std::complex<RVHP> get_tree(momentum_configuration<RVHP>& mc, const std::vector<int>& ind){return d_tree_ptr->get_value(mc,ind);};
	//! returns the cut part
	virtual std::complex<R> get_rational(momentum_configuration<R>& mc, const std::vector<int>& ind){std::complex<R> ret=_rational_part->get_value(mc,ind); _accuracy=_rational_part->get_accuracy();return ret;};
	virtual std::complex<RHP> get_rational(momentum_configuration<RHP>& mc, const std::vector<int>& ind){std::complex<RHP> ret=_rational_part->get_value(mc,ind); _accuracy=_rational_part->get_accuracy();return ret;};
	virtual std::complex<RVHP> get_rational(momentum_configuration<RVHP>& mc, const std::vector<int>& ind){std::complex<RVHP> ret=_rational_part->get_value(mc,ind); _accuracy=_rational_part->get_accuracy();return ret;};
	//! returns the rational part
	virtual SeriesC<R> get_cut(momentum_configuration<R>& mc, const std::vector<int>& ind){SeriesC<R> ret=_cut_part->get_value(mc,ind); _accuracy=_cut_part->get_accuracy();return ret;};
	virtual SeriesC<RHP> get_cut(momentum_configuration<RHP>& mc, const std::vector<int>& ind){SeriesC<RHP> ret=_cut_part->get_value(mc,ind); _accuracy=_cut_part->get_accuracy();return ret;};
	virtual SeriesC<RVHP> get_cut(momentum_configuration<RVHP>& mc, const std::vector<int>& ind){SeriesC<RVHP> ret=_cut_part->get_value(mc,ind); _accuracy=_cut_part->get_accuracy();return ret;};
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> get_tree(momentum_configuration<RGMP>& mc, const std::vector<int>& ind){return d_tree_ptr->get_value(mc,ind);};
	virtual SeriesC<RGMP> get_cut(momentum_configuration<RGMP>& mc, const std::vector<int>& ind){SeriesC<RGMP> ret=_cut_part->get_value(mc,ind); _accuracy=_cut_part->get_accuracy();return ret;};
	virtual std::complex<RGMP> get_rational(momentum_configuration<RGMP>& mc, const std::vector<int>& ind){std::complex<RGMP> ret=_rational_part->get_value(mc,ind); _accuracy=_rational_part->get_accuracy();return ret;};
#endif
	
//! returns the tree
	virtual std::complex<R> get_tree(const eval_param<R>& ep){return d_tree_ptr->get_value(ep);};
	virtual std::complex<RHP> get_tree(const eval_param<RHP>& ep){return d_tree_ptr->get_value(ep);};
	virtual std::complex<RVHP> get_tree(const eval_param<RVHP>& ep){return d_tree_ptr->get_value(ep);};
	//! returns the cut part
	virtual std::complex<R> get_rational(const eval_param<R>& ep){std::complex<R> ret=_rational_part->get_value(ep); 
			_accuracy=_rational_part->get_accuracy();
        		_conjugate_rat_part=conj(ret); 
			return ret;};
	virtual std::complex<RHP> get_rational(const eval_param<RHP>& ep){std::complex<RHP> ret=_rational_part->get_value(ep); 
			_accuracy=_rational_part->get_accuracy();
        		//_conjugate_rat_part=conj(to_double(ret)); 
        		_conjugate_rat_part_HP=conj(ret); 
			return ret;};
	virtual std::complex<RVHP> get_rational(const eval_param<RVHP>& ep){std::complex<RVHP> ret=_rational_part->get_value(ep); 
			_accuracy=_rational_part->get_accuracy();
        		_conjugate_rat_part_VHP=conj(ret); 
			return ret;};
	//! returns the rational part
	virtual SeriesC<R> get_cut(const eval_param<R>& ep){SeriesC<R> ret=_cut_part->get_value(ep);_accuracy=_cut_part->get_accuracy();
        _conjugate_cut_part=_cut_part->get_conjugate_cut_part();   return ret;};
	virtual SeriesC<RHP> get_cut(const eval_param<RHP>& ep){SeriesC<RHP> ret=_cut_part->get_value(ep); _accuracy=_cut_part->get_accuracy();
        _conjugate_cut_part_HP=_cut_part->get_conjugate_cut_part_HP();   return ret;};
	virtual SeriesC<RVHP> get_cut(const eval_param<RVHP>& ep){SeriesC<RVHP> ret=_cut_part->get_value(ep); _accuracy=_cut_part->get_accuracy();
        _conjugate_cut_part_VHP=_cut_part->get_conjugate_cut_part_VHP();   return ret;};

#ifdef BH_USE_GMP
	virtual std::complex<RGMP> get_tree(const eval_param<RGMP>& ep){return d_tree_ptr->get_value(ep);};
	virtual std::complex<RGMP> get_rational(const eval_param<RGMP>& ep){std::complex<RGMP> ret=_rational_part->get_value(ep); _accuracy=_rational_part->get_accuracy();
        _conjugate_rat_part_GMP=conj(ret); return ret;};
#endif

	//! Constructor
	One_Loop_Helicity_Amplitude(const process& pro, color_structure cs, Rational_factory<Rational_base>* RRF=Rational_factory<Rational_base>::default_rational_factory() ) ;
	One_Loop_Helicity_Amplitude(const process& pro, color_structure cs, Rational_factory<Rational_base>* RRF, cut_part_factory<Cut_Part_base>* CPF) ;
	virtual ~One_Loop_Helicity_Amplitude();

	virtual void dry_run(const std::vector<int>& ind){_cut_part->dry_run(ind);_rational_part->dry_run(ind);}
	void set_mu(int index){_cut_part->set_mu(index);};
	void set_mu_HP(int index){_cut_part->set_mu_HP(index);};
	void set_mu_VHP(int index){_cut_part->set_mu_VHP(index);};
#ifdef BH_USE_GMP
	void set_mu_GMP(int index){_cut_part->set_mu_GMP(index);};
	void set_mu(R val,const RHP& val_HP,const RVHP& val_VHP,const RGMP& val_GMP){_cut_part->set_mu(val,val_HP,val_VHP,val_GMP);};
#else
	void set_mu(R val,RHP val_HP,RVHP val_VHP){_cut_part->set_mu(val,val_HP,val_VHP);};
#endif
	void set_mu(const multi_precision_constant& multi){_cut_part->set_mu(multi);};
	void set_scheme(scheme sc){_scheme=sc;};
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return _accuracy;};
	SeriesC<R> get_conjugate_amplitude(){return (_conjugate_cut_part+_conjugate_rat_part);};
	SeriesC<RHP> get_conjugate_amplitude_HP(){return (_conjugate_cut_part_HP+_conjugate_rat_part_HP);};
	SeriesC<RVHP> get_conjugate_amplitude_VHP(){return (_conjugate_cut_part_VHP+_conjugate_rat_part_VHP);};
#ifdef BH_USE_GMP
	SeriesC<RGMP> get_conjugate_amplitude_GMP(){return (_conjugate_cut_part_GMP+_conjugate_rat_part_GMP);};
#endif
private :
	void init(const process& pro, color_structure cs,Rational_factory<Rational_base>* RRF,cut_part_factory<Cut_Part_base>* CPF);
	template <class T> SeriesC<T> eval_fn(momentum_configuration<T>&,const std::vector<int>&);
};


typedef One_Loop_Helicity_Amplitude OneLoopHelAmpl;


#if 0
class OneLoopAmplitude : public Amplitude {
	std::vector<raw_bubble> r_bubbles;
	std::vector<raw_triangle> r_triangles;
	std::vector<raw_box> r_boxes;
	std::map<int,raw_part*> code_to_rpart;
	std::vector<One_Loop_Helicity_Amplitude<Cut_Part,Rec_Rational>*> hel_amplitudes;

public :
	//! constructor
	OneLoopAmplitude(const std::vector<particle*>& p,const  std::vector<particle_ID>& possible_props);
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
	One_Loop_Helicity_Amplitude<Cut_Part,Rec_Rational>* get_hel_ampl(size_t i);
	virtual ~OneLoopAmplitude();
};
#endif
}


#endif /*ONELOOPHELAMPL_H_*/
