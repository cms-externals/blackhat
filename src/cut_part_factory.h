/*!\file cut_part_factory.h
\brief Header for the cut_part factories momentum configurations
*/

#include "cut_part.h"

#ifndef CUT_PART_FACTORY_H_
#define CUT_PART_FACTORY_H_

namespace BH {

class IR_checked_Cut_Part;

//namespace Skeleton {
//	class FHZ_factory_wS;
//	class Darren_factory_wS;
//}

//! abstract class for cut_part factories
template <class cut_part_type> class cut_part_factory {
public:
	typedef cut_part_type Cut_Part_Type;
	virtual Cut_Part_Type* new_cut_part(const process&,color_structure)=0;
	static cut_part_factory<cut_part_type>* s_default_cut_part_factory(const process& pro); // the cut part factory for the standard case
	virtual ~cut_part_factory(){};
};

#ifdef SWIG
%template(CutPartFactory) cut_part_factory<Cut_Part_base>;
#endif

template <> cut_part_factory<Cut_Part_base>* cut_part_factory<Cut_Part_base>::s_default_cut_part_factory(const process& pro);

class Known_cut_part_factory : public cut_part_factory<Cut_Part_base> {
public:
	virtual Cut_Part_base* new_cut_part(const process&,color_structure);
	static Known_cut_part_factory* s_default_KCPF;
	virtual ~Known_cut_part_factory(){};
};

#ifdef USE_OLD_CPF

//! concrete class for cut_part factories
template <class CUTD_FACTORY> class custom_cut_part_factory : public cut_part_factory<Cut_Part_base> {
public:
	virtual Cut_Part_base* new_cut_part(const process&,color_structure);
};

//typedef custom_cut_part_factory<Darren_factory_new> QCD_cut_part_factory;
//typedef custom_cut_part_factory<Darren_factory> Darren_cut_part_factory;
//typedef custom_cut_part_factory<Darren_factory_new> Darren_new_cut_part_factory;
//typedef custom_cut_part_factory<FHZ_factory> FHZ_cut_part_factory;
//typedef custom_cut_part_factory<Skeleton::FHZ_factory_wS> FHZ_cut_part_factory_wS;
//typedef custom_cut_part_factory<Skeleton::Darren_factory_wS> Darren_cut_part_factory_wS;

//! concrete class for cut_part factories
template <class CUTD_FACTORY> class custom_IRC_cut_part_factory : public cut_part_factory<Cut_Part_base> {
public:
	virtual Cut_Part_base* new_cut_part(const process&,color_structure);
};

typedef custom_IRC_cut_part_factory<Darren_factory_new> IRC_cut_part_factory ;
typedef custom_IRC_cut_part_factory<Darren_factory_new> Darren_new_IRC_cut_part_factory ;
typedef custom_IRC_cut_part_factory<FHZ_factory> FHZ_IRC_cut_part_factory ;

#endif


}
#endif /* CUT_PART_FACTORY_H_ */
