/*
 * worker_utils.cpp
 *
 *  Created on: 31 Jul 2009
 *      Author: daniel
 */

#include "worker_utils.h"
#include <istream>
#include <cassert>
#include "process.h"
#include <map>
#include <memory>
#include "eval_param.h"
#include "BH_debug.h"
using namespace std;

namespace BH {

namespace worker {

typedef map<long,particle*> pdg_map;

template <class T> struct do_delete : public std::unary_function<T*,void> {
        void operator()(T* ptr) { delete ptr;}
};





pdg_map* get_map(){
        pdg_map* pdg_to_particle_map=new pdg_map;
        pdg_to_particle_map->insert(pair<long,particle*>(1,&quark));
        pdg_to_particle_map->insert(pair<long,particle*>(11,&lepton));
        pdg_to_particle_map->insert(pair<long,particle*>(21,&gluon));
        pdg_to_particle_map->insert(pair<long,particle*>(8,&photon));
        pdg_to_particle_map->insert(pair<long,particle*>(1000,&gluino));
        //pdg_to_particle_map->insert(pair<long,particle*>(-1,&scalar));
        pdg_to_particle_map->insert(pair<long,particle*>(-1,&gluon_massive_scalar));
        pdg_to_particle_map->insert(pair<long,particle*>(-2,&quark_massive));
        pdg_to_particle_map->insert(pair<long,particle*>(-3,&gluino_massive));
        pdg_to_particle_map->insert(pair<long,particle*>(-4,&scalar_massive));
        pdg_to_particle_map->insert(pair<long,particle*>(-5,&scalar));
        pdg_to_particle_map->insert(pair<long,particle*>(-6,&gluon_massive));

        return pdg_to_particle_map;
}


void write(const process& PRO,ostream& os){

         os << "P " << PRO.n() << " ";

         for (int i=1;i<=PRO.n();i++){
                 os << PRO.p(i).type()->pdg_code()<< " ";
                 os << PRO.p(i).helicity()<< " ";
                 os << PRO.p(i).flavor()<< " ";
                 os << PRO.p(i).is_anti_particle()<< " ";
                 if ( PRO.p(i).mass_label() == 0){
                	 os << "ZM ";
                 } else
                 if ( PRO.p(i).mass_label() == rat_ext_mass.label() ){
                	 os << "REM ";
                 } else
                 if ( PRO.p(i).mass_label() == top_ext_mass.label() ){
                	 os << "RET ";
                 } else
                     if ( PRO.p(i).mass_label() == top_mass.label() ){
                    	 os << "TQ ";
                     } else
                         {
                        	 os << "Z ";
                         }
         }
         os << PRO.pcode() << " ";
}

particle_ID get_Particle_ID(long pdg_code,short hel,short flavor,bool ap){
        static auto_ptr<pdg_map> pdg_to_particle_map(get_map());
        particle* par=pdg_to_particle_map->operator[](pdg_code);
        if (! par){
        	throw BHerror("Unsupported particle");
        }
        return particle_ID(par,hel,flavor,ap);
}

bool read_process_from_stream(process& ret,istream& is){
        string title;
        is >> title;
        assert(title=="P");
        int nbr;
        is >> nbr;
        assert(nbr >= 3 );
        assert(nbr < 10 );
        vector<particle_ID> particles;
        long pdg_code;
        bool ap;
        short hel;
        short flavor;
        string mass_lab;
         for (int i=1;i<=nbr;i++){
                 is >> pdg_code;
                 is >> hel;
                 is >> flavor;
                 is >> ap;
                 is >> mass_lab;
                 if ( mass_lab == "ZM"){
                     particles.push_back(get_Particle_ID(pdg_code,hel,flavor,ap));
                 }
                 else if ( mass_lab == "REM"){
                     particles.push_back(get_Particle_ID(pdg_code,hel,flavor,ap));
                 }
                 else if ( mass_lab == "RET"){
                     particles.push_back(get_Particle_ID(pdg_code,hel,flavor,ap));
                 }
                 else if ( mass_lab == "TQ"){
                     particles.push_back(get_Particle_ID(pdg_code,hel,flavor,ap));
                 } else {
                     particles.push_back(get_Particle_ID(pdg_code,hel,flavor,ap));
                 }

         }
         long pcode;
         is >> pcode;
         ret=process(particles,pcode);
         BH_DEBUG_MESSAGE3("Read process ",ret," from file.");
         return true;
}


process read_process_from_stream(istream& is){
		process ret;
		read_process_from_stream(ret,is);
		return ret;
}




} /* worker */
} /* BH */
