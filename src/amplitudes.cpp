/*!\file amplitudes.cpp
\brief implementation of the amplitudes: OneLoopHelAmpl, OneLoopAmplitude
*/
#include "amplitudes.h"
#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "integrals.h"
#include "BH_A0.h"
#include "OneLoopHelAmpl.h"
#include "IR_checked.h"
#include "cut_part_factory.h"
#include "scheme.h"
#include "process_utils.h"
#include "process_utils_massive.h"
#include "mode_dependent_typedefs.h"  // for TREE_TYPE and TREE_FACTORY_TYPE
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "settings.h"
#include <cassert>
#include "BH_debug.h"

#if BH_USE_GMP
#include "gmp_r.h"
#endif
#define _VERBOSE 0
#define _VERBOSE_ZER0 0 //Used to switch on the reasons why we set an amplitude to zero

using namespace std;
namespace BH {


ostream& operator<<(ostream& s, HelAmpl& p){
	return s<<"A("<< p.d_process << ")";
}


void HelAmpl::print(){
	cout<<"A("<< d_process << ")"<<endl;
}

//C HelAmpl::eval(mom_conf mc){
//	_MESSAGE("Trying to evaluate:");
//	cout<<"A("<< pro << ")"<<endl;
//	_MESSAGE(" with ::eval(), but I don't know how to do it... Returned 0.");
//	return 0.;
//}

bool HelAmpl::is_zero() const {
	_MESSAGE("Trying to evaluate is_zero() for ");
	cout<<"A("<< d_process << ")"<<endl;
	_MESSAGE("but I don't know this amplitude... Returned false. ");
	return false;
}

bool HelAmpl::is_split_helicity(){
	size_t nbr_changes=0;
	for (int i=1;i<=d_process.n();i++){
		if (d_process.p(i).helicity()!=d_process.p(i%d_process.n()+1).helicity() ) {
			nbr_changes++;
//			_PRINT(i);
		}
	}
	return (nbr_changes==2);
}

size_t count_ph(const process& p, const particle_ID& t){
	size_t c=0;
	for (size_t i=1;i<=p.n();i++) {if (p.p(i).is_a(t)) c++;}
	return c;
}

size_t count_ph(const process& p, const particle& t){
	size_t c=0;
	for (size_t i=1;i<=p.n();i++) {if (p.p(i).is_a(t)) c++;}
	return c;
}

size_t count_masses(const process& p){
	size_t c=0;
	for (size_t i=1;i<=p.n();i++) {
		if(p.p(i).mass_label()>0) c++;
	}
	return c;
}

// A function for finding how many non-color ordered particles there are
int count_unordered(const process& p){
	int c=0;
	for (size_t i=1;i<=p.n();i++) {if (!p.p(i).ordered(particle::ord_type)) c++;}
	return c;
}

int count_unordered(const process& p, int level){
	int c=0;
	for (size_t i=1;i<=p.n();i++) {if (p.p(i).ordered(level)) c++;}
	return c;
}


bool is_zero_massive(const process& pro){
#if _VERBOSE_ZER0==1
	_MESSAGE2("Checking Massive amp=",pro);
#endif

	// Does something with the flavours of gluons
	process temp_pro=pro;

	// after fixing the flavors, there should be no negative flavor gluons, if there are some, it means fix_flavors failed and returned the original process
	for(size_t i=1;i<=pro.n();i++){
		if((temp_pro.p(i).is_a(gluon))&&(temp_pro.p(i).flavor()<0)){
#if _VERBOSE_ZER0==1
			_MESSAGE("Negative flavour gluons remain");
#endif
			return true;
		}
	}

	// Perform a count of all the different types of particle we could have
	size_t n_gluons_massive=count_ph(pro,gluon_massive_scalar)+count_ph(pro,gluon_massive);
	size_t n_scalars_massless=count_ph(pro,scalar);
	size_t n_scalars_massive=count_ph(pro,scalar_massive);
	size_t n_gluons=count_ph(pro,gluon);
	size_t n_gluinos=count_ph(pro,gluino);
	size_t n_massive_quarks=count_ph(pro,quark_massive);
	size_t n_massive_gluinos=count_ph(pro,gluino_massive);
	size_t n_quarks=count_ph(pro,quark);
	size_t n_leptons=count_ph(pro,lepton);
	size_t n_photons=count_ph(pro,photon);
	size_t n_higgs=count_ph(pro,higgs);
	size_t n_tot_quarks=n_quarks+n_massive_quarks+n_leptons;//We treat leptons as quarks for testing purposes
	size_t n_tot_gluinos=n_gluinos+n_massive_gluinos;
	// and whether they are even or not
	size_t n_quarks_mod=n_quarks%2;
	size_t n_gluinos_mod=n_gluinos%2;
	size_t n_quarks_massive_mod=n_massive_quarks%2;
	size_t n_gluinos_massive_mod=n_massive_gluinos%2;
	size_t n_tot_quarks_mod=n_tot_quarks%2;
	size_t n_tot_gluinos_mod=n_tot_gluinos%2;
	size_t n_gluons_massive_mod=n_gluons_massive%2;
//    size_t n_gluons_hel_massive_mod=n_gluons_hel_massive%2;
	size_t n_scalars_massless_mod=n_scalars_massless%2;
	size_t n_scalars_massive_mod=n_scalars_massive%2;

	//We must have an even number of quarks and/or gluinos of any flavour and mass
	if(n_tot_quarks_mod!=0||n_tot_gluinos_mod!=0){
#if _VERBOSE_ZER0==1
		_MESSAGE("Odd number of quarks");
#endif
		return true;
	}


	// Find all the different fermion flavours in the amplitude
	vector<int> flavours;
	vector<int> flavours_gluino;
	for(int i=1;i<=temp_pro.n();i++)
	{
		if((temp_pro.p(i).is_a(quark))||(temp_pro.p(i).is_a(quark_massive))){
			flavours.push_back(temp_pro.p(i).flavor());
		}
		if((temp_pro.p(i).is_a(gluino))||(temp_pro.p(i).is_a(gluino_massive))){
			//We mark all fermion flavours by adding 10000 to them, hopefully no actual 10,000 flavour lines
			// appear in other parts of the code
			flavours.push_back(temp_pro.p(i).flavor()+10000);
		}
	}
	sort(flavours.begin(),flavours.end());
	flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

#if _VERBOSE_ZER0==1
	if(flavours.size()>0){
		cout << "list of flavours=";
		_vector_cout(flavours);
	}
#endif

	// For each flavour make sure we have a particle and an anti-particle exist
	//  also if both entries are massless then they need to conserve helicity
	// First construct vectors to store the checked results of the different tests
	vector<int> test_part, test_hel, test_mass, test_nmbr;
	for(int itst=0;itst<flavours.size();itst++){
		test_part.push_back(0);
		test_hel.push_back(0);
		test_mass.push_back(0);
		test_nmbr.push_back(0);
	}

	// We also want to make sure that all same flavours are next to each other
	int first=-100;
	bool blastfound=false;
	int last_flavour=0;
	int non_adjacent=0;
	for(int i=1;i<=temp_pro.n();i++)
	{
		if(temp_pro.p(i).is_a(quark)||temp_pro.p(i).is_a(quark_massive)){
			int current_flavour=temp_pro.p(i).flavor();
			//Count the number of times we dont find a match
			if(last_flavour!=current_flavour){
				//We never find a match on the first time around so we want to store this
				if(first==-100){
					first=current_flavour;
				}
				non_adjacent++;
			}
			last_flavour=current_flavour;

			//Find the location in test for this flavour
			int ilook=0;
			while((flavours[ilook]!=temp_pro.p(i).flavor())&&(ilook++<flavours.size()));

			// Record that we have another quark
			test_nmbr[ilook]+=1;

			// Record the properties of the quark we have
			if(temp_pro.p(i).is_anti_particle()){
				test_part[ilook]+=1;
			}
			else{
				test_part[ilook]+=2;
			}
			if(temp_pro.p(i).mass_label()>0){
				test_mass[ilook]+=1;
			}
			else{
				test_mass[ilook]+=0;
			}
			test_hel[ilook]+=temp_pro.p(i).helicity();
		} else if(temp_pro.p(i).is_a(gluino)||temp_pro.p(i).is_a(gluino_massive)){
			int current_flavour=temp_pro.p(i).flavor()+10000;
			//Count the number of times we dont find a match
			if(last_flavour!=current_flavour){
				//We never find a match on the first time around so we want to store this
				if(first==-100){
					first=current_flavour;
				}
				non_adjacent++;
			}
			last_flavour=current_flavour;

			//Find the location in test for this flavour
			int ilook=0;
			while((flavours[ilook]!=temp_pro.p(i).flavor()+10000)&&(ilook++<flavours.size()));

			// Record that we have another fermion
			test_nmbr[ilook]+=1;

			// Record the properties of the fermion we have
			if(temp_pro.p(i).is_anti_particle()){
				test_part[ilook]+=1;
			}
			else{
				test_part[ilook]+=2;
			}
			if(temp_pro.p(i).mass_label()>0){
				test_mass[ilook]+=1;
			}
			else{
				test_mass[ilook]+=0;
			}
			test_hel[ilook]+=temp_pro.p(i).helicity();
		}
	}
	//We want to find the flavour of the next quark after the last quark (i.e. check how we wrap around) we stored this in first
	if(first!=last_flavour){
		non_adjacent++;
	}
//_MESSAGE6("process=",temp_pro," non_adjacent=",non_adjacent," flavours=",flavours.size());
	// If we have more than one missmatched pair we have a vanishing result, we also take into account the fact that we always get
	// an additional 1 added to the count because the first check always fails.
	if(non_adjacent>flavours.size()+2){
		return true;
	}
//#if _VERBOSE_ZER0==1
//	_MESSAGE("We have the following four test vectors");
//	if(flavours.size()>0){
//		_vector_cout(test_part);
//		_vector_cout(test_mass);
//		_vector_cout(test_hel);
//		_vector_cout(test_nmbr);
//  }
//#endif

	// Below we will use this information to check that all test are passed

	//We cannot have more than two quarks/gluinos of any flavour
	if(test_nmbr.size()>0){
		if(*max_element(test_nmbr.begin(),test_nmbr.end())>2){
#if _VERBOSE_ZER0==1
			_MESSAGE("More than two quarks of the same flavour");
#endif
			return true;
		}
	}

	// All quark flavours have 2 or fewer entries at this point. If there were flavour changing
	// particles then we would need to match them up. This would imply an even number of
	// mismatched flavours. We would still require that there was a quark-anti-quark pair for each pair of
	// mismatched flavours. Below we compute these quantities but we need only the nmbr_missmatched_flavours
	int test_part_tot=0;
	int cross_part_tot=0;
	int nmbr_missmatched_flavours=0;
	for(int itst=0;itst<test_part.size();itst++){
		// If there is only a single quark of that particular flavour then its "partner" assuming
		//  a scalar is present can be of a different flavour (this no longer the case, scalars preserve flavour)
		if(test_nmbr[itst]==1){
			nmbr_missmatched_flavours++;
			cross_part_tot+=test_part[itst];
		}
		else{
			test_part_tot+=test_part[itst];
			// If both quarks are massless then we must have helicity conservation
			if(test_mass[itst]==0&&test_hel[itst]!=0){
#if _VERBOSE_ZER0==1
				_MESSAGE("Mismatched helicities for massless quarks of the same flavour");
#endif
				return true;
			}
		}
	}

	// We must have no mismatched flavours
	if(nmbr_missmatched_flavours>0){
#if _VERBOSE_ZER0==1
		_MESSAGE("Mismatched flavour quarks");
#endif
		return true;
	}
	// We must have a quark anti-quark pair for each flavour
	for(int isch=0;isch<test_part.size();isch++){
		if(test_part[isch]!=3){
#if _VERBOSE_ZER0==1
				_MESSAGE("Mismatched number of particles and anti-particles for a flavour");
#endif
				return true;
		}
	}


	// Tests involving scalar amplitudes
	if(n_gluons_massive>0||n_scalars_massive>0){
		// Rules for scalars decompose into two groups depending upon whether we have an even or odd number of scalars
		if(n_gluons_massive_mod==1||n_scalars_massive_mod==1){
			// If we have an odd number of massive gluons and only gluons then we have zero
			if((pro.n()-n_gluons_massive)==n_gluons||(pro.n()-n_scalars_massive)==n_gluons){
#if _VERBOSE_ZER0==1
				_MESSAGE("Odd nmb of scalars and only gluons");
#endif
				return true;
			}
            
            //Otherwise We must have a massive-massless quark flavour if we have an odd number of scalars
            bool not_found_missmatch=true;
            for(int j=0;j<test_mass.size();j++){
                if(test_mass[j]==1){
                    not_found_missmatch=false;
                    break;
                }
            }
            if(not_found_missmatch){
#if _VERBOSE_ZER0==1
                _MESSAGE("Odd number of scalars but no massive-massless same flavour quark pair.");
#endif
                return true;
            }
		}
		else{
			// For an even number of scalars there must not be any massless-massive quark or gluinos pairs
			if(n_quarks_massive_mod!=0||n_gluinos_massive_mod!=0||n_gluinos_mod!=0||n_quarks_mod!=0){
#if _VERBOSE_ZER0==1
				_MESSAGE("Even number of scalars but massive-massless quark/gluino pairs present.");
#endif
				return true;
			}
		}
		// Tests on the three point vertex with a scalar
		if(temp_pro.n()==3){
            //If we have a single scalar then we get a non zero result only when we have opposite helicity gluons attached
			int hel=0;
            int R_cnt=0;
			/***************************/
            bool bhiggs(false);
            /***************************/
			for(int icnt=0;icnt<3;icnt++){
				hel+=pro.p(icnt+1).helicity();
                if(pro.p(icnt+1).helicity()==0&&pro.p(icnt+1).type()==gluon_massive){
                    R_cnt++;
                }
				
                /***************************/
                if(pro.p(icnt+1).type()==higgs){bhiggs=true;}
                /***************************/
			}
			if(hel!=0&&R_cnt==1&&!bhiggs&&n_gluons==2){// If we have opposite helicity legs then we have 0+1-1=0 otherwise we must have two helicities the same
#if _VERBOSE_ZER0==1
				_MESSAGE("Single scalar and not opposite helicity gluon 3-pt function");
#endif
				return true;
			}
            
			// Only s-s-g and Qb-q-s or Lb-G-s or ph-s-s (as well as gM-g-s) three point vertices are non-zero
			if((!(((n_gluons_massive==2||n_scalars_massive==2)&&n_gluons==1)||((n_gluons_massive==2||n_scalars_massive==2)&&n_higgs==1)
					||((n_gluons_massive==1||n_scalars_massive==1)&&((n_massive_quarks==1&&n_quarks==1)
								||(n_massive_gluinos==1&&n_gluinos==1)))))&&(n_gluons_massive==1||n_scalars_massive==1)){
//			if(((n_massive_quarks!=1)||(n_massive_gluinos!=1))&&(n_gluons_massive!=2)){
#if _VERBOSE_ZER0==1
				_MESSAGE("Three point scalar with wrong number of massive quarks");
#endif
				return true;
			}
		}
	}
    
//    if(n_gluons_hel_massive>0){
//        if(n_gluons_hel_massive_mod==1&&n_gluons_massive==0){
//#if _VERBOSE_ZER0==1
//            _MESSAGE("Odd nmb of massive gluons");
//#endif
//            return true;
//        }
//		
//
//    }

	// Tests involving leptons
	if(n_leptons>0){
		// NOTE the test below are temporary
		// We must have at least two leptons in the amplitude (i.e. we cannot split the llb pair)
		if(n_leptons%2==1){
#if _VERBOSE_ZER0==1
			_MESSAGE("Only a single lepton in the amplitude, need two");
#endif
			return true;
		}

		//We cannot have just leptons and scalars
		if(temp_pro.n()==4&&(n_gluons_massive==2||n_scalars_massive==2)){
#if _VERBOSE_ZER0==1
			_MESSAGE("Only two leptons and two scalars in the amplitude");
#endif
			return true;
		}
		//We cannot have just leptons and gluinos
		if((n_quarks+n_massive_quarks==0)&&(n_tot_gluinos>0)){
#if _VERBOSE_ZER0==1
			_MESSAGE("Only leptons and massive gluinos");
#endif
			return true;
		}
	}

	//Tests involving photons
	if(n_photons>0){
		//We must have at least 2 quarks gluinos or leptons in the amplitude or this is zero
		if(n_tot_quarks<2){
#if _VERBOSE_ZER0==1
			_MESSAGE("No quarks or leptons with a photon in the amplitude");
#endif
			return true;
		}
	}

	// If we have all but one helicity the same then we have zero as long as they do not contain more than 4 quarks some of which are massive
	//  or a higgs boson and a pair of scalars
	int matchhel=0;
	for(int i=1;i<pro.n();i++){
		if(pro.p(i).helicity()==pro.p(i+1).helicity()){
			matchhel++;
		}
	}
	if((matchhel==pro.n()-1)&&(!((n_tot_quarks+n_tot_gluinos>3)&&(n_massive_quarks+n_massive_gluinos>0)||(n_higgs>0)))){
#if _VERBOSE_ZER0==1
		_MESSAGE("Helicity mismatch");
#endif
		return true;
	}

//	//If we have a three point vertex and a higgs particle then some cases are zero
	if(n_higgs>0){
		// First see what type of higgs it is, we have seperate rules for the full higgs instead of the scalar parts.
		bool full_higgs=false;
		for(int i=1;i<=pro.n();i++){
			if(pro.p(i).flavor()==0&&pro.p(i).is_a(higgs)){
				full_higgs=true;
			}
		}
				
		if(pro.n()==3){//Is this a 3-pt vertex
			if(!(n_gluons==2||(n_gluons_massive==2||n_scalars_massive==2))){// Assume all three point vertices including quarks are zero, non-zero only for gluons and scalars
#if _VERBOSE_ZER0==1
				_MESSAGE("Higgs and non-gluon 3-pt function");
#endif
				return true;
			}
		}
		
		// If it is the full higgs H then
		if(full_higgs){
			if(pro.n()==3){//Is this a 3-pt vertex
				// Non zero only if both gluons have the same helicity, or there are scalars
				int hel=0;
				for(int icnt=0;icnt<3;icnt++){
					hel+=pro.p(icnt+1).helicity();
				}
				
				if(hel==0&&(n_gluons_massive==0||n_scalars_massive==0)){// If we have opposite helicity legs then we have 0+1-1=0 otherwise we must have two helicities the same, or a scalar and a gluon
#if _VERBOSE_ZER0==1
			_MESSAGE("Higgs and opposite helicity gluon 3-pt function");
#endif
					return true;
				}

			}
		}
		
		
//		if(pro.n()==3){//Is this a 3-pt vertex
//			if(!(n_gluons==2||n_gluons_massive==2)){// Assume all three point vertices including quarks are zero, non-zero only for gluons and scalars
//#if _VERBOSE_ZER0==1
//		_MESSAGE("Higgs and non-gluon 3-pt function");
//#endif
//				return true;
//			}
//            // Non zero only if both gluons have the same helicity, or there are scalars
//            int hel=0;
//            int phtype=1;//Assume we have a ph, change to a phd if we find one
//            for(int icnt=0;icnt<3;icnt++){
//                if(pro.p(icnt+1).is_a(phd)){phtype=-1;}
//                hel+=pro.p(icnt+1).helicity();
//            }
//            if(hel==0&&n_gluons_massive==0){// If we have opposite helicity legs then we have 0+1-1=0 otherwise we must have two helicities the same, or a scalar and a gluon
//#if _VERBOSE_ZER0==1
//		_MESSAGE("Higgs and opposite helicity gluon 3-pt function");
//#endif
//                return true;
//            }
//            
//            if((phtype==1&&hel>0)||(phtype==-1&&hel<0)){// If we have a positive helcity gluon and a scalar or two positive helicity gluons a phi amplitude is zero and the reverse for a phi dagger amplitude
//#if _VERBOSE_ZER0==1
//                _MESSAGE("Higgs and wrong helicity gluon 3-pt function");
//#endif
//                return true;
//            }
//            
//			//If there are only two other massless scalars then this is zero
//			if(n_scalars_massless==2||n_scalars_massive==2){
//				return true;
//			}
//		}
//
//
//		//If we have one or fewer negative helicity particles for a phi then it is zero
//		// or conversely one or fewer positive helicity legs for a phid. This translates
//		// to either a positive or negative helicity count n-1,n-3,...,-(n-3) and -(n-1)
//		// being the only possible outcomes unless we have a scalar present in which case
//		// we do not have zero
//		int phtype=1;//Assume we have a ph, change to a phd if we find one
//		int hel_sum=0;
//		for(int i=1;i<=pro.n();i++){
//			if(pro.p(i).is_a(phd)){phtype=-1;}
//			hel_sum+=pro.p(i).helicity();
//		}
//		int bound(pro.n()-3);//as pro.n() returns a size t (which is unsigned) we get a size_t for pro.n()-3, so here we implicitly create a signed integer
//		if(n_gluons_massive==0&&n_massive_quarks==0&&n_massive_gluinos==0){
//			if(phtype>0&&(hel_sum>=bound)){
//#if _VERBOSE_ZER0==1
//		_MESSAGE2("Higgs ph with not enough negative legs : ",pro);
//#endif
//			return true;
//		}
//		if(phtype<0&&(hel_sum<=-bound)){
//#if _VERBOSE_ZER0==1
//			_MESSAGE2("Higgs phd and not enough positive legs : ",pro);
//#endif
//			return true;
//			}
//		}
//		
//		// If we have two massive quarks a higgs and the remainder gluons then the gluon needs to be of the right helicity type for the higgs scalar
//		if(n_massive_quarks==2&&n_gluons+n_massive_quarks+1==pro.n()){
//			// We sum all the gluon helicities
//			int glue_hel=0;
//			for(int i=1;i<=pro.n();i++){
//				if(pro.p(i).is_a(gluon)){glue_hel+=pro.p(i).helicity();}
//			}
//			if(phtype>0){
//				// We have a ph therefore the all plus gluon amplitude is zero
//				if(glue_hel==n_gluons){
//#if _VERBOSE_ZER0==1
//					_MESSAGE2("Higgs ph and all positive helicity gluons : ",pro);
//#endif
//					return true;
//				}
//			}
//			else{
//				// We have a phs therefore we need a positive helicity gluon
//				if(glue_hel==-n_gluons){
//#if _VERBOSE_ZER0==1
//					_MESSAGE2("Higgs phd and all negative helicity gluons : ",pro);
//#endif
//					return true;
//				}
//			}
//		}

	}
	
	
	// Tests involving massless scalar amplitudes
	if(n_scalars_massless>0){
		// Rules for scalars decompose into two groups depending upon whether we have an even or odd number of scalars
		if(n_scalars_massless_mod==1){
			// If we have an odd number of scalars and only gluons then we have zero
			if((pro.n()-n_scalars_massless)==n_gluons){
#if _VERBOSE_ZER0==1
				_MESSAGE("Odd nmb of scalars and only gluons");
#endif
				return true;
			}
			//We must have a massive-massless quark flavour if we have an odd number of scalars
			bool not_found_missmatch=true;
			for(int j=0;j<test_mass.size();j++){
				if(test_mass[j]==1){
					not_found_missmatch=false;
					break;
				}
			}
			if(not_found_missmatch){
#if _VERBOSE_ZER0==1
				_MESSAGE("Odd number of scalars but no massive-massless same flavour quark pair.");
#endif
				return true;
			}
		}
		else{
			// For an even number of scalars there must not be any massless-massive quark or gluinos pairs
			if(n_quarks_massive_mod!=0||n_gluinos_massive_mod!=0||n_gluinos_mod!=0||n_quarks_mod!=0){
#if _VERBOSE_ZER0==1
				_MESSAGE("Even number of scalars massless but massive-massless quark/gluino pairs present.");
#endif
				return true;
			}
		}
		// Tests on the three point vertex with a scalar
		if(temp_pro.n()==3){
			// Only s-s-g and Qb-q-s or Lb-G-s or ph-s-s three point vertices are non-zero
			if(!((n_scalars_massless==2&&n_gluons==1)||(n_scalars_massless==2&&n_higgs==1)
				 ||(n_scalars_massless==1&&((n_massive_quarks==1&&n_quarks==1)
								   ||(n_massive_gluinos==1&&n_gluinos==1))))){
				//			if(((n_massive_quarks!=1)||(n_massive_gluinos!=1))&&(n_scalars!=2)){
#if _VERBOSE_ZER0==1
				_MESSAGE("Three point scalar massless with wrong number of massive quarks");
#endif
				return true;
			}
		}
	}	
	
	// Tests involving massless scalar amplitudes
	if(n_scalars_massive>0){
		// Rules for scalars decompose into two groups depending upon whether we have an even or odd number of scalars
		if(n_scalars_massive_mod==1){
			// If we have an odd number of scalars and only gluons then we have zero
			if((pro.n()-n_scalars_massive)==n_gluons){
#if _VERBOSE_ZER0==1
				_MESSAGE("Odd nmb of scalars and only gluons");
#endif
				return true;
			}
			//We must have a massive-massless quark flavour if we have an odd number of scalars
			bool not_found_missmatch=true;
			for(int j=0;j<test_mass.size();j++){
				if(test_mass[j]==1){
					not_found_missmatch=false;
					break;
				}
			}
			if(not_found_missmatch){
#if _VERBOSE_ZER0==1
				_MESSAGE("Odd number of scalars but no massive-massless same flavour quark pair.");
#endif
				return true;
			}
		}
		else{
			// For an even number of scalars there must not be any massless-massive quark or gluinos pairs
			if(n_quarks_massive_mod!=0||n_gluinos_massive_mod!=0||n_gluinos_mod!=0||n_quarks_mod!=0){
#if _VERBOSE_ZER0==1
				_MESSAGE("Even number of scalars massless but massive-massless quark/gluino pairs present.");
#endif
				return true;
			}
		}
		// Tests on the three point vertex with a scalar
		if(temp_pro.n()==3){
			// Only s-s-g and Qb-q-s or Lb-G-s or ph-s-s three point vertices are non-zero
			if(!((n_scalars_massive==2&&n_gluons==1)||(n_scalars_massive==2&&n_higgs==1)
				 ||(n_scalars_massive==1&&((n_massive_quarks==1&&n_quarks==1)
											||(n_massive_gluinos==1&&n_gluinos==1))))){
				//			if(((n_massive_quarks!=1)||(n_massive_gluinos!=1))&&(n_scalars!=2)){
#if _VERBOSE_ZER0==1
				_MESSAGE("Three point scalar massless with wrong number of massive quarks");
#endif
				return true;
			}
		}
	}
	
	// If there are no scalar type particles
	if(n_gluons_massive+n_scalars_massive+n_scalars_massive+n_scalars_massless==0){
		//We must have an even number of massless-massive quarks or gluinos if no scalars are present and there are massive quarks or gluinos
		if(((n_quarks+n_gluinos)%2!=0&&(n_massive_quarks+n_massive_gluinos)%2!=0)&&(n_massive_quarks+n_massive_gluinos>0)){
#if _VERBOSE_ZER0==1
			_MESSAGE("Uneven number of massless-massive quark pairs which are not allowed without a scalar_massive/scalar_massless/gluon_massive_scalar");
#endif
			return true;
		}
	}
	

#if _VERBOSE_ZER0==1
	_MESSAGE("Not Zero");
#endif

	return false;
}

bool is_zero_massless(const process& pro){
	
	BH_DEBUG_MESSAGE2("Checking Massless amp=",pro);
	
	long pc=pro.pcode();

	size_t n_gluons=(pc)%10;
	size_t n_quarks=(pc/10)%10;
	size_t n_leptons=(pc/100)%10;
	size_t n_photons=(pc/100000)%10;
	size_t n_gluinos=(pc/1000000)%10;
	size_t n_scalars_massless=(pc/1000000000)%10;


	switch (n_quarks){
		case  1:case  3:case  5: case  7: case 9:  return true;   // odd number of quarks
	}
	switch (n_leptons){
		case  1:case  3:case  5: case  7: case 9:  return true;   // odd number of leptons
	}
	switch (n_gluinos){
		case  1:case  3:case  5: case  7: case 9:  return true;   // odd number of gluinos
	}


	switch (pc/10){
		case 200020:  return true;
		case 210000:  return true;
		case 210020:  return true;
		case 210022:  return true;// set to zero amplitudes with multiple QED couplings
		case 400020:  return true;
		case  10022:  return true;
	}


	// Does something with the flavours of gluons
	process temp_pro=pro;

	// after fixing the flavors, there should be no negative flavor gluons, if there are some, it means fix_flavors failed and returned the original process
	for(size_t i=1;i<=pro.n();i++){
		if((temp_pro.p(i).flavor()<0)&&(temp_pro.p(i).is_a(gluon))){
#if _VERBOSE_ZER0==1
			_MESSAGE("Negative flavour gluons remain");
#endif
			return true;
		}
	}


	size_t n_tot_quarks=n_quarks;//We treat leptons as quarks for testing purposes
	// and weather they are even or not
	size_t n_quarks_mod=n_quarks%2;
	size_t n_tot_quarks_mod=n_tot_quarks%2;
	if(n_gluinos%2 ==1){
#if _VERBOSE_ZER0==1
		_MESSAGE("Odd number of quarks");
#endif
		return true;
		}
	//We must have an even number of quarks of any flavour and mass
	if(n_tot_quarks_mod>0){
#if _VERBOSE_ZER0==1
		_MESSAGE("Odd number of quarks");
#endif
		return true;
	}
//new
	//new adjacent flavor constraint
	//starting from a given flavored quark or gluino we walk around the process.
	//we may only encounter a second flavored quark, if intermediate flavors are saturated
	//we count gluino flavors negative to avoid conflicts with quarks
	if(n_quarks>2||n_gluinos>2||(n_quarks+n_gluinos+n_leptons)>2){
		std::vector<int> non_saturated_flav; non_saturated_flav.push_back(0);
		int flav_loc(0);
		//ad we walk around the process the vector on_saturated_flav should grow, 
		//as we encounter quarks and gluinos, but it should have size one when we 
		//reach the last particle. if and only if not some flavor lines must have crossed.
		for(int i=1;i<pro.n()+1;i++){
			//_MESSAGE(non_saturated_flav);
			if(pro.p(i).is_a(gluon)) continue;
			else if(pro.p(i).is_a(quark)){
				flav_loc=pro.p(i).flavor();
				if(flav_loc==non_saturated_flav.back()) non_saturated_flav.pop_back();
				else non_saturated_flav.push_back(flav_loc);
			}
			else if(pro.p(i).is_a(gluino)){
				flav_loc=-pro.p(i).flavor();// we keep track of gluinos using negative flavor entries
				if(flav_loc==non_saturated_flav.back()) non_saturated_flav.pop_back();
				else non_saturated_flav.push_back(flav_loc);
			}
			else if(pro.p(i).is_a(gluino)){
				flav_loc=-pro.p(i).flavor();// we keep track of gluinos using negative flavor entries
				if(flav_loc==non_saturated_flav.back()) non_saturated_flav.pop_back();
				else non_saturated_flav.push_back(flav_loc);
			}
			else if(pro.p(i).is_a(lepton)){
				flav_loc=1000+pro.p(i).flavor();// we keep track of leptons by shifting their flavor by 1000
				if(flav_loc==non_saturated_flav.back()) non_saturated_flav.pop_back();
				else non_saturated_flav.push_back(flav_loc);
			}
		}
		if(non_saturated_flav.size()>1){
#if _VERBOSE_ZER0==1
		_MESSAGE2(pro,": zero because of crossing flavor-lines");
#endif
			return true;
		}		
	}
//also checks that leptons do no sit between gluinos (G l l G forbidden)
//also checks that photons do no sit between gluinos (G y G forbidden)
//for this we check if there is an even number of gluinos 
//an leptons and quarks to either side of the leptons
	if(n_quarks==2&&n_gluinos>1&&(n_leptons)>1){
		size_t count_gluinos(0);	
		int nbr=pro.n();
		int pos_start(0);
		//start from quark
		while(!pro.p((pos_start++)%nbr+1).is_a(quark));
		//move on to first lepton
		while(!pro.p((pos_start++)%nbr+1).is_a(lepton));

		int pos_loc(pos_start%nbr);
		while(!pro.p((pos_loc)%nbr+1).is_a(quark)){
			if(pro.p(pos_loc%nbr+1).is_a(gluino)) count_gluinos++;
			pos_loc++;
		}
		if(count_gluinos%2==1){
			//_MESSAGE("lepton gluino conflict"); 
			return true;
		}
		else{
			int pos_loc(pos_start%nbr+nbr);
			while(!(pro.p((pos_loc)%nbr+1).is_a(quark))){
				if(pro.p(pos_loc%nbr+1).is_a(gluino)) count_gluinos++;
				pos_loc--;
			}
			if(count_gluinos%2==1){
				//_MESSAGE("lepton gluino conflict"); 
				return true;
			}
		}
	}
	
	
	
/* removed because of more general above function
	// if there is two pairs of quraks, the quarks of the same flavor have to be adjacent

	if (n_quarks == 4){
		vector<particle_ID> all_quarks;
		remove_copy_if(pro.begin(),pro.end(),back_inserter(all_quarks),not1(is_of_type(quark)));
		if ( all_quarks[0].flavor()==all_quarks[2].flavor() && all_quarks[0].flavor()!=all_quarks[1].flavor() ) {
		#if _VERBOSE_ZER0==1
		_MESSAGE("if there is two pairs of quraks, the quarks of the same flavor have to be adjacent");
		#endif

			return true;
		}
	}

	if (n_quarks == 2 && n_gluinos == 2){
		process intern=replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1);
		vector<particle_ID> all_quarks;

		remove_copy_if(intern.begin(),intern.end(),back_inserter(all_quarks),not1(is_of_type(quark)));
		if ( all_quarks[0].flavor()==all_quarks[2].flavor() && all_quarks[0].flavor()!=all_quarks[1].flavor() ) {
		#if _VERBOSE_ZER0==1
		_MESSAGE("if there is a pair of quarks and a pair of gluinos, the fermion line should not cross");
		#endif

			return true;
		}
	}
*/
	
	
	
	

	// Find all the different fermion flavours in the amplitude
	vector<int> flavours;
	for(int i=1;i<=temp_pro.n();i++)
	{
		if((temp_pro.p(i).is_a(quark))||(temp_pro.p(i).is_a(quark_massive))){
			flavours.push_back(temp_pro.p(i).flavor());
		}
	}
	sort(flavours.begin(),flavours.end());
	flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());
#if _VERBOSE_ZER0==1
	if(flavours.size()>0){
		cout << "list of flavours=";
		_vector_cout(flavours);
	}
#endif

	// For each flavour make sure we have a particle and an anti-particle
	//  also if both entries are massless then they need to conserve helicity
	// First construct a vector to store the checked results of the particle-anti-particle test
	vector<int> test_part, test_hel, test_nmbr, test_lepton;
	for(int itst=0;itst<flavours.size();itst++){
		test_part.push_back(0);
		test_hel.push_back(0);
		test_nmbr.push_back(0);
	}

	for(int i=1;i<=temp_pro.n();i++)
	{
		if((temp_pro.p(i).is_a(quark))||(temp_pro.p(i).is_a(quark_massive))){
			//Find the location in test for this flavour
			int ilook=0;
			while((flavours[ilook]!=temp_pro.p(i).flavor())&&(ilook++<flavours.size()));

			// Record that we have another quark
			test_nmbr[ilook]+=1;

			// Record the properties of the quark we have
			if(temp_pro.p(i).is_anti_particle()){
				test_part[ilook]+=1;
			}
			else{
				test_part[ilook]+=2;
			}
			test_hel[ilook]+=temp_pro.p(i).helicity();
		}
	}
#if _VERBOSE_ZER0==1
	_MESSAGE("We have the following four test vectors");
	if(flavours.size()>0){
		_vector_cout(test_part);
		_vector_cout(test_hel);
		_vector_cout(test_nmbr);
  }
#endif

	// Below we will use this information to check that all test are passed

	//We cannot have more than two quarks of any flavour
	if(test_nmbr.size()>0){
		if(*max_element(test_nmbr.begin(),test_nmbr.end())>2){
#if _VERBOSE_ZER0==1
			_MESSAGE("More than two quarks of the same flavour");
#endif
			return true;
		}
	}


	// Unless we have matching numbers of quark and anti-quarks (independant of flavour) then we have zero
	int test_part_tot=0;
	int cross_part_tot=0;
	int nmbr_missmatched_flavours=0;
	for(int itst=0;itst<test_part.size();itst++){
		// If there is only a single quark of that particluar flavour then we have a missmatched flavour
		if(test_nmbr[itst]==1){
			nmbr_missmatched_flavours++;
			cross_part_tot+=test_part[itst];
		}
		else{
			test_part_tot+=test_part[itst];
		}
	}
	//might be too spcific if we want Ws around
	//one should then
	for(int itst=0;itst<test_part.size();itst++){
		if(test_part[itst]!=3){
#if _VERBOSE_ZER0==1
		_MESSAGE2("Missmatched massless quarks and anti-quarks of flavor: ",flavours[itst]);
#endif
			return true;
		}
	}
	
	// In the completely massless case with no flavour changing leptons we must have no mismatched
	//  quarks
	if(nmbr_missmatched_flavours>0){
#if _VERBOSE_ZER0==1
		_MESSAGE("Missmatched massless quarks");
#endif
		return true;
	}

	// We must have quark helicity conservation
	int helsum=0;
	for(int itst=0;itst<test_hel.size();itst++){
		helsum+=test_hel[itst];
	}
	if(helsum!=0){
#if _VERBOSE_ZER0==1
				_MESSAGE("Missmatched helicities for massless quarks");
#endif
		return true;
	}

	// We either must have a quark anti-quark pair in each flavour
	if(test_part_tot!=test_part.size()*3){
#if _VERBOSE_ZER0==1
		_MESSAGE("Missmatched number of particles and anti-particles, no scalars");
#endif
		return true;
	}

	// Tests involving leptons
	if(n_leptons>0){
		// NOTE the test below are temporay
		// We must have at least two leptons in the amplitude (i.e. we cannot split the llb pair)
		if(n_leptons%2==1){
#if _VERBOSE_ZER0==1
			_MESSAGE("Only a single lepton in the amplitude, need two");
#endif
			return true;
		}
//		// We must have a flavour change in the other quarks, i.e. we require at least two unmatched flavours in the amplitude + an extra one for the leptons themselves
//		if(n_quarks>0&&nmbr_missmatched_flavours<2){
//#if _VERBOSE_ZER0==1
//			_MESSAGE("No flavour change in the quark with leptons present");
//#endif
//			return true;
//		}

		//We cannot have lepton and gluon amplitudes without a quark present
		if(n_quarks==0&&n_gluons>0){
#if _VERBOSE_ZER0==1
			_MESSAGE("Cannot have lepton and gluon amplitudes without a quark present");
#endif
			return true;
		}
	}

	//Tests involving photons
	if(n_photons>0){
		//If we only have photons and gluons in the amplitude then it must be zero
		if(n_tot_quarks==0&&n_leptons==0){
#if _VERBOSE_ZER0==1
			_MESSAGE("Photons and gluons only");
#endif
			return true;
		}
	}

	// If we have all but one helicity the same then we have zero
	int matchhelm=0;
	int matchhelp=0;
	int matchhel0=0;
	for(int i=1;i<=pro.n();i++){
		switch(pro.p(i).helicity()){
		case 1:
			matchhelp++;
			break;
		case -1:
			matchhelm++;
			break;
		case 0:
			matchhel0++;
			break;
		}
	}
//	_MESSAGE9(pro," : matchhelp=",matchhelp,", matchhelm=",matchhelm,", matchhel0=",matchhel0," pro.n()=",pro.n());

	if(((matchhelp>(pro.n()-1))||(matchhelm>(pro.n()-1)))&&(pro.n()==3)&&(n_gluons==3)){
#if _VERBOSE_ZER0==1
		_MESSAGE("Helicity mismatch");
#endif
       
		return true;
	}

	if(((matchhelp>(pro.n()-2))||(matchhelm>(pro.n()-2)))&&(pro.n()>3)){
#if _VERBOSE_ZER0==1
		_MESSAGE("Helicity mismatch");
#endif
		//allow for Tr(G^3) effective vertex amplitudes: need to allow more amplitudes
        if(settings::general::s_use_g3_coupling){return false;}
       
		return true;
	}

#if _VERBOSE_ZER0==1
		_MESSAGE("Testing gluinos");
#endif
	if ( ((pc/1000000)%10) != 0){  // if there are gluinos
	vector<particle_ID> all_gluinos;
	remove_copy_if(pro.begin(),pro.end(),back_inserter(all_gluinos),not1(is_of_type(gluino)));
	if (! all_gluinos.empty()){
		vector<int> all_gluino_flavors=find_all_flavors(pro,gluino);
		for (int fl=0;fl< all_gluino_flavors.size();fl++){
			vector<particle_ID> all_G;
			remove_copy_if(pro.begin(),pro.end(),back_inserter(all_G),not1(is_of_type_and_flavor(gluino,all_gluino_flavors[fl])));
			if (count_if(all_G.begin(),all_G.end(),is_of_type_and_pap(gluino,true))!=count_if(all_G.begin(),all_G.end(),is_of_type_and_pap(gluino,false))){
				return true;
			}
			if (count_if(all_G.begin(),all_G.end(),has_negative_helicity())!=count_if(all_G.begin(),all_G.end(),has_positive_helicity())){
				return true;
			}
		}

	}
	}

#if _VERBOSE_ZER0==1
		_MESSAGE("Testing quarks");
#endif
	if ( ((pc/10)%10) != 0) {  // if there are massless quarks
		vector<particle_ID> all_quarks;
		remove_copy_if(pro.begin(),pro.end(),back_inserter(all_quarks),not1(is_of_type(quark)));
		if (! all_quarks.empty()){
			vector<int> all_quark_flavors=find_all_flavors(pro,quark);
			for (int fl=0;fl< all_quark_flavors.size();fl++){
				vector<particle_ID> all_Q;
				remove_copy_if(pro.begin(),pro.end(),back_inserter(all_Q),not1(is_of_type_and_flavor_mod(quark,all_quark_flavors[fl],10)));
				//remove_copy_if(pro.begin(),pro.end(),back_inserter(all_Q),not1(is_of_type_and_flavor(quark,all_quark_flavors[fl])));

				if (count_if(all_Q.begin(),all_Q.end(),is_of_type_and_pap(quark,true))!=count_if(all_Q.begin(),all_Q.end(),is_of_type_and_pap(quark,false))){
					return true;
				}
				if (count_if(all_Q.begin(),all_Q.end(),has_negative_helicity())!=count_if(all_Q.begin(),all_Q.end(),has_positive_helicity())){
					return true;
				}
			}
		}
	}
#if _VERBOSE_ZER0==1
		_MESSAGE("Testing leptons");
#endif

	if ( ((pc/100)%10) != 0) {  // if there are leptons
		vector<particle_ID> all_leptons;
		remove_copy_if(pro.begin(),pro.end(),back_inserter(all_leptons),not1(is_of_type(lepton)));
		if (! all_leptons.empty()){
			vector<int> all_lepton_flavors=find_all_flavors(pro,lepton);
			for (int fl=0;fl< all_lepton_flavors.size();fl++){
				vector<particle_ID> all_Q;
				remove_copy_if(pro.begin(),pro.end(),back_inserter(all_Q),not1(is_of_type_and_flavor(lepton,all_lepton_flavors[fl])));

				if (count_if(all_Q.begin(),all_Q.end(),is_of_type_and_pap(lepton,true))!=count_if(all_Q.begin(),all_Q.end(),is_of_type_and_pap(lepton,false))){
					return true;
				}
				if (count_if(all_Q.begin(),all_Q.end(),has_negative_helicity())!=count_if(all_Q.begin(),all_Q.end(),has_positive_helicity())){
					return true;
				}
			}
		}
	}


	if(n_scalars_massless>0){
		// Check we have an even number of scalars
		if(n_scalars_massless%2==1){//Our massless scalars do not convert gluons to quarks etc.
			return true;
		}
		
		// If the remaining particles are gluons then we must have at least one different helicity, for amplitudes with more than 3 legs
		if(pro.n()-2==n_gluons&&pro.n()>3){
			int totalhel(0);
			for(int isum=0;isum<pro.n();isum++){
				totalhel+=pro[isum].helicity();
			}
			if(abs(totalhel)==n_gluons){
				return true;
			}
		}
	}


    /*********************
     * checks photons and quarks are aligned 
     *
     * 1) photons need quarks of identical flavor to couple to
     * 2) photons should not be trapped between wrong flavor quarks
     *
     * warning: interrete quark flavors MODULO 10
     */
    
    if(n_quarks>1 && n_photons >1){
    
        vector<int> flavor_ph;
        vector<int> flavor_q;

        for(size_t i=0;i<pro.n();i++){
            if(pro[i].is_a(quark) && pro[i].is_anti_particle() && count(flavor_q.begin(), flavor_q.end(), pro[i].flavor()%10)==0) flavor_q.push_back(pro[i].flavor()%10);
            if(pro[i].is_a(photon) && count(flavor_ph.begin(), flavor_ph.end(), pro[i].flavor()%10)==0 ) flavor_ph.push_back(pro[i].flavor());
        }
  
 
        size_t n_fl_q(flavor_q.size());
        size_t n_fl_ph(flavor_ph.size());
        /* if fewer quark flavors than photon flavors tree has to vanish
         */
        if(n_fl_q<n_fl_ph){
                 BH_DEBUG_MESSAGE3("Zero tree: ",pro," Photon quark flavor missmatch.");
		return true;
	}

        for(size_t i=0;i<n_fl_ph;i++){
        
            bool match_flavor(false);
            size_t pos(0), flavor_loc(flavor_ph[i]);
            while( pos!=n_fl_q && flavor_q[pos]!=flavor_loc ) pos++;
            
            /* did not match flavor if we have pos == n_fl_q
             */
            if(pos==n_fl_q){
                 BH_DEBUG_MESSAGE3("Zero tree: ",pro," Photon quark flavor missmatch.");
                return true;
            }
        }
    }
     


    /********************
     * end photons quark checks
     */
        



#if _VERBOSE_ZER0==1
	_MESSAGE("Not Zero");
#endif
	return false;
}


TreeHelAmpl::TreeHelAmpl(const process& pro): HelAmpl(pro)  {
	TREE_FACTORY_TYPE TF;
	d_tree_ptr=TF.new_tree(pro);
}

complex<R> TreeHelAmpl::eval(momentum_configuration<R>& mc,const vector<int>& ind){
	return d_tree_ptr->eval(mc,ind);
}
complex<RHP> TreeHelAmpl::eval(momentum_configuration<RHP>& mc,const vector<int>& ind){
	return d_tree_ptr->eval(mc,ind);
}
complex<RVHP> TreeHelAmpl::eval(momentum_configuration<RVHP>& mc,const vector<int>& ind){
	return d_tree_ptr->eval(mc,ind);
}

#if BH_USE_GMP
complex<RGMP> TreeHelAmpl::eval(momentum_configuration<RGMP>& mc,const vector<int>& ind){
	return d_tree_ptr->eval(mc,ind);
}
#endif

complex<R> TreeHelAmpl::eval(const eval_param<R>& ep){
	return d_tree_ptr->eval(ep);
}
complex<RHP> TreeHelAmpl::eval(const eval_param<RHP>& ep){
	return d_tree_ptr->eval(ep);
}
complex<RVHP> TreeHelAmpl::eval(const eval_param<RVHP>& ep){
	return d_tree_ptr->eval(ep);
}

size_t count_massless_scalars(const process& p){
	size_t c=0;
	for (size_t i=1;i<=p.n();i++) {
		if(p.p(i).type()==scalar) c++;
	}
	return c;
}
	
bool Tree_is_zero(const process& pro){
    // For now we split the rules up into massless cases and massive cases
    if(count_masses(pro)>0){
        return is_zero_massive(pro);
    }
    else{
        if(count_massless_scalars(pro)>0){
            return is_zero_massive(pro);
        }
        else{
            return is_zero_massless(pro);
        }
    }
}

bool TreeHelAmpl::is_zero() const
{
		return Tree_is_zero(d_process);
}

TreeHelAmpl::TreeHelAmpl(const TreeHelAmpl& T): HelAmpl(T.get_process()){
	TREE_FACTORY_TYPE TF;
	d_tree_ptr=TF.new_tree(T.get_process());
}

TreeHelAmpl TreeHelAmpl::operator=(const TreeHelAmpl& T){
	delete d_tree_ptr;
	HelAmpl::operator=(T);
	TREE_FACTORY_TYPE TF;
	d_tree_ptr=TF.new_tree(T.get_process());
	return *this;
}

TreeHelAmpl::~TreeHelAmpl(){ delete d_tree_ptr; }




}
