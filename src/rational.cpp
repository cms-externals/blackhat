/*!\file rational.cpp
\implementation of the rational recursion relations using Darren's ratrecurs.cpp
*/
#include "BH_utilities.h"
#include "rational.h"
#include "BH_A0.h"
#include "cut_part_factory.h"
#ifndef BH_PUBLIC
#include "partial_order.h"
#include "ordering.h"
#include "rec_tree.h"
#endif
#include "BH_debug.h"



#define  _VERBOSE 0

#define _INCLUDE_MASSIVE_HEL_GLUON 0 //If we want to include the gpM and gmM particles then set this to 1

#if _VERBOSE
//Definitions needed for debugging purposes
#include "known_rational.h"
#include "rec_rational.h"
#include "inf_rat.h"
#endif


using namespace std;

namespace BH {

#ifndef BH_PUBLIC

void gen_all_partitions(const process& pro,vector<part>& leglist,int amp_i,int amp_j, const vector<int>& indlst)
{
	//Now set up the vector of partitions
	size_t s, t, pos, leni, lenj, n=pro.n();
	if(amp_i<amp_j){
		for(s=amp_i+1;(s>amp_i)&&(s<=amp_j);s==n?s=1:s++){
			for(amp_j==n?t=1:t=amp_j+1;(t>amp_j)||(t<=amp_i);t==n?t=1:t++){
				part v(2);
				lenj=(t-s+n)%n;
				leni=n-lenj;
				//Only store a partition if it is of length 2 or more
				if(lenj<2||leni<2)continue;
				// Insert all the legs on the left hand side
				for(size_t i1=0;i1<lenj;i1++){
					pos=((i1+s)==n)?n:(i1+s)%n;
					v.add(1,plabel(pro.p(indlst[pos-1]+1),indlst[pos-1]));
				}
				// Insert all the legs on the right hand side
				for(size_t i1=0;i1<leni;i1++){
					pos=((i1+t)==n)?n:(i1+t)%n;
					v.add(2,plabel(pro.p(indlst[pos-1]+1),indlst[pos-1]));
				}
				v.set_code(10*s+t);
				leglist.push_back(v);
			}
		}
	}
	else{
		for(s=amp_j+1;(s>amp_j)&&(s<=amp_i);s==n?s=1:s++){
			for(amp_i==n?t=1:t=amp_i+1;(t>amp_i)||(t<=amp_j);t==n?t=1:t++){
				part v(2);
				leni=(t-s+n)%n;
				lenj=n-leni;
				//Only store a partition if it is of length 2 or more
				if(lenj<2||leni<2)continue;
				// Insert all the legs on the left hand side
				for(size_t i1=0;i1<leni;i1++){
					pos=((i1+s)==n)?n:(i1+s)%n;
					v.add(2,plabel(pro.p(indlst[pos-1]+1),indlst[pos-1]));
				}
				// Insert all the legs on the right hand side
				for(size_t i1=0;i1<lenj;i1++){
					pos=((i1+t)==n)?n:(i1+t)%n;
					v.add(1,plabel(pro.p(indlst[pos-1]+1),indlst[pos-1]));
				}
				v.set_code(10*t+s);
				leglist.push_back(v);
			}
		}
	}
}

void construct_partitions(const process& pro,const vector<particle_ID>& possible_props,vector<part>& leglist,int amp_i,int amp_j){
	// Set up an ordered list to "label" each leg as we do not have an ind at this point
	vector<int> input;
	for(int i=0;i<pro.n();input.push_back(i++));

	// Generate the ind lists of all the possible orderings we can have
	ordering_constraint pro_ord=BH::ordering::get_ordering_constraint(pro);

	// Convert all entries to zero based rather than one based this will save doing a-1
	// all over the place in the evaluation code
	vector<vector<int> >::iterator list_itr;
	vector<int>::iterator ord_list_itr;
	for(list_itr=pro_ord.strong.begin();list_itr!=pro_ord.strong.end();list_itr++){
		for(ord_list_itr=list_itr->begin();ord_list_itr!=list_itr->end();(*ord_list_itr++)--);
	}
	for(list_itr=pro_ord.weak.begin();list_itr!=pro_ord.weak.end();list_itr++){
		for(ord_list_itr=list_itr->begin();ord_list_itr!=list_itr->end();(*ord_list_itr++)--);
	}
    BH_DEBUG_MESSAGE2("Ordering cnstr (ZB) : ",pro_ord);

	vector<vector<int> > list_ord=partially_ordered_vectors(input,pro_ord);

    BH_DEBUG(
	cout << "Ordered lists:" << endl;
	for(int iii=0;iii<list_ord.size();iii++){
		_vector_cout(list_ord[iii]);
        }; 
    )

	//Now compute all possible splits from these and store all the possible partitions in leglist
	// the order each time round depends upon the ordering of the external legs, the positions of i and j will
	//  also depend upon this ordering.
	vector<part> init_leglist;
	for(list_itr=list_ord.begin();list_itr!=list_ord.end();list_itr++){
		// Find the absolute position in the particular permutation of the shifted legs
		//  as it is on these legs that we will shift
		size_t rot_amp_i=0, rot_amp_j=0;
		for(;rot_amp_i<list_itr->size()&&(*list_itr)[rot_amp_i]!=(amp_i-1);rot_amp_i++);
		for(;rot_amp_j<list_itr->size()&&(*list_itr)[rot_amp_j]!=(amp_j-1);rot_amp_j++);
//		cout << "list_ord ";
//		_vector_cout(*list_itr);
		BH_DEBUG_MESSAGE4("New i=",rot_amp_i," New j=",rot_amp_j);
		gen_all_partitions(pro,init_leglist,rot_amp_i+1,rot_amp_j+1,*list_itr);
	}

	//leglist will have duplicate entries from splits that could come from different orderings, remove them
	//  we do this by sorting the entries and removing the non-unique ones. To do this
	//  we need a comparison function which sorts the legs relative to the weakly ordered legs

	// We want to compare in an unordered way the weakly ordered terms only, so find which softly ordered
	//  particles are in the trees and store them in a set
	set<int> order_with_respect_to;
	for(int cor=1;cor<3;cor++){
		for(int i=0;i<init_leglist[0].c(cor).size();i++){
//			if(init_leglist[0].c(cor)[i].is_a(lepton)||init_leglist[0].c(cor)[i].is_a(photon)||init_leglist[0].c(cor)[i].is_a(higgs)){
			if(init_leglist[0].c(cor)[i].get_ordered()>0){
				order_with_respect_to.insert(init_leglist[0].c(cor)[i].ind());
			}
		}
	}

//	_MESSAGE("Before sort");
//	for(vector<part>::iterator it=init_leglist.begin();it!=init_leglist.end();it++){
//		_MESSAGE(*it);
//	}

	sort(init_leglist.begin(),init_leglist.end(),unordered_compare_part(order_with_respect_to));
	vector<part>::iterator leglist_new_end=unique(init_leglist.begin(),init_leglist.end(),unordered_equal_part(order_with_respect_to));

//	_MESSAGE("After sort");
//	for(vector<part>::iterator it=init_leglist.begin();it!=leglist_new_end;it++){
//		_MESSAGE(*it);
//	}
//	_vector_cout(possible_props);

	// Now add the intermediate legs by summing over the possible propagator legs and their helicities
	int no_possible_props=possible_props.size();
	for(vector<part>::iterator leglist_itr=init_leglist.begin();leglist_itr!=leglist_new_end;leglist_itr++){
		for(size_t i2=0;i2<no_possible_props;i2++){
			vector<plabel> left=leglist_itr->c(1);
			vector<plabel> right=leglist_itr->c(2);

			//Insert the propagator legs |+><-| etc we do this by inverting the
			//  possible_props in the right hand vertex
			// Also we attach a particle label of -1 to these legs unless it is
			//  a weakly ordered particle (photon, lepton etc. ) in which case it is marked as -2 so that the
			//  sort and reduction of dulicate entries ahead will work.
			int insert_index=-1;
			if(possible_props[i2].is_a(lepton)||possible_props[i2].is_a(photon)||possible_props[i2].is_a(higgs)){
				insert_index=-2;
			}
			left.push_back(plabel(possible_props[i2],insert_index));
			right.insert(right.begin(),plabel(possible_props[i2].conjugate(),insert_index));
			part new_process_to_insert(2);
			new_process_to_insert.set_code(leglist_itr->get_code());
			new_process_to_insert.set_corner(1,left);
			new_process_to_insert.set_corner(2,right);
			leglist.push_back(new_process_to_insert);
		}
	}
    
//    _MESSAGE("After intermediate leg insertion:");
//    for(vector<part>::iterator it=leglist.begin();it!=leglist.end();it++){
//    		_MESSAGE(*it);
//    }


	//Finally remove any duplicate trees that contain unordered photons that
	// were introduced by a photon possible propagator
	order_with_respect_to.insert(-2);
	sort(leglist.begin(),leglist.end(),unordered_compare_part_prop(order_with_respect_to));
	leglist_new_end=unique(leglist.begin(),leglist.end(),unordered_equal_part_prop(order_with_respect_to));
	leglist.erase(leglist_new_end,leglist.end());
    
   BH_DEBUG_MESSAGE("Final result:");
   BH_DEBUG( for(vector<part>::iterator it=leglist.begin();it!=leglist.end();it++){
        _MESSAGE(*it);
   })

}

vector<particle_ID> possible_propagators(const process& pro){
	vector<particle_ID> result;

	// Code for deciding weather to include a scalar propagator or not in an all quark amplitude
	vector<int> flavours;
	bool found_gluon_massive=false;
   	bool found_massive=false;
	for(int i=1;i<=pro.n();i++)
	{
		if(pro.p(i).is_a(quark)){
			flavours.push_back(pro.p(i).flavor());
		}
		if(pro.p(i).is_a(quark_massive)){
			flavours.push_back(pro.p(i).flavor());
			found_massive=true;
		}
		if(pro.p(i).is_a(gluino)){
			//We mark all glunio flavours by adding 10000 to them, hopefully no actual 10,000 flavour lines
			// appear in other parts of the code
			flavours.push_back(pro.p(i).flavor()+10000);
		}
		if(pro.p(i).is_a(gluino_massive)){
			//We mark all glunio flavours by adding 10000 to them, hopefully no actual 10,000 flavour lines
			// appear in other parts of the code
			flavours.push_back(pro.p(i).flavor()+10000);
			found_massive=true;
		}
		// If there is a massive gluon present then we must include massive and massless fermions of the same flavour
		if(pro.p(i).is_a(scalar_massive)||pro.p(i).is_a(gluon_massive)){
			found_gluon_massive=true;
		}
	}
	sort(flavours.begin(),flavours.end());
	flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

	// Only if we have a pair of quark-massive-quarks or gluino-massive gluino (or a gluino/quark pair) vertices can we have an internal scalar
	if(flavours.size()>1){
		vector<int> test_mass;
		for(int itst=0;itst<flavours.size();itst++){
			test_mass.push_back(0);
		}

		for(int i=1;i<=pro.n();i++)
		{
			if(pro.p(i).is_a(quark)||pro.p(i).is_a(quark_massive)){
				//Find the location in test for this flavour
				int ilook=0;
				while((flavours[ilook]!=pro.p(i).flavor())&&(ilook++<flavours.size()));

				if(pro.p(i).mass_label()>0){
					test_mass[ilook]+=1;
				}
				else{
					test_mass[ilook]+=0;
				}
			} else if(pro.p(i).is_a(gluino)||pro.p(i).is_a(gluino_massive)){
				//Find the location in test for this flavour
				int ilook=0;
				while((flavours[ilook]!=pro.p(i).flavor()+10000)&&(ilook++<flavours.size()));

				if(pro.p(i).mass_label()>0){
					test_mass[ilook]+=1;
				}
				else{
					test_mass[ilook]+=0;
				}
			}
		}

		if(find(test_mass.begin(),test_mass.end(),1)!=test_mass.end()){
#if _INCLUDE_MASSIVE_HEL_GLUON==0
			result.push_back(particle_ID(scalar_massive,0,1,false));
#else
			result.push_back(particle_ID(gluon_massive,0,1,false));	
#endif
		}
	}


	//Now the rest of the code
	for(int i=1;i<=pro.n();i++){
		  if( pro.p(i).is_a(quark) ) {
			  result.push_back(particle_ID(quark, 1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(quark,-1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(quark, 1,pro.p(i).flavor(),false));
			  result.push_back(particle_ID(quark,-1,pro.p(i).flavor(),false));

			  //We add in a possible massive quark if we have a massless, for the masses we use all the masses we found
			  if(found_gluon_massive||found_massive){
				  result.push_back(particle_ID(quark_massive, 1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(quark_massive,-1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(quark_massive, 1,pro.p(i).flavor(),false));
				  result.push_back(particle_ID(quark_massive,-1,pro.p(i).flavor(),false));
			  }

			  //also need to have gluons in case of quark-only amplitudes
			  result.push_back(particle_ID(gluon, 1,1,false));
			  result.push_back(particle_ID(gluon,-1,1,false));

		  }
		  if( pro.p(i).is_a(lepton) ) {
			  result.push_back(particle_ID(lepton, 1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(lepton,-1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(lepton, 1,pro.p(i).flavor(),false));
			  result.push_back(particle_ID(lepton,-1,pro.p(i).flavor(),false));

				result.push_back(yp);
				result.push_back(ym);
			}
		  if( pro.p(i).is_a(quark_massive) ) {
			  result.push_back(particle_ID(quark_massive, 1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(quark_massive,-1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(quark_massive, 1,pro.p(i).flavor(),false));
			  result.push_back(particle_ID(quark_massive,-1,pro.p(i).flavor(),false));

			  // We need to insert the correct non dynamic mass depending upon wheather we have normal massive quarks or not
			  if(pro.p(i).mass_label()==rat_ext_mass.label()){// If the massive quark is for the rat_ext then we need a massless quark
			  if(found_gluon_massive){
				  result.push_back(particle_ID(quark, 1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(quark,-1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(quark, 1,pro.p(i).flavor(),false));
				  result.push_back(particle_ID(quark,-1,pro.p(i).flavor(),false));
			  }
		  }
			  else{ // We need to insert the appropriate rat_ext massive quark or not

			  }
		  }
		  if( pro.p(i).is_a(gluon) ) {
			  result.push_back(particle_ID(gluon, 1,pro.p(i).flavor(),false));
			  result.push_back(particle_ID(gluon,-1,pro.p(i).flavor(),false));
		  }
		  if( pro.p(i).is_a(gluon_massive) ) {
			  result.push_back(particle_ID(gluon_massive,0,pro.p(i).flavor(),false));
#if _INCLUDE_MASSIVE_HEL_GLUON==1
              // When we include these terms they need to be inserted in as well
              result.push_back(particle_ID(gluon_massive,1,pro.p(i).flavor(),false));
              result.push_back(particle_ID(gluon_massive,-1,pro.p(i).flavor(),false));
              result.push_back(particle_ID(gluon, 1,pro.p(i).flavor(),false));
              result.push_back(particle_ID(gluon,-1,pro.p(i).flavor(),false));
#endif
		  }
		  if( pro.p(i).is_a(gluino) ) {
			  result.push_back(particle_ID(gluino, 1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(gluino,-1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(gluino, 1,pro.p(i).flavor(),false));
			  result.push_back(particle_ID(gluino,-1,pro.p(i).flavor(),false));

			  //We add in a possible massive quark if we have a massless, for the masses we use all the masses we found
			  if(found_gluon_massive||found_massive){
				  result.push_back(particle_ID(gluino_massive, 1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(gluino_massive,-1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(gluino_massive, 1,pro.p(i).flavor(),false));
				  result.push_back(particle_ID(gluino_massive,-1,pro.p(i).flavor(),false));
			  }

//also need to have gluons in case of quark-only amplitudes
			  result.push_back(particle_ID(gluon, 1,1,false));
			  result.push_back(particle_ID(gluon,-1,1,false));

		  }
		  if( pro.p(i).is_a(gluino_massive) ) {
			  result.push_back(particle_ID(gluino_massive, 1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(gluino_massive,-1,pro.p(i).flavor(),true));
			  result.push_back(particle_ID(gluino_massive, 1,pro.p(i).flavor(),false));
			  result.push_back(particle_ID(gluino_massive,-1,pro.p(i).flavor(),false));

			  // We need to insert the correct non dynamic mass depending upon wheather we have normal massive quarks or not
			  if(pro.p(i).mass_label()==rat_ext_mass.label()){// If the massive quark is for the rat_ext then we need a massless quark
			  if(found_gluon_massive){
				  result.push_back(particle_ID(gluino, 1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(gluino,-1,pro.p(i).flavor(),true));
				  result.push_back(particle_ID(gluino, 1,pro.p(i).flavor(),false));
				  result.push_back(particle_ID(gluino,-1,pro.p(i).flavor(),false));
				  }
			  }else{ // We need to insert the appropriate rat_ext massive gluino or not

			  }
		  }
		if( pro.p(i).is_a(scalar) ) {
			result.push_back(particle_ID(scalar, 0,1,false));
		}
		if( pro.p(i).is_a(scalar_massive) ) {
			result.push_back(particle_ID(scalar_massive, 0,pro.p(i).flavor(),false));
		}
        if( pro.p(i).is_a(gluon_massive) && pro.p(i).helicity()==0) {
			result.push_back(particle_ID(gluon_massive,0,pro.p(i).flavor(),false));
#if _INCLUDE_MASSIVE_HEL_GLUON==1
            result.push_back(particle_ID(gluon_massive, 1,pro.p(i).flavor(),false));
            result.push_back(particle_ID(gluon_massive,-1,pro.p(i).flavor(),false));
            result.push_back(particle_ID(gluon, 1,pro.p(i).flavor(),false));
            result.push_back(particle_ID(gluon,-1,pro.p(i).flavor(),false));
#endif
        }
	
#if _INCLUDE_MASSIVE_HEL_GLUON==1
		if( pro.p(i).is_a(gluon_massive)  && pro.p(i).helicity()!=0) {
			result.push_back(particle_ID(gluon_massive,0,pro.p(i).flavor(),false));
            result.push_back(particle_ID(gluon_massive, 1,pro.p(i).flavor(),false));
            result.push_back(particle_ID(gluon_massive,-1,pro.p(i).flavor(),false));
            result.push_back(particle_ID(gluon, 1,pro.p(i).flavor(),false));
            result.push_back(particle_ID(gluon,-1,pro.p(i).flavor(),false));
        }
#endif

	}
	sort(result.begin(),result.end());
	result.erase(unique(result.begin(), result.end()),result.end());
	return result;
}


#endif
/*
 *
 *
 * Code for Rec_Pair
 *
 *
 *
 */


// Functions for doing the different types of i and j shifts
template <class T> size_t massless_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr)
{
	// Create the shifted momenta
	size_t indi=ind[i];
	size_t indj=ind[j];
	complex<T> z=-Psqr/(mc.spab(indi,P,indj));
	ind[i]=mc.insert(mc.L(indi),mc.Lt(indi)-z*mc.Lt(indj));
	ind[j]=mc.insert(mc.L(indj)+z*mc.L(indi),mc.Lt(indj));

	// Now return the on-shell propagator momentum
	int pij=mc.insert(mc.L(indi),mc.Lt(indj));
	return mc.insert(mc.mom(P)+z*mc.mom(pij));
}

template <class T> size_t massive_i_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr)
{
	size_t indi=ind[i];
	size_t indj=ind[j];

	// Construct the masslessly projected momenta for the shift, we assume only a single mass
	size_t indalti=mc.insert(mc.mom(indi)-T(0.5)*(mc.mom(*(ind.end()-1)).E()/mc.sp(indi,indj))*mc.mom(indj));

	// Compute the solution to the pole in z
	complex<T> z=-Psqr/(mc.spab(indalti,P,indj));

	// Set up the shift i and j legs, with i of mass m^2
	size_t shiftalt=mc.insert(mc.L(indalti),mc.Lt(indj));
	ind[i]=mc.insert(mc.mom(indi)-z*mc.mom(shiftalt),_mt_massive);
	ind[j]=mc.insert(mc.L(indj)+z*mc.L(indalti),mc.Lt(indj));

	// Return the intermediate propagator momentum
	return mc.insert(mc.mom(P)+z*mc.mom(shiftalt));
}

template <class T> size_t massive_j_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr)
{
	size_t indi=ind[i];
	size_t indj=ind[j];

	// Construct the masslessly projected momenta for the shift, we assume only a single mass
	size_t indaltj=mc.insert(mc.mom(indj)-T(0.5)*(mc.mom(*(ind.end()-1)).E()/mc.sp(indi,indj))*mc.mom(indi));

	// Compute the solution to the pole in z
	complex<T> z=-Psqr/(mc.spab(indi,P,indaltj));

	// Set up the shift i and j legs, and j has mass m^2
	size_t shiftalt=mc.insert(mc.L(indi),mc.Lt(indaltj));
	ind[i]=mc.insert(mc.L(indi),mc.Lt(indi)-z*mc.Lt(indaltj));
	ind[j]=mc.insert(mc.mom(indj)+z*mc.mom(shiftalt),_mt_massive);

	// Return the intermediate propagator momentum
	return mc.insert(mc[P]+z*mc[shiftalt]);
}

template <class T> size_t massive_prop_massive_i_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr)
{
	size_t indi=ind[i];
	size_t indj=ind[j];

	// Construct the masslessly projected momenta for the shift, we assume only a single mass
	size_t indalti=mc.insert(mc.mom(indi)-T(0.5)*(mc.mom(*(ind.end()-1)).E()/mc.sp(indi,indj))*mc.mom(indj));

	// Compute the solution to the pole in z
	complex<T> z=-Psqr/(mc.spab(indalti,P,indj));

	// Set up the shift i and j legs, with i of mass m^2
	size_t shiftalt=mc.insert(mc.L(indalti),mc.Lt(indj));
	ind[i]=mc.insert(mc.mom(indi)-z*mc.mom(shiftalt),_mt_massive);
	ind[j]=mc.insert(mc.L(indj)+z*mc.L(indalti),mc.Lt(indj));

	// Return the intermediate propagator momentum
	return mc.insert(mc.mom(P)+z*mc.mom(shiftalt),_mt_massive);
}

template <class T> size_t massive_prop_massive_j_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr)
{
	size_t indi=ind[i];
	size_t indj=ind[j];

	// Construct the masslessly projected momenta for the shift, we assume only a single mass
	size_t indaltj=mc.insert(mc.mom(indj)-T(0.5)*(mc.mom(*(ind.end()-1)).E()/mc.sp(indi,indj))*mc.mom(indi));

	// Compute the solution to the pole in z
	complex<T> z=-Psqr/(mc.spab(indi,P,indaltj));

	// Set up the shift i and j legs, and j has mass m^2
	size_t shiftalt=mc.insert(mc.L(indi),mc.Lt(indaltj));
	ind[i]=mc.insert(mc.L(indi),mc.Lt(indi)-z*mc.Lt(indaltj));
	ind[j]=mc.insert(mc.mom(indj)+z*mc.mom(shiftalt),_mt_massive);

	// Return the intermediate propagator momentum
	return mc.insert(mc[P]+z*mc[shiftalt],_mt_massive);
}

template <class T> size_t massive_ij_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr)
{
	//	 If any of the shifted legs is massive then we cannot use its momenta directly instead we peform a massless projection
	//	  in the direction of the other leg. If this other leg is also massive then we perform a double massless projection,
	//	  in the same way as K1flat and K2flat in the cut code, i.e. we define
	//	  hati=i-z\eta and hatj=j+z\eta with \eta=<i_flat|\gamma^{\mu}|j_flat>
	//	  with
	//	  i_flat=(i-(m_i/gamma)j)/(1-m_i m_j/gamma^2), j_flat=(j-(m_j/gamma)i)/(1-m_i m_j/gamma^2)
	//	  with gamma=i.j+\sqrt((i.j)^2-m_i m_j) and the m_i and m_j are the masses of i and j respectively.
	//	 In the massless case these just degenerate into the usual choices

	//Assumes a single mass
	size_t indi=ind[i];
	size_t indj=ind[j];
	complex<T> mi=mc.mom(*(ind.end()-1)).E();
	complex<T> mj=mc.mom(*(ind.end()-1)).E();
	complex<T> gamma,ij=mc.sp(indi,indj);
	if(ij.real()<T(0.)){
		gamma=ij-sqrt(ij*ij-mi*mj);
	}
	else{
		gamma=ij+sqrt(ij*ij-mi*mj);
	}
	complex<T> den=T(1.)/(T(1.)-mi*mj/pow(gamma,2));
	// Construct the masslessly projected momenta
	size_t indalti=mc.insert((mc.mom(indi)-(mi/gamma)*mc.mom(indj))*den);
	size_t indaltj=mc.insert((mc.mom(indj)-(mj/gamma)*mc.mom(indi))*den);

	// Compute the solution to the pole in z
	complex<T> z=-Psqr/(mc.spab(indalti,P,indaltj));

	// Set up the shift i and j legs
	size_t shiftalt=mc.insert(mc.L(indalti),mc.Lt(indaltj));
	ind[i]=mc.insert(mc.mom(indi)-z*mc.mom(shiftalt),_mt_massive);
	ind[j]=mc.insert(mc.mom(indj)+z*mc.mom(shiftalt),_mt_massive);

	// Return the intermediate propagator momentum
	return mc.insert(mc.mom(P)+z*mc.mom(shiftalt));
}

// Functions for doing the different types of i and j shifts
template <class T> void massless_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j)
{
	momentum<complex<T> > ij_mom=PfLLt(ep.p(j)->Lt(),ep.p(i)->L());
	complex<T> z=-Psqr/(T(2)*(P*ij_mom));
    
	// Insert the shifted legs
	shifti=Cmom<T>(ep.p(i)->L(),ep.p(i)->Lt()-z*ep.p(j)->Lt());
	shiftj=Cmom<T>(ep.p(j)->L()+z*ep.p(i)->L(),ep.p(j)->Lt());
    
	// Now compute the on-shell propagator momentum
	Phat=P+z*ij_mom;
    
    // Set up the reference momenta
    ref_i=ep.ref();
    ref_j=ep.ref();
}

template <class T> void massive_i_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j)
{
    Cmom<T> alti=(*ep.p(i))-(eval_param<T>::mass2(imass)/(T(2)*(ep.p(i)->P()*ep.p(j)->P())))*(*ep.p(j));
    lambda<T> lpart(alti.L());
    momentum<complex<T> > ij_mom=PfLLt(ep.p(j)->Lt(),alti.L());
//	lambda<T> lpart(ep.p(i)->Sm()*ep.p(j)->Lt());
//	momentum<complex<T> > ij_mom=PfLLt(ep.p(j)->Lt(),lpart);
	complex<T> z=-Psqr/(T(2)*(P*ij_mom));
    
	// Insert the shifted legs
	shifti=Cmom<T>(ep.p(i)->P()-z*ij_mom,_mt_massive);
    //	shiftj=Cmom<T>(ep.p(j)->L()+z*alti.L(),ep.p(j)->Lt());
	shiftj=Cmom<T>(ep.p(j)->L()+z*lpart,ep.p(j)->Lt());
	    
	// Now compute the on-shell propagator momentum
	Phat=P+z*ij_mom;
    
    // Set up the reference momenta
    ref_i=ep.p(j);
    ref_j=ep.ref();
}

template <class T> void massive_j_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j)
{
    Cmom<T> altj=(*ep.p(j))-(eval_param<T>::mass2(jmass)/(T(2)*(ep.p(i)->P()*ep.p(j)->P())))*(*ep.p(i));
    lambdat<T> ltpart(altj.Lt());
    momentum<complex<T> > ij_mom=PfLLt(altj.Lt(),ep.p(i)->L());
//	lambdat<T> ltpart(ep.p(j)->Sm()*ep.p(i)->L());
//	momentum<complex<T> > ij_mom=PfLLt(ltpart,ep.p(i)->L());
	complex<T> z=-Psqr/(T(2)*(P*ij_mom));
    
	// Insert the shifted legs
    //	shifti=Cmom<T>(ep.p(i)->L(),ep.p(i)->Lt()-z*altj.Lt());
	shifti=Cmom<T>(ep.p(i)->L(),ep.p(i)->Lt()-z*ltpart);
	shiftj=Cmom<T>(ep.p(j)->P()+z*ij_mom,_mt_massive);
   
	// Now compute the on-shell propagator momentum
	Phat=P+z*ij_mom;

    // Set up the reference momenta
    ref_i=ep.ref();
    ref_j=ep.p(i);
}

template <class T> void massive_prop_massive_i_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j)
{
	Cmom<T> alti=(*ep.p(i))-(eval_param<T>::mass2(imass)/(T(2)*(ep.p(i)->P()*ep.p(j)->P())))*(*ep.p(j));
    lambda<T> lpart(alti.L());
	momentum<complex<T> > ij_mom=PfLLt(ep.p(j)->Lt(),alti.L());
//	lambda<T> lpart(ep.p(i)->Sm()*ep.p(j)->Lt());
//    momentum<complex<T> > ij_mom=PfLLt(ep.p(j)->Lt(),lpart);
	complex<T> z=-Psqr/(T(2)*(P*ij_mom));
    
	// Insert the shifted legs
	shifti=Cmom<T>(ep.p(i)->P()-z*ij_mom,_mt_massive);
	// shiftj=Cmom<T>(ep.p(j)->L()+z*alti.L(),ep.p(j)->Lt());
    shiftj=Cmom<T>(ep.p(j)->L()+z*lpart,ep.p(j)->Lt());
    
	// Now compute the on-shell propagator momentum
	Phat=Cmom<T>(P+z*ij_mom,_mt_massive);
    
    // Set up the reference momenta
    ref_i=ep.p(j);
    ref_j=ep.ref();
}

template <class T> void massive_prop_massive_j_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j)
{
	Cmom<T> altj=(*ep.p(j))-(eval_param<T>::mass2(jmass)/(T(2)*(ep.p(i)->P()*ep.p(j)->P())))*(*ep.p(i));
    lambdat<T> ltpart(altj.Lt());
	momentum<complex<T> > ij_mom=PfLLt(altj.Lt(),ep.p(i)->L());
//	lambdat<T> ltpart(ep.p(j)->Sm()*ep.p(i)->L());
//    momentum<complex<T> > ij_mom=PfLLt(ltpart,ep.p(i)->L());
	complex<T> z=-Psqr/(T(2)*(P*ij_mom));
    
	// Insert the shifted legs
	// shifti=Cmom<T>(ep.p(i)->L(),ep.p(i)->Lt()-z*altj.Lt());
    shifti=Cmom<T>(ep.p(i)->L(),ep.p(i)->Lt()-z*ltpart);
	shiftj=Cmom<T>(ep.p(j)->P()+z*ij_mom,_mt_massive);
    
	// Now compute the on-shell propagator momentum
	Phat=Cmom<T>(P+z*ij_mom,_mt_massive);
    
    // Set up the reference momenta
    ref_i=ep.ref();
    ref_j=ep.p(i);
}

template <class T> void massive_ij_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j)
{
//	//	 If any of the shifted legs is massive then we cannot use its momenta directly instead we perform a massless projection
//	//	  in the direction of the other leg. If this other leg is also massive then we perform a double massless projection,
//	//	  in the same way as K1flat and K2flat in the cut code, i.e. we define
//	//	  hati=i-z\eta and hatj=j+z\eta with \eta=<i_flat|\gamma^{\mu}|j_flat>
//	//	  with
//	//	  i_flat=(i-(m_i/gamma)j)/(1-m_i m_j/gamma^2), j_flat=(j-(m_j/gamma)i)/(1-m_i m_j/gamma^2)
//	//	  with gamma=i.j+\sqrt((i.j)^2-m_i m_j) and the m_i and m_j are the masses of i and j respectively.
//	//	 In the massless case these just degenerate into the usual choices
//
//	//Assumes a single mass
//	size_t indi=ind[i];
//	size_t indj=ind[j];
//	complex<T> mi=mc.mom(*(ind.end()-1)).E();
//	complex<T> mj=mc.mom(*(ind.end()-1)).E();
//	complex<T> gamma,ij=mc.sp(indi,indj);
//	if(ij.real()<T(0.)){
//		gamma=ij-sqrt(ij*ij-mi*mj);
//	}
//	else{
//		gamma=ij+sqrt(ij*ij-mi*mj);
//	}
//	complex<T> den=T(1.)/(T(1.)-mi*mj/pow(gamma,2));
//	// Construct the masslessly projected momenta
//	size_t indalti=mc.insert((mc.mom(indi)-(mi/gamma)*mc.mom(indj))*den);
//	size_t indaltj=mc.insert((mc.mom(indj)-(mj/gamma)*mc.mom(indi))*den);
//
//	// Compute the solution to the pole in z
//	complex<T> z=-Psqr/(mc.spab(indalti,P,indaltj));
//
//	// Set up the shift i and j legs
//	size_t shiftalt=mc.insert(mc.L(indalti),mc.Lt(indaltj));
//	ind[i]=mc.insert(mc.mom(indi)-z*mc.mom(shiftalt),_mt_massive);
//	ind[j]=mc.insert(mc.mom(indj)+z*mc.mom(shiftalt),_mt_massive);
//
//	// Return the intermediate propagator momentum
//	return mc.insert(mc.mom(P)+z*mc.mom(shiftalt));
}

Rec_Pair::Rec_Pair(Rec_BB* left, Rec_BB* right, part& partition, int ii, int jj, Const_Fact_Fn* fact) : _part(partition), i(ii), j(jj),
		_epl(partition.c(1).size()), _epr(partition.c(2).size()),
		_epl_HP(partition.c(1).size()), _epr_HP(partition.c(2).size()),
		_epl_VHP(partition.c(1).size()), _epr_VHP(partition.c(2).size())
#if BH_USE_GMP
 ,_epl_GMP(partition.c(1).size()), _epr_GMP(partition.c(2).size())
#endif
{
	daughters.push_back(left);
	daughters.push_back(right);
	daughters.push_back(fact);

	//Set up the vectors for the indices that will be passed
	for(int i=0;i<partition.c(1).size();i++){
		indshiftl.push_back(0);
	}
	for(int i=0;i<partition.c(2).size();i++){
		indshiftr.push_back(0);
	}

	//Store the maximum size of the indexes as we will need these often
	maxl=indshiftl.size();
	maxr=indshiftr.size();
	max=maxl+maxr-2;

	//Find the location of j in indshiftl
	shifted_ind_j=0;
	while(shifted_ind_j<maxl-1&&get_part().c(1)[shifted_ind_j].ind()!=jj){shifted_ind_j++;}

	//Find the location of i in indshiftr
	shifted_ind_i=1;
	while(shifted_ind_i<maxr&&get_part().c(2)[shifted_ind_i].ind()!=ii){shifted_ind_i++;}

	// Store the locations of the shifted legs in the eval_params now
	_epr.set(0,&_PHat);
	_epr_HP.set(0,&_PHat_HP);
	_epr_VHP.set(0,&_PHat_VHP);
	
	_epl.set(maxl-1,&_mPHat);
	_epl_HP.set(maxl-1,&_mPHat_HP);
	_epl_VHP.set(maxl-1,&_mPHat_VHP);
    
#if BH_USE_GMP
    _epr_GMP.set(0,&_PHat_GMP);
    _epl_GMP.set(maxl-1,&_mPHat_GMP);
#endif
}

Rec_Pair_massive::Rec_Pair_massive(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, int imass, int jmass, int massive_shift, Const_Fact_Fn* fact) : Rec_Pair(left,right,partition,ii,jj,fact), _imass(imass), _jmass(jmass), d_massive_shift(massive_shift)
{
	switch(massive_shift){
	case 0: // Both are massless
		shift_ij=&massless_shift_ij<R>;
		shift_ij_HP=&massless_shift_ij<RHP>;
		shift_ij_VHP=&massless_shift_ij<RVHP>;
		shift_ij_ep=&massless_shift_ij_ep<R>;
		shift_ij_ep_HP=&massless_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massless_shift_ij_ep<RVHP>;
		break;
	case 1: // i is massive and j is massless
		shift_ij=&massive_i_shift_ij<R>;
		shift_ij_HP=&massive_i_shift_ij<RHP>;
		shift_ij_VHP=&massive_i_shift_ij<RVHP>;
		shift_ij_ep=&massive_i_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_i_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_i_shift_ij_ep<RVHP>;
		break;
	case 2: // j is massive and i is massless
		shift_ij=&massive_j_shift_ij<R>;
		shift_ij_HP=&massive_j_shift_ij<RHP>;
		shift_ij_VHP=&massive_j_shift_ij<RVHP>;
		shift_ij_ep=&massive_j_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_j_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_j_shift_ij_ep<RVHP>;
		break;
	case 3: // They are both massive
		shift_ij=&massive_ij_shift_ij<R>;
		shift_ij_HP=&massive_ij_shift_ij<RHP>;
		shift_ij_VHP=&massive_ij_shift_ij<RVHP>;
		shift_ij_ep=&massive_ij_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_ij_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_ij_shift_ij_ep<RVHP>;
		break;
	}

#ifdef BH_USE_GMP
	switch (massive_shift){
	case 0: // Both are massless
		shift_ij_GMP=&massless_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massless_shift_ij_ep<RGMP>;
		break;
	case 1: // i is massive and j is massless
		shift_ij_GMP=&massive_i_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massive_i_shift_ij_ep<RGMP>;
		break;
	case 2: // j is massive and i is massless
		shift_ij_GMP=&massive_j_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massive_j_shift_ij_ep<RGMP>;
		break;
	case 3: // They are both massive
		shift_ij_GMP=&massive_ij_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massive_ij_shift_ij_ep<RGMP>;
		break;
	}
#endif


	//We must also add the extra entries at the end of the vectors for the ref vector and mass
	indshiftl.push_back(0);
	indshiftl.push_back(0);
	indshiftr.push_back(0);
	indshiftr.push_back(0);
}

Rec_Pair_massive_prop::Rec_Pair_massive_prop(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, size_t in_mass_leg, int imass, int jmass, int massive_shift, Const_Fact_Fn* fact) :  Rec_Pair_massive(left,right,partition,ii,jj,imass,jmass,massive_shift,fact), _mass_leg(in_mass_leg)
{
	//We use modified shift functions here so change them to the correct functions
	switch(massive_shift){
	case 1: // i is massive and j is massless
		shift_ij=&massive_prop_massive_i_shift_ij<R>;
		shift_ij_HP=&massive_prop_massive_i_shift_ij<RHP>;
		shift_ij_VHP=&massive_prop_massive_i_shift_ij<RVHP>;
		shift_ij_ep=&massive_prop_massive_i_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_prop_massive_i_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_prop_massive_i_shift_ij_ep<RVHP>;
		break;
	case 2: // j is massive and i is massless
		shift_ij=&massive_prop_massive_j_shift_ij<R>;
		shift_ij_HP=&massive_prop_massive_j_shift_ij<RHP>;
		shift_ij_VHP=&massive_prop_massive_j_shift_ij<RVHP>;
		shift_ij_ep=&massive_prop_massive_j_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_prop_massive_j_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_prop_massive_j_shift_ij_ep<RVHP>;
		break;
	}

}


Rec_Pair_massive_prop_massless_shift::Rec_Pair_massive_prop_massless_shift(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, size_t in_mass_leg, Const_Fact_Fn* fact)  :  Rec_Pair(left,right,partition,ii,jj,fact), _mass_leg(in_mass_leg)
{
	//We must also add the extra entries at the end of the vectors for the ref vector and mass
	indshiftl.push_back(0);
	indshiftl.push_back(0);
	indshiftr.push_back(0);
	indshiftr.push_back(0);
}

Rec_Pair_massive_unshifted::Rec_Pair_massive_unshifted(Rec_BB* left,Rec_BB* right,part& partition,int ii,int jj, Const_Fact_Fn* fact)  :  Rec_Pair(left,right,partition,ii,jj,fact)
{
	//We must also add the extra entries at the end of the vectors for the ref vector and mass
	indshiftl.push_back(0);
	indshiftl.push_back(0);
	indshiftr.push_back(0);
	indshiftr.push_back(0);
}


template <class T> complex<T>  Rec_Pair::Rec_Pair_eval(momentum_configuration<T>& mc,const vector<int>& ind){
#if _VERBOSE
	_MESSAGE3("-------------------------Begin Rec_Pair_eval ",get_part().get_code()," -----------------------------------");
#endif
	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom=mc.mom(ind[c1[0].ind()]);
	indshiftl[0]=ind[c1[0].ind()];
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		Psum_mom+=mc.mom(ind[c1[psiter].ind()]);
		indshiftl[psiter]=ind[c1[psiter].ind()];
	}

	// Create the shifted momenta
	size_t indi=ind[i];
	size_t indj=ind[j];
	momentum<complex<T> > ij_mom=PfLLt(mc.Lt(indj),mc.L(indi));
	complex<T> Psqr=Psum_mom.square();
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));

	// Now compute the on-shell propagator momentum
	Cmom<T> Phat_mom(Psum_mom+z*ij_mom);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	indshiftr[0]=mc.insert(Phat_mom);
	indshiftl[maxl-1]=mc.insert(Phat_mom.L(),-Phat_mom.Lt());
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		indshiftr[kr]=ind[c2[kr].ind()];
	}

	// Insert the shifted legs
	indshiftr[shifted_ind_i]=mc.insert(mc.L(indi),mc.Lt(indi)-z*mc.Lt(indj));
	indshiftl[shifted_ind_j]=mc.insert(mc.L(indj)+z*mc.L(indi),mc.Lt(indj));

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
	_MESSAGE7(tempres," : From ",left()->eval(mc,indshiftl)," * (",T(1.)/Psqr,") * ",right()->eval(mc,indshiftr));
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");
	return tempres;
#else
	return Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
#endif

}

template <class T> complex<T>  Rec_Pair_massive::Rec_Pair_eval(momentum_configuration<T>& mc,const vector<int>& ind){
#if _VERBOSE
	_MESSAGE3("---------------------------Begin Rec_Pair_massive ",get_part().get_code()," ---------------------------------");
#endif
	// Compute the shifted legs, shiftBA gives an [i<j shift
	momentum<complex<T> > Psum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	for(size_t psiter=0;psiter<get_part().c(1).size()-1;psiter++){
		Psum_mom+=mc.mom(ind[get_part().c(1)[psiter].ind()]);
	}
	size_t P=mc.insert(Psum_mom,_mt_massive);
	complex<T> Psqr=mc.m2(P);

	vector<int> shifted_legs=ind;
	size_t Phat=get_shifted_ij(mc,shifted_legs,P,Psqr);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	size_t kl=0;
	for(;kl<get_part().c(1).size()-1;kl++){
		indshiftl[kl]=shifted_legs[get_part().c(1)[kl].ind()];
	}
	indshiftl[kl]=mc.insert(-mc[Phat]);
	indshiftr[0]=Phat;
	size_t kr=1;
	for(;kr<get_part().c(2).size();kr++){
		indshiftr[kr]=shifted_legs[get_part().c(2)[kr].ind()];
	}

	// In the case of massive amplitudes the mass is sent passed as an extra momentum
	//  as are any reference vectors needed these need to be tagged on to any amplitude
	//  we have produced.
	indshiftl[maxl]=ind[max];
	indshiftr[maxr]=ind[max];
	indshiftl[maxl+1]=ind[max+1];
	indshiftr[maxr+1]=ind[max+1];

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
	_MESSAGE7(tempres," : From ",left()->eval(mc,indshiftl)," * (",T(1.)/Psqr,") * ",right()->eval(mc,indshiftr));
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval massive ",get_part().get_code()," ----------------------------------");
	return tempres;
#else
	return complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr;
#endif

}

template <class T> complex<T>  Rec_Pair_massive_prop::Rec_Pair_eval(momentum_configuration<T>& mc,const vector<int>& ind){
#if _VERBOSE
	_MESSAGE3("--------------------------Begin Rec_Pair_massive_prop ",get_part().get_code()," ----------------------------------");
#endif
	// Compute the shifted legs, shiftBA gives an [i<j shift
	momentum<complex<T> > Psum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	for(size_t psiter=0;psiter<get_part().c(1).size()-1;psiter++){
		Psum_mom+=mc.mom(ind[get_part().c(1)[psiter].ind()]);
	}
	size_t P=mc.insert(Psum_mom,_mt_massive);
	// The mass is stored in the E component of the last element passed down in the ind
	vector<int>::const_iterator endit=ind.end()-1;
	complex<T> Psqr=mc.m2(P)-mc.mom(*endit).E();
	vector<int> shifted_legs=ind;
	size_t Phat=get_shifted_ij(mc,shifted_legs,P,Psqr);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	size_t kl=0;
	for(;kl<get_part().c(1).size()-1;kl++){
		indshiftl[kl]=shifted_legs[get_part().c(1)[kl].ind()];
	}
	indshiftl[kl]=mc.insert(-mc[Phat]);
	indshiftr[0]=Phat;
	size_t kr=1;
	const std::vector<plabel>& c2=get_part().c(2);
	for(;kr<get_part().c(2).size();kr++){
		indshiftr[kr]=shifted_legs[c2[kr].ind()];
	}

	// In the case of massive amplitudes the mass is sent passed as an extra momentum
	//  as are any reference vectors needed these need to be tagged on to any amplitude
	//  we have produced.
	indshiftl[maxl]=ind[max];
	indshiftr[maxr]=ind[max];
	indshiftl[maxl+1]=ind[max+1];
	indshiftr[maxr+1]=ind[max+1];

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
	_MESSAGE7(tempres," : From ",left()->eval(mc,indshiftl)," * (",T(1.)/Psqr,") * ",right()->eval(mc,indshiftr));
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval massive_prop ",get_part().get_code()," ----------------------------------");
	return tempres;
#else
	return complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr;
#endif
}

template <class T> complex<T>  Rec_Pair_massive_prop_massless_shift::Rec_Pair_eval(momentum_configuration<T>& mc,const vector<int>& ind){
#if _VERBOSE
	_MESSAGE3("-------------------------Begin Rec_Pair_massive_prop_massless_shift ",get_part().get_code()," -----------------------------------");
#endif
	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom=mc.mom(ind[c1[0].ind()]);
	indshiftl[0]=ind[c1[0].ind()];
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		Psum_mom+=mc.mom(ind[c1[psiter].ind()]);
		indshiftl[psiter]=ind[c1[psiter].ind()];
	}

	// Create the shifted momenta
	size_t indi=ind[i];
	size_t indj=ind[j];
	momentum<complex<T> > ij_mom=PfLLt(mc.Lt(indj),mc.L(indi));
	complex<T> Psqr=Psum_mom.square()-mc.mom(ind[max+1]).E();
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));

	// Now compute the on-shell propagator momentum
	momentum<complex<T> > Phat_mom(Psum_mom+z*ij_mom);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	indshiftr[0]=mc.insert(Phat_mom,_mt_massive);
	indshiftl[maxl-1]=mc.insert(-Phat_mom,_mt_massive);
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		indshiftr[kr]=ind[c2[kr].ind()];
	}

	// Insert the shifted legs
	indshiftr[shifted_ind_i]=mc.insert(mc.L(indi),mc.Lt(indi)-z*mc.Lt(indj));
	indshiftl[shifted_ind_j]=mc.insert(mc.L(indj)+z*mc.L(indi),mc.Lt(indj));

	// In the case of massive amplitudes the mass is sent passed as an extra momentum
	//  as are any reference vectors needed these need to be tagged on to any amplitude
	//  we have produced.
	indshiftl[maxl]=ind[max];
	indshiftr[maxr]=ind[max];
	indshiftl[maxl+1]=ind[max+1];
	indshiftr[maxr+1]=ind[max+1];

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
	_MESSAGE7(tempres," : From ",left()->eval(mc,indshiftl)," * (",T(1.)/Psqr,") * ",right()->eval(mc,indshiftr));
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");
	return tempres;
#else
	return Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
#endif

}

template <class T> complex<T>  Rec_Pair_massive_unshifted::Rec_Pair_eval(momentum_configuration<T>& mc,const vector<int>& ind){
#if _VERBOSE
	_MESSAGE3("-------------------------Begin Rec_Pair_massive_unshifted ",get_part().get_code()," -----------------------------------");
#endif
	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom=mc.mom(ind[c1[0].ind()]);
	indshiftl[0]=ind[c1[0].ind()];
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		Psum_mom+=mc.mom(ind[c1[psiter].ind()]);
		indshiftl[psiter]=ind[c1[psiter].ind()];
	}

	// Create the shifted momenta
	size_t indi=ind[i];
	size_t indj=ind[j];
	momentum<complex<T> > ij_mom=PfLLt(mc.Lt(indj),mc.L(indi));
	complex<T> Psqr=Psum_mom.square();
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));

	// Now compute the on-shell propagator momentum
	// Now compute the on-shell propagator momentum
	Cmom<T> Phat_mom(Psum_mom+z*ij_mom);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	indshiftr[0]=mc.insert(Phat_mom);
	indshiftl[maxl-1]=mc.insert(Phat_mom.L(),-Phat_mom.Lt());
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		indshiftr[kr]=ind[c2[kr].ind()];
	}

	// Insert the shifted legs
	indshiftr[shifted_ind_i]=mc.insert(mc.L(indi),mc.Lt(indi)-z*mc.Lt(indj));
	indshiftl[shifted_ind_j]=mc.insert(mc.L(indj)+z*mc.L(indi),mc.Lt(indj));

	// In the case of massive amplitudes the mass is sent passed as an extra momentum
	//  as are any reference vectors needed these need to be tagged on to any amplitude
	//  we have produced.
	indshiftl[maxl]=ind[max];
	indshiftr[maxr]=ind[max];
	indshiftl[maxl+1]=ind[max+1];
	indshiftr[maxr+1]=ind[max+1];

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
	_MESSAGE7(tempres," : From ",left()->eval(mc,indshiftl)," * (",T(1.)/Psqr,") * ",right()->eval(mc,indshiftr));
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");
	return tempres;
#else
	return Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(mc,indshiftl)*right()->eval(mc,indshiftr)*(fact()->eval(mc,indshiftr)))/Psqr);
#endif

}

#ifndef BH_PUBLIC


ostream& operator<<(ostream& s, Rec_BB& RBB){
	s <<"Unspecified Rec_BB containing "  << RBB.nbr_daughters() << " daughters: " << endl;
	for (int i=0;i<RBB.nbr_daughters();i++){
		s << *(RBB.get_daughter(i+1));
	}
	return s<< "end of unspecified Rec_BB " << endl;
}
ostream& operator<<(ostream& s, Rec_Pair& RP){
	s <<"Rec_Pair: " << endl;
	s << "left:" << *(RP.get_daughter(1));
	return s << "right:" << *(RP.get_daughter(2));
}

ostream& operator<<(ostream& s, Rec_Tree& RT){
	return s<<"Tree for process: "<< RT.get_process() << ";";
}
ostream& operator<<(ostream& s, Known_Rec_Tree& RT){
	return s<<"Known tree for process: "<< RT.get_process() << ";";
}
ostream& operator<<(ostream& s, Unknown_Rec_Tree& RT){
	return s<<"Unknown tree for process: "<< RT.get_process() << ", using ["<< RT.get_i()+1 << ","<< RT.get_j()+1<< "> shift";
}

#endif

Const_Fact_Fn unity_real(1);
Const_Fact_Fn* Const_Fact_Fn::unity=&unity_real;


template complex<R>  Rec_Pair_massive_prop::Rec_Pair_eval(momentum_configuration<R>& mc,const vector<int>& ind);
template complex<RHP>  Rec_Pair_massive_prop::Rec_Pair_eval(momentum_configuration<RHP>& mc,const vector<int>& ind);
template complex<RVHP>  Rec_Pair_massive_prop::Rec_Pair_eval(momentum_configuration<RVHP>& mc,const vector<int>& ind);
#if BH_USE_GMP
template complex<RGMP>  Rec_Pair_massive_prop::Rec_Pair_eval(momentum_configuration<RGMP>& mc,const vector<int>& ind);
#endif


#if BH_USE_GMP
std::complex<RGMP> Const_Fact_Fn::eval_and_increase_precision_if_needed(){
	if (RGMP::get_current_precision()> d_GMP_precision){
		if (d_den==0){
			if (d_num==0){
				//no 'analytical value known', simply increase the precision
				_value_GMP=CGMP(RGMP(_value_GMP.real(),RGMP::get_current_precision()),RGMP(_value_GMP.imag(),RGMP::get_current_precision()));
			} else {
				_value_GMP=RGMP(d_num);
			}
		} else {
			_value_GMP=RGMP(d_num)/RGMP(d_den);
		}
		d_GMP_precision=RGMP::get_current_precision();
	}

	return _value_GMP;
}
#endif

}
