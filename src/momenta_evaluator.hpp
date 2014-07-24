/*
 * momenta_evaluator.hpp
 *
 *  Created on: 10 Aug 2009
 *      Author: daniel
 */

#ifndef MOMENTA_EVALUATOR_HPP_
#define MOMENTA_EVALUATOR_HPP_

#include "index_vector.h"
#include "momenta_evaluator.h"
#include "eval_param.h"
#include "coeff_param.h"
#include "BH_debug.h"

using namespace std;

namespace BH {

namespace cut {

namespace Darren {

template <int TPOINTSBUB,int YPOINTS> template <class T,int CPOINTS,class cutDbase> void  bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::compute_momenta(momentum_configuration<T>& mc,const std::vector<int>& ind, MomentaGridType (&result)[2], MomentaGridType (&resultm)[2],coeffparam<T,CPOINTS>& tp,const cutDbase& cb){
	// Construct K1 by summing over all the legs of c(1) using ind to get their location in the momconf
	momentum<complex<T> > K1sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	for(size_t k1iter=1;k1iter<=cb.corner_size(1);k1iter++){
		K1sum_mom+=mc.mom(ind[cb.corner_ind(1,k1iter)-1]);
	}
	momentum<std::complex<T> > K1=K1sum_mom;
	complex<T> S1=K1.square();

	//Construct chi so that all the coefficients of l^{\mu} are of order 1
	// this is achieved by scaling the vector (1,0,1,0) (we want to avoid a chi proportional to one of the external legs)
	// so that gamma=S1

	momentum<complex<T> > chi_init(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,-1));
	// Now rescale to the actual chi
	Cmom<T>  chic=Cmom<T>(S1/(T(2)*(chi_init*K1))*chi_init);
	complex<T> gammab= S1;
	Cmom<T> K1flatbc=Cmom<T>(K1-chic.P());

	const complex<T>* eval_pts;
	bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_eval_pts(eval_pts);

	const complex<T>* ypoint;
	bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_yeval_pts(ypoint);

	for(int j=0;j<YPOINTS;j++){
		for(int i=0;i<TPOINTSBUB;i++){
			//Construct the cut-momenta
			((result[0])[j])[i]=mc.insert((ypoint[j]/eval_pts[i])*K1flatbc.Lt()+chic.Lt(),eval_pts[i]*K1flatbc.L()+(T(1)-ypoint[j])*chic.L());
			((result[1])[j])[i]=mc.insert((ypoint[j]-T(1))*K1flatbc.Lt()+eval_pts[i]*chic.Lt(),K1flatbc.L()-(ypoint[j]/eval_pts[i])*chic.L());
			((resultm[0])[j])[i]=mc.insert(-mc.Lt(((result[0])[j])[i]),mc.L(((result[0])[j])[i]));
			((resultm[1])[j])[i]=mc.insert(-mc.Lt(((result[1])[j])[i]),mc.L(((result[1])[j])[i]));
		}
	}

	tp.K1=K1;
	tp.S1=S1;
	tp.chic=chic;
	tp.gammab=S1;
	tp.K1flatbc=K1flatbc;

}

template <int TPOINTSBUB,int YPOINTS> template <class T,int CPOINTS,class cutDbase> void  bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::compute_momenta(const eval_param<T>& ep, Cmom<T> (&momenta)[2][YPOINTS][TPOINTSBUB],Cmom<T> (&momenta_m)[2][YPOINTS][TPOINTSBUB] ,coeffparam<T,CPOINTS>& tp,const cutDbase& cb){
	// Construct K1 by summing over all the legs of c(1) using ind to get their location in the momconf
	momentum<complex<T> > K1sum_mom(ep.p(cb.corner_ind(1,1)-1)->P());
	for(size_t k1iter=2;k1iter<=cb.corner_size(1);k1iter++){
		K1sum_mom+=ep.p(cb.corner_ind(1,k1iter)-1)->P();
	}
	momentum<std::complex<T> > K1=K1sum_mom;
	complex<T> S1=K1.square();

	//Construct chi so that all the coefficients of l^{\mu} are of order 1
	// this is achieved by scaling the vector (1,0,1,0) (we want to avoid a chi proportional to one of the external legs)
	// so that gamma=S1
	//	momentum<complex<T> > chi_init(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,-1));
	//	tp.chic=Cmom<T>(tp.S1/(T(2)*(chi_init*tp.K1))*chi_init);

	momentum<complex<T> > chi_init(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,1));
	// Now rescale to the actual chi
	Cmom<T>  chic=Cmom<T>(S1/(T(2)*(chi_init*K1))*chi_init);
	complex<T> gammab= S1;
	Cmom<T> K1flatbc=Cmom<T>(K1-chic.P());

	const complex<T>* eval_pts;
	bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_eval_pts(eval_pts);

	const complex<T>* ypoint;
	bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_yeval_pts(ypoint);

	// Using the eval_param, for now we have to put the incoming momenta into the left or right side


	//We check the error of our computation by looking at the 1/t component of the bubble numerator structure.
	// Due to our parameterization this should always have a coefficient of zero.
    for(int j=0;j<YPOINTS;j++){
    	for(int i=0;i<TPOINTSBUB;i++){
    		//Construct the cut-momenta
    		((momenta[0])[j])[i]=Cmom<T>((ypoint[j]/eval_pts[i])*K1flatbc.Lt()+chic.Lt(),eval_pts[i]*K1flatbc.L()+(T(1)-ypoint[j])*chic.L());
    		((momenta[1])[j])[i]=Cmom<T>((ypoint[j]-T(1))*K1flatbc.Lt()+eval_pts[i]*chic.Lt(),K1flatbc.L()-(ypoint[j]/eval_pts[i])*chic.L());
    		((momenta_m[0])[j])[i]=Cmom<T>(-((momenta[0])[j])[i].Lt(),((momenta[0])[j])[i].L());
    		((momenta_m[1])[j])[i]=Cmom<T>(-((momenta[1])[j])[i].Lt(),((momenta[1])[j])[i].L());
    		BH_DEBUG_MESSAGE4("ypoint[",j,"]: ",ypoint[j]);
    		BH_DEBUG_MESSAGE4("eval_pts[",i,"]: ",eval_pts[i]);
    	}
    }
	tp.K1=K1;
	tp.S1=S1;
	tp.chic=chic;
	tp.gammab=S1;
	tp.K1flatbc=K1flatbc;

}

//In cut_Darren.hpp
	template <class cutDbase> void Do_GenIndicesBub(cutDbase& cut,const vector<int>& ind, vector<int*>& pmoml, vector<int*>& pmmoml, std::vector<int> (&indlst)[2]);

template <class MomentumEvaluator,class CutType,int CPOINTS> template <class T> void Normal_Corner_Tree_Strategy<MomentumEvaluator,CutType,CPOINTS>::fill_trees(
		BH::momentum_configuration<T>& mc,
		std::vector<int> const& ind,
	    std::vector<std::complex<T> >& trees_result_1,
	    std::vector<std::complex<T> >& trees_result_2,
	    coeffparam<T,CPOINTS>& tp,CutType& self
){

	const complex<T>* eval_pts;
    MomentumEvaluator::get_eval_pts(eval_pts);

	const complex<T>* ypoint;
	MomentumEvaluator::get_yeval_pts(ypoint);

	typename MomentumEvaluator::MomentaGridType momenta[2];
	typename MomentumEvaluator::MomentaGridType momentam[2];

	this->template compute_momenta<T,CPOINTS,CutType>(mc,ind,momenta,momentam, tp,self);


	//We check the error of our computation by looking at the 1/t component of the bubble numerator structure.
	// Due to our parameterization this should always have a coefficient of zero.
	std::vector<int*> pmoml(2),pmmoml(2);
	std::vector<int> indlst[2];
	indlst[0].assign(self.corner_size(1)+2,0);
	indlst[1].assign(self.corner_size(2)+2,0);
    Do_GenIndicesBub(self,ind, pmoml, pmmoml,indlst);


    for(int j=0;j<MomentumEvaluator::Ypoints;j++){
            for(int i=0;i<MomentumEvaluator::Tpointsbub;i++){
        			//Construct the cut-momenta
        			*(pmoml[0])=momenta[0][j][i];
        			*(pmoml[1])=momenta[1][j][i];
        			*(pmmoml[0])=momentam[0][j][i];
        			*(pmmoml[1])=momentam[1][j][i];
        			//Construct the two-particle cut at the momenta above
        			int current_index=j*MomentumEvaluator::Tpointsbub+i;
                    trees_result_1[current_index]=self.eval_tree(1,mc,indlst[0]);
                    trees_result_2[current_index]=self.eval_tree(2,mc,indlst[1]);
        }
    }

};

template <class MomentumEvaluator,class CutType,int CPOINTS> template <class T> void Normal_Corner_Tree_Strategy<MomentumEvaluator,CutType,CPOINTS>::fill_trees(
		BH::momentum_configuration<T>& mc,
		std::vector<int> const& ind,
		std::vector<std::complex<T> >& trees_result_1,
		std::vector<std::complex<T> >& trees_result_2,
		std::vector<std::complex<T> >& trees_result_3,
		coeffparam<T,CPOINTS>& tp,CutType& self,int k1leg,int k2leg,int k3leg,int masslessleg_type,int reverse
){
	BH_DEBUG_PRINT(masslessleg_type);

	int momenta[3*MomentumEvaluator::Tpointstri];
	int momentam[3*MomentumEvaluator::Tpointstri];

	this->template compute_momenta<T,CutType>(mc,ind,momenta,momentam, tp,self,k1leg,k2leg,k3leg,masslessleg_type);

	std::vector<int> indlst[3];
	indlst[0].assign(self.corner_size(1)+2,0);
	indlst[1].assign(self.corner_size(2)+2,0);
	indlst[2].assign(self.corner_size(3)+2,0);


    size_t mm;
    for(mm=1;mm<=self.corner_size(1);mm++){
            (indlst[0])[mm]=ind[self.corner_ind(1,mm)-1];
    }
    for(mm=1;mm<=self.corner_size(2);mm++){
            (indlst[1])[mm]=ind[self.corner_ind(2,mm)-1];
    }
    for(mm=1;mm<=self.corner_size(3);mm++){
            (indlst[2])[mm]=ind[self.corner_ind(3,mm)-1];
    }

//    static_cast<triangle_Darren<CutType> >(self);

    if(reverse==1){
	for(int icirc=0;icirc<MomentumEvaluator::Tpointstri;icirc++){
		(indlst[k1leg-1]).back()=momenta[3*icirc+1];
		(indlst[k3leg-1]).back()=momenta[3*icirc+2];
		(indlst[k2leg-1]).back()=momenta[3*icirc+0];

		(indlst[k1leg-1]).front()=momentam[3*icirc+0];
		(indlst[k3leg-1]).front()=momentam[3*icirc+1];
		(indlst[k2leg-1]).front()=momentam[3*icirc+2];

		BH_DEBUG(
				if (icirc==0){
					for (int j=0;j<indlst[0].size();j++){
						_MESSAGE6 ("corner ",0," icirc ",icirc," : ",mc.p(indlst[0][j]).E());
					}
					for (int j=0;j<indlst[1].size();j++){
						_MESSAGE6 ("corner ",0," icirc ",icirc," : ",mc.p(indlst[1][j]).E());
					}
					for (int j=0;j<indlst[2].size();j++){
						_MESSAGE6 ("corner ",0," icirc ",icirc," : ",mc.p(indlst[2][j]).E());
					}
				}
		)


		BH_DEBUG(


			for (int cor=0;cor<3;cor++){
				int n=indlst[cor].size()-2;
				momentum<complex<T> > pcheck;

				for (int j=0;j<indlst[cor].size();j++){
					pcheck+=mc.mom(indlst[cor][j]);
					_MESSAGE6 ("corner ",cor," j ",j," : ",mc.p(indlst[cor][j]).E());
				}
				if (abs(pcheck.E())>1e-12){
					_MESSAGE("pcheck failed!");
					_MESSAGE("===============");
//					_PRINT(p0);
//					_PRINT(p1);
//					_PRINT(p2);
					_PRINT(mc.mom(indlst[cor][0]));
					_PRINT(mc.mom(indlst[cor][n+1]));
//					_PRINT(true_offset[corner]);
//					_PRINT(corner_map[corner]);
					_MESSAGE("===============");
				}
				if (abs(mc.mom(indlst[cor][0])*mc.mom(indlst[cor][0])) > 1e-12){
					_MESSAGE2("not massless: ",mc.mom(indlst[cor][0]));
				}
				if (abs(mc.mom(indlst[cor][n+1])*mc.mom(indlst[cor][n+1])) > 1e-12){
					_MESSAGE2("not massless: ",mc.mom(indlst[cor][n+1]));
				}
			}


		)

		int current_index=icirc;

		trees_result_1[current_index]=self.eval_tree(1,mc,indlst[0]);
		trees_result_2[current_index]=self.eval_tree(2,mc,indlst[1]);
		trees_result_3[current_index]=self.eval_tree(3,mc,indlst[2]);

	}
//	BH_DEBUG(
//			for (int j=0;j<indlst[0].size();j++){
//				_MESSAGE6 ("corner ",1," icirc: ",0," : ",mc.p(indlst[0][j]));
//			}
//		for (int j=0;j<indlst[1].size();j++){
//			_MESSAGE6 ("corner ",2," icirc: ",0," : ",mc.p(indlst[1][j]));
//		}
//		for (int j=0;j<indlst[2].size();j++){
//			_MESSAGE6 ("corner ",3," icirc: ",0," : ",mc.p(indlst[2][j]));
//		}
//	)
    }else{
    	for(int icirc=0;icirc<MomentumEvaluator::Tpointstri;icirc++){
    		(indlst[k1leg-1]).back()=momenta[3*icirc+0];
    		(indlst[k3leg-1]).back()=momenta[3*icirc+1];
    		(indlst[k2leg-1]).back()=momenta[3*icirc+2];

    		(indlst[k1leg-1]).front()=momentam[3*icirc+1];
    		(indlst[k3leg-1]).front()=momentam[3*icirc+2];
    		(indlst[k2leg-1]).front()=momentam[3*icirc+0];


    		int current_index=icirc;

    		trees_result_1[current_index]=self.eval_tree(1,mc,indlst[0]);
    		trees_result_2[current_index]=self.eval_tree(2,mc,indlst[1]);
    		trees_result_3[current_index]=self.eval_tree(3,mc,indlst[2]);
    	}
    }

};


template <class MomentumEvaluator,class CutType,int CPOINTS> template <class T> void Normal_Corner_Tree_Strategy<MomentumEvaluator,CutType,CPOINTS>::fill_trees(
		const eval_param<T>& ep,
        eval_param<T>**& epc,                                                                                                                                                        
		std::vector<std::complex<T> >& trees_result_1,
		std::vector<std::complex<T> >& trees_result_2,
		coeffparam<T,CPOINTS>& tp,CutType& self
){
	Cmom<T> momenta[2][MomentumEvaluator::Ypoints][MomentumEvaluator::Tpointsbub];
 	Cmom<T> momenta_m[2][MomentumEvaluator::Ypoints][MomentumEvaluator::Tpointsbub];

	this->template compute_momenta<T,CPOINTS,CutType>(ep,momenta,momenta_m, tp,self);

	for(int mm=0;mm<self.corner_size(1);mm++){
		epc[0]->set(mm+1,ep.p(self.corner_ind(1,mm+1)-1));
	}
	for(int mm=0;mm<self.corner_size(2);mm++){
		epc[1]->set(mm+1,ep.p(self.corner_ind(2,mm+1)-1));
	}
	
	for(int j=0;j<MomentumEvaluator::Ypoints;j++){
		for(int i=0;i<MomentumEvaluator::Tpointsbub;i++){
			//Construct the cut-momenta
			*(epc[1]->back())=&momenta[0][j][i];   // l[0]
			*(epc[0]->back())=&momenta[1][j][i];  // l[1]
			*(epc[0]->begin())=&momenta_m[0][j][i];  // ml[0]
			*(epc[1]->begin())=&momenta_m[1][j][i];  // ml[1]

			//Construct the two-particle cut at the momenta above
			int current_index=j*MomentumEvaluator::Tpointsbub+i;
			trees_result_1[current_index]=self.eval_tree(1,*epc[0]);
			trees_result_2[current_index]=self.eval_tree(2,*epc[1]);

	    }
	}
};





template <class MomentumEvaluator,class CutType,int CPOINTS> template <class T> void Skeleton_Corner_Tree_Strategy<MomentumEvaluator,CutType,CPOINTS>::fill_trees(
		BH::momentum_configuration<T>& mc,
		std::vector<int> const& ind,
		std::vector<std::complex<T> >& trees_result_1,
		std::vector<std::complex<T> >& trees_result_2,
		coeffparam<T,CPOINTS>& tp,CutType& self
){
	_MESSAGE("USING SKELETON");
	const complex<T>* eval_pts;
	MomentumEvaluator::get_eval_pts(eval_pts);

	const complex<T>* ypoint;
	MomentumEvaluator::get_yeval_pts(ypoint);

	typename MomentumEvaluator::MomentaGridType momenta[2];
	typename MomentumEvaluator::MomentaGridType momentam[2];


	//needed to pass the coeffparam to the rest of the code. It should be removes eventually
	this->template compute_momenta<T,CPOINTS,CutType>(mc,ind,momenta,momentam, tp,self);


	Index_Vector IV(ind);
	int pos;
	typename CutType::TargetSkeletonType* BS = self.d_skeleton_user.get_skeleton_without_construction(IV.get_permutation_code(),self,pos);

	if ( BS == 0 ){
		_MESSAGE("The skeleton should already be there !");
		_PRINT(IV);
		CutType::TargetSkeletonType::FactoryType::s_default_SF->print_state();
	}
	BS->prepare(mc);

// this is only needed if there are more than one case per permutation
//	//these are the indices of the trees in the tree manager of the skeleton
//	int tree_index_1 = self.get_tree_pos(IV.get_permutation_code(),Skeleton::cached_corner_tree_info::not_reversed,1);
//	int tree_index_2 = self.get_tree_pos(IV.get_permutation_code(),Skeleton::cached_corner_tree_info::not_reversed,1);
	//these are the indices of the trees in the tree manager of the skeleton
	int tree_index_1 = self.get_tree_pos(IV.get_permutation_code(),1);
	int tree_index_2 = self.get_tree_pos(IV.get_permutation_code(),2);


	for(int j=0;j<MomentumEvaluator::Ypoints;j++){
		for(int i=0;i<MomentumEvaluator::Tpointsbub;i++){
			int current_index=j*MomentumEvaluator::Tpointsbub+i;

			trees_result_1[current_index]=BS->d_ctm.d_computed_tree_values[0][tree_index_1][current_index];
			trees_result_2[current_index]=BS->d_ctm.d_computed_tree_values[1][tree_index_2][current_index];
		}
	}

};

template <class MomentumEvaluator,class CutType,int CPOINTS> template <class T> void Skeleton_Corner_Tree_Strategy<MomentumEvaluator,CutType,CPOINTS>::fill_trees(
		BH::momentum_configuration<T>& mc,
		std::vector<int> const& ind,
		std::vector<std::complex<T> >& trees_result_1,
		std::vector<std::complex<T> >& trees_result_2,
		std::vector<std::complex<T> >& trees_result_3,
		coeffparam<T,CPOINTS>& tp,CutType& self,int k1leg,int k2leg,int k3leg,int masslessleg_type,int reverse
){
//	typename MomentumEvaluator::MomentaGridType momenta[3];
//	typename MomentumEvaluator::MomentaGridType momentam[3];

	int momenta[3*MomentumEvaluator::Tpointstri];
	int momentam[3*MomentumEvaluator::Tpointstri];


	Index_Vector IV(ind);
	int pos;
	typename CutType::TargetSkeletonType* BS = self.d_skeleton_user.get_skeleton_without_construction(IV.get_permutation_code(),self,pos);
	int corner_map=self.d_skeleton_user.get_corner_map(pos);

	BH_DEBUG_PRINT(self);
	BH_DEBUG_MESSAGE2("USING SKELETON: ",*BS);
	BH_DEBUG_PRINT(masslessleg_type);


	int local_k1leg;
	int local_k2leg;
	int local_k3leg;

	switch(corner_map){
	case 123: {
		local_k1leg=BS->d_k1leg;
		local_k2leg=BS->d_k2leg;
		local_k3leg=BS->d_k3leg;
	} break;
	case 213: {
//		local_k1leg=BS->d_k3leg;
//		local_k2leg=BS->d_k1leg;
//		local_k3leg=BS->d_k2leg;

		local_k1leg=BS->d_k2leg;
		local_k2leg=BS->d_k1leg;
		local_k3leg=BS->d_k3leg;
	} break;
	case 132: {
		local_k1leg=BS->d_k1leg;
		local_k2leg=BS->d_k3leg;
		local_k3leg=BS->d_k2leg;
	} break;
	default: {
		local_k1leg=BS->d_k1leg;
		local_k2leg=BS->d_k2leg;
		local_k3leg=BS->d_k3leg;
	} break;
	}

	BH_DEBUG(
			_PRINT(typeid(CutType).name());
			BH_DEBUG_MESSAGE( "%%%%  corner map: ");
			_PRINT(corner_map);
			_PRINT(BS->d_k1leg);
			_PRINT(BS->d_k2leg);
			_PRINT(BS->d_k3leg);
			_PRINT(local_k1leg);
			_PRINT(local_k2leg);
			_PRINT(local_k3leg);
			_PRINT(k1leg);
			_PRINT(k2leg);
			_PRINT(k3leg);

			momentum<complex<T> > K[3];
			for (int corner =1; corner <=3;corner++){
				for (int ii=1;ii<=self.corner_size(corner);ii++){
					K[corner-1]+=mc.mom(self.corner_ind(corner,ii));
				}
			}
			_PRINT(K[0]);
			_PRINT(K[1]);
			_PRINT(K[2]);
			_PRINT(BS->get_external_momentum(1));
			_PRINT(BS->get_external_momentum(2));
			_PRINT(BS->get_external_momentum(3));

			_MESSAGE("%%%%");

	)


	// use local k1leg, but use the masslessleg from the cut because the skeleton does not know about helicities
	this->template compute_momenta<T,CutType>(mc,ind,momenta,momentam, tp,self, local_k1leg, local_k2leg,local_k3leg, masslessleg_type);



	if ( BS == 0 ){
		_MESSAGE("The skeleton should already be there !");
		_PRINT(IV);
		CutType::TargetSkeletonType::FactoryType::s_default_SF->print_state();
		throw;
	}
// the skeleton should already have been prepared
//	BS->prepare(mc);

// this is only needed if there are more than one case per permutation
//	//these are the indices of the trees in the tree manager of the skeleton
//	int tree_index_1 = self.get_tree_pos(IV.get_permutation_code(),Skeleton::cached_corner_tree_info::not_reversed,1);
//	int tree_index_2 = self.get_tree_pos(IV.get_permutation_code(),Skeleton::cached_corner_tree_info::not_reversed,1);
	//these are the indices of the trees in the tree manager of the skeleton
	int tree_index_1 = self.get_tree_pos(IV.get_permutation_code(),1);
	int tree_index_2 = self.get_tree_pos(IV.get_permutation_code(),2);
	int tree_index_3 = self.get_tree_pos(IV.get_permutation_code(),3);

	for(int icirc=0;icirc<MomentumEvaluator::Tpointstri;icirc++){

		int current_index=icirc;

		trees_result_1[current_index]=BS->d_ctm.d_computed_tree_values[0][tree_index_1][current_index];
		trees_result_2[current_index]=BS->d_ctm.d_computed_tree_values[1][tree_index_2][current_index];
		trees_result_3[current_index]=BS->d_ctm.d_computed_tree_values[2][tree_index_3][current_index];
	}

};

template <class MomentumEvaluator,class CutType,int CPOINTS> template <class T> void Skeleton_Corner_Tree_Strategy<MomentumEvaluator,CutType,CPOINTS>::fill_trees(
		const eval_param<T>& ep,
        eval_param<T>**& epc,
		std::vector<std::complex<T> >& trees_result_1,
		std::vector<std::complex<T> >& trees_result_2,
		coeffparam<T,CPOINTS>& tp,CutType& self
){
	BH_DEBUG_MESSAGE("USING SKELETON");
	const complex<T>* eval_pts;
	MomentumEvaluator::get_eval_pts(eval_pts);

	const complex<T>* ypoint;
	MomentumEvaluator::get_yeval_pts(ypoint);

	Cmom<T> momenta[2][MomentumEvaluator::Ypoints][MomentumEvaluator::Tpointsbub];
	Cmom<T> momenta_m[2][MomentumEvaluator::Ypoints][MomentumEvaluator::Tpointsbub];


	//needed to pass the coeffparam to the rest of the code. It should be removes eventually
	this->template compute_momenta<T,CPOINTS,CutType>(ep,momenta,momenta_m, tp,self);

	const eval_param_with_permutation<T>* epwp=dynamic_cast<const eval_param_with_permutation<T>*>(&ep);
	long perm_code = epwp->get_permutation_code();
	int pos;
	typename CutType::TargetSkeletonType* BS = self.d_skeleton_user.get_skeleton_without_construction(perm_code,self,pos);

	if ( BS == 0 ){
		_MESSAGE("The skeleton should already be there !");
		_PRINT(perm_code);
		CutType::TargetSkeletonType::FactoryType::s_default_SF->print_state();
	}
	BS->prepare(*(epwp->get_mc()));

// this is only needed if there are more than one case per permutation
//	//these are the indices of the trees in the tree manager of the skeleton
//	int tree_index_1 = self.get_tree_pos(IV.get_permutation_code(),Skeleton::cached_corner_tree_info::not_reversed,1);
//	int tree_index_2 = self.get_tree_pos(IV.get_permutation_code(),Skeleton::cached_corner_tree_info::not_reversed,1);
	//these are the indices of the trees in the tree manager of the skeleton
	int tree_index_1 = self.get_tree_pos(perm_code,1);
	int tree_index_2 = self.get_tree_pos(perm_code,2);


	for(int j=0;j<MomentumEvaluator::Ypoints;j++){
		for(int i=0;i<MomentumEvaluator::Tpointsbub;i++){
			int current_index=j*MomentumEvaluator::Tpointsbub+i;

			trees_result_1[current_index]=BS->d_ctm.d_computed_tree_values[0][tree_index_1][current_index];
			trees_result_2[current_index]=BS->d_ctm.d_computed_tree_values[1][tree_index_2][current_index];
		}
	}

};


template <class T> inline RVHP convertToRVHP(const T&);
#if BH_USE_GMP
template <> inline RVHP convertToRVHP(const RGMP& g){ return g.to_double(); };  // probably it is sufficient
#endif
template <class T> inline RVHP convertToRVHP(const T& t){ return RVHP(t); };

template <int YPOINTS,int TPOINTSBUB,class MomentumEvaluator>  template <class T> complex<T> General_Bubble_Combiner<YPOINTS,TPOINTSBUB,MomentumEvaluator>::combine(std::complex<T> (&ampl)[YPOINTS],std::complex<T> (&ampl_err)[YPOINTS],std::complex<T>& trisuby,std::complex<T>& amp_error,const RVHP& tri_sub_acc,RVHP& accuracy){
	
	// Extract the coefficients of y
	const complex<T>* ymatrixpoint;
	MomentumEvaluator::get_y_matrix_eval_pts(ymatrixpoint);
	//complex<T> result_y[YPOINTS]={complex<T>(0,0)};
	complex<T> result(0,0);

	for(int i=0;i<YPOINTS;i++){
		for(int j=0;j<YPOINTS;j++){
			result+=ymatrixpoint[j+i*(YPOINTS)]*ampl[j];
			//result_y[i]+=ymatrixpoint[j+i*(YPOINTS)]*ampl[j];
			//_MESSAGE8("i. ",i,", j. ",j,"==",ymatrixpoint[j+i*(YPOINTS)]," * ",ampl[j]);
		}
	}
	
	//_PRINT(result);
//	for(int i=0;i<YPOINTS;i++){
//		_MESSAGE4("y:",i," = ",result_y[i]);
//	}

	//For the error estimate we only care about the ty^{YPOINTS-1} contribution
	complex<T> err_ycoeff(0,0);
	for(int j=0;j<YPOINTS;j++){
		err_ycoeff+=ymatrixpoint[j+(YPOINTS-1)*YPOINTS]*ampl_err[j];
	}

	// Store the accuracy
	RVHP comp_acc=convertToRVHP(T(1)/abs(amp_error-err_ycoeff));

	if(comp_acc<tri_sub_acc){
		accuracy=comp_acc;
	}
	else{
		accuracy=tri_sub_acc;
	}

	//_PRINT(trisuby);

	return complex<T>(0,-1)*(result-trisuby)/T(TPOINTSBUB);
}

template <int YPOINTS,int TPOINTSBUB,class MomentumEvaluator>  template <class T> complex<T> Normal_Bubble_Combiner<YPOINTS,TPOINTSBUB,MomentumEvaluator>::combine(std::complex<T> (&ampl)[YPOINTS],std::complex<T> (&ampl_err)[YPOINTS],std::complex<T>& trisuby,std::complex<T>& amp_error,const RVHP& tri_sub_acc,RVHP& accuracy){

	// Store the accuracy
	RVHP comp_acc=convertToRVHP<T>(T(1)/abs(amp_error-ampl_err[0]));

	if(comp_acc<tri_sub_acc){
		accuracy=comp_acc;
	}
	else{
		accuracy=tri_sub_acc;
	}

//	_MESSAGE6("bubbles bits : ",ampl[0],", ",T(3)*ampl[1],", ",trisuby);
//	_MESSAGE6("BUB : ",complex<T>(0,-1)*(ampl[0]+T(3)*ampl[1]),"+",complex<T>(0,1)*trisuby,"/",T(16));
	BH_DEBUG_PRINT(ampl[0]);
	BH_DEBUG_PRINT(ampl[1]);

	complex<T> result=(complex<T>(0,-1)*(ampl[0]+T(3)*ampl[1])+trisuby)/T(16);

	return result;
}



template <int CPOINTS, int TPOINTSTRI> template <class T,class cutDbase> void triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::compute_momenta(
		momentum_configuration<T>& mc,
		const std::vector<int>& ind ,
		int (&result)[3*TPOINTSTRI],
		int (&resultm)[3*TPOINTSTRI],
		coeffparam<T,CPOINTS>& tp2,    
		const cutDbase& cb,
		int k1leg,int k2leg,int k3leg,int masslessleg_type){
	momentum<complex<T> > K1sum_mom(mc.mom(ind[cb.corner_ind(k1leg,1)-1]));
	for(size_t kiter=2;kiter<=cb.corner_size(k1leg);kiter++){
		K1sum_mom+=mc.mom(ind[cb.corner_ind(k1leg,kiter)-1]);
	}
	momentum<complex<T> > K2sum_mom(-mc.mom(ind[cb.corner_ind(k2leg,1)-1]));
	for(size_t kiter=2;kiter<=cb.corner_size(k2leg);kiter++){
		K2sum_mom-=mc.mom(ind[cb.corner_ind(k2leg,kiter)-1]);
	}
	tp2.K1=K1sum_mom;
	tp2.K2=K2sum_mom;
	tp2.S1=K1sum_mom.square();
	tp2.S2=K2sum_mom.square();

	complex<T> K1K2=K1sum_mom*K2sum_mom;
	tp2.alp1=complex<T>(1,0);

	BH_DEBUG_PRINT(K1sum_mom);
	BH_DEBUG_PRINT(K2sum_mom);

    if(masslessleg_type!=0){
    	BH_DEBUG_MESSAGE("masslessleg_type!=0");
    	//K2 is massless so create spinors for the K2 leg
    	Cmom<T> K2m(K2sum_mom);
    	// For gamma we have
    	complex<T> gammap(T(2)*K1K2);
    	tp2.f1=complex<T>(1,0);
    	tp2.f2=gammap/tp2.S1;
    	tp2.alp2=complex<T>(0,0);
    	tp2.gamma=tp2.S1;

    	complex<T> gamfac(tp2.S1/gammap);
    	tp2.K1flatc=Cmom<T>(K1sum_mom-gamfac*K2sum_mom);
    	if(masslessleg_type>0){
        	BH_DEBUG_MESSAGE("masslessleg_type>0");
    		// Rescale the spinors of K2 correctly
       		tp2.K2flatc=Cmom<T>(K2m.Lt(),gamfac*K2m.L());

       		//Set up the rescaling factors for l1 and l2
    		tp2.vec_alp[0]=tp2.alp1;
    		tp2.vec_alp[1]=tp2.alp2;
    		tp2.vec_alp[2]=complex<T>(-1,0);
    		tp2.vec_alp[3]=complex<T>(0,0);
    		tp2.vec_alp[4]=complex<T>(0,0);
    		tp2.vec_alp[5]=(T(1)-gammap/tp2.S1);

    		tp2.vec1c=Cmom<T>(tp2.K2flatc.Lt(),tp2.K1flatc.L());
    		tp2.vec2c=Cmom<T>(tp2.K1flatc.Lt(),tp2.K2flatc.L());
    	}
    	else{
        	BH_DEBUG_MESSAGE("masslessleg_type<0");
    		// Rescale the spinors of K2 correctly
    		tp2.K2flatc=Cmom<T>(gamfac*K2m.Lt(),K2m.L());

       		//Set up the rescaling factors for l1 and l2
    		tp2.vec_alp[0]=tp2.alp2;
    		tp2.vec_alp[1]=tp2.alp1;
    		tp2.vec_alp[2]=complex<T>(0,0);
    		tp2.vec_alp[3]=complex<T>(-1,0);
    		tp2.vec_alp[4]=(T(1)-gammap/tp2.S1);
    		tp2.vec_alp[5]=complex<T>(0,0);

    		tp2.vec1c=Cmom<T>(tp2.K1flatc.Lt(),tp2.K2flatc.L());
    		tp2.vec2c=Cmom<T>(tp2.K2flatc.Lt(),tp2.K1flatc.L());
    	}
    }
    else{
    	BH_DEBUG_MESSAGE("masslessleg_type==0");
    	// Compute the gamma of K1flatp and K2flatp so we can define f1 and f2.
    	complex<T> gammap=(K1K2+sqrt(pow(K1K2,2)-tp2.S1*tp2.S2));
    	tp2.f1=(tp2.S1*tp2.S2-pow(gammap,2))/(tp2.S2*(tp2.S1-gammap));
    	tp2.f2=(tp2.S1*tp2.S2-pow(gammap,2))/(tp2.S1*(tp2.S2-gammap));
    	tp2.alp2=complex<T>(1,0);
    	// Now convert gammap into gamma for the rest of the computation
    	tp2.gamma=gammap/(tp2.f1*tp2.f2);
    	tp2.vec_alp[0]=tp2.alp1;
    	tp2.vec_alp[1]=tp2.alp2;
 	    //Set up the rescaling factors for l1 and l2
    	tp2.vec_alp[3]=(T(1)-(gammap*gammap-tp2.S1*tp2.S2)/(tp2.S2*(gammap-tp2.S1)));
    	tp2.vec_alp[2]=(T(1)-(gammap*gammap-tp2.S1*tp2.S2)/(gammap*(gammap-tp2.S2)));
    	tp2.vec_alp[5]=(T(1)-(gammap*gammap-tp2.S1*tp2.S2)/(gammap*(gammap-tp2.S1)));
    	tp2.vec_alp[4]=(T(1)-(gammap*gammap-tp2.S1*tp2.S2)/(tp2.S1*(gammap-tp2.S2)));

    	//Create the unscaled K1flat and K2flat
    	complex<T> denalpha=T(1)/(pow(gammap,2)-tp2.S1*tp2.S2);
    	Cmom<T> K1flatcI(denalpha*gammap*(gammap*K1sum_mom-tp2.S1*K2sum_mom));
    	Cmom<T> K2flatcI(denalpha*gammap*(gammap*K2sum_mom-tp2.S2*K1sum_mom));

    	//Now rescale one of the spinors to increase numerical stability
    	tp2.K1flatc=Cmom<T>(K1flatcI.Lt(),(T(1)/tp2.f1)*K1flatcI.L());
    	tp2.K2flatc=Cmom<T>((T(1)/tp2.f2)*K2flatcI.Lt(),K2flatcI.L());

    	// Now set up the vectors for the last two terms in the l basis
    	tp2.vec1c=Cmom<T>(tp2.K1flatc.Lt(),tp2.K2flatc.L());
 	  	tp2.vec2c=Cmom<T>(tp2.K2flatc.Lt(),tp2.K1flatc.L());

    }

	BH_DEBUG_MESSAGE("in momentum_evaluator");
	BH_DEBUG_PRINT(tp2.vec1c);
	BH_DEBUG_PRINT(tp2.vec2c);
	BH_DEBUG_PRINT(tp2.K1flatc);
	BH_DEBUG_PRINT(tp2.K2flatc);


	// Get the points to evaluate at the correct precision
	const complex<T>* eval_pts;
	box_Darren_eval_points<TPOINTSTRI>::get_eval_pts(eval_pts);

	BH_DEBUG_PRINT(tp2.vec_alp[1]);
	BH_DEBUG_PRINT(tp2.vec2c.Lt());
	BH_DEBUG_PRINT(eval_pts[0]);
	BH_DEBUG_PRINT(tp2.vec1c.Lt());
	BH_DEBUG_PRINT(tp2.vec1c.L());
	BH_DEBUG_PRINT(eval_pts[0]);
	BH_DEBUG_PRINT(tp2.vec_alp[0]);
	BH_DEBUG_PRINT(tp2.vec2c.L());

	for(int icirc=0; icirc<TPOINTSTRI; icirc++){
		//Construct the cut-momenta simply by recomputing the triple cut at each pole
		//  this is the slowest method as it does not reuse any information
		result[3*icirc+0]=mc.insert(tp2.vec_alp[1]*tp2.vec2c.Lt()*(T(1)/eval_pts[icirc])+tp2.vec1c.Lt(),tp2.vec1c.L()*eval_pts[icirc]+tp2.vec_alp[0]*tp2.vec2c.L());
		result[3*icirc+1]=mc.insert(tp2.vec_alp[2]*tp2.vec2c.Lt()*(T(1)/eval_pts[icirc])+tp2.vec1c.Lt(),tp2.vec1c.L()*eval_pts[icirc]+tp2.vec_alp[3]*tp2.vec2c.L());
		result[3*icirc+2]=mc.insert(tp2.vec_alp[4]*tp2.vec2c.Lt()*(T(1)/eval_pts[icirc])+tp2.vec1c.Lt(),tp2.vec1c.L()*eval_pts[icirc]+tp2.vec_alp[5]*tp2.vec2c.L());
		resultm[3*icirc+0]=mc.insert(-mc.Lt(result[3*icirc+0]),mc.L(result[3*icirc+0]));
		resultm[3*icirc+1]=mc.insert(-mc.Lt(result[3*icirc+1]),mc.L(result[3*icirc+1]));
		resultm[3*icirc+2]=mc.insert(-mc.Lt(result[3*icirc+2]),mc.L(result[3*icirc+2]));


	}
	BH_DEBUG_PRINT(mc.mom(result[0]));

};






} /* Darren */

} /* cut */

} /* BH */



#endif /* MOMENTA_EVALUATOR_HPP_ */
