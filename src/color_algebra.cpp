/*
 * color_algebra.cpp
 *
 */


#ifndef COLOR_ALGEBRA_H_
#define COLOR_ALGEBRA_H_



#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "color_algebra.h"
#include "settings.h"

using namespace std;
namespace BH {

#if BH_USE_GMP
#define REPEAT_ARG(X) X,X,X,X
#else
#define REPEAT_ARG(X) X,X,X
#endif

color_index::color_index(std::string ct,int lb) : label(lb) {
	if(ct.compare("adjoint")==0){
		kind=adjoint;
	}
	else if(ct.compare("fundamental")==0){
		kind=fundamental;
	}
	else{
		std::cout<<" index type for: ("<<ct<<","<<lb<<") declared as 'unknown', possibilities are 'adjoint' or 'fundamental'"<<endl;
		kind=unknown;
	}
}

color_index::color_index(){ label=-1; kind=unknown; };

std::ostream& operator<<(std::ostream& s, const color_index& ci){
	if(ci.kind==ci.fundamental)
		return s<<"i"<<ci.label;
	else if(ci.kind==ci.adjoint)
		return s<<"a"<<ci.label;
	// unknown!
	else
		return s<<"unknown"<<ci.label;
}


adjoint_index::adjoint_index(int lb) {
	label=lb;
	kind=adjoint;
}

fundamental_index::fundamental_index(int lb) {
	label=lb;
	kind=fundamental;
}

void print_scs_const_error(std::vector<color_index*> cis){
	int cissize=cis.size();
	std::cout<<"color_string construction failed: inconsistent use of "<<endl
		<<"color indices in: {";
	for(int kk=0;kk<(cissize-1);kk++){
		std::cout<<cis[kk];	std::cout<<",";
	}
	std::cout<<cis[cissize-1]; std::cout<<"}\n";
}

color_string::color_string(){ adjoint_degree=0; fundamental_degree=0; };

color_string::color_string(std::vector<color_index*> cis): adjoint_degree(0), fundamental_degree(0)
{
	int cissize=cis.size();
//	std::cout<<"cissize: "<<cissize<<endl;
	if ( cissize>0 ) {
		if(cis[0]->kind==color_index::fundamental){
			if(cis[cissize-1]->kind==color_index::fundamental){
				fundamental_degree++;
				fundamental_color_labels.push_back(cis[0]->label);
				fundamental_color_indices.push_back(*cis[0]);
				fundamental_color_labels.push_back(cis[cissize-1]->label);
				fundamental_color_indices.push_back(*cis[cissize-1]);
				for(int kk=1;kk<(cissize-1);kk++){
					if(cis[kk]->kind==color_index::adjoint){
						adjoint_degree++;
						adjoint_color_labels.push_back(cis[kk]->label);
						adjoint_color_indices.push_back(*cis[kk]);
					}
					else{
						print_scs_const_error(cis);
						break;
					}
				}
			}
			else{
				print_scs_const_error(cis);
				std::cout<<"returned identity!"<<endl;
			}
		}
		else{
			for(int kk=0;kk<cissize;kk++){
				if(cis[kk]->kind==color_index::adjoint){
					adjoint_degree++;
					adjoint_color_labels.push_back(cis[kk]->label);
					adjoint_color_indices.push_back(*cis[kk]);
				}
				else{
					print_scs_const_error(cis);
					break;
				}
			}
		}
	}
}

//returns true if it order changes, false otherwise.
bool color_string::rotate_trace(){

	if(fundamental_degree==0&&adjoint_degree>1){
		std::vector<size_t>::iterator min_iter;
		min_iter=min_element(adjoint_color_labels.begin(),adjoint_color_labels.end());
		size_t offset=min_iter-adjoint_color_labels.begin();
		
		if(offset>0){	
			rotate(adjoint_color_labels.begin(),min_iter,adjoint_color_labels.end());
			rotate(adjoint_color_indices.begin(),adjoint_color_indices.begin()+offset,adjoint_color_indices.end());
		
			return true;
		}
	}
	return false;
}

color_string::color_string(const color_string& scs){
	adjoint_degree=scs.adjoint_degree;
	fundamental_degree=scs.fundamental_degree;
	adjoint_color_indices=scs.adjoint_color_indices;
	fundamental_color_indices=scs.fundamental_color_indices;
	adjoint_color_labels=scs.adjoint_color_labels;
	fundamental_color_labels=scs.fundamental_color_labels;
}



color_string& color_string::operator=(const color_string& scs){
	adjoint_degree=scs.adjoint_degree;
	fundamental_degree=scs.fundamental_degree;
	adjoint_color_indices=scs.adjoint_color_indices;
	fundamental_color_indices=scs.fundamental_color_indices;
	adjoint_color_labels=scs.adjoint_color_labels;
	fundamental_color_labels=scs.fundamental_color_labels;
	return *this;
}

color_string::color_string(const color_index& ci1, const color_index& ci2,
				const color_index& ci3,const color_index& ci4,const color_index& ci5,
				const color_index& ci6,const color_index& ci7,const color_index& ci8){
	vector<color_index> vci;
		vci.push_back(ci1);	vci.push_back(ci2);	vci.push_back(ci3);	vci.push_back(ci4);
		vci.push_back(ci5);	vci.push_back(ci6);	vci.push_back(ci7);	vci.push_back(ci8);
	vector<color_index*> vcip;
	for(int kk=0;kk<8;kk++){
		if(vci[kk].label!=-1){
			vcip.push_back(&vci[kk]);
		}
		else
			break;
	}
	color_string local_scs(vcip);

	this->adjoint_degree=local_scs.adjoint_degree;
	this->fundamental_degree=local_scs.fundamental_degree;
	this->adjoint_color_indices=local_scs.adjoint_color_indices;
	this->fundamental_color_indices=local_scs.fundamental_color_indices;
	this->adjoint_color_labels=local_scs.adjoint_color_labels;
	this->fundamental_color_labels=local_scs.fundamental_color_labels;
}

std::ostream& operator<<(std::ostream& s, const color_string& cs){
	if(cs.adjoint_degree==0 && cs.fundamental_degree==0){
		s<<"Id()";
	}
	else if(cs.adjoint_degree==0){
		s<<"delta";
		s<<"_{"<<cs.fundamental_color_indices[0];
		s<<"}^{"<<cs.fundamental_color_indices[1];
		s<<"}";
	}
	else if(cs.fundamental_degree==0){
		s<<"Tr(T^{";
		for(int kk=0;kk<(cs.adjoint_degree-1);kk++){
			s<<cs.adjoint_color_indices[kk]<<",";
		}	
		s<<cs.adjoint_color_indices[cs.adjoint_degree-1]<<"})";
	}
	else{
		s<<"(T^{";
		for(int kk=0;kk<(cs.adjoint_degree-1);kk++){
			s<<cs.adjoint_color_indices[kk]<<",";
		}	
		s<<cs.adjoint_color_indices[cs.adjoint_degree-1]<<"})";
		s<<"_{"<<cs.fundamental_color_indices[0];
		s<<"}^{"<<cs.fundamental_color_indices[1];
		s<<"}";
	}
	return s;	
}


int color_string::conjugate(){
	//return minus sign when odd number of adjoint matrices	
	int sign(1);
	//CHECK IF WE NEED THIS SIGN: if(adjoint_degree%2==1){sign=-1;};
	//reverse adjoint indices and exchange fundamental indices	
	std::vector<size_t>::iterator begin,end;
	std::vector<color_index>::iterator begin_ind,end_ind;

	begin=adjoint_color_labels.begin();
	end=adjoint_color_labels.end();
	begin_ind=adjoint_color_indices.begin();
	end_ind=adjoint_color_indices.end();

	std::reverse(begin,end);
	std::reverse(begin_ind,end_ind);

	begin=fundamental_color_labels.begin();
	end=fundamental_color_labels.end();
	begin_ind=fundamental_color_indices.begin();
	end_ind=fundamental_color_indices.end();

	std::reverse(begin,end);
	std::reverse(begin_ind,end_ind);

	return(sign);
}

bool operator==(const color_string& cs1, const color_string& cs2){
	bool compare;
	compare=(cs1.fundamental_color_labels==cs2.fundamental_color_labels&&
		cs1.adjoint_color_labels==cs2.adjoint_color_labels);
	return compare;
}

//specified ordering
bool operator<(const color_string& cs1, const color_string& cs2){
	if(cs1==cs2){return false;}
	else if(cs1.fundamental_degree<cs2.fundamental_degree){return true;}
	else if(cs1.fundamental_degree>cs2.fundamental_degree){return false;}
	else if(cs1.fundamental_degree>0&&cs1.fundamental_degree==cs2.fundamental_degree){
		if(cs1.fundamental_color_labels[0]<cs2.fundamental_color_labels[0]){return true;}
		else if(cs1.fundamental_color_labels[0]>cs2.fundamental_color_labels[0]){return false;}
		else if(cs1.fundamental_color_labels[0]==cs2.fundamental_color_labels[0]){
			if(cs1.fundamental_color_labels[1]<cs2.fundamental_color_labels[1]){return true;}
			else if(cs1.fundamental_color_labels[1]>cs2.fundamental_color_labels[1]){return false;}
		}
	}
	else if(cs1.adjoint_degree<cs2.adjoint_degree){return true;}
	else if(cs1.adjoint_degree>cs2.adjoint_degree){return false;}
	else if(cs1.adjoint_degree>0&&cs1.adjoint_degree==cs2.adjoint_degree){
		for(int k=0;k<cs1.adjoint_degree;k++){
			if(cs1.adjoint_color_labels[k]<cs2.adjoint_color_labels[k]){return true;}
			else if(cs1.adjoint_color_labels[k]>cs2.adjoint_color_labels[k]){return false;}
		}
	}
	else{return false;}
}

bool compare_cs(color_string* cs1, color_string* cs2){
	return ((*cs1)<(*cs2));
}

color_tensor* color_string::simplify(){
// to be done: color = color //. {TrC[a__]:>TrC@@RotateLeft[{a},Position[{a},Sort[{a}][[1]]][[1,1]]-1]};
// notice: The Id color_sting is treated as one and _not_ Tr(one)=Nc

if(this->rotate_trace()){
	color_constant one(1,0);
	color_tensor ct2(one,*this);
	return( new color_tensor(ct2));
}

if(fundamental_degree>0){
  	//color = color /. IndexDeltaF[i_,i_]-> Nc;
	if(fundamental_color_labels[0]==fundamental_color_labels[1]&&adjoint_color_indices.size()==1){
		std::vector<color_index*> vci2_pointers;
		color_constant zero(0,0);
		color_string cs2(vci2_pointers);
		color_tensor ct2(zero,cs2);
		return( new color_tensor(ct2));
	}
  	//color = color /. TraceF[x_,{i_,i_}] :> TrC@@x;
	if(fundamental_color_labels[0]==fundamental_color_labels[1]&&adjoint_color_indices.size()>0){
		std::vector<color_index*> vci2_pointers;
		for(int kk=0;kk<adjoint_color_indices.size();kk++){vci2_pointers.push_back(&(adjoint_color_indices[kk]));}
		color_constant one(1,0);
		color_string cs2(vci2_pointers);
		color_tensor ct2(one,cs2);
		return( new color_tensor(ct2));
	}
  	//color = color /. {TraceF[{},{i1_,i1_}] -> Nc};
	if(fundamental_color_labels[0]==fundamental_color_labels[1]&&adjoint_color_indices.size()==0){
		std::vector<color_index*> vci2_pointers;
		color_constant Nc(1,1);
		color_string cs2(vci2_pointers);
		color_tensor ct2(Nc,cs2);
		return( new color_tensor(ct2));
	}
  	//color = color //. {TraceF[{x___, a_, y___, a_,z___},{i1_,i2_}] -> (TraceF[{x, z},{i1,i2}]*TrC[y] - 1/Nc*TraceF[{x, y, z},{i1,i2}])};
	if(adjoint_degree>2){
		for(int kk=0;kk<adjoint_degree-1;kk++){
			for(int ll=kk+1;ll<adjoint_degree;ll++){
				if(adjoint_color_labels[kk]==adjoint_color_labels[ll]){
					std::vector<color_index*> vci2_1_pointers, vci2_2_pointers, vci3_pointers;
					color_constant cc2_1(1,0);
					color_constant cc2_2(1,0);
					color_constant minus_over_Nc(-1,-1);

					vci2_1_pointers.push_back(&(fundamental_color_indices[0]));
					for(int i=0;i<kk;i++){vci2_1_pointers.push_back(&(adjoint_color_indices[i]));}
					for(int i=ll+1;i<adjoint_degree;i++){vci2_1_pointers.push_back(&(adjoint_color_indices[i]));}
					vci2_1_pointers.push_back(&(fundamental_color_indices[1]));
					switch(ll-kk){
						//y has no indices
						case 1: {color_constant Nc(1,1);cc2_2=Nc;} break;
						//multiple indices in y
						default: {for(int i=kk+1;i<ll;i++){vci2_2_pointers.push_back(&(adjoint_color_indices[i]));};} break;
					}
					vci3_pointers.push_back(&(fundamental_color_indices[0]));
					for(int i=0;i<kk;i++){vci3_pointers.push_back(&(adjoint_color_indices[i]));}
					for(int i=kk+1;i<ll;i++){vci3_pointers.push_back(&(adjoint_color_indices[i]));}
					for(int i=ll+1;i<adjoint_degree;i++){vci3_pointers.push_back(&(adjoint_color_indices[i]));}
					vci3_pointers.push_back(&(fundamental_color_indices[1]));

					color_string cs2_1(vci2_1_pointers), cs2_2(vci2_2_pointers), cs3(vci3_pointers);
					color_tensor ct2_1(cc2_1,cs2_1), ct2_2(cc2_2,cs2_2), ct3(minus_over_Nc,cs3);
					return( new color_tensor(ct2_1*ct2_2+ct3));
				}	
			}
		}
	}
}
else if(fundamental_degree==0){
  	//color = color /. TrC[a_]-> 0;
	if(adjoint_degree==1){
		std::vector<color_index*> vci2_pointers;
		color_constant zero(0,0);
		color_string cs2(vci2_pointers);
		color_tensor ct2(zero,cs2);
		return( new color_tensor(ct2));
	}
  	//color = color /. {TrC[a_, a_] -> (Nc^2 - 1)};
	if(adjoint_degree==2){
		if(adjoint_color_labels[0]==adjoint_color_labels[1]){
		std::vector<color_index*> vci2_pointers;
		color_constant Nc(1,1);
	       	color_constant mone(-1,0);
//		color_constant cc=Nc*Nc+one;
		color_constant cc=((Nc*Nc)+mone);
		color_string cs2(vci2_pointers);
		color_tensor ct2(cc,cs2);
		return( new color_tensor(ct2));
		}
	}
  	//color = color //. {TrC[x___, a_, y___, a_,z___] -> (TrC[x, z]*TrC[y] - 1/Nc*TrC[x, y, z])};
	if(adjoint_degree>2){
		for(int kk=0;kk<adjoint_degree-1;kk++){
			for(int ll=kk+1;ll<adjoint_degree;ll++){
				if(adjoint_color_labels[kk]==adjoint_color_labels[ll]){
					std::vector<color_index*> vci2_1_pointers, vci2_2_pointers, vci3_pointers;
					color_constant cc2_1(1,0);
					color_constant cc2_2(1,0);
					color_constant minus_over_Nc(-1,-1);
					for(int i=0;i<kk;i++){vci2_1_pointers.push_back(&(adjoint_color_indices[i]));}
					for(int i=ll+1;i<adjoint_degree;i++){vci2_1_pointers.push_back(&(adjoint_color_indices[i]));}
					switch(ll-kk){
						//y has no indices
						case 1: {color_constant Nc(1,1);cc2_2=Nc; break;}
						//y is only one index
						case 2: {color_constant zero(0,0);cc2_2=zero;break;}
						//multiple indices in y
						default: {
							//y has all indeces
							if(kk==0&&ll==adjoint_degree-1){color_constant Nc(1,1);cc2_2=Nc;}	 
							for(int i=kk+1;i<ll;i++){vci2_2_pointers.push_back(&(adjoint_color_indices[i]));};break;}
					}
					for(int i=0;i<kk;i++){vci3_pointers.push_back(&(adjoint_color_indices[i]));}
					for(int i=kk+1;i<ll;i++){vci3_pointers.push_back(&(adjoint_color_indices[i]));}
					for(int i=ll+1;i<adjoint_degree;i++){vci3_pointers.push_back(&(adjoint_color_indices[i]));}

					color_string cs2_1(vci2_1_pointers), cs2_2(vci2_2_pointers), cs3(vci3_pointers);
					color_tensor ct2_1(cc2_1,cs2_1), ct2_2(cc2_2,cs2_2), ct3(minus_over_Nc,cs3);
					return( new color_tensor(ct2_1*ct2_2+ct3));
				}	
			}
		}
	}
}
	//default: no simplification possible, so return zero-pointer
	return 0;
}

color_tensor* color_string::simplify(const color_string& cs1){
// notice: The Id color_sting is treated as one and _not_ Tr(one)=Nc

// reduce produxct of two identity color strings to one identity color string
	if(fundamental_degree==0&&
		adjoint_degree==0&&
		cs1.fundamental_degree==0&&
		cs1.adjoint_degree==0){
	std::vector<color_index*> vci2_pointers;
	color_constant one(1,0);

	color_string cs2(vci2_pointers);
	color_tensor ct2(one,cs2);

	return( new color_tensor(ct2));
	}

//color = color /. TraceF[x_,{i_,j_}]*TraceF[y_,{j_,k_}] :> TraceF[Join[x,y],{i,k}];
	if(fundamental_degree>0&&cs1.fundamental_degree>0){
	
		if(cs1.fundamental_color_labels[0]==fundamental_color_labels[1]){
			std::vector<color_index> vci2;
			
			vci2.push_back(fundamental_color_indices[0]);
			for(int kk=0;kk<adjoint_color_indices.size();kk++){
				vci2.push_back(adjoint_color_indices[kk]);
			}	
			for(int kk=0;kk<cs1.adjoint_color_indices.size();kk++){
				vci2.push_back(cs1.adjoint_color_indices[kk]);
			}	
			vci2.push_back(cs1.fundamental_color_indices[1]);	
		std::vector<color_index*> vci2_pointers;
		for(int kk=0;kk<vci2.size();kk++){vci2_pointers.push_back(&vci2[kk]);}
		
		color_constant cc(1,0);
		color_string cs2(vci2_pointers);
		return( new color_tensor(cc,cs2));
		}
		if(cs1.fundamental_color_labels[1]==fundamental_color_labels[0]){
			std::vector<color_index> vci2;
			vci2.push_back(cs1.fundamental_color_indices[0]);
			for(int kk=0;kk<cs1.adjoint_color_indices.size();kk++){
				vci2.push_back(cs1.adjoint_color_indices[kk]);
			}	
			for(int kk=0;kk<adjoint_color_indices.size();kk++){
				vci2.push_back(adjoint_color_indices[kk]);
			}	
			vci2.push_back(fundamental_color_indices[1]);	

		std::vector<color_index*> vci2_pointers;
		for(int kk=0;kk<vci2.size();kk++){vci2_pointers.push_back(&(vci2[kk]));}
		
		color_constant cc(1,0);
		color_string cs2(vci2_pointers);
		return( new color_tensor(cc,cs2));
		}}
//color = color //. {TraceF[{w___, a_, x___},{i1_,i2_}]*TraceF[{y___, a_, z___},{i3_,i4_}] -> TraceF[{w,z},{i1,i4}]TraceF[{y,x},{i3,i2}]-1/Nc*TraceF[{w,x},{i1,i2}]*TraceF[{y, z},{i3,i4}]};
// here we use notation: cc2 * vci2_1 * vci2_2 + cc3* vci3_1 * vci3_2 
		size_t kk_max(cs1.adjoint_color_indices.size()),ll_max(adjoint_color_indices.size());
	for(int kk=0;kk<kk_max;kk++){
		for(int ll=0;ll<ll_max;ll++){
			if(cs1.adjoint_color_labels[kk]==adjoint_color_labels[ll]){
		std::vector<color_index> vci2_1,vci2_2,vci3_1,vci3_2;
		color_constant one(1,0),cc3(-1,-1);
	if(fundamental_degree>0&&cs1.fundamental_degree>0){
		vci2_1.push_back(fundamental_color_indices[0]);	
		vci3_1.push_back(fundamental_color_indices[0]);	
		vci2_2.push_back(cs1.fundamental_color_indices[0]);	
		vci3_2.push_back(cs1.fundamental_color_indices[0]);	
			for(int i=0;i<ll;i++)	vci2_1.push_back(adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci2_1.push_back(cs1.adjoint_color_indices[i]);
			
			for(int i=0;i<kk;i++)	vci2_2.push_back(cs1.adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci2_2.push_back(adjoint_color_indices[i]);
			
			for(int i=0;i<ll;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			
			for(int i=0;i<kk;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
		vci2_1.push_back(cs1.fundamental_color_indices[1]);	
		vci3_1.push_back(fundamental_color_indices[1]);	
		vci2_2.push_back(fundamental_color_indices[1]);	
		vci3_2.push_back(cs1.fundamental_color_indices[1]);	
	}
 //color = color //. {TraceF[{w___, a_, x___},{i1_,i2_}]*TrC[y___, a_, z___] -> (TraceF[{w, z, y,x},{i1,i2}] - 1/Nc*TraceF[{w,x},{i1,i2}]*TrC[z, y])};
	if(fundamental_degree>0&&cs1.fundamental_degree==0&&cs1.adjoint_degree>1){
		vci2_1.push_back(fundamental_color_indices[0]);
		vci3_1.push_back(fundamental_color_indices[0]);	
//		vci2_2.push_back(cs1.fundamental_color_indices[0]);	
//		vci3_2.push_back(cs1.fundamental_color_indices[0]);	
			for(int i=0;i<ll;i++)	vci2_1.push_back(adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci2_1.push_back(cs1.adjoint_color_indices[i]);
			for(int i=0;i<kk;i++)	vci2_1.push_back(cs1.adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci2_1.push_back(adjoint_color_indices[i]);
			
			for(int i=0;i<ll;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			
			for(int i=0;i<kk;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
		vci2_1.push_back(fundamental_color_indices[1]);	
		vci3_1.push_back(fundamental_color_indices[1]);	
//		vci2_2.push_back(fundamental_color_indices[1]);	
//		vci3_2.push_back(fundamental_color_indices[1]);	
	}
//color = color //. {TrC[y___, a_, z___] TraceF[{w___, a_, x___},{i1_,i2_}]-> (TraceF[{w, z, y,x},{i1,i2}] - 1/Nc*TrC[z, y]*TraceF[{w,x},{i1,i2}])};
	if(fundamental_degree==0&&adjoint_degree>1&&cs1.fundamental_degree>0){
//		vci2_1.push_back(fundamental_color_indices[0]);
//		vci3_1.push_back(fundamental_color_indices[0]);	
		vci2_2.push_back(cs1.fundamental_color_indices[0]);	
		vci3_2.push_back(cs1.fundamental_color_indices[0]);	
			for(int i=0;i<kk;i++)	vci2_2.push_back(cs1.adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci2_2.push_back(adjoint_color_indices[i]);
			for(int i=0;i<ll;i++)	vci2_2.push_back(adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci2_2.push_back(cs1.adjoint_color_indices[i]);
			
			for(int i=0;i<ll;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			
			for(int i=0;i<kk;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
//		vci2_1.push_back(fundamental_color_indices[1]);	
//		vci3_1.push_back(fundamental_color_indices[1]);	
		vci2_2.push_back(cs1.fundamental_color_indices[1]);	
		vci3_2.push_back(cs1.fundamental_color_indices[1]);	
	}
// color = color //. TrC[w___, a_, x___]*TrC[y___, a_, z___] -> (TrC[x, w, z, y] - 1/Nc*TrC[x, w]*TrC[z, y])
	if(fundamental_degree==0&&adjoint_degree>1&&cs1.fundamental_degree==0&&cs1.adjoint_degree>1){
			for(int i=0;i<ll;i++)	vci2_1.push_back(adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci2_1.push_back(cs1.adjoint_color_indices[i]);
			for(int i=0;i<kk;i++)	vci2_1.push_back(cs1.adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci2_1.push_back(adjoint_color_indices[i]);
			
			for(int i=0;i<ll;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			for(int i=ll+1;i<ll_max;i++)	vci3_1.push_back(adjoint_color_indices[i]);
			
			for(int i=0;i<kk;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
			for(int i=kk+1;i<kk_max;i++)	vci3_2.push_back(cs1.adjoint_color_indices[i]);
	}

		std::vector<color_index*> vci2_1_pointers, vci2_2_pointers, vci3_1_pointers, vci3_2_pointers;
		for(int i=0;i<vci2_1.size();i++){vci2_1_pointers.push_back(&(vci2_1[i]));}
		for(int i=0;i<vci2_2.size();i++){vci2_2_pointers.push_back(&(vci2_2[i]));}
		for(int i=0;i<vci3_1.size();i++){vci3_1_pointers.push_back(&(vci3_1[i]));}
		for(int i=0;i<vci3_2.size();i++){vci3_2_pointers.push_back(&(vci3_2[i]));}

		color_string cs2_1(vci2_1_pointers),cs2_2(vci2_2_pointers),cs3_1(vci3_1_pointers),cs3_2(vci3_2_pointers);
		color_tensor ct2_1(one,cs2_1), ct2_2(one,cs2_2), ct3_1(cc3,cs3_1), ct3_2(one,cs3_2);

		return( new color_tensor(ct2_1*ct2_2+ct3_1*ct3_2));
			};};
		}
				
	return 0;
	}


color_constant::color_constant(int num,int power){
	if(power<0){
		for(int kk=1;kk<std::abs(power);kk++){
			negative_coeffs.push_back(multi_precision_constant(0));
			negative_coeffs_fracs.push_back(multi_precision_fraction(0));
		}
		negative_coeffs.push_back(multi_precision_constant(num));
		negative_coeffs_fracs.push_back(multi_precision_fraction(num));
	}
	else{
		for(int kk=0;kk<power;kk++){
			positive_coeffs.push_back(multi_precision_constant(0));
			positive_coeffs_fracs.push_back(multi_precision_fraction(0));
		}
		positive_coeffs.push_back(multi_precision_constant(num));
		positive_coeffs_fracs.push_back(multi_precision_fraction(num));
	}
}


color_constant::color_constant(int num,int den,int power){
	if(power<0){
		for(int kk=1;kk<std::abs(power);kk++){
			negative_coeffs.push_back(multi_precision_constant(0));
			negative_coeffs_fracs.push_back(multi_precision_fraction(0));
		}
		negative_coeffs.push_back(multi_precision_constant(REPEAT_ARG(R(num)/R(den))));
		negative_coeffs_fracs.push_back(multi_precision_fraction(num,den));
	}
	else{
		for(int kk=0;kk<power;kk++){
			positive_coeffs.push_back(multi_precision_constant(0));
			positive_coeffs_fracs.push_back(multi_precision_fraction(0));
		}
		positive_coeffs.push_back(multi_precision_constant(REPEAT_ARG(R(num)/R(den))));
		positive_coeffs_fracs.push_back(multi_precision_fraction(num,den));
	}
}

bool color_constant::is_zero(){
	bool result;
	if(positive_coeffs.size()==1&&
	   negative_coeffs.size()==0){
		result=(R(positive_coeffs[0])==0.&&positive_coeffs_fracs[0].m_num==0);
		return result;}
	
	size_t zero_coeffs(0);
	for(int kk=0;kk<positive_coeffs_fracs.size();kk++){
		if(R(positive_coeffs[kk])==0.&&
			positive_coeffs_fracs[kk].m_num==0){zero_coeffs++;}
	}
	for(int kk=0;kk<negative_coeffs_fracs.size();kk++){
		if(R(negative_coeffs[kk])==0.&&
			negative_coeffs_fracs[kk].m_num==0){zero_coeffs++;}
	}
	result=(zero_coeffs==positive_coeffs_fracs.size()+negative_coeffs_fracs.size());
	return result;
}


void color_constant::project_to_Nc_powers(int power_max, int power_min){

    if(power_max<power_min){
        _MESSAGE("bad color constant projection"); throw;
    }

	int neg_coeff,pos_coeff;
	neg_coeff=negative_coeffs.size();pos_coeff=positive_coeffs.size();

    if(power_max>-1){
        if(power_min>-1){
                negative_coeffs.clear();
                negative_coeffs_fracs.clear();
                //int i_max=min(power_min,pos_coeff-1);
                int i_max=min(power_min,pos_coeff);
                for(size_t i=0; i<i_max;i++){
#if BH_USE_GMP
         	        positive_coeffs[i]=multi_precision_constant(0,0,0,0);
#else
         	        positive_coeffs[i]=multi_precision_constant(0,0,0);
#endif
         	        positive_coeffs_fracs[i]=multi_precision_fraction(0);
                }
                for(size_t i=power_max+1;i<pos_coeff;i++){
#if BH_USE_GMP
         	        positive_coeffs[i]=multi_precision_constant(0,0,0,0);
#else
         	        positive_coeffs[i]=multi_precision_constant(0,0,0);
#endif
    	            positive_coeffs_fracs[i]=multi_precision_fraction(0);
                }
        }
        else{
            for(size_t i=-power_min;i<neg_coeff;i++){
#if BH_USE_GMP
         	    negative_coeffs[i]=multi_precision_constant(0,0,0,0);
#else
         	    negative_coeffs[i]=multi_precision_constant(0,0,0);
#endif
    	    	negative_coeffs_fracs[i]=multi_precision_fraction(0);
            } 
            for(size_t i=power_max+1;i<pos_coeff;i++){
#if BH_USE_GMP
                positive_coeffs[i]=multi_precision_constant(0,0,0,0);
#else
                positive_coeffs[i]=multi_precision_constant(0,0,0);
#endif
       	        positive_coeffs_fracs[i]=multi_precision_fraction(0);
            }
        }
    }
    else{
        positive_coeffs.clear();
        positive_coeffs_fracs.clear();
#if BH_USE_GMP
        positive_coeffs.push_back(multi_precision_constant(0,0,0,0));
#else
        positive_coeffs.push_back(multi_precision_constant(0,0,0));
#endif
    	positive_coeffs_fracs.push_back(multi_precision_fraction(0));
        
        // power_min<0
        int i_max=min(-power_max-1,neg_coeff);
        for(size_t i=0;i<i_max;i++){
#if BH_USE_GMP
            negative_coeffs[i]=multi_precision_constant(0,0,0,0);
#else
            negative_coeffs[i]=multi_precision_constant(0,0,0);
#endif
    	    negative_coeffs_fracs[i]=multi_precision_fraction(0);
        }
        for(size_t i=-power_min;i<neg_coeff;i++){
#if BH_USE_GMP
            negative_coeffs[i]=multi_precision_constant(0,0,0,0);
#else
            negative_coeffs[i]=multi_precision_constant(0,0,0);
#endif
    	    negative_coeffs_fracs[i]=multi_precision_fraction(0);
        }
    }


    vector<multi_precision_constant >::reverse_iterator it=negative_coeffs.rbegin();
    vector<multi_precision_fraction >::reverse_iterator it_f=negative_coeffs_fracs.rbegin();
    while( it!=negative_coeffs.rend() && R(*it)==0. && R(*it_f)==0. ){
        negative_coeffs.pop_back(); 
        negative_coeffs_fracs.pop_back(); 
        it++;it_f++;
    }

    it=positive_coeffs.rbegin();
    it_f=positive_coeffs_fracs.rbegin();
    while( it!=positive_coeffs.rend() &&  R(*it)==0. && R(*it_f)==0.){
        positive_coeffs.pop_back(); 
        positive_coeffs_fracs.pop_back(); 
        it++;it_f++;
    }

    if(positive_coeffs.size()==0 && positive_coeffs_fracs.size()==0){
#if BH_USE_GMP
        positive_coeffs.push_back(multi_precision_constant(0,0,0,0));
#else
        positive_coeffs.push_back(multi_precision_constant(0,0,0));
#endif
    	positive_coeffs_fracs.push_back(multi_precision_fraction(0));
    }

    return;
}





double color_constant::eval(){
	double result(0);
	double Nc=BH::settings::BH_interface_settings::s_nc;
	int negp=this->negative_coeffs.size();
	int posp=this->positive_coeffs.size();
	for(int i=negp;i>0;i--){
			result+=R(this->negative_coeffs_fracs[i-1])/std::pow(Nc,i);
	}
	for(int i=0;i<posp;i++){
			result+=R(this->positive_coeffs_fracs[i])*std::pow(Nc,i);
	}
	return result;
};


std::ostream& operator<<(std::ostream& s, const color_constant& cc){
	int negp=cc.negative_coeffs.size();
	int posp=cc.positive_coeffs.size();
	s<<"(";
	// negative powers
	for(int kk=negp;kk>1;kk--){
		if(cc.negative_coeffs_fracs[kk-1].m_num!=0){
			s<< R(cc.negative_coeffs_fracs[kk-1])<<"/Nc^"<<kk;
			for(int jj=(kk-1);jj>0;jj--){
				if(cc.negative_coeffs_fracs[jj-1].m_num!=0){
					s<<" + ";
					break;
				}
			}
		}
	}
	if(negp>0){
		if(cc.negative_coeffs_fracs[0].m_num!=0){
			s<<R(cc.negative_coeffs_fracs[0])<<"/Nc^1";
		}
	}
	if(posp>0 && negp>0){
		for(int jj=0;jj<posp;jj++){
			if(cc.positive_coeffs_fracs[jj].m_num!=0){
				s<<" + ";
				break;
			}
		}
	}
	// positive powers
	for(int kk=0;kk<posp;kk++){
		if(cc.positive_coeffs_fracs[kk].m_num!=0){
			if(kk==0){s<<R(cc.positive_coeffs_fracs[kk])<<"*Nc^0";}
					else {s<<R(cc.positive_coeffs_fracs[kk])<<"*Nc^"<<kk;}
			for(int jj=(kk+1);jj<posp;jj++){
				if(cc.positive_coeffs_fracs[jj].m_num!=0){
					s<<" + ";
					break;
				}
			}
		}
	}
	if(negp==0&&posp==1&&cc.positive_coeffs_fracs[0].m_num==0){
		s<<"0*Nc^0";};
	s<<")";
	if(negp==0&&posp==0)
		s<<"empty color_constant";
return s;
}

const color_constant operator+(const color_constant& cc1,const color_constant& cc2){
	color_constant cc3;
	int np1,np2,np3;
	int pp1,pp2,pp3;
	np1=cc1.negative_coeffs_fracs.size();
	np2=cc2.negative_coeffs_fracs.size();
	pp1=cc1.positive_coeffs_fracs.size();
	pp2=cc2.positive_coeffs_fracs.size();
	if(pp1>pp2)
		//bug Fernando :pp3=np1;
		pp3=pp1;
	else
		pp3=pp2;
	if(np1>np2)
		np3=np1;
	else
		np3=np2;
	// positive powers for sum
	for(int kk=0;kk<pp3;kk++){
		if(kk<pp1 && kk<pp2){
			cc3.positive_coeffs.push_back(multi_precision_constant(
						REPEAT_ARG(cc1.positive_coeffs[kk]+cc2.positive_coeffs[kk])
						));
			cc3.positive_coeffs_fracs.push_back(cc1.positive_coeffs_fracs[kk]+cc2.positive_coeffs_fracs[kk]);
		}
		else if(kk<pp1){
			cc3.positive_coeffs.push_back(cc1.positive_coeffs[kk]);
			cc3.positive_coeffs_fracs.push_back(cc1.positive_coeffs_fracs[kk]);
		}
		else{
			cc3.positive_coeffs.push_back(cc2.positive_coeffs[kk]);
			cc3.positive_coeffs_fracs.push_back(cc2.positive_coeffs_fracs[kk]);
		}	
	}
	// negative powers for sum
	for(int kk=0;kk<np3;kk++){
		if(kk<np1 && kk<np2){
			cc3.negative_coeffs.push_back(multi_precision_constant(
						REPEAT_ARG(cc1.negative_coeffs[kk]+cc2.negative_coeffs[kk])
						));
			cc3.negative_coeffs_fracs.push_back(cc1.negative_coeffs_fracs[kk]+cc2.negative_coeffs_fracs[kk]);
		}
		else if(kk<np1){
			cc3.negative_coeffs.push_back(cc1.negative_coeffs[kk]);
			cc3.negative_coeffs_fracs.push_back(cc1.negative_coeffs_fracs[kk]);
		}
		else{
			cc3.negative_coeffs.push_back(cc2.negative_coeffs[kk]);
			cc3.negative_coeffs_fracs.push_back(cc2.negative_coeffs_fracs[kk]);
		}	
	}
	return cc3;
}

const color_constant operator-(const color_constant& cc1,const color_constant& cc2){
	color_constant cc3;
	int np1,np2,np3;
	int pp1,pp2,pp3;
	np1=cc1.negative_coeffs_fracs.size();
	np2=cc2.negative_coeffs_fracs.size();
	pp1=cc1.positive_coeffs_fracs.size();
	pp2=cc2.positive_coeffs_fracs.size();
	if(pp1>pp2)
		pp3=pp1;
	else
		pp3=pp2;
	if(np1>np2)
		np3=np1;
	else
		np3=np2;
	// positive powers for sum
	for(int kk=0;kk<pp3;kk++){
		if(kk<pp1 && kk<pp2){
			cc3.positive_coeffs.push_back(multi_precision_constant(
						REPEAT_ARG(cc1.positive_coeffs[kk]-cc2.positive_coeffs[kk])
						));
			cc3.positive_coeffs_fracs.push_back(cc1.positive_coeffs_fracs[kk]-cc2.positive_coeffs_fracs[kk]);
		}
		else if(kk<pp1){
			cc3.positive_coeffs.push_back(cc1.positive_coeffs[kk]);
			cc3.positive_coeffs_fracs.push_back(cc1.positive_coeffs_fracs[kk]);
		}
		else{
			cc3.positive_coeffs.push_back(multi_precision_constant(REPEAT_ARG(-cc2.positive_coeffs[kk])));
			cc3.positive_coeffs_fracs.push_back(-cc2.positive_coeffs_fracs[kk]);
		}	
	}
	// negative powers for sum
	for(int kk=0;kk<np3;kk++){
		if(kk<np1 && kk<np2){
			cc3.negative_coeffs.push_back(multi_precision_constant(
						REPEAT_ARG(cc1.negative_coeffs[kk]-cc2.negative_coeffs[kk])
						));
			cc3.negative_coeffs_fracs.push_back(cc1.negative_coeffs_fracs[kk]-cc2.negative_coeffs_fracs[kk]);
		}
		else if(kk<np1){
			cc3.negative_coeffs.push_back(cc1.negative_coeffs[kk]);
			cc3.negative_coeffs_fracs.push_back(cc1.negative_coeffs_fracs[kk]);
		}
		else{
			cc3.negative_coeffs.push_back(multi_precision_constant(REPEAT_ARG(-cc2.negative_coeffs[kk])));
			cc3.negative_coeffs_fracs.push_back(-cc2.negative_coeffs_fracs[kk]);
		}	
	}
	return cc3;
}


color_constant& color_constant::operator=(const color_constant& cc1){
	positive_coeffs=cc1.positive_coeffs;
	positive_coeffs_fracs=cc1.positive_coeffs_fracs;
	negative_coeffs=cc1.negative_coeffs;
	negative_coeffs_fracs=cc1.negative_coeffs_fracs;
	return *this;
}


color_constant& color_constant::operator+=(const color_constant& cc1){
	*this=*this+cc1;
	return *this;
}

const color_constant operator*(const color_constant& cc1,const color_constant& cc2){
	color_constant cc3;
	vector<color_constant> vcc;
	int np1,np2;
	int pp1,pp2;
	np1=cc1.negative_coeffs_fracs.size();
	np2=cc2.negative_coeffs_fracs.size();
	pp1=cc1.positive_coeffs_fracs.size();
	pp2=cc2.positive_coeffs_fracs.size();
	for(int kk=-np1;kk<pp1;kk++){
		color_constant clocal;
		if(kk<0){
			for(int jj=-np2;jj<pp2;jj++){
				if(jj<0){
					color_constant ccl2(cc1.negative_coeffs_fracs[-kk-1].m_num*cc2.negative_coeffs_fracs[-jj-1].m_num,
							cc1.negative_coeffs_fracs[-kk-1].m_den*cc2.negative_coeffs_fracs[-jj-1].m_den,
						jj+kk);
//					std::cout<<"h1 "<<kk<<" "<<jj<<endl;
//					std::cout<<cc1.negative_coeffs_fracs[-kk-1].m_num<<" "<<cc2.negative_coeffs_fracs[-jj-1].m_num<<endl;
//					ccl2.print_constant();std::cout<<endl;
					clocal+=ccl2;
				}
				else{
					color_constant ccl2(cc1.negative_coeffs_fracs[-kk-1].m_num*cc2.positive_coeffs_fracs[jj].m_num,
							cc1.negative_coeffs_fracs[-kk-1].m_den*cc2.positive_coeffs_fracs[jj].m_den,
						jj+kk);
//					std::cout<<"h2 "<<kk<<" "<<jj<<endl;
//					std::cout<<cc1.negative_coeffs_fracs[-kk-1].m_num<<" "<<cc2.positive_coeffs_fracs[jj].m_num<<endl;
//					ccl2.print_constant();std::cout<<endl;
					clocal+=ccl2;
				}
			}
		}
		else{
			for(int jj=-np2;jj<pp2;jj++){
				if(jj<0){
					color_constant ccl2(cc1.positive_coeffs_fracs[kk].m_num*cc2.negative_coeffs_fracs[-jj-1].m_num,
							cc1.positive_coeffs_fracs[kk].m_den*cc2.negative_coeffs_fracs[-jj-1].m_den,
						jj+kk);
//					std::cout<<"h3 "<<kk<<" "<<jj<<endl;
//					std::cout<<cc1.positive_coeffs_fracs[kk].m_num<<" "<<cc2.negative_coeffs_fracs[-jj-1].m_num<<endl;
//					ccl2.print_constant();std::cout<<endl;
					clocal+=ccl2;
				}
				else{
					color_constant ccl2(cc1.positive_coeffs_fracs[kk].m_num*cc2.positive_coeffs_fracs[jj].m_num,
							cc1.positive_coeffs_fracs[kk].m_den*cc2.positive_coeffs_fracs[jj].m_den,
						jj+kk);
//					std::cout<<"h4 "<<kk<<" "<<jj<<endl;
//					std::cout<<cc1.positive_coeffs_fracs[kk].m_num<<" "<<cc2.positive_coeffs_fracs[jj].m_num<<endl;
//					ccl2.print_constant();std::cout<<endl;
					clocal+=ccl2;
				}
			}
		}
		cc3+=clocal;
	}
	return cc3;
}



const color_constant operator*(int i,const color_constant& cc2){
	color_constant cc1(i,0);
	color_constant cc3=cc1*cc2;
	return cc3;
	}

bool operator==(const color_constant& cc1, const color_constant& cc2){
	color_constant cc3=cc1-cc2;
	return cc3.is_zero();
};


single_color_tensor& single_color_tensor::operator=(const single_color_tensor& sct){
	color_coeff=sct.color_coeff;
	color_string_vector=sct.color_string_vector;
	return *this;
};


single_color_tensor::single_color_tensor(const color_constant& coeff, 
		color_index ci1,
		color_index ci2,
		color_index ci3,
		color_index ci4,
		color_index ci5,
		color_index ci6,
		color_index ci7,
		color_index ci8,
		color_index ci9,
		color_index ci10,
		color_index ci11): color_coeff(coeff){
if(color_coeff.is_zero()){
	color_string cs;
	color_constant cc(0,0);
	color_coeff=cc;
	color_string_vector.push_back(new color_string(cs));
}
else{	
	color_coeff=coeff;
	vector<color_index> vci;
		vci.push_back(ci1);
		vci.push_back(ci2);
		vci.push_back(ci3);
		vci.push_back(ci4);
		vci.push_back(ci5);
		vci.push_back(ci6);
		vci.push_back(ci7);
		vci.push_back(ci8);
		vci.push_back(ci9);
		vci.push_back(ci10);
		vci.push_back(ci11);
	vector<color_index> vcip;
	for(int kk=0;kk<11;kk++){
		if(vci[kk].label!=-1){
			vcip.push_back(vci[kk]);
		}
		else
			break;
	}
	for(int kk=0;kk<vcip.size();kk++){
		vector<color_index*> cs;
		if(vcip[kk].kind==color_index::fundamental){
			cs.push_back(&vcip[kk]);kk++;
			while(vcip[kk].kind==color_index::adjoint){
				cs.push_back(&vcip[kk]);kk++;
			}
			cs.push_back(&vcip[kk]);}
		else {	
			cs.push_back(&vcip[kk]);kk++;
			while(vcip[kk].kind==color_index::adjoint){
				cs.push_back(&vcip[kk]);kk++;
			}
		}
	color_string_vector.push_back(new color_string(cs));
	}
}
}

single_color_tensor::single_color_tensor(const color_constant& coeff, std::vector<color_string* > csv): color_coeff(coeff){
	if(color_coeff.is_zero()){
		color_string cs;
		color_constant cc(0,0);
		color_coeff=cc;
		color_string_vector.push_back(new color_string(cs));
	}
	else{
		color_string_vector=csv;
	}
};

single_color_tensor::single_color_tensor(const color_constant& coeff,const color_string& csv): color_coeff(coeff){
	if(color_coeff.is_zero()){
		color_string cs;
		color_constant cc(0,0);
		color_coeff=cc;
		color_string_vector.push_back(new color_string(cs));
	}
	else{
		color_string_vector.push_back(new color_string(csv));
	};
}
	
	



single_color_tensor::single_color_tensor(const single_color_tensor& sct): color_coeff(sct.color_coeff), color_string_vector(sct.color_string_vector){};

const single_color_tensor operator*(const single_color_tensor& cc1,const single_color_tensor& cc2){
	color_constant clocal;
	clocal=(cc1.color_coeff*cc2.color_coeff);
	std::vector<color_string* > color_string_vector;
	size_t cc1_lngth,cc2_lngth;
		cc1_lngth=cc1.color_string_vector.size();
		cc2_lngth=cc2.color_string_vector.size();
	for(int kk=0;kk<cc1_lngth;kk++){
		color_string_vector.push_back(cc1.color_string_vector[kk]);
	}
	for(int kk=0;kk<cc2_lngth;kk++){
		color_string_vector.push_back(cc2.color_string_vector[kk]);
	}
	single_color_tensor sct(clocal,color_string_vector);
	return sct;	
}

const single_color_tensor operator*(const color_constant& cc,const single_color_tensor& cc1){
	single_color_tensor sct(cc*cc1.color_coeff,cc1.color_string_vector);
	return sct;	
}

color_tensor* single_color_tensor::simplify(){
	size_t kk_max(color_string_vector.size());

	//drop Id() expressions in single_color_tensors
	if(kk_max>1&&color_string_vector[0]->adjoint_degree==0&&color_string_vector[0]->fundamental_degree==0){
		std::vector<color_string* >  v_cs;
		for(int i=1;i<kk_max;i++){v_cs.push_back(color_string_vector[i]);}
		return new color_tensor(color_coeff,v_cs);
	}

	//general simplifucations and contractions
	for(int kk=0;kk<kk_max;kk++){
		color_tensor* ct_pointer=color_string_vector[kk]->simplify();
		if(ct_pointer!=0){
			std::vector<color_string* >  v_cs;
			for(int i=0;i<kk;i++){v_cs.push_back(color_string_vector[i]);}	
			for(int i=kk+1;i<kk_max;i++){v_cs.push_back(color_string_vector[i]);}	
			color_tensor ct(color_coeff, v_cs);
			return( new color_tensor(ct*(*ct_pointer)));
		}
	}
	for(int kk=0;kk<kk_max-1;kk++){
		for(int ll=kk+1;ll<kk_max;ll++){
			color_tensor* ct_pointer=color_string_vector[kk]->simplify(*(color_string_vector[ll]));
			if(ct_pointer!=0){
			std::vector<color_string* >  v_cs;
			for(int i=0;i<kk;i++){v_cs.push_back(color_string_vector[i]);}	
			for(int i=kk+1;i<ll;i++){v_cs.push_back(color_string_vector[i]);}	
			for(int i=ll+1;i<kk_max;i++){v_cs.push_back(color_string_vector[i]);}	
			color_tensor ct(color_coeff, v_cs);
			return( new color_tensor(ct*(*ct_pointer)));
			}
		}
	}
	return 0;
}


void single_color_tensor::conjugate(){
	int sign(1);
	for(int kk=0;kk<color_string_vector.size();kk++){
	sign=sign*(color_string_vector[kk]->conjugate());
	}
	color_coeff=sign*color_coeff;
	}

bool single_color_tensor::is_zero(){
	return color_coeff.is_zero();
}


//order: 
//(1) order of present form of single_color_tensors
//(2) fewer color_string
//(3) same number of color_strings: ask for order of color_strings
//(4) identical color_strings: return false
bool operator<(const single_color_tensor&sct1, const single_color_tensor&sct2){
	if(sct1.color_string_vector.size()<sct2.color_string_vector.size()){return true;}
	else if(sct1.color_string_vector.size()>sct2.color_string_vector.size()){return false;}
	else if(sct1.color_string_vector.size()==sct2.color_string_vector.size()){
		if(sct1.color_string_vector.size()==0){return false;}
		else{
			for(int k=0;k<sct1.color_string_vector.size();k++){
				if(*(sct1.color_string_vector[k])<*(sct2.color_string_vector[k])){return true;}
				if(*(sct2.color_string_vector[k])<*(sct1.color_string_vector[k])){return false;}
			}
			return false;
		}
		
		return false;}
}
bool compare_sct(single_color_tensor* sct1, single_color_tensor* sct2){
	return(*sct1<*sct2);
}

//equal: 
//(1) prefactor does not matter, only color strings
//(2) color_string_vector is not ordered into cannonical form;
//	single_color_tensors identical after reordering, are not recognized as identical

bool operator==(const single_color_tensor&sct1, const single_color_tensor&sct2){
	if(sct1.color_string_vector.size()==sct2.color_string_vector.size()){
		if(sct1.color_string_vector.size()==0){return true;}
		else {
			size_t k(0) ;
			while(*(sct1.color_string_vector[k])==*(sct2.color_string_vector[k])){
//				_PRINT(*(sct1.color_string_vector[k]));
//				_PRINT(*(sct2.color_string_vector[k]));
				k++;
				if(k==sct1.color_string_vector.size()){return true;};
			};
		return false;
		}
	}
	else {return false;}
}

//checks if single_color_tensor is sorted
bool single_color_tensor::sortedQ(){
	size_t size(color_string_vector.size());
	if(size<2){return true;}
	else{	
		for(int k=0;k<size-1;k++){
			if(!(*(color_string_vector[k])<*(color_string_vector[k+1]))){return false;}
		}
	}
	// new here
	return true;
}

//primitive sorting aglorithm using ordering of color_strings
//retruns true if expression changed, false otherwise
bool single_color_tensor::sort(){
	size_t size(color_string_vector.size());
	if(sortedQ()){return false;}
	else if(size>1){
		std::sort(color_string_vector.begin(),color_string_vector.end(),compare_cs);
		return true;
	}
	return false;
}

color_constant* single_color_tensor::project_to_color_constant(){
	return (&color_coeff); 
};

color_constant* single_color_tensor::project_to_single_color_tensor(const single_color_tensor& sct ){
	if((*this)==sct){
		return &(this->color_coeff);
	}
	else {return 0;}
};


std::ostream& operator<<(std::ostream& s, const single_color_tensor& sct){
	s<<sct.color_coeff<<"*";
	for(int kk=0;kk<sct.color_string_vector.size();kk++){
	s<<*(sct.color_string_vector[kk]);
	}
	return s;
}


color_tensor& color_tensor::operator=(const color_tensor& ct){
	single_color_tensors=ct.single_color_tensors;
	return *this;
};



color_tensor::color_tensor(const single_color_tensor& sct){
//	single_color_tensor local_sct(sct);
//	single_color_tensors.push_back(&local_sct);
	single_color_tensors.push_back(new single_color_tensor(sct));
	}

color_tensor::color_tensor(std::vector<single_color_tensor* > v_sct){
	for(int kk=0;kk<v_sct.size();kk++){
	single_color_tensors.push_back(new single_color_tensor(*v_sct[kk]));
	};
	}


color_tensor::color_tensor(const color_constant & coeff, std::vector<color_string* > csv){
	single_color_tensors.push_back(new single_color_tensor(coeff, csv));
	};

color_tensor::color_tensor(const color_constant & coeff, const color_string& csv){
	single_color_tensors.push_back(new single_color_tensor(coeff, csv));
	};


color_tensor::color_tensor(const color_constant& coeff, 
				color_index ci1,
				color_index ci2,
				color_index ci3,
				color_index ci4,
				color_index ci5,
				color_index ci6,
				color_index ci7,
				color_index ci8,
				color_index ci9,
				color_index ci10,
				color_index ci11){
	single_color_tensors.push_back(new single_color_tensor(coeff, ci1, ci2, ci3, ci4, ci5, ci6, ci7, ci8, ci9, ci10, ci11));
	}

const color_tensor operator+(const color_tensor& ct1,const color_tensor& ct2){
	color_tensor ct(ct1);
	for(int kk=0;kk<ct2.single_color_tensors.size();kk++){
		ct.single_color_tensors.push_back(ct2.single_color_tensors[kk]);
	}
	return ct;	
	}

const color_tensor operator-(const color_tensor& ct1,const color_tensor& ct2){
	color_tensor ct(ct1);
	color_constant minus(-1,0);
	for(int kk=0;kk<ct2.single_color_tensors.size();kk++){
		single_color_tensor sct=minus*(*(ct2.single_color_tensors[kk]));
		ct.single_color_tensors.push_back(&sct);
	}
	return ct;	
	}

const color_tensor operator*(const color_constant& cc,const color_tensor& ct1){
	std::vector<single_color_tensor*> v_sct;
	for(int kk=0;kk<ct1.single_color_tensors.size();kk++){
		v_sct.push_back(new single_color_tensor(cc*(*ct1.single_color_tensors[kk])));
	}
	color_tensor ct(v_sct);
	return ct;	
	}


const color_tensor operator*(const color_tensor& ct1,const color_tensor& ct2){
	std::vector<single_color_tensor*> v_sct;
	for(int kk=0;kk<ct1.single_color_tensors.size();kk++){
		for(int ll=0;ll<ct2.single_color_tensors.size();ll++){
		single_color_tensor sct((*ct1.single_color_tensors[kk])*(*ct2.single_color_tensors[ll]));
		v_sct.push_back(new single_color_tensor((*ct1.single_color_tensors[kk])*(*ct2.single_color_tensors[ll])));
		}
	}
	color_tensor ct(v_sct);
	return ct;	
	}


bool color_tensor::single_simplify(){
	int kk_max(single_color_tensors.size());
	//_PRINT(kk_max);
	std::vector<single_color_tensor* >  v_sct;
	for(int kk=0;kk<kk_max;kk++){
		if(!(single_color_tensors[kk]->is_zero())){
			v_sct.push_back(single_color_tensors[kk]);
		}
	}
	if(v_sct.size()==0){v_sct.push_back(single_color_tensors[0]);};

	single_color_tensors=v_sct;
	kk_max=single_color_tensors.size();
	for(int kk=0;kk<kk_max;kk++){
		color_tensor* ct_pointer=single_color_tensors[kk]->simplify();
		if(ct_pointer!=0){
			v_sct.clear();
			for(int i=0;i<kk;i++){v_sct.push_back(single_color_tensors[i]);}	
			for(int i=0;i<ct_pointer->single_color_tensors.size();i++){v_sct.push_back(ct_pointer->single_color_tensors[i]);}	
			for(int i=kk+1;i<kk_max;i++){v_sct.push_back(single_color_tensors[i]);}	
			
			
			single_color_tensors=v_sct;
			return true;
		}
	}
	return false;
};

void color_tensor::sort(){
	int size(single_color_tensors.size());
	for(int k=0;k<size;k++){
		single_color_tensors[k]->sort();
	}
#if 0
	//section to check weak ordering
	cout<<"-------------------------------------------------"<<endl;
	for(int i=0;i<single_color_tensors.size();i++){
		for(int j=0;j<single_color_tensors.size();j++){
		cout<<i<<":"<<j<<": "<<(*(single_color_tensors[i])<*(single_color_tensors[j]))<<"="<<(*(single_color_tensors[j])<*(single_color_tensors[i]))<<"=?"<<(*(single_color_tensors[j])==*(single_color_tensors[i]))<<endl;
		_PRINT(*(single_color_tensors[i]));_PRINT(*(single_color_tensors[j]));
		}
	}
	cout<<"-------------------------------------------------"<<endl;
#endif
	if(size>1){
		std::stable_sort(single_color_tensors.begin(),single_color_tensors.end(),compare_sct);
	}
	return;
}

//returns true if two single_color_tensors could be summed, false otherwise.
bool color_tensor::shorten(){
	int size(single_color_tensors.size());
	if(size<2){return false;};
	for(int k=0;k<size-1;k++){
		if(*(single_color_tensors[k])==*(single_color_tensors[k+1])){
			std::vector<single_color_tensor* > scts;
			color_constant cc=(single_color_tensors[k]->color_coeff)+(single_color_tensors[k+1]->color_coeff);
//			for(int l=0;l<size-1;l++){
			for(int l=0;l<size;l++){
				if(l!=k&&l!=k+1){scts.push_back(single_color_tensors[l]);}
				else if(l==k){
					scts.push_back(new single_color_tensor(cc,single_color_tensors[k]->color_string_vector));
				}
			}
			single_color_tensors=scts;
			return true;
		}
	}
	return false;
}

bool color_tensor::is_zero(){
	for(int i=0;i<single_color_tensors.size();i++){
		if(!single_color_tensors[i]->is_zero()){
			return false; 
		}
	}
	return true;
}


void color_tensor::conjugate(){
	size_t size(single_color_tensors.size());
	if(size>0){
		for(int k=0;k<size;k++){
			single_color_tensors[k]->conjugate();
		}
	}
}



void color_tensor::simplify(){
	//int kk(0);
	while(this->single_simplify()){
	//	kk++;
	//	std::cout<<kk<<std::endl;
	//	_PRINT(*this);
	};
	//_MESSAGE("end single_simplify");
	//_PRINT(*this);
	//_MESSAGE("start sort");
	this->sort();
	//_PRINT(*this);
	while(this->shorten()){
	//	_MESSAGE("shorten");
	//	kk++;
	//	std::cout<<kk<<std::endl;
	//	_PRINT(*this);
	};
	//_MESSAGE("simplified");
	//_PRINT(*this);
	return;
};


color_constant* color_tensor::project_to_color_constant(){
	return single_color_tensors[0]->project_to_color_constant(); 
};



color_constant* color_tensor::project_to_color_tensor(const color_tensor& ct){
	std::vector<color_constant*> cc_vec;
	color_constant cc(0,0);
	if(ct.single_color_tensors.size()>1){
		_MESSAGE("Trying to project on non-primitive color-tensor. Can only project on color_tensors that contain only one single_color_tensor."); throw;
	}
	single_color_tensor* sct=ct.single_color_tensors[0];
	for(int i=0;i<single_color_tensors.size();i++){
		cc_vec.push_back(single_color_tensors[i]->project_to_single_color_tensor(*sct));
	}
	for(int i=0;i<cc_vec.size();i++){
		if(cc_vec[i]!=0){cc=cc+(*(cc_vec[i]));};
	}
	return new color_constant(cc);
};


std::ostream& operator<<(std::ostream& s, const color_tensor& ct){
	if(ct.single_color_tensors.size()>0){
		for(int kk=0;kk<ct.single_color_tensors.size()-1;kk++){
		      	s<<*(ct.single_color_tensors[kk])<<" + "<<endl;
		}
		s<<*(ct.single_color_tensors[ct.single_color_tensors.size()-1]);
	}
	else {s<<"empty tensor";};
	return s;
	
}

}

#endif /* COLOR_ALGEBRA_H_ */


