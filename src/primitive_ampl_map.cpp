/*
 * primitive_ampl_map.cpp
 *
 *  Created on: 15 Jan 2011
 *      Author: 
 */
#include <algorithm>
#include "BH_utilities.h"
#include "primitive_ampl_map.h"
#include "helcode.h"

using namespace std;

namespace BH {

//modified in the sense that we do not disinguish between gluinos and quarks
long compute_pcode(const vector<plabel>& p){
  	int nbr_q=0;
  	int nbr_l=0;
  	int nbr_g=0;
  	int nbr_higgs=0;
  	int nbr_g_massive=0;
  	int nbr_q_massive=0;
  	int nbr_photon=0;
  	int nbr_gluino=0;
  	int nbr_gluino_massive=0;
	int nbr_sc_massless=0;
	int nbr_sc_massive=0;

  	for (size_t i=0;i<p.size();i++) {
		switch (p[i].type()->pdg_code()){
			case 1: nbr_q++; break; // quark
			case 11: nbr_l++; break; // lepton
			case 21: nbr_g++; break; // gluon
			case 25: nbr_higgs++; break; // higgs
			case -1: nbr_g_massive++; break; // massive gluon
			case -2: nbr_q_massive++; break; // massive quark
			case 8: nbr_photon++; break; // photon
			case 1000: nbr_gluino++; break; // gluinos
			case -3: nbr_gluino_massive++; break; // massive gluinos
			case -4: nbr_sc_massive++; break; // massive scalars
			case -5: nbr_sc_massless++; break; // massless scalars
  	  	}
  	}

  	return
  		nbr_g
  		+10*(nbr_q+nbr_gluino)
  		+100*nbr_l
  		+1000*nbr_g_massive
  		+10000*(nbr_q_massive+nbr_gluino_massive)
		+100000*nbr_photon
		+100000000*nbr_higgs
		+1000000000*nbr_sc_massless
		+10000000000*nbr_sc_massive;
}



// flip process and indices to take adventage of flip symmetry and of LT analytics
// needed to call correct analytics
// needed to speed up computation

void flip_cs_at(size_t pos, string& cs){
    if(cs==""||cs=="glue"||cs=="nf") return;
    else if(cs[0]=='n'){
        //nf-cases
        if(cs[pos+2]=='L') cs[pos+2]='R';
        else               cs[pos+2]='L';
    }
    else{
        if(cs[pos]=='L') cs[pos]='R';
        else             cs[pos]='L';
    }
    return;}

void flip_cs(string& cs){
    if(cs=="glue"||cs=="nf") return;
    
    size_t offset(0), n(cs.size());
    if(cs[0]=='n') offset=2; //nf...-cases
            
    for(size_t i=offset;i<n;i++){
        if(cs[i]=='L') cs[i]='R';
        else           cs[i]='L';
    }
    return;
}    

void flip_qb_to_q(vector<plabel>& plabels, string & cs){

    /*
    cout<<cs<<" ";
    for(size_t i=0;i<plabels.size();i++){
        cout<<plabels[i]<<" ";
    }
    cout<<endl;
    */
    flip_cs_at(0,cs);
    size_t offset(0);
    while(!plabels[offset].is_a(quark)) offset++;
    plabels[offset].set_is_antiparticle(false); 
    offset++;
    while(!plabels[offset].is_a(quark)) offset++;
    plabels[offset].set_is_antiparticle(true);


    /*
    cout<<cs<<" ";
    for(size_t i=0;i<plabels.size();i++){
        cout<<plabels[i]<<" ";
    }
    cout<<endl;
    */

    return;
}

bool is_Ltype_cs(string cs){
    if(cs=="glue"||cs=="nf") return true;
    if(cs[0]=='n') return cs[2]=='L';
    return cs[0]=='L';
}



/* 
 * arrange (anti-)quarks and (anti-)gluinos with flavor bigger than 
 * one to appeare in specified order
 * assumse ascending gap-less flavor labels: flavors {1,2,3,... }
 */
void arrange_quarks_and_cs(vector<plabel>& plabels, string& cs){
    //external colored fermions
    size_t extq((compute_pcode(plabels)%100)/20);

    vector<bool> q_is_anti(extq,true); //vector to specify which cs/flavors have to be flipped
    vector<size_t> q_set_flavor(extq,1); //vector to specify which cs/flavors have to be flipped
    size_t flavor(0),flavor_count(1);
    for(size_t i=0;i<plabels.size();i++){
        flavor=plabels[i].flavor();
        if(flavor>1){
            if(q_is_anti[flavor-1]){
                flavor_count++;
                q_is_anti[flavor-1]=false;
                q_set_flavor[flavor-1]=flavor_count;
                plabels[i].set_flavor(flavor_count);
                if(!plabels[i].is_anti_particle()) {
                    plabels[i].set_is_antiparticle(true);   
                    if(cs!="") flip_cs_at(flavor-1,cs);
                }
            }
            else{
                plabels[i].set_flavor(q_set_flavor[flavor-1]);
                if(plabels[i].is_anti_particle()){
                    plabels[i].set_is_antiparticle(false);
                }
            }
        }
    }
    if(cs!=""){
        string cs_new(cs);
        size_t n(cs.size()-q_set_flavor.size());
        for(size_t i=1; i<q_set_flavor.size();i++){
           cs[n+i]=cs_new[n+q_set_flavor[i]-1]; 
        }
    }
    return;
}


void rot_qm_pro_ind(vector<plabel>& plabels, double & sign, string & cs, short length_X=0, bool qbm_first=true){
    //rotates pro and ind: {... qm ...} to {qm .... }
    //rotates pro and ind: {... qbm ...} to {qbm ... }
    //for case of (qm .... qbp) return ( qbm ... qp ) with R/L...T -> L/R ... T
    //length_X number of non-qcd particles to be ignored
    if(qbm_first && 
            plabels[0].is_a(quark)&&
            plabels[0].helicity()==-1&&
            plabels[0].is_anti_particle()) return;
    else if(!qbm_first &&
            plabels[0].is_a(quark)&&
            plabels[0].helicity()==1&&
            plabels[0].is_anti_particle()) return;

    size_t offset(0), qm_offset(0), qp_offset(0);
    while(!plabels[offset].is_a(quark)) offset++;
    if(plabels[offset].helicity()==-1){qm_offset=offset;}
    else {qp_offset=offset;};
    offset++;
    while(!plabels[offset].is_a(quark)) offset++;
    if(plabels[offset].helicity()==-1){qm_offset=offset;}
    else {qp_offset=offset;};

    size_t n=plabels.size()-length_X;
    
    // for leptons and photons flipping photon emission to 
    // other side of quark line gives a sign
    if(qbm_first){
        rotate(plabels.begin(),plabels.begin()+qm_offset,plabels.end()-length_X);
        if(length_X>0 && qm_offset>qp_offset) sign*=-1;

        if(!plabels[0].is_anti_particle()){
            qp_offset=(n+qp_offset-qm_offset)%n;
            plabels[0].set_is_antiparticle(true);
            plabels[qp_offset].set_is_antiparticle(false);
            flip_cs_at(0,cs); //flip first L/R-turner
        }
    }
    else{
        rotate(plabels.begin(),plabels.begin()+qp_offset,plabels.end()-length_X);
        if(length_X>0 && qp_offset>qp_offset) sign*=-1;
 
        if(!plabels[0].is_anti_particle()){
            qm_offset=(n+qm_offset-qp_offset)%n;
            plabels[0].set_is_antiparticle(true);
            plabels[qp_offset].set_is_antiparticle(false);
            flip_cs_at(0,cs); //flip first L/R-turner
        }
    }

    return;
}

// expect quark in first position
void which_quark_hel_first_tree_X(vector<plabel>& plabels, double & sign, short first_quark_hel, short length_X=0){

    if(first_quark_hel!=plabels[0].helicity()){
    
        size_t pos_q2(1);
        string cs="";

        while(!plabels[pos_q2].is_a(quark)) pos_q2++;
            plabels[0].set_is_antiparticle(true);
            plabels[pos_q2].set_is_antiparticle(false);
            
            rotate(plabels.begin(),plabels.begin()+pos_q2,plabels.end()-length_X);
            sign*=-1;  // sign from lepton flip to other side of quark line
           
            reverse(plabels.begin()+1,plabels.end()-length_X);
            if((plabels.size()-length_X)%2==1) sign*=-1;  

            arrange_quarks_and_cs(plabels,cs);

    } 
}

//assume lepton positions in the end of plabels
void sort_leptons(vector<plabel>& plabels, double & sign){
   
    /*
     * sort leptons to have lepton with lowest index first
     * convert anti-leptons to leptons to have lepton come first
     *
     */

    size_t n(plabels.size());
    if(plabels[n-2].ind()>plabels[n-1].ind()){
        reverse(plabels.end()-2,plabels.end());
        sign*=-1;
    }
    if(plabels[n-2].is_anti_particle()){
        plabels[n-1].set_is_antiparticle(true);
        plabels[n-2].set_is_antiparticle(false);
    }
}



void rot_qm_pro_ind_tree_X(vector<plabel>& plabels, double & sign, short length_X=0){

    size_t pos_q1(0), pos_q2(0);
    size_t x,y,n(plabels.size()-length_X);
    while(!plabels[pos_q1].is_a(quark)) pos_q1++;
    pos_q2=pos_q1+1;
    while(!plabels[pos_q2].is_a(quark)) pos_q2++;
    x=pos_q2-pos_q1-1;
    y=n-2-x;
    if(pos_q1!=0){
        rotate(plabels.begin(),plabels.begin()+pos_q1,plabels.end()-length_X);
        pos_q2=pos_q2-pos_q1;
    }
    //cases q(i) x  q(j) y X: 
    //  - i<j
    //  - i>j
    //  - x>y
    //  - x<y
            
    if(x<y){
        if(plabels[0].ind()<plabels[pos_q2].ind()){
            reverse(plabels.begin()+1,plabels.end()-length_X);
            if(length_X==2){
                if(n%2==1) sign*=(-1); //signs from reverse and 
                //lepton flip to other side of quark line
                //and lepton exchange
            }
            else if(length_X==1){
                if(n%2==1) sign*=(-1); 
            }
            else if(length_X==0){
                if(n%2==1) sign*=(-1); 
            }

        }
        else if(plabels[0].ind()>plabels[pos_q2].ind()){
            plabels[0].set_is_antiparticle(true);
            plabels[pos_q2].set_is_antiparticle(false);
            rotate(plabels.begin(),plabels.begin()+pos_q2,plabels.end()-length_X);
            if(length_X>0) sign*=-1;  // sign from lepton flip to other side of quark line
        }
    }
    else if(x>=y){
        if(plabels[0].ind()>plabels[pos_q2].ind()){
            plabels[0].set_is_antiparticle(true);
            plabels[pos_q2].set_is_antiparticle(false);
            reverse(plabels.begin(),plabels.end()-length_X);
            if(length_X==2){
                if(n%2==1) sign*=(-1);  //signs from revese and 
                sign*=-1;               //reverse of leptons
            }
            else if(length_X==1){
                if(n%2==0) sign*=(-1); //sign from reverse
            }
            else if(length_X==0){
                if(n%2==1) sign*=(-1); //sign from reverse
            }

            if(y>0){
                rotate(plabels.begin(),plabels.begin()+y,plabels.end()-length_X);
            }
        }
    }

}




void rot_qm_pro_ind_X(vector<plabel>& plabels, double & sign, string & cs, short length_X=0, bool qm_first=true){
   
    /*
     * for qm_first=true (qm_first=false puts qp into first position) 
     * rotates pro and ind: {... qm ...} to {qm .... }
     * rotates pro and ind: {... qbm ...} to {qbm ... }
     * for case of (qbm .... qbp) return ( qm ... qp ) with L...T -> R ... T
     *
     */
    
    size_t n(plabels.size());

    if(qm_first
            &&plabels[0].is_a(quark)
            && plabels[0].helicity()==-1
            && !plabels[0].is_anti_particle())  return;
    else if(!qm_first
            &&plabels[0].is_a(quark)
            && plabels[0].helicity()==1
            && !plabels[0].is_anti_particle())  return;

    //find positions of quarks
    size_t offset(0), qm_offset(0), qp_offset(0);
    while(!plabels[offset].is_a(quark)) offset++;
    if(plabels[offset].helicity()==-1){qm_offset=offset;}
    else {qp_offset=offset;};
    offset++;
    while(!plabels[offset].is_a(quark)) offset++;
    if(plabels[offset].helicity()==-1){qm_offset=offset;}
    else {qp_offset=offset;};

    //reverse process if necessary
    //leave X part in same place
    if((qm_first && (qm_offset>qp_offset))
            ||(!qm_first && (qp_offset>qm_offset))){
        reverse(plabels.begin(),plabels.end()-length_X);
        if(length_X>0) reverse(plabels.end()-length_X,plabels.end());
        flip_cs(cs);    
        qp_offset=(n-1)-length_X-qp_offset;
        qm_offset=(n-1)-length_X-qm_offset;

        if(n%2==1) sign*=(-1);
    }
    
    //if first quark is anti-quark turn into quark
    //flip color first structure label R/L -> L/R
    if((qp_offset<qm_offset)
            && plabels[qp_offset].is_anti_particle()){
        plabels[qp_offset].set_is_antiparticle(false);
        plabels[qm_offset].set_is_antiparticle(true);
        flip_cs_at(0,cs);
    }
    else if((qm_offset<qp_offset)
            && plabels[qm_offset].is_anti_particle()){
        plabels[qm_offset].set_is_antiparticle(false);
        plabels[qp_offset].set_is_antiparticle(true);
        flip_cs_at(0,cs);
    }

    //rotate qm/qp into first position
    //leave X part in same place
    if(qm_first) offset=qm_offset;
    else offset=qp_offset;
    rotate(plabels.begin(),plabels.begin()+offset,plabels.end()-length_X);
    return;
} 



void flip_pro_ind(vector<plabel>& plabels, double & sign, string & cs, short length_X=0){
    //flips pro and ind: {12345...n} to {1 n ... 3 2} with (-1)^n
    //flips color structures
  
    if((plabels.size())%2==1) sign*=(-1);
 
    // reorder leptons
    // flip leptons/photons to other side of quark line
    // no net sign for leptons
    if(length_X==1) sign*=-1; 
    
    reverse(plabels.begin()+1,plabels.end()-length_X);
    flip_cs(cs);
}


//adapted from regular helcode to take plabels as arguments
int helcode_g(const vector<plabel>& p){
    int res=0;int base=1;
    for (size_t i=0;i<p.size();i++){
     	switch (p[i].helicity()) {
	    case 1: res+=base;  break;
	}
	base*=2;
    }
    return res;
}


/* 1) we rotate particla 1 into first position
 * 2) then, with the following function, we rotate into canonicial helicity configurations in order to:
 * 	a) reduce the memory usage
 * 	b) adapt to the available analytics
 */

void cannonical_gluon_hel_configs(vector<plabel>& plabels){

	size_t n(plabels.size());
	int hel_sum(0);
	for(size_t i=0;i<n;i++) hel_sum+=plabels[i].helicity();
	
	//if(hel_sum<0) hel_sum+=n;
	//else hel_sum-=n; //greater equal zero

	switch(n){
		case 6:{
		   	if(hel_sum>0){
				//_PRINT(hel_sum);
				for(size_t i=0;i<5;i++){
					//_PRINT(helcode_g(plabels));
					switch(helcode_g(plabels)){
						//target helcodes
						case 60: case 58: case 54: return; //MHV 
						default: {
            						rotate(plabels.begin(),plabels.begin()+1,plabels.end());
							 } break;		   
					}
				}
				//_PRINT(helcode_g(plabels));
				return;
			}
		     	else{
				//_PRINT(hel_sum);
				for(size_t i=0;i<5;i++){
					//_PRINT(helcode_g(plabels));
					switch(helcode_g(plabels)){
						//target helcodes
						case 3: case 5: case 9:			//MHVb
						case 7: case 11: case 21:  		//NMHV +++--- , ++-+-- , +-+-+-
						case 56: case 52: case 42:  return;  	//NMHV ---+++ , --+-++ , -+-+-+
						default: {
            						rotate(plabels.begin(),plabels.begin()+1,plabels.end());
							 } break;		   
					}
				}
				//_PRINT(helcode_g(plabels));
				return;
			}	
		       } break;
		default: return;
	}

	return;
}


void canonical_pro(vector<plabel>& plabels, double & sign, short & conjQ, string & cs){

    long code(compute_pcode(plabels));

    switch(code/10){
        case 0: { 
                /* map all primitive amplitudes to
                 *
                 * m/p(1) m/p(i) ... m/p(j) with i<j
                 *
                 * */

	    //map to minimal set of amplitudes for given momentum asingment
            size_t offset(0), n(plabels.size());
            while(plabels[offset].ind()!=1) offset++;
            if(offset!=0) rotate(plabels.begin(),plabels.begin()+offset,plabels.end());
            if(plabels[1].ind()>plabels[n-1].ind()){
                reverse(plabels.begin()+1,plabels.end());
                if(n%2==1) sign*=(-1);
            }
       
	    /*map to minimal set of helicity amplitudes to reduce memory usage
	     * and allow complete usage of analytics
	     */
	    cannonical_gluon_hel_configs(plabels);
	    


	    //need to check this
	    if(code>5) conjugateQ(plabels, sign, conjQ, cs);
          
        } break;
        case 2: case 4: case 6: case 8: {
                /* map all primitive amplitudes to
                 *
                 * 1) (qbm ... qp ... ) cs: (L/R,...)
                 * 2) (R,...) -> (L,...)
                 *
                 */

            //needed for analytics (in order to call right ampls)
            if(code%10+code/10<6){
                rot_qm_pro_ind(plabels,sign,cs);
                if (not (cs==""|| is_Ltype_cs(cs))) flip_pro_ind(plabels,sign,cs);
                arrange_quarks_and_cs(plabels,cs);
            }
            //needed for numerics (in order to reduce computation)
            if(code%10+code/10>5) conjugateQ(plabels, sign, conjQ, cs);
  

        } break;
        case 22: case 24: case 26: case 28: {
                /* map all primitive amplitudes to
                 *
                 *
                 *
                 */
            short n_leptons(2);
            sort_leptons(plabels,sign);
            if(plabels.back().helicity()==1){
                bool qm_first(false);
                rot_qm_pro_ind_X(plabels,sign,cs,n_leptons,qm_first);
            }
            else{
                bool qm_first(true);
                rot_qm_pro_ind_X(plabels,sign,cs,n_leptons,qm_first);
            }
            sort_leptons(plabels,sign);
            if (!is_Ltype_cs(cs)) flip_pro_ind(plabels,sign,cs,n_leptons);
            arrange_quarks_and_cs(plabels,cs);
            
            conjugateQ(plabels, sign, conjQ, cs);
                                                
        } break;
        case 10002: case 10004: case 10006: case 10008: {
                /* map all primitive amplitudes to
                 *
                 *     qp .... qbm ... yp
                 * or:
                 *     qm .... qbp ... ym
                 *
                 */
           
#if 1                                                           
            short n_photons(1);
            if(plabels.back().helicity()==1){
                bool qm_first(false);
                rot_qm_pro_ind_X(plabels,sign,cs,n_photons,qm_first);
            }
            else{
                bool qm_first(true);
                rot_qm_pro_ind_X(plabels,sign,cs,n_photons,qm_first);
            }
            if (!is_Ltype_cs(cs)) flip_pro_ind(plabels,sign,cs,n_photons);
            arrange_quarks_and_cs(plabels,cs);
            
            conjugateQ(plabels, sign, conjQ, cs);
           
#endif                                                            
        } break;
    default: return;
    }

};


void conjugateQ(vector<plabel>& plabels, double & sign, short & conjQ, string & cs){

    size_t code(compute_pcode(plabels));
    size_t n=plabels.size();
    short in_fermions(0);

    //need this sign although it cancels in interferece terms
    //a conjugate amplitude might multiply a non-conjugate one
    for(size_t i=0;i<plabels.size();i++){
        if(plabels[i].ind()<3 && plabels[i].is_fermion()) in_fermions++;
    }
   
    //sign for conjugation depending on sign of energy component of momentum
    switch(in_fermions){
        case 2: in_fermions=1;break;
        case 1: in_fermions=-1;break;
        default: in_fermions=1;break;
    }

    //_MESSAGE2("in_fermions",in_fermions);
    //sign depending on quark-line dependent power of I
    if(((code/100+code%100/10)/2)%2==0) in_fermions*=-1;


    switch(code/10){
        case 0: { 
                /* gluons only
		 * conjugate if first gluon has positive helicity and beyond 5-points
                 *
                 * */
           	if(code==6){
			int hel_sum(0); 
		    	for(size_t i=0;i<n;i++) hel_sum+=plabels[i].helicity();
			hel_sum=n-abs(hel_sum);	
			// helsum>4 condition switches off conj of MHV ampls, for which we have analystics
	            	if(plabels[0].helicity()==1 && hel_sum>4 ) conjQ=in_fermions;
	            	//if(plabels[0].helicity()==1) conjQ=in_fermions;
		}

        } break;
        case 2: case 4: case 6: case 8: {
                /* 
		 * quarks and gluons only
		 * map all primitive amplitudes to
                 *
                 */

             if(plabels[0].helicity()==-1) conjQ=in_fermions;


        } break;
        case 22: case 24: case 26: case 28: {
                /* expecting input:
                 *
                 *  (1) qp .... qbm ... lp(i) lbm(j) ; with i<j
                 *      or
                 *  (2) qm .... qbp ... lm(i) lbp(j) ; with i<j
                 * 
                 * for case (2) set flags for complex conjugation
                 */
            
             if(plabels[0].helicity()==-1) conjQ=in_fermions;

        } break;
        case 10002: case 10004: case 10006: case 10008: {
                /* expecting input:
                 *
                 *  (1) qp .... qbm ... yp
                 *      or
                 *  (2) qm .... qbp ... ym
                 * 
                 * for case (2) set flags for complex conjugation
                 */

            if(plabels.back().helicity()==-1) conjQ=in_fermions;
            //if(plabels.back().helicity()==1) conjQ=in_fermions;
     } break;
    }


};


void conjugateQ_tree(vector<plabel>& plabels, double & sign, short & conjQ){

    size_t code(compute_pcode(plabels));
    size_t n=plabels.size();
    short in_fermions(0);

    //need this sign although it cancels in interferece terms
    //a conjugate amplitude might multiply a non-conjugate one
    for(size_t i=0;i<plabels.size();i++){
        if(plabels[i].ind()<3 && plabels[i].is_fermion()) in_fermions++;
    }
   
    //sign for conjugation depending on sign of energy component of momentum
    switch(in_fermions){
        case 2: in_fermions=1;break;
        case 1: in_fermions=-1;break;
        default: in_fermions=1;break;
    }

    //_MESSAGE2("in_fermions",in_fermions);
    //sign depending on quark-line dependent power of I
    if(((code/100+code%100/10)/2)%2==0) in_fermions*=-1;


    switch(code/10){
        case 0: { 
                /* conjugate if first gluon has positive helicity and beyond 5-points
                 *
                 * */
            
            if(code>5 && plabels[0].helicity()==1) conjQ=in_fermions;
          
        } break;
        case 2: case 4: case 6: case 8: {
                /* map all primitive amplitudes to
                 *
                 *  
                 *
                 */

            if(plabels[0].helicity()==1) conjQ=in_fermions;

        } break;
        case 22: case 24: case 26: case 28: {
                /* map all primitive amplitudes to
                 *
                 * depending on momentum labels in the lepton pair
                 * set flags for complex conjugation
                 *
                 * assume input is canonical process qm...qbp lm lbp
                 */


            size_t n_leptons(2), hel(0);
            if(plabels[n-2].helicity()==1){
                conjQ=in_fermions;
                hel=1;
                which_quark_hel_first_tree_X(plabels,sign,hel,n_leptons);
            }
            else{
                hel=-1;
                which_quark_hel_first_tree_X(plabels,sign,hel,n_leptons);
            }
                
                                            
        } break;
        case 10002: case 10004: case 10006: case 10008: {
                /* map all primitive amplitudes to
                 *
                 *     qp .... qbm ... yp
                 * or:
                 *     qm .... qbp ... ym
                 *
                 */

            size_t n_photons(1), hel(0);
            if(plabels[n-1].helicity()==1){
                conjQ=in_fermions;
                hel=1;
                which_quark_hel_first_tree_X(plabels,sign,hel,n_photons);
            }
            else{
                hel=-1;
                which_quark_hel_first_tree_X(plabels,sign,hel,n_photons);
            }

        } break;
    }

};




void canonical_pro_tree(vector<plabel>& plabels, double & sign, short & conjQ){

    long code(compute_pcode(plabels));
    string cs=""; //dummy color_structure to be able reuse functions for virtual parts

    switch(code/10){
        case 0: { 
                /* map all primitive amplitudes to
                 *
                 * m/p(1) m/p(i) ... m/p(j) with i<j
                 *
                 * */

            size_t offset(0), n(plabels.size());
            while(plabels[offset].ind()!=1) offset++;
            if(offset!=0) rotate(plabels.begin(),plabels.begin()+offset,plabels.end());
            if(plabels[1].ind()>plabels[n-1].ind()){
                reverse(plabels.begin()+1,plabels.end());
                if(n%2==1) sign*=(-1);
            }
            //conjugateQ_tree(plabels, sign, conjQ);
          
        } break;
        case 2: case 4: case 6: case 8: {
                /* map all primitive amplitudes to
                 *
                 * 1) (qbm ... qp ... ) cs: (L/R,...)
                 * 2) (R,...) -> (L,...)
                 *
                 */
            rot_qm_pro_ind_tree_X(plabels,sign);
            arrange_quarks_and_cs(plabels,cs);
            //conjugateQ_tree(plabels, sign, conjQ);
  

        } break;
        case 22: case 24: case 26: case 28: {
                /* map all primitive amplitudes to
                 *
                 *
                 *
                 */
            short n_leptons(2);
            rot_qm_pro_ind_tree_X(plabels,sign,n_leptons);
            sort_leptons(plabels,sign); 
            arrange_quarks_and_cs(plabels,cs);
            conjugateQ_tree(plabels, sign, conjQ);
            
        } break;
        case 10002: case 10004: case 10006: case 10008: {
                /* map all primitive amplitudes to
                 *
                 *     qp .... qbm ... yp
                 * or:
                 *     qm .... qbp ... ym
                 *
                 */
            short n_photons(1);
            rot_qm_pro_ind_tree_X(plabels,sign,n_photons);
            arrange_quarks_and_cs(plabels,cs);
            conjugateQ_tree(plabels, sign, conjQ);
        } break;
    }

};



}
