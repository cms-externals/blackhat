#include "particles.h"
#include "process.h"
#include "BH_utilities.h"
#include <algorithm>
#include "BH_debug.h"

using namespace std;

namespace BH {

size_t nbr_of_flavors(const process& pro,const particle& type);


int helcode_phi_2q2Q(const process& pro){
  int factor = 0;
  int result = 0;
  int base  = 6;
  int power = 1;
  int nlegs = pro.n();
  vector<int> flavors;


  for(int i=1;i<=nlegs;i++){
	  if( (*pro.p(i).type()) == quark ) {
		  vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
		  if ( pos == flavors.end() ){
			  flavors.push_back(pro.p(i).flavor());
		  }
	  }

	  if( (*pro.p(i).type()) == quark ) {
		  if( pro.p(i).flavor() == flavors[0] ) {
			  if (pro.p(i).helicity() == -1 ){factor = 1;} else {factor = 2;}
		  }
		  else if( pro.p(i).flavor() == flavors[1] ) {
			  if (pro.p(i).helicity() == -1 ){factor = 4;} else {factor = 5;}
		  }
		  else {_WARNING("ERROR: wrong ptype to helcode_2q2Q (flavor not found)");}
	  }
	  else if( (*pro.p(i).type()) == gluon || (*pro.p(i).type()) == photon ) {
		  if      (pro.p(i).helicity() == -1 ) {factor = 0;} else {factor = 3;}
	  }

	  else if (pro.p(i).is_a(higgs)){factor = 14;}	  

	  else {
		  _WARNING("ERROR2: wrong ptype to helcode_2q2Q (type not found) ");
		  _PRINT(*pro.p(i).type());
	 }
	  result += power*factor;
	  power *= base;
  }

  return(result);
}




  int helcode_Ng1ph(const process& pro){ // n gluons and 1 phi 
  // std::cout << "ENTERING helcode_Ng1ph" <<  std::endl;
  int factor = 0;
  int result = 0;
  int base  = 4;
  int power = 1;
  for (int i=1;i<=pro.n();i++){
        if  (pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if (pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==false){factor = 1;}
    	else if (pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==true){factor = 2;}	
	else if (pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 3;}
	else {std::cout << "ERROR: wrong ptype to helcode_Ng1ph for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    result += power*factor;
  //   std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
    power *= base;
    }
  return(result);
}


int helcode_g(const process& p){
    int res=0;int base=1;
    for (size_t i=1;i<=p.n();i++){
     	switch (p.p(i).helicity()) {
	    case 1: res+=base;  break;
	}
	base*=2;
    }
    return res;
}

int helcode_2q(const process& pro){
  // std::cout << "ENTERING helcode_2q" <<  std::endl;
  int factor = 0;
  int result = 0;
  int base  = 4;
  int power = 1;
  for (int i=1;i<=pro.n();i++){
    if  (pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    else if (pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 ){factor = 1;}
    else if (pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 ){factor = 2;}
    else if (pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 3;}
    else {std::cout << "ERROR: wrong ptype to helcode_2q for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    result += power*factor;
  //   std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
    power *= base;
    }
  return(result);
}


int helcode_4q(const process& pro){
  //std::cout << "ENTERING helcode_4q" <<  std::endl;
  int factor = 0;
  int result = 0;
  int base  = 6;
  int power = 1;
  for (int i=1;i<=pro.n();i++){
    if      (pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    else if (pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    else if (pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()== false ){factor = 2;}
    else if (pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== true ){factor = 3;}
    else if (pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()== true ){factor = 4;}
    else if (pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 5;}
    else {std::cout << "ERROR: wrong ptype to helcode_4q" << std::endl;}
    result += power*factor;
 //    std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
    power *= base;
    }
  return(result);
}


int helcode_2q2Q(const process& pro){
  int factor = 0;
  int result = 0;
  int base  = 6;
  int power = 1;
  int nlegs = pro.n();
  vector<int> flavors;

  for(int i=1;i<=nlegs;i++){
	  if( (*pro.p(i).type()) == quark ) {
		  vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
		  if ( pos == flavors.end() ){
			  flavors.push_back(pro.p(i).flavor());
		  }
	  }

	  if( (*pro.p(i).type()) == quark ) {
		  if( pro.p(i).flavor() == flavors[0] ) {
			  if (pro.p(i).helicity() == -1 ){factor = 1;} else {factor = 2;}
		  }
		  else if( pro.p(i).flavor() == flavors[1] ) {
			  if (pro.p(i).helicity() == -1 ){factor = 4;} else {factor = 5;}
		  }
		  else {_WARNING("ERROR: wrong ptype to helcode_2q2Q (flavor not found)");}
	  }
	  else if( (*pro.p(i).type()) == gluon || (*pro.p(i).type()) == photon ) {
		  if      (pro.p(i).helicity() == -1 ) {factor = 0;} else {factor = 3;}
	  }
	  else {
		  _WARNING("ERROR: wrong ptype to helcode_2q2Q (type not found) ");
	 }
	  result += power*factor;
	  power *= base;
  }

  return(result);
}




// This switch only works for 6 quarks of 3 different flavor with NO gluons
int helcode_2q1_2q2_2q3(const process & pro){
   int factor = 0;
  int result = 0;
  int base  = 6;
  int power = 1;
  int nlegs = pro.n();
  vector<int> flavors;
  //std::cout << pro << std::endl;
  for(int i=1;i<=nlegs;i++){
  if( (*pro.p(i).type()) == quark ) {
   vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),
                                                        pro.p(i).flavor() ) ;
    if ( pos == flavors.end() ){
		  flavors.push_back(pro.p(i).flavor());
		  }
  }
  //  _PRINT(i);
  //  _PRINT(*pro.p(i).type()==quark);
  // _PRINT( pro.p(i).flavor() == flavors[0]);
  //_PRINT( pro.p(i).flavor() == flavors[1]);
  // _PRINT( pro.p(i).flavor() == flavors[2]);
  if( (*pro.p(i).type()) == quark ) {
     if    ( pro.p(i).flavor() == flavors[0] ) {
        if (pro.p(i).helicity() == -1 ){factor = 0;} else {factor = 1;}
     }
     else if( pro.p(i).flavor() == flavors[1] ) {
        if (pro.p(i).helicity() == -1 ){factor = 2;} else {factor = 3;}
     }
     else if( pro.p(i).flavor() == flavors[2] ) {
        if (pro.p(i).helicity() == -1 ){factor = 4;} else {factor = 5;}
     }
     else {_WARNING("ERROR: wrong ptype to helcode_2q1_2q2_2q3 (flavor not found)");}
  }
 else {_WARNING("ERROR: wrong ptype to helcode_2q1_2q2_2q3 (type not found) ");
 }
	  result += power*factor;
	  power *= base;
	  // _PRINT(result);
	  // _PRINT(power);
  }
  return(result);

}


// This switch only works for 6 quarks of 2 different flavor with NO gluons
int helcode_4q1_2q2(const process & pro){
  int factor = 0;
  int result = 0;
  int base  = 6;
  int power = 1;
  int nlegs = pro.n();
  vector<int> flavors;
  //std::cout << pro << std::endl;
  for(int i=1;i<=nlegs;i++){
  if( (*pro.p(i).type()) == quark ) {
   vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),
                                                        pro.p(i).flavor() ) ;
    if ( pos == flavors.end() ){
		  flavors.push_back(pro.p(i).flavor());
		  }
  }
  //  _PRINT(i);
  //  _PRINT(*pro.p(i).type()==quark);
  // _PRINT( pro.p(i).flavor() == flavors[0]);
  //_PRINT( pro.p(i).flavor() == flavors[1]);
  // _PRINT( pro.p(i).flavor() == flavors[2]);
  if( (*pro.p(i).type()) == quark ) {
     if    ( pro.p(i).flavor() == flavors[0] ) {
        if (pro.p(i).helicity() == -1 ){
	   if(pro.p(i).is_anti_particle()== false){factor = 0;}
           else {factor = 1;}
	}
        else if (pro.p(i).helicity() == 1 ){
           if(pro.p(i).is_anti_particle()== false){factor = 2;}
           else {factor = 3;}
	}
        else {_WARNING("ERROR: wrong helicity to helcode_4q1_2q2 ");
        }
     }
     else if( pro.p(i).flavor() == flavors[1] ) {
        if (pro.p(i).helicity() == -1 ){factor = 4;} else {factor = 5;}
     }
  }
  else {_WARNING("ERROR: wrong ptype to helcode_4q1_2q2 (type not found) "); }
	  result += power*factor;
	  power *= base;
	  // _PRINT(result);
	  // _PRINT(power);
  }
  return(result);

}



// This switch is a special switch for testing 6 quarks of different flavor with no gluons
int helcode_2q1_2q2_2q3Special(const process & pro){
  /*
switch value               process
   1:     {Q[m, 1], Q[p, 1], Q[m, 2], Q[p, 2], Q[m, 3], Q[p, 3]}
   2:     {Q[m, 1], Q[p, 1], Q[p, 2], Q[m, 2], Q[m, 3], Q[p, 3]}
   3:     {Q[m, 1], Q[p, 1], Q[m, 2], Q[p, 2], Q[p, 3], Q[m, 3]}
   4:     {Q[p, 1], Q[m, 1], Q[m, 2], Q[p, 2], Q[p, 3], Q[m, 3]}
   5:     {Q[p, 1], Q[m, 1], Q[m, 2], Q[p, 2], Q[m, 3], Q[p, 3]}
   6:     {Q[m, 1], Q[p, 1], Q[p, 2], Q[m, 2], Q[p, 3], Q[m, 3]}

   7:     {Q[m, 1], Q[m, 2], Q[p, 2], Q[p, 1], Q[p, 3], Q[m, 3]}
   8:     {Q[m, 1], Q[p, 2], Q[m, 2], Q[p, 1], Q[m, 3], Q[p, 3]}
   9:     {Q[m, 1], Q[m, 2], Q[p, 2], Q[p, 1], Q[m, 3], Q[p, 3]}
   10:    {Q[m, 1], Q[p, 2], Q[m, 2], Q[p, 1], Q[p, 3], Q[m, 3]}

*/
  int result = 0;
  int nlegs = pro.n();
  vector<int> flavors;
  for(int i=1;i<=nlegs;i++){
     vector<int>::iterator pos= find (
                       flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
		  if ( pos == flavors.end() ){
			  flavors.push_back(pro.p(i).flavor());
		  }
  }
  //     1:     {Q[m, 1], Q[p, 1], Q[m, 2], Q[p, 2], Q[m, 3], Q[p, 3]}
  if( pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[0] && pro.p(2).helicity() ==  1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() == -1 &&
      pro.p(4).flavor() == flavors[1] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() == -1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() ==  1 )
    {result = 1;}
  //   2:     {Q[m, 1], Q[p, 1], Q[p, 2], Q[m, 2], Q[m, 3], Q[p, 3]}
  else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[0] && pro.p(2).helicity() ==  1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() ==  1 &&
      pro.p(4).flavor() == flavors[1] && pro.p(4).helicity() == -1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() == -1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() ==  1 )
    {result = 2;}
  //  3:     {Q[m, 1], Q[p, 1], Q[m, 2], Q[p, 2], Q[p, 3], Q[m, 3]}
   else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[0] && pro.p(2).helicity() ==  1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() == -1 &&
      pro.p(4).flavor() == flavors[1] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() ==  1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() == -1 )
    {result = 3;}
  //   4:     {Q[p, 1], Q[m, 1], Q[m, 2], Q[p, 2], Q[p, 3], Q[m, 3]}
   else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() ==  1 &&
      pro.p(2).flavor() == flavors[0] && pro.p(2).helicity() == -1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() == -1 &&
      pro.p(4).flavor() == flavors[1] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() ==  1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() == -1 )
    {result = 4;}

  //   5:     {Q[p, 1], Q[m, 1], Q[m, 2], Q[p, 2], Q[m, 3], Q[p, 3]}
     else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() ==  1 &&
      pro.p(2).flavor() == flavors[0] && pro.p(2).helicity() == -1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() == -1 &&
      pro.p(4).flavor() == flavors[1] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() == -1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() ==  1 )
    {result = 5;}

  //   6:     {Q[m, 1], Q[p, 1], Q[p, 2], Q[m, 2], Q[p, 3], Q[m, 3]}
     else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[0] && pro.p(2).helicity() ==  1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() ==  1 &&
      pro.p(4).flavor() == flavors[1] && pro.p(4).helicity() == -1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() ==  1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() == -1 )
    {result = 6;}

  //   7:     {Q[m, 1], Q[m, 2], Q[p, 2], Q[p, 1], Q[p, 3], Q[m, 3]}
     else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[1] && pro.p(2).helicity() == -1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() ==  1 &&
      pro.p(4).flavor() == flavors[0] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() ==  1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() == -1 )
    {result = 7;}

  //   8:     {Q[m, 1], Q[p, 2], Q[m, 2], Q[p, 1], Q[m, 3], Q[p, 3]}
      else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[1] && pro.p(2).helicity() ==  1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() == -1 &&
      pro.p(4).flavor() == flavors[0] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() == -1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() ==  1 )
    {result = 8;}
  //   9:     {Q[m, 1], Q[m, 2], Q[p, 2], Q[p, 1], Q[m, 3], Q[p, 3]}
      else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[1] && pro.p(2).helicity() == -1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() ==  1 &&
      pro.p(4).flavor() == flavors[0] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() == -1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() ==  1 )
    {result = 8;}
  //   10:    {Q[m, 1], Q[p, 2], Q[m, 2], Q[p, 1], Q[p, 3], Q[m, 3]}

      else if(
      pro.p(1).flavor() == flavors[0] && pro.p(1).helicity() == -1 &&
      pro.p(2).flavor() == flavors[1] && pro.p(2).helicity() ==  1 &&
      pro.p(3).flavor() == flavors[1] && pro.p(3).helicity() == -1 &&
      pro.p(4).flavor() == flavors[0] && pro.p(4).helicity() ==  1 &&
      pro.p(5).flavor() == flavors[2] && pro.p(5).helicity() ==  1 &&
      pro.p(6).flavor() == flavors[2] && pro.p(6).helicity() == -1 )
    {result = 8;}
   else {
	  _WARNING("ERROR: wrong ptype to helcode_2q2Q (type not found) ");
	 }
  return(result);
}


/*
int helcode_2q_2q2_2q3(const process& pro){
  int factor = 0;
  int result = 0;
  int base  = 8;
  int power = 1;
  int nlegs = pro.n();
  vector<int> flavors;

  for(int i=1;i<=nlegs;i++){
	  if( (*pro.p(i).type()) == quark ) {
		  vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
		  if ( pos == flavors.end() ){
			  flavors.push_back(pro.p(i).flavor());
		  }
	  }
	  if( (*pro.p(i).type()) == quark ) {
		  if( pro.p(i).flavor() == flavors[0] ) {
			  if (pro.p(i).helicity() == -1 ){factor = 1;} else {factor = 2;}
		  }
		  else if( pro.p(i).flavor() == flavors[1] ) {
			  if (pro.p(i).helicity() == -1 ){factor = 4;} else {factor = 5;}
		  }
		  else if( pro.p(i).flavor() == flavors[2] ) {
			  if (pro.p(i).helicity() == -1 ){factor = 6;} else {factor = 7;}
		  else {_WARNING("ERROR: wrong ptype to helcode_2q2Q (flavor not found)");}
	  }
	  else if( (*pro.p(i).type()) == gluon ) {
		  if      (pro.p(i).helicity() == -1 ) {factor = 0;} else {factor = 3;}
	  }
	  else {
		  _WARNING("ERROR: wrong ptype to helcode_2q2Q (type not found) ");
	 }
	  result += power*factor;
	  power *= base;
  }

  return(result);
}

*/

int helcode_2q2l(const process& pro){
   int factor = 0;
  int result = 0;
  int base  = 6;
  int power = 1;
  int nlegs = pro.n();
    for (int i=1;i<=nlegs;i++)
      {
    	if  (pro.p(i).is_a(quark) ) {
   			if (pro.p(i).helicity() == -1 ){factor = 1;} else {factor = 2;}
    	}
    	else if (pro.p(i).is_a( gluon) ) {
  	  	  if      (pro.p(i).helicity() == -1 ) {factor = 0;} else {factor = 3;}
    	}
    	else if( (pro.p(i).is_a( lepton ) )) {
  	  	  if      (pro.p(i).helicity() == -1 ) {factor = 4;} else {factor = 5;}
    	}
       else {_WARNING("ERROR: wrong ptype to helcode_2q2e (type not found) ");}
     result += power*factor;
//     _PRINT(result);
     power *= base;
      }

//  std::cout << "helcode_2q2l: " << result << " for " << pro  << std::endl;
  return(result);
}

int helcode_2q1y(const process& pro){
   int factor = 0;
  int result = 0;
  int base  = 6;
  int power = 1;
  int nlegs = pro.n();
    for (int i=1;i<=nlegs;i++)
      {
    	if  (pro.p(i).is_a(quark) ) {
   			if (pro.p(i).helicity() == -1 ){factor = 1;} else {factor = 2;}
    	}
    	else if (pro.p(i).is_a( gluon) ) {
  	  	  if      (pro.p(i).helicity() == -1 ) {factor = 0;} else {factor = 3;}
    	}
    	else if( (pro.p(i).is_a( photon ) )) {
  	  	  if      (pro.p(i).helicity() == -1 ) {factor = 4;} else {factor = 5;}
    	}
       else {_WARNING("ERROR: wrong ptype to helcode_2q1y (type not found) ");}
     result += power*factor;
//     _PRINT(result);
     power *= base;
      }

//  std::cout << "helcode_2q2l: " << result << " for " << pro  << std::endl;
  return(result);
}


// assumes llb ordering
long helcode_2q2l2Q(const process& pro){
  int factor = 0;
  long result = 0;
  long base  = 8;
  long power = 1;
  int nlegs = pro.n();
  vector<int> flavors;
  int flavor_1=-100;
  int flavor_2=-100;
  int flavor_close,flavor_far;
  size_t pos_q,pos_q2,pos_qb,pos_qb2,pos_l=0,pos_lb=0;
  for(int i=1;i<=nlegs;i++){
	  if( (*pro.p(i).type()) == quark ) {
		  if (flavor_1==-100){
			  flavor_1=pro.p(i).flavor();
		  }
		  else {
			  if  (flavor_1 == pro.p(i).flavor()){
			  }
			  else {
				  if (flavor_2==-100){
					  flavor_2=pro.p(i).flavor();
				  } else {
					  if  (flavor_2 == pro.p(i).flavor()){
					  }
				  } ;

			  }
		  } ;
	  }
	  if( (*pro.p(i).type()) == lepton ) {
		  if (pro.p(i).is_anti_particle()){
				  pos_lb=i;
			  } else {
				   pos_l=i;
			  }
	  }
  }
//	  size_t pos_q_left=0,pos_q_right=0;
if (pos_lb==0 || pos_l==0) return 0;

	  int offset=0;
	while ( ! (pro.p((pos_lb+offset)%pro.n()+1).is_a(quark)  )){++offset;}
	size_t pos_q_right=(pos_lb+offset)%pro.n()+1;
	offset =-2;
	while ( ! (pro.p((pro.n()+pos_l+offset)%pro.n()+1).is_a(quark)  )){--offset;}
	size_t pos_q_left=(pro.n()+pos_l+offset)%pro.n()+1;

		  if ( (pro.p(pos_q_left).flavor() ==  pro.p(pos_q_right).flavor())  && (pro.p(pos_q_right).is_anti_particle()!= pro.p(pos_q_left).is_anti_particle())){
			  if (pro.p(pos_q_left).flavor() == flavor_1){
				  flavor_close=flavor_1;
				  flavor_far=flavor_2;
			  }
			  else {
			flavor_close=flavor_2;
			flavor_far=flavor_1;
		}

	} else {  // qq QQ ll configuration, the photon couples to the qq fermion line
		flavor_close=pro.p(pos_q_right).flavor();
		flavor_far=pro.p(pos_q_left).flavor();
	}

	  for(int i=1;i<=nlegs;i++){

		  if( (*pro.p(i).type()) == quark ) {
			  if( pro.p(i).flavor() == flavor_close ) {
				  if (pro.p(i).helicity() == -1 ){factor = 0;} else {factor = 1;}
			  }
			  else if( pro.p(i).flavor() == flavor_far ) {
				  if (pro.p(i).helicity() == -1 ){factor = 2;} else {factor = 3;}
			  }
		  }
		  else if( (*pro.p(i).type()) == gluon ) {
			  if (pro.p(i).helicity() == -1 ) {factor = 4;} else {factor = 5;}
		  }
		  else if( (*pro.p(i).type()) == lepton ) {
			  if (pro.p(i).helicity() == -1 ) {factor = 6;} else {factor = 7;}
		  }
	  result += power*factor;
	  power *= base;
  }
//_MESSAGE4("process: ",pro," result: ",result);

return(result);
}

long helcode_2q2l2G(const process& pro){
  int factor = 0;
  long result = 0;
  long base  = 8;
  long power = 1;
  int nlegs = pro.n();
  vector<int> flavors;
  int flavor_1=-100;
  int flavor_2=-100;
  int flavor_close,flavor_far;

  for(int i=1;i<=nlegs;i++){

	  if( (*pro.p(i).type()) == quark ) {
			  if (pro.p(i).helicity() == -1 ){factor = 0;} else {factor = 1;}
	  }
	  else if( (*pro.p(i).type()) == gluino ) {
			  if (pro.p(i).helicity() == -1 ){factor = 2;} else {factor = 3;}
	  }
	  else if( (*pro.p(i).type()) == gluon ) {
		  if (pro.p(i).helicity() == -1 ) {factor = 4;} else {factor = 5;}
	  }
	  else if( (*pro.p(i).type()) == lepton ) {
		  if (pro.p(i).helicity() == -1 ) {factor = 6;} else {factor = 7;}
	  }
  result += power*factor;
  power *= base;
  }
//_MESSAGE4("process: ",pro," result: ",result);

return(result);
}

long helcode_2q2G1y(const process& pro){
  int factor = 0;
  long result = 0;
  long base  = 8;
  long power = 1;
  int nlegs = pro.n();
  vector<int> flavors;
  int flavor_1=-100;
  int flavor_2=-100;
  int flavor_close,flavor_far;

  for(int i=1;i<=nlegs;i++){

	  if( (*pro.p(i).type()) == quark ) {
			  if (pro.p(i).helicity() == -1 ){factor = 0;} else {factor = 1;}
	  }
	  else if( (*pro.p(i).type()) == gluino ) {
			  if (pro.p(i).helicity() == -1 ){factor = 2;} else {factor = 3;}
	  }
	  else if( (*pro.p(i).type()) == gluon ) {
		  if (pro.p(i).helicity() == -1 ) {factor = 4;} else {factor = 5;}
	  }
	  else if( (*pro.p(i).type()) == photon ) {
		  if (pro.p(i).helicity() == -1 ) {factor = 6;} else {factor = 7;}
	  }
  result += power*factor;
  power *= base;
  }
  BH_DEBUG_MESSAGE4("process: ",pro," helcode: ",result);

return(result);
}


long helcode_2q2l2Q_zvi(const process& pro){
  int factor = 0;
  long result = 0;
  long base  = 8;
  long power = 1;
  int nlegs = pro.n();
  vector<int> flavors;
  int flavor_1=-100;
  int flavor_2=-100;
  int flavor_close,flavor_far;

  for(int i=1;i<=nlegs;i++){

	  if( (*pro.p(i).type()) == quark ) {
		  if (pro.p(i).flavor() == 1){
			  if (pro.p(i).helicity() == -1 ){factor = 0;} else {factor = 1;}
		  } else {
			  if (pro.p(i).helicity() == -1 ){factor = 2;} else {factor = 3;}
		  }
	  }

	  else if( (*pro.p(i).type()) == gluon ) {
		  if (pro.p(i).helicity() == -1 ) {factor = 4;} else {factor = 5;}
	  }
	  else if( (*pro.p(i).type()) == lepton ) {
		  if (pro.p(i).helicity() == -1 ) {factor = 6;} else {factor = 7;}
	  }
  result += power*factor;
  power *= base;
  }
//_MESSAGE4("process: ",pro," result: ",result);

return(result);
}

long helcode_phi_1q(const process& pro){
    // std::cout << "ENTERING helcode_phi_1q" <<  std::endl;
    long factor = 0;
    long result = 0;
    long base  = 0x10;
    long power = 1;

    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
    	if((pro.p(i).type()->stat()==particle::fermion)){
    		//Only add it to the list if it is not a lepton
    		if(pro.p(i).is_not_a(lepton)){
    			flavours.push_back(pro.p(i).flavor());
    		}
    	}
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}
    	else if(pro.p(i).is_a(gluon_massive) && pro.p(i).helicity()== 0){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 6;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 8;}
    	else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()!=0){factor = 9;}
    	else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()!=0){factor = 0xA;}
        else if(pro.p(i).is_a(gluon_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0xB;}
    	else if(pro.p(i).is_a(gluon_massive) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 0xC;}
		else if(pro.p(i).is_a(scalar_massive) && pro.p(i).helicity()== 0){factor = 0xD;}
        else if(pro.p(i).is_a(higgs) && pro.p(i).flavor()==0){factor = 0xE;}
    	else {std::cout << "ERROR: wrong ptype to helcode_phi_1q for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}

long helcode_phi_SM_1q(const process& pro){
    // std::cout << "ENTERING helcode_phi_1q" <<  std::endl;
    long factor = 0;
    long result = 0;
    long base  = 0x10;
    long power = 1;
    
    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
        if((pro.p(i).type()->stat()==particle::fermion)){
            //Only add it to the list if it is not a lepton
            if(pro.p(i).is_not_a(lepton)){
                flavours.push_back(pro.p(i).flavor());
            }
        }
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());
    
    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
        if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
        else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
        else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}
        else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}
        else if(pro.p(i).is_a(scalar_massive) && pro.p(i).helicity()== 0){factor = 4;}
        else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}
        else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 6;}
        else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}
        else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 8;}
        else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()!=0){factor = 9;}
        else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()!=0){factor = 0xA;}  
		else if(pro.p(i).is_a(higgs) && pro.p(i).flavor()==0){factor = 0xB;}
        else {std::cout << "ERROR: wrong ptype to helcode_phi_1q for " << pro <<  std::endl;throw;}
        result += power*factor;
        power *= base;
    }
    
    return(result);
}

    
long helcode_phi_2q(const process& pro){
    // std::cout << "ENTERING helcode_phi_2q" <<  std::endl;
    long factor = 0;
    long result = 0;
    long base  = 0x100;
    long power = 1;

    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
    	if((pro.p(i).type()->stat()==particle::fermion)){
    		//Only add it to the list if it is not a lepton
    		if(pro.p(i).is_not_a(lepton)){
    			flavours.push_back(pro.p(i).flavor());
    		}
    	}
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 6;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[1]){factor = 8;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[1]){factor = 9;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xA;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xB;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xC;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xD;}
    	else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()!=0){factor = 0xE;}
    	else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()!=0){factor = 0xF;}
		else if(pro.p(i).is_a(higgs) && pro.p(i).flavor()==0){factor = 0x10;}
    	else {std::cout << "ERROR: wrong ptype to helcode_phi_2q for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}

long helcode_phi_2q_full(const process& pro){
    // std::cout << "ENTERING helcode_phi_2q" <<  std::endl;
    long factor = 0;
    long result = 0;
    long base  = 0x100;
    long power = 1;

    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
    	if((pro.p(i).type()->stat()==particle::fermion)){
    		//Only add it to the list if it is not a lepton
    		if(pro.p(i).is_not_a(lepton)){
    			flavours.push_back(pro.p(i).flavor());
    		}
    	}
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}
    	else if(pro.p(i).is_a(gluon_massive) && pro.p(i).helicity()== 0){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 6;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 8;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[1]){factor = 9;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[1]){factor = 0xA;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xB;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xC;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xD;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xE;}
    	else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()!=0){factor = 0xF;}
    	else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()!=0){factor = 0x10;}
		else if(pro.p(i).is_a(higgs) && pro.p(i).flavor()==0){factor = 0x11;}
    	else {std::cout << "ERROR: wrong ptype to helcode_phi_2q for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}


int helcode_2qs_massive(const process& pro){
  int factor = 0;
  int result = 0;
  int base  = 0x10;
  int power = 1;

  // We count down rather than up so that the number is in the correct order of the process as read left to right
  for (int i=pro.n();i>0;i--){
	  if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
	  else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
	  else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 ){factor = 2;}
	  else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 ){factor = 3;}
	  else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
	  else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false ){factor = 5;}
	  else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false ){factor = 6;}
	  else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true ){factor = 7;}
	  else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true ){factor = 8;}
	  else {std::cout << "ERROR: wrong ptype to helcode_2qs_massive for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
	  result += power*factor;
	  //   std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
	  power *= base;
  }

  return(result);
}

int helcode_2Q2qs_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 0x10;
    int power = 1;

    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
    	if((pro.p(i).type()->stat()==particle::fermion)){
    		flavours.push_back(pro.p(i).flavor());
    	}
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}//q-
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}//q+
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}//Q-
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 6;}//Q+
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}//Qb-
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 8;}//Qb+
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[1]){factor = 9;}//q2-
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[1]){factor = 0xA;}//q2+
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xB;}//Q2-
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xC;}//Q2+
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xD;}//Qb2-
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xE;}//Qb2+
    	else {std::cout << "ERROR: wrong ptype to helcode_2qs_massive for " << pro <<  std::endl;/*throw;*/}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}

int helcode_2Q1g1y_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 0x10;
    int power = 1;

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== -1){factor = 2;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== +1){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false){factor = 5;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false){factor = 6;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true){factor = 7;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true){factor = 8;}
     	else if(pro.p(i).is_a(lepton) && pro.p(i).helicity()== -1){factor = 9;}
     	else if(pro.p(i).is_a(lepton) && pro.p(i).helicity()== +1){factor = 0xA;}
     	else {std::cout << "ERROR: wrong ptype to helcode_2qs_massive for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}

int helcode_2Q2qs_lepton_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 20;
    int power = 1;

    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
    	if((pro.p(i).type()->stat()==particle::fermion)){
    		//Only add it to the list if it is not a lepton
    		if(pro.p(i).is_not_a(lepton)){
    			flavours.push_back(pro.p(i).flavor());
    		}
    	}
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 6;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 8;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[1]){factor = 9;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[1]){factor = 0xA;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xB;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xC;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xD;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xE;}
    	else if(pro.p(i).is_a(lepton) && pro.p(i).helicity()== -1){factor = 0xF;}
    	else if(pro.p(i).is_a(lepton) && pro.p(i).helicity()== +1){factor = 0x10;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== -1){factor = 0x11;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== +1){factor = 0x12;}
    	else {std::cout << "ERROR: wrong ptype to helcode_2qs_massive for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}


int helcode_2Ls_massive(const process& pro){
  int factor = 0;
  int result = 0;
  int base  = 0x10;
  int power = 1;

  // We count down rather than up so that the number is in the correct order of the process as read left to right
  for (int i=pro.n();i>0;i--){
	  if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
	  else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
	  else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== -1 ){factor = 2;}
	  else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== +1 ){factor = 3;}
	  else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
	  else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false ){factor = 5;}
	  else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false ){factor = 6;}
	  else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true ){factor = 7;}
	  else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true ){factor = 8;}
	  else {std::cout << "ERROR: wrong ptype to helcode_2Ls_massive for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
	  result += power*factor;
	  //   std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
	  power *= base;
  }

  return(result);
}

// We dont care about flavour here
int helcode_2L2qs_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 0x10;
    int power = 1;

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 ){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 ){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false ){factor = 5;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false ){factor = 6;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true ){factor = 7;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true ){factor = 8;}
    	else {std::cout << "ERROR: wrong ptype to helcode_2L2Gs_massive for " << pro <<  std::endl;/*throw;*/}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}

// We dont care about flavour here
int helcode_2L2qs1y_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 0x10;
    int power = 1;

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(photon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 ){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 ){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false ){factor = 5;}//Lm
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false ){factor = 6;}//Lp
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true ){factor = 7;}//Lbm
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true ){factor = 8;}//Lbp
    	else {std::cout << "ERROR: wrong ptype to helcode_2L2Gs_massive for " << pro <<  std::endl;/*throw;*/}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}


int helcode_2L2Gs_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 0x10;
    int power = 1;

    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
    	if((pro.p(i).type()->stat()==particle::fermion)){
    		flavours.push_back(pro.p(i).flavor());
    	}
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 6;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 8;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[1]){factor = 9;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[1]){factor = 0xA;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xB;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xC;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xD;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xE;}
    	else {std::cout << "ERROR: wrong ptype to helcode_2L2Gs_massive for " << pro <<  std::endl;/*throw;*/}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}


int helcode_2L2Gs_lepton_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 20;
    int power = 1;

    //Find all the different fermion flavors
    vector<int> flavours;
    for(int i=1;i<=pro.n();i++)
    {
    	if((pro.p(i).type()->stat()==particle::fermion)){
    		//Only add it to the list if it is not a lepton
    		if(pro.p(i).is_not_a(lepton)){
    			flavours.push_back(pro.p(i).flavor());
    		}
    	}
    }
    sort(flavours.begin(),flavours.end());
    flavours.erase(unique(flavours.begin(),flavours.end()),flavours.end());

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[0] ){factor = 2;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[0] ){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 5;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[0] ){factor = 6;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 7;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[0] ){factor = 8;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== -1 && pro.p(i).flavor()==flavours[1]){factor = 9;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== +1 && pro.p(i).flavor()==flavours[1]){factor = 0xA;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xB;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()==flavours[1] ){factor = 0xC;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xD;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()==flavours[1] ){factor = 0xE;}
    	else if(pro.p(i).is_a(lepton) && pro.p(i).helicity()== -1){factor = 0xF;}
    	else if(pro.p(i).is_a(lepton) && pro.p(i).helicity()== +1){factor = 0x10;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== -1){factor = 0x11;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== +1){factor = 0x12;}
    	else {std::cout << "ERROR: wrong ptype to helcode_2L2Gs_massive_lepton for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}

int helcode_1Q1q1G1L_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 0x10;
    int power = 1;

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false){factor = 5;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false){factor = 6;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true){factor = 7;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true){factor = 8;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== -1){factor = 9;}
    	else if(pro.p(i).is_a(gluino) && pro.p(i).helicity()== +1){factor = 0xA;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false){factor = 0xB;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false){factor = 0xC;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true){factor = 0xD;}
    	else if(pro.p(i).is_a(gluino_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true){factor = 0xE;}
    	else {std::cout << "ERROR: wrong ptype to helcode_2qs_massive for " << pro <<  std::endl;/*throw;*/}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}


// We dont care about flavour here
int helcode_2Q2qs1y_massive(const process& pro){
    // std::cout << "ENTERING helcode_2q" <<  std::endl;
    int factor = 0;
    int result = 0;
    int base  = 0x10;
    int power = 1;

    // We count down rather than up so that the number is in the correct order of the process as read left to right
    for (int i=pro.n();i>0;i--){
    	if(pro.p(i).is_a(photon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
    	else if(pro.p(i).is_a(photon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 ){factor = 2;}
    	else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 ){factor = 3;}
    	else if(pro.p(i).is_a(scalar_massive)){factor = 4;}
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==false ){factor = 5;}//Lm
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==false ){factor = 6;}//Lp
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()==true ){factor = 7;}//Lbm
    	else if(pro.p(i).is_a(quark_massive) && pro.p(i).helicity()== +1 && pro.p(i).is_anti_particle()==true ){factor = 8;}//Lbp
    	else {std::cout << "ERROR: wrong ptype to helcode_2L2Gs_massive for " << pro <<  std::endl;/*throw;*/}
    	result += power*factor;
    	power *= base;
    }

    return(result);
}

int helcode_2s(const process& pro){
	int factor = 0;
	int result = 0;
	int base  = 4;
	int power = 1;
	for(int i=1;i<=pro.n();i++){
		if((pro.p(i).is_a(gluon))&&(pro.p(i).helicity()==-1)){factor = 0;}
		else if(pro.p(i).is_a(scalar_massive)){factor = 2;}
		else if((pro.p(i).is_a(gluon))&&(pro.p(i).helicity()==1)){factor = 3;}
		else {std::cout << "ERROR: wrong ptype " << pro.p(i) << " to helcode_2s" << std::endl;}
		result += power*factor;
		//cout<< i << " " << factor << " " << power << " " << base <<endl;
		power *= base;
    }
//	_MESSAGE4("helcode_2s : for process ",pro," have code ",result);
	return(result);
}


long helcode(const process& pro){
	switch (pro.pcode()/10){
	case 0: return helcode_g(pro);
	case 2: return helcode_2q(pro);
	case 4: {
		switch (nbr_of_flavors(pro,quark)){
		case 1: return helcode_4q(pro);
		case 2: return helcode_2q2Q(pro);
		}
	}
	case 22: return helcode_2q2l(pro);
	case 10002: return helcode_2q1y(pro);
	case 200022: return helcode_2q2l2G(pro);
	case 210002: return helcode_2q2G1y(pro);

}
}
	
int helcode_2qs_massless(const process& pro){
	int factor = 0;
	int result = 0;
	int base  = 0x10;
	int power = 1;
	
	// We count down rather than up so that the number is in the correct order of the process as read left to right
	for (int i=pro.n();i>0;i--){
		if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== -1 && pro.p(i).is_anti_particle()== false  ){factor = 0;}
		else if(pro.p(i).is_a(gluon) && pro.p(i).helicity()== 1 && pro.p(i).is_anti_particle()== false  ){factor = 1;}
		else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== -1 ){factor = 2;}
		else if(pro.p(i).is_a(quark) && pro.p(i).helicity()== +1 ){factor = 3;}
		else if(pro.p(i).is_a(scalar)){factor = 4;}
		else if(pro.p(i).is_a(scalar_massive)){factor = 5;}
		else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==false && pro.p(i).flavor()!=0){factor = 9;}
    	else if(pro.p(i).is_a(higgs) && pro.p(i).is_anti_particle()==true && pro.p(i).flavor()!=0){factor = 0xA;}
		else if(pro.p(i).is_a(higgs) && pro.p(i).flavor()==0){factor = 0xB;}
		else {std::cout << "ERROR: wrong ptype to helcode_2qs_massless for " << pro <<  std::endl;throw BHerror("Wrong ptype");}
		result += power*factor;
		//   std::cout << i << " " << factor << " " << power << " " << base <<std::endl;
		power *= base;
	}
	
	return(result);
}
    
long helcode_2Gsc(const process& pro){
    int factor = 0;
    long result = 0;
    int base  = 0x10;
    int power = 1;
	
	// We count down rather than up so that the number is in the correct order of the process as read left to right
    for(int i=pro.n();i>0;i--){
        if((pro.p(i).is_a(gluon))&&(pro.p(i).helicity()==-1)){factor = 0;}
        else if((pro.p(i).is_a(gluon))&&(pro.p(i).helicity()==1)){factor = 1;}
        else if(pro.p(i).is_a(gluon_massive)&&(pro.p(i).helicity()==0)){factor = 2;}
        else if((pro.p(i).is_a(gluon_massive))&&(pro.p(i).helicity()==-1)){factor = 3;}
        else if(pro.p(i).is_a(gluon_massive)&&(pro.p(i).helicity()==1)){factor = 4;}
        else {std::cout << "ERROR: wrong ptype " << pro.p(i) << " to helcode_2Gsc" << std::endl;}
        result += power*factor;
        //cout<< i << " " << factor << " " << power << " " << base <<endl;
        power *= base;
    }
    //	_MESSAGE4("helcode_2s : for process ",pro," have code ",result);
    return(result);
}
	
}
