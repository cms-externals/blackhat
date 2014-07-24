/*
 * matrix elements for lm lbp -> 2q5g2l
 *
 *  Created on: Oct 21, 2008
 *
 */


#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "BH_interface_impl.h"




#define _VERBOSE 0
#define single_hel 0
#define Zboson 0 //set 0 ... for Ws only in order to use VSM with prefactor 
#define threads 3 //[0,3,4,6] = optimized splitup of computation load for [arbitrary,3,4,6]; splits with helicity

// various switches: photonZW:  0 for photon
// 				1 for photon+Z
// 				2 for Z->neutrinos only
// 				3 for W only
#define _PHOTON_ONLY 1     // 0 photon only; 1 for all other


using namespace std;

namespace BH {

int tree_count(0);




vector<int> index9(vector<int> & indin, int num){
	  int index(num);
	  int n=indin.size()-2;
	  vector<int> ind;
	  for(int i=0;i<n;i++) {
	       int j=index/int(std::pow(10.,n-1-i));
	       index-=j*int(std::pow(10.,n-1-i));
	       ind.push_back(indin[j-1]);}
	       ind.push_back(indin[n]);
	       ind.push_back(indin[n+1]);
	  return ind;
}


process q5gq(process & pro, int num){
	tree_count++;
	ph_type h1,h2,h3,h4,h5,h6,h7,h8,h9;
        int n=pro.n()-2;
        int index(num);
        vector<int> ind;
        for(int i=0;i<n;i++) {
                int j=index/int(std::pow(10.,n-1-i));
                index-=j*int(std::pow(10.,n-1-i));
                ind.push_back(j);}
	h1=pro.p(ind[1-1]);
	h2=pro.p(ind[2-1]);
	h3=pro.p(ind[3-1]);
	h4=pro.p(ind[4-1]);
	h5=pro.p(ind[5-1]);
	h6=pro.p(ind[6-1]);
	h7=pro.p(ind[7-1]);
	h8=pro.p(8);
	h9=pro.p(9);
	return process(h1,h2,h3,h4,h5,h6,h7,h8,h9);	
}


Squared_ME* A_loop_2q_5g_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, int color, int tree_color, const vector<ph_type> _ph_type,QCDorder lo_or_nlo, Squared_ME* SM=0){


size_t T_1234567_up;size_t L_color_1234567_up;
size_t T_1234657_up;size_t L_color_1234657_up;
size_t T_1235467_up;size_t L_color_1235467_up;
size_t T_1235647_up;size_t L_color_1235647_up;
size_t T_1236457_up;size_t L_color_1236457_up;
size_t T_1236547_up;size_t L_color_1236547_up;
size_t T_1243567_up;size_t L_color_1243567_up;
size_t T_1243657_up;size_t L_color_1243657_up;
size_t T_1245367_up;size_t L_color_1245367_up;
size_t T_1245637_up;size_t L_color_1245637_up;
size_t T_1246357_up;size_t L_color_1246357_up;
size_t T_1246537_up;size_t L_color_1246537_up;
size_t T_1253467_up;size_t L_color_1253467_up;
size_t T_1253647_up;size_t L_color_1253647_up;
size_t T_1254367_up;size_t L_color_1254367_up;
size_t T_1254637_up;size_t L_color_1254637_up;
size_t T_1256347_up;size_t L_color_1256347_up;
size_t T_1256437_up;size_t L_color_1256437_up;
size_t T_1263457_up;size_t L_color_1263457_up;
size_t T_1263547_up;size_t L_color_1263547_up;
size_t T_1264357_up;size_t L_color_1264357_up;
size_t T_1264537_up;size_t L_color_1264537_up;
size_t T_1265347_up;size_t L_color_1265347_up;
size_t T_1265437_up;size_t L_color_1265437_up;
size_t T_1324567_up;size_t L_color_1324567_up;
size_t T_1324657_up;size_t L_color_1324657_up;
size_t T_1325467_up;size_t L_color_1325467_up;
size_t T_1325647_up;size_t L_color_1325647_up;
size_t T_1326457_up;size_t L_color_1326457_up;
size_t T_1326547_up;size_t L_color_1326547_up;
size_t T_1342567_up;size_t L_color_1342567_up;
size_t T_1342657_up;size_t L_color_1342657_up;
size_t T_1345267_up;size_t L_color_1345267_up;
size_t T_1345627_up;size_t L_color_1345627_up;
size_t T_1346257_up;size_t L_color_1346257_up;
size_t T_1346527_up;size_t L_color_1346527_up;
size_t T_1352467_up;size_t L_color_1352467_up;
size_t T_1352647_up;size_t L_color_1352647_up;
size_t T_1354267_up;size_t L_color_1354267_up;
size_t T_1354627_up;size_t L_color_1354627_up;
size_t T_1356247_up;size_t L_color_1356247_up;
size_t T_1356427_up;size_t L_color_1356427_up;
size_t T_1362457_up;size_t L_color_1362457_up;
size_t T_1362547_up;size_t L_color_1362547_up;
size_t T_1364257_up;size_t L_color_1364257_up;
size_t T_1364527_up;size_t L_color_1364527_up;
size_t T_1365247_up;size_t L_color_1365247_up;
size_t T_1365427_up;size_t L_color_1365427_up;
size_t T_1423567_up;size_t L_color_1423567_up;
size_t T_1423657_up;size_t L_color_1423657_up;
size_t T_1425367_up;size_t L_color_1425367_up;
size_t T_1425637_up;size_t L_color_1425637_up;
size_t T_1426357_up;size_t L_color_1426357_up;
size_t T_1426537_up;size_t L_color_1426537_up;
size_t T_1432567_up;size_t L_color_1432567_up;
size_t T_1432657_up;size_t L_color_1432657_up;
size_t T_1435267_up;size_t L_color_1435267_up;
size_t T_1435627_up;size_t L_color_1435627_up;
size_t T_1436257_up;size_t L_color_1436257_up;
size_t T_1436527_up;size_t L_color_1436527_up;
size_t T_1452367_up;size_t L_color_1452367_up;
size_t T_1452637_up;size_t L_color_1452637_up;
size_t T_1453267_up;size_t L_color_1453267_up;
size_t T_1453627_up;size_t L_color_1453627_up;
size_t T_1456237_up;size_t L_color_1456237_up;
size_t T_1456327_up;size_t L_color_1456327_up;
size_t T_1462357_up;size_t L_color_1462357_up;
size_t T_1462537_up;size_t L_color_1462537_up;
size_t T_1463257_up;size_t L_color_1463257_up;
size_t T_1463527_up;size_t L_color_1463527_up;
size_t T_1465237_up;size_t L_color_1465237_up;
size_t T_1465327_up;size_t L_color_1465327_up;
size_t T_1523467_up;size_t L_color_1523467_up;
size_t T_1523647_up;size_t L_color_1523647_up;
size_t T_1524367_up;size_t L_color_1524367_up;
size_t T_1524637_up;size_t L_color_1524637_up;
size_t T_1526347_up;size_t L_color_1526347_up;
size_t T_1526437_up;size_t L_color_1526437_up;
size_t T_1532467_up;size_t L_color_1532467_up;
size_t T_1532647_up;size_t L_color_1532647_up;
size_t T_1534267_up;size_t L_color_1534267_up;
size_t T_1534627_up;size_t L_color_1534627_up;
size_t T_1536247_up;size_t L_color_1536247_up;
size_t T_1536427_up;size_t L_color_1536427_up;
size_t T_1542367_up;size_t L_color_1542367_up;
size_t T_1542637_up;size_t L_color_1542637_up;
size_t T_1543267_up;size_t L_color_1543267_up;
size_t T_1543627_up;size_t L_color_1543627_up;
size_t T_1546237_up;size_t L_color_1546237_up;
size_t T_1546327_up;size_t L_color_1546327_up;
size_t T_1562347_up;size_t L_color_1562347_up;
size_t T_1562437_up;size_t L_color_1562437_up;
size_t T_1563247_up;size_t L_color_1563247_up;
size_t T_1563427_up;size_t L_color_1563427_up;
size_t T_1564237_up;size_t L_color_1564237_up;
size_t T_1564327_up;size_t L_color_1564327_up;
size_t T_1623457_up;size_t L_color_1623457_up;
size_t T_1623547_up;size_t L_color_1623547_up;
size_t T_1624357_up;size_t L_color_1624357_up;
size_t T_1624537_up;size_t L_color_1624537_up;
size_t T_1625347_up;size_t L_color_1625347_up;
size_t T_1625437_up;size_t L_color_1625437_up;
size_t T_1632457_up;size_t L_color_1632457_up;
size_t T_1632547_up;size_t L_color_1632547_up;
size_t T_1634257_up;size_t L_color_1634257_up;
size_t T_1634527_up;size_t L_color_1634527_up;
size_t T_1635247_up;size_t L_color_1635247_up;
size_t T_1635427_up;size_t L_color_1635427_up;
size_t T_1642357_up;size_t L_color_1642357_up;
size_t T_1642537_up;size_t L_color_1642537_up;
size_t T_1643257_up;size_t L_color_1643257_up;
size_t T_1643527_up;size_t L_color_1643527_up;
size_t T_1645237_up;size_t L_color_1645237_up;
size_t T_1645327_up;size_t L_color_1645327_up;
size_t T_1652347_up;size_t L_color_1652347_up;
size_t T_1652437_up;size_t L_color_1652437_up;
size_t T_1653247_up;size_t L_color_1653247_up;
size_t T_1653427_up;size_t L_color_1653427_up;
size_t T_1654237_up;size_t L_color_1654237_up;
size_t T_1654327_up;size_t L_color_1654327_up;

#if Zboson==1
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,ind[7],ind[8],_ph_type);
//-------------------------------------------
#endif

double Nc(n_c);
multi_precision_fraction one(1);

for(int l1 = 1; l1 < 6;++l1){
	for(int l2 = 1; l2 < 6; ++l2){
		for(int l3 = 1; l3 < 6; ++l3){
			for(int l4 = 1; l4 < 6; ++l4){
				for(int l5 = 1; l5 < 6; ++l5){
			if((l1!=l2)&&(l1!=l3)&&(l1!=l4)&&(l1!=l5)&&(l2!=l3)&&(l2!=l4)&&(l2!=l5)&&(l3!=l4)&&(l3!=l5)&&(l4!=l5)){					


int i1=ind.at(0);
int i2=ind.at(l1);
int i3=ind.at(l2);
int i4=ind.at(l3);
int i5=ind.at(l4);
int i6=ind.at(l5);
int i7=ind.at(6);
int i8=ind.at(7);
int i9=ind.at(8);

ph_type h1=pro.p(1);
ph_type h2=pro.p(l1+1);
ph_type h3=pro.p(l2+1);
ph_type h4=pro.p(l3+1);
ph_type h5=pro.p(l4+1);
ph_type h6=pro.p(l5+1);
ph_type h7=pro.p(7);
ph_type h8=pro.p(8);
ph_type h9=pro.p(9);

process PRO(h1,h2,h3,h4,h5,h6,h7,h8,h9);
vector<int> IND;
IND.push_back(i1);
IND.push_back(i2);
IND.push_back(i3);
IND.push_back(i4);
IND.push_back(i5);
IND.push_back(i6);
IND.push_back(i7);
IND.push_back(i8);
IND.push_back(i9);




#if Zboson==1
T_1234567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1234567),index9(IND,1234567),_prop_hel_fn),index9(IND,1234567));
SM->add_tree(cached_cross_term_md(T_1234567_up,T_1234567_up,one,4*std::pow(-1 + std::pow(Nc,2),5)/std::pow(Nc,4))); 

if(tree_color==0){ 
//T_1234567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1234567),index9(IND,1234567),_prop_hel_fn),index9(IND,1234567));
T_1234657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1234657),index9(IND,1234657),_prop_hel_fn),index9(IND,1234657));
T_1235467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1235467),index9(IND,1235467),_prop_hel_fn),index9(IND,1235467));
T_1235647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1235647),index9(IND,1235647),_prop_hel_fn),index9(IND,1235647));
T_1236457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1236457),index9(IND,1236457),_prop_hel_fn),index9(IND,1236457));
T_1236547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1236547),index9(IND,1236547),_prop_hel_fn),index9(IND,1236547));
T_1243567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1243567),index9(IND,1243567),_prop_hel_fn),index9(IND,1243567));
T_1243657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1243657),index9(IND,1243657),_prop_hel_fn),index9(IND,1243657));
T_1245367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1245367),index9(IND,1245367),_prop_hel_fn),index9(IND,1245367));
T_1245637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1245637),index9(IND,1245637),_prop_hel_fn),index9(IND,1245637));
T_1246357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1246357),index9(IND,1246357),_prop_hel_fn),index9(IND,1246357));
T_1246537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1246537),index9(IND,1246537),_prop_hel_fn),index9(IND,1246537));
T_1253467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1253467),index9(IND,1253467),_prop_hel_fn),index9(IND,1253467));
T_1253647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1253647),index9(IND,1253647),_prop_hel_fn),index9(IND,1253647));
T_1254367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1254367),index9(IND,1254367),_prop_hel_fn),index9(IND,1254367));
T_1254637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1254637),index9(IND,1254637),_prop_hel_fn),index9(IND,1254637));
T_1256347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1256347),index9(IND,1256347),_prop_hel_fn),index9(IND,1256347));
T_1256437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1256437),index9(IND,1256437),_prop_hel_fn),index9(IND,1256437));
T_1263457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1263457),index9(IND,1263457),_prop_hel_fn),index9(IND,1263457));
T_1263547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1263547),index9(IND,1263547),_prop_hel_fn),index9(IND,1263547));
T_1264357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1264357),index9(IND,1264357),_prop_hel_fn),index9(IND,1264357));
T_1264537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1264537),index9(IND,1264537),_prop_hel_fn),index9(IND,1264537));
T_1265347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1265347),index9(IND,1265347),_prop_hel_fn),index9(IND,1265347));
T_1265437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1265437),index9(IND,1265437),_prop_hel_fn),index9(IND,1265437));
T_1324567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1324567),index9(IND,1324567),_prop_hel_fn),index9(IND,1324567));
T_1324657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1324657),index9(IND,1324657),_prop_hel_fn),index9(IND,1324657));
T_1325467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1325467),index9(IND,1325467),_prop_hel_fn),index9(IND,1325467));
T_1325647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1325647),index9(IND,1325647),_prop_hel_fn),index9(IND,1325647));
T_1326457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1326457),index9(IND,1326457),_prop_hel_fn),index9(IND,1326457));
T_1326547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1326547),index9(IND,1326547),_prop_hel_fn),index9(IND,1326547));
T_1342567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1342567),index9(IND,1342567),_prop_hel_fn),index9(IND,1342567));
T_1342657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1342657),index9(IND,1342657),_prop_hel_fn),index9(IND,1342657));
T_1345267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1345267),index9(IND,1345267),_prop_hel_fn),index9(IND,1345267));
T_1345627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1345627),index9(IND,1345627),_prop_hel_fn),index9(IND,1345627));
T_1346257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1346257),index9(IND,1346257),_prop_hel_fn),index9(IND,1346257));
T_1346527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1346527),index9(IND,1346527),_prop_hel_fn),index9(IND,1346527));
T_1352467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1352467),index9(IND,1352467),_prop_hel_fn),index9(IND,1352467));
T_1352647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1352647),index9(IND,1352647),_prop_hel_fn),index9(IND,1352647));
T_1354267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1354267),index9(IND,1354267),_prop_hel_fn),index9(IND,1354267));
T_1354627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1354627),index9(IND,1354627),_prop_hel_fn),index9(IND,1354627));
T_1356247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1356247),index9(IND,1356247),_prop_hel_fn),index9(IND,1356247));
T_1356427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1356427),index9(IND,1356427),_prop_hel_fn),index9(IND,1356427));
T_1362457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1362457),index9(IND,1362457),_prop_hel_fn),index9(IND,1362457));
T_1362547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1362547),index9(IND,1362547),_prop_hel_fn),index9(IND,1362547));
T_1364257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1364257),index9(IND,1364257),_prop_hel_fn),index9(IND,1364257));
T_1364527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1364527),index9(IND,1364527),_prop_hel_fn),index9(IND,1364527));
T_1365247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1365247),index9(IND,1365247),_prop_hel_fn),index9(IND,1365247));
T_1365427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1365427),index9(IND,1365427),_prop_hel_fn),index9(IND,1365427));
T_1423567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1423567),index9(IND,1423567),_prop_hel_fn),index9(IND,1423567));
T_1423657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1423657),index9(IND,1423657),_prop_hel_fn),index9(IND,1423657));
T_1425367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1425367),index9(IND,1425367),_prop_hel_fn),index9(IND,1425367));
T_1425637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1425637),index9(IND,1425637),_prop_hel_fn),index9(IND,1425637));
T_1426357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1426357),index9(IND,1426357),_prop_hel_fn),index9(IND,1426357));
T_1426537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1426537),index9(IND,1426537),_prop_hel_fn),index9(IND,1426537));
T_1432567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1432567),index9(IND,1432567),_prop_hel_fn),index9(IND,1432567));
T_1432657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1432657),index9(IND,1432657),_prop_hel_fn),index9(IND,1432657));
T_1435267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1435267),index9(IND,1435267),_prop_hel_fn),index9(IND,1435267));
T_1435627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1435627),index9(IND,1435627),_prop_hel_fn),index9(IND,1435627));
T_1436257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1436257),index9(IND,1436257),_prop_hel_fn),index9(IND,1436257));
T_1436527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1436527),index9(IND,1436527),_prop_hel_fn),index9(IND,1436527));
T_1452367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1452367),index9(IND,1452367),_prop_hel_fn),index9(IND,1452367));
T_1452637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1452637),index9(IND,1452637),_prop_hel_fn),index9(IND,1452637));
T_1453267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1453267),index9(IND,1453267),_prop_hel_fn),index9(IND,1453267));
T_1453627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1453627),index9(IND,1453627),_prop_hel_fn),index9(IND,1453627));
T_1456237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1456237),index9(IND,1456237),_prop_hel_fn),index9(IND,1456237));
T_1456327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1456327),index9(IND,1456327),_prop_hel_fn),index9(IND,1456327));
T_1462357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1462357),index9(IND,1462357),_prop_hel_fn),index9(IND,1462357));
T_1462537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1462537),index9(IND,1462537),_prop_hel_fn),index9(IND,1462537));
T_1463257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1463257),index9(IND,1463257),_prop_hel_fn),index9(IND,1463257));
T_1463527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1463527),index9(IND,1463527),_prop_hel_fn),index9(IND,1463527));
T_1465237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1465237),index9(IND,1465237),_prop_hel_fn),index9(IND,1465237));
T_1465327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1465327),index9(IND,1465327),_prop_hel_fn),index9(IND,1465327));
T_1523467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1523467),index9(IND,1523467),_prop_hel_fn),index9(IND,1523467));
T_1523647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1523647),index9(IND,1523647),_prop_hel_fn),index9(IND,1523647));
T_1524367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1524367),index9(IND,1524367),_prop_hel_fn),index9(IND,1524367));
T_1524637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1524637),index9(IND,1524637),_prop_hel_fn),index9(IND,1524637));
T_1526347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1526347),index9(IND,1526347),_prop_hel_fn),index9(IND,1526347));
T_1526437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1526437),index9(IND,1526437),_prop_hel_fn),index9(IND,1526437));
T_1532467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1532467),index9(IND,1532467),_prop_hel_fn),index9(IND,1532467));
T_1532647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1532647),index9(IND,1532647),_prop_hel_fn),index9(IND,1532647));
T_1534267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1534267),index9(IND,1534267),_prop_hel_fn),index9(IND,1534267));
T_1534627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1534627),index9(IND,1534627),_prop_hel_fn),index9(IND,1534627));
T_1536247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1536247),index9(IND,1536247),_prop_hel_fn),index9(IND,1536247));
T_1536427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1536427),index9(IND,1536427),_prop_hel_fn),index9(IND,1536427));
T_1542367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1542367),index9(IND,1542367),_prop_hel_fn),index9(IND,1542367));
T_1542637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1542637),index9(IND,1542637),_prop_hel_fn),index9(IND,1542637));
T_1543267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1543267),index9(IND,1543267),_prop_hel_fn),index9(IND,1543267));
T_1543627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1543627),index9(IND,1543627),_prop_hel_fn),index9(IND,1543627));
T_1546237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1546237),index9(IND,1546237),_prop_hel_fn),index9(IND,1546237));
T_1546327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1546327),index9(IND,1546327),_prop_hel_fn),index9(IND,1546327));
T_1562347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1562347),index9(IND,1562347),_prop_hel_fn),index9(IND,1562347));
T_1562437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1562437),index9(IND,1562437),_prop_hel_fn),index9(IND,1562437));
T_1563247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1563247),index9(IND,1563247),_prop_hel_fn),index9(IND,1563247));
T_1563427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1563427),index9(IND,1563427),_prop_hel_fn),index9(IND,1563427));
T_1564237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1564237),index9(IND,1564237),_prop_hel_fn),index9(IND,1564237));
T_1564327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1564327),index9(IND,1564327),_prop_hel_fn),index9(IND,1564327));
T_1623457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1623457),index9(IND,1623457),_prop_hel_fn),index9(IND,1623457));
T_1623547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1623547),index9(IND,1623547),_prop_hel_fn),index9(IND,1623547));
T_1624357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1624357),index9(IND,1624357),_prop_hel_fn),index9(IND,1624357));
T_1624537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1624537),index9(IND,1624537),_prop_hel_fn),index9(IND,1624537));
T_1625347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1625347),index9(IND,1625347),_prop_hel_fn),index9(IND,1625347));
T_1625437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1625437),index9(IND,1625437),_prop_hel_fn),index9(IND,1625437));
T_1632457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1632457),index9(IND,1632457),_prop_hel_fn),index9(IND,1632457));
T_1632547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1632547),index9(IND,1632547),_prop_hel_fn),index9(IND,1632547));
T_1634257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1634257),index9(IND,1634257),_prop_hel_fn),index9(IND,1634257));
T_1634527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1634527),index9(IND,1634527),_prop_hel_fn),index9(IND,1634527));
T_1635247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1635247),index9(IND,1635247),_prop_hel_fn),index9(IND,1635247));
T_1635427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1635427),index9(IND,1635427),_prop_hel_fn),index9(IND,1635427));
T_1642357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1642357),index9(IND,1642357),_prop_hel_fn),index9(IND,1642357));
T_1642537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1642537),index9(IND,1642537),_prop_hel_fn),index9(IND,1642537));
T_1643257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1643257),index9(IND,1643257),_prop_hel_fn),index9(IND,1643257));
T_1643527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1643527),index9(IND,1643527),_prop_hel_fn),index9(IND,1643527));
T_1645237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1645237),index9(IND,1645237),_prop_hel_fn),index9(IND,1645237));
T_1645327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1645327),index9(IND,1645327),_prop_hel_fn),index9(IND,1645327));
T_1652347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1652347),index9(IND,1652347),_prop_hel_fn),index9(IND,1652347));
T_1652437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1652437),index9(IND,1652437),_prop_hel_fn),index9(IND,1652437));
T_1653247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1653247),index9(IND,1653247),_prop_hel_fn),index9(IND,1653247));
T_1653427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1653427),index9(IND,1653427),_prop_hel_fn),index9(IND,1653427));
T_1654237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1654237),index9(IND,1654237),_prop_hel_fn),index9(IND,1654237));
T_1654327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1654327),index9(IND,1654327),_prop_hel_fn),index9(IND,1654327));
}
#endif

#if Zboson==0
T_1234567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1234567),index9(IND,1234567)),index9(IND,1234567));
SM->add_tree(cached_cross_term_md(T_1234567_up,T_1234567_up,one,4*std::pow(-1 + std::pow(Nc,2),5)/std::pow(Nc,4))); 

if(tree_color==0){ 
//T_1234567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1234567),index9(IND,1234567)),index9(IND,1234567));
T_1234657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1234657),index9(IND,1234657)),index9(IND,1234657));
T_1235467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1235467),index9(IND,1235467)),index9(IND,1235467));
T_1235647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1235647),index9(IND,1235647)),index9(IND,1235647));
T_1236457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1236457),index9(IND,1236457)),index9(IND,1236457));
T_1236547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1236547),index9(IND,1236547)),index9(IND,1236547));
T_1243567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1243567),index9(IND,1243567)),index9(IND,1243567));
T_1243657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1243657),index9(IND,1243657)),index9(IND,1243657));
T_1245367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1245367),index9(IND,1245367)),index9(IND,1245367));
T_1245637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1245637),index9(IND,1245637)),index9(IND,1245637));
T_1246357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1246357),index9(IND,1246357)),index9(IND,1246357));
T_1246537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1246537),index9(IND,1246537)),index9(IND,1246537));
T_1253467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1253467),index9(IND,1253467)),index9(IND,1253467));
T_1253647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1253647),index9(IND,1253647)),index9(IND,1253647));
T_1254367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1254367),index9(IND,1254367)),index9(IND,1254367));
T_1254637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1254637),index9(IND,1254637)),index9(IND,1254637));
T_1256347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1256347),index9(IND,1256347)),index9(IND,1256347));
T_1256437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1256437),index9(IND,1256437)),index9(IND,1256437));
T_1263457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1263457),index9(IND,1263457)),index9(IND,1263457));
T_1263547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1263547),index9(IND,1263547)),index9(IND,1263547));
T_1264357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1264357),index9(IND,1264357)),index9(IND,1264357));
T_1264537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1264537),index9(IND,1264537)),index9(IND,1264537));
T_1265347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1265347),index9(IND,1265347)),index9(IND,1265347));
T_1265437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1265437),index9(IND,1265437)),index9(IND,1265437));
T_1324567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1324567),index9(IND,1324567)),index9(IND,1324567));
T_1324657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1324657),index9(IND,1324657)),index9(IND,1324657));
T_1325467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1325467),index9(IND,1325467)),index9(IND,1325467));
T_1325647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1325647),index9(IND,1325647)),index9(IND,1325647));
T_1326457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1326457),index9(IND,1326457)),index9(IND,1326457));
T_1326547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1326547),index9(IND,1326547)),index9(IND,1326547));
T_1342567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1342567),index9(IND,1342567)),index9(IND,1342567));
T_1342657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1342657),index9(IND,1342657)),index9(IND,1342657));
T_1345267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1345267),index9(IND,1345267)),index9(IND,1345267));
T_1345627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1345627),index9(IND,1345627)),index9(IND,1345627));
T_1346257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1346257),index9(IND,1346257)),index9(IND,1346257));
T_1346527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1346527),index9(IND,1346527)),index9(IND,1346527));
T_1352467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1352467),index9(IND,1352467)),index9(IND,1352467));
T_1352647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1352647),index9(IND,1352647)),index9(IND,1352647));
T_1354267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1354267),index9(IND,1354267)),index9(IND,1354267));
T_1354627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1354627),index9(IND,1354627)),index9(IND,1354627));
T_1356247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1356247),index9(IND,1356247)),index9(IND,1356247));
T_1356427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1356427),index9(IND,1356427)),index9(IND,1356427));
T_1362457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1362457),index9(IND,1362457)),index9(IND,1362457));
T_1362547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1362547),index9(IND,1362547)),index9(IND,1362547));
T_1364257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1364257),index9(IND,1364257)),index9(IND,1364257));
T_1364527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1364527),index9(IND,1364527)),index9(IND,1364527));
T_1365247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1365247),index9(IND,1365247)),index9(IND,1365247));
T_1365427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1365427),index9(IND,1365427)),index9(IND,1365427));
T_1423567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1423567),index9(IND,1423567)),index9(IND,1423567));
T_1423657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1423657),index9(IND,1423657)),index9(IND,1423657));
T_1425367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1425367),index9(IND,1425367)),index9(IND,1425367));
T_1425637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1425637),index9(IND,1425637)),index9(IND,1425637));
T_1426357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1426357),index9(IND,1426357)),index9(IND,1426357));
T_1426537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1426537),index9(IND,1426537)),index9(IND,1426537));
T_1432567_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1432567),index9(IND,1432567)),index9(IND,1432567));
T_1432657_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1432657),index9(IND,1432657)),index9(IND,1432657));
T_1435267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1435267),index9(IND,1435267)),index9(IND,1435267));
T_1435627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1435627),index9(IND,1435627)),index9(IND,1435627));
T_1436257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1436257),index9(IND,1436257)),index9(IND,1436257));
T_1436527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1436527),index9(IND,1436527)),index9(IND,1436527));
T_1452367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1452367),index9(IND,1452367)),index9(IND,1452367));
T_1452637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1452637),index9(IND,1452637)),index9(IND,1452637));
T_1453267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1453267),index9(IND,1453267)),index9(IND,1453267));
T_1453627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1453627),index9(IND,1453627)),index9(IND,1453627));
T_1456237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1456237),index9(IND,1456237)),index9(IND,1456237));
T_1456327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1456327),index9(IND,1456327)),index9(IND,1456327));
T_1462357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1462357),index9(IND,1462357)),index9(IND,1462357));
T_1462537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1462537),index9(IND,1462537)),index9(IND,1462537));
T_1463257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1463257),index9(IND,1463257)),index9(IND,1463257));
T_1463527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1463527),index9(IND,1463527)),index9(IND,1463527));
T_1465237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1465237),index9(IND,1465237)),index9(IND,1465237));
T_1465327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1465327),index9(IND,1465327)),index9(IND,1465327));
T_1523467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1523467),index9(IND,1523467)),index9(IND,1523467));
T_1523647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1523647),index9(IND,1523647)),index9(IND,1523647));
T_1524367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1524367),index9(IND,1524367)),index9(IND,1524367));
T_1524637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1524637),index9(IND,1524637)),index9(IND,1524637));
T_1526347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1526347),index9(IND,1526347)),index9(IND,1526347));
T_1526437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1526437),index9(IND,1526437)),index9(IND,1526437));
T_1532467_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1532467),index9(IND,1532467)),index9(IND,1532467));
T_1532647_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1532647),index9(IND,1532647)),index9(IND,1532647));
T_1534267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1534267),index9(IND,1534267)),index9(IND,1534267));
T_1534627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1534627),index9(IND,1534627)),index9(IND,1534627));
T_1536247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1536247),index9(IND,1536247)),index9(IND,1536247));
T_1536427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1536427),index9(IND,1536427)),index9(IND,1536427));
T_1542367_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1542367),index9(IND,1542367)),index9(IND,1542367));
T_1542637_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1542637),index9(IND,1542637)),index9(IND,1542637));
T_1543267_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1543267),index9(IND,1543267)),index9(IND,1543267));
T_1543627_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1543627),index9(IND,1543627)),index9(IND,1543627));
T_1546237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1546237),index9(IND,1546237)),index9(IND,1546237));
T_1546327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1546327),index9(IND,1546327)),index9(IND,1546327));
T_1562347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1562347),index9(IND,1562347)),index9(IND,1562347));
T_1562437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1562437),index9(IND,1562437)),index9(IND,1562437));
T_1563247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1563247),index9(IND,1563247)),index9(IND,1563247));
T_1563427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1563427),index9(IND,1563427)),index9(IND,1563427));
T_1564237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1564237),index9(IND,1564237)),index9(IND,1564237));
T_1564327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1564327),index9(IND,1564327)),index9(IND,1564327));
T_1623457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1623457),index9(IND,1623457)),index9(IND,1623457));
T_1623547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1623547),index9(IND,1623547)),index9(IND,1623547));
T_1624357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1624357),index9(IND,1624357)),index9(IND,1624357));
T_1624537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1624537),index9(IND,1624537)),index9(IND,1624537));
T_1625347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1625347),index9(IND,1625347)),index9(IND,1625347));
T_1625437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1625437),index9(IND,1625437)),index9(IND,1625437));
T_1632457_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1632457),index9(IND,1632457)),index9(IND,1632457));
T_1632547_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1632547),index9(IND,1632547)),index9(IND,1632547));
T_1634257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1634257),index9(IND,1634257)),index9(IND,1634257));
T_1634527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1634527),index9(IND,1634527)),index9(IND,1634527));
T_1635247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1635247),index9(IND,1635247)),index9(IND,1635247));
T_1635427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1635427),index9(IND,1635427)),index9(IND,1635427));
T_1642357_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1642357),index9(IND,1642357)),index9(IND,1642357));
T_1642537_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1642537),index9(IND,1642537)),index9(IND,1642537));
T_1643257_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1643257),index9(IND,1643257)),index9(IND,1643257));
T_1643527_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1643527),index9(IND,1643527)),index9(IND,1643527));
T_1645237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1645237),index9(IND,1645237)),index9(IND,1645237));
T_1645327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1645327),index9(IND,1645327)),index9(IND,1645327));
T_1652347_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1652347),index9(IND,1652347)),index9(IND,1652347));
T_1652437_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1652437),index9(IND,1652437)),index9(IND,1652437));
T_1653247_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1653247),index9(IND,1653247)),index9(IND,1653247));
T_1653427_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1653427),index9(IND,1653427)),index9(IND,1653427));
T_1654237_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1654237),index9(IND,1654237)),index9(IND,1654237));
T_1654327_up=SM->add(new CTree_with_prefactor(q5gq(PRO,1654327),index9(IND,1654327)),index9(IND,1654327));
}
#endif


if(tree_color==0){ 
//SM->add_tree(cached_cross_term_md(T_1234567_up,T_1234567_up,one,4*std::pow(-1 + std::pow(Nc,2),5)/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1234657_up,T_1234567_up,one,4*-(std::pow(-1 + std::pow(Nc,2),4)/std::pow(Nc,4)))); 
SM->add_tree(cached_cross_term_md(T_1235467_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 3*std::pow(Nc,2) + 3*std::pow(Nc,4) - std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1235647_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1236457_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1236547_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,2) - std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1243567_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 3*std::pow(Nc,2) + 3*std::pow(Nc,4) - std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1243657_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1245367_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1245637_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1246357_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1246537_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1253467_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1253647_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1254367_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,2) - std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1254637_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1256347_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1256437_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1263457_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1263547_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1264357_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1264537_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1265347_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1265437_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - 3*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1324567_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 3*std::pow(Nc,2) + 3*std::pow(Nc,4) - std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1324657_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1325467_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1325647_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1326457_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1326547_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1342567_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1342657_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1345267_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1345627_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1346257_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1346527_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1352467_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1352647_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1354267_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1354627_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1356247_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1356427_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1362457_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1362547_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1364257_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1364527_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1365247_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1365427_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1423567_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1423657_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1425367_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1425637_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1426357_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1426537_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1432567_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,2) - std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1432657_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1435267_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1435627_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1436257_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1436527_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1452367_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - 2*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1452637_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1453267_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1453627_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1456237_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1456327_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1462357_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1462537_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1463257_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1463527_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1465237_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1465327_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1523467_up,T_1234567_up,one,4*((1 - std::pow(Nc,2))*(-1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1523647_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1524367_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1524637_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1526347_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1526437_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1532467_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1532647_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1534267_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1534627_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1536247_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1536427_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1542367_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1542637_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1543267_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - 3*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1543627_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1546237_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1546327_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1562347_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1562437_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1563247_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1563427_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - 4*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1564237_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - 4*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1564327_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 5*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1623457_up,T_1234567_up,one,4*(-1 + std::pow(Nc,2))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1623547_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1624357_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1624537_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1625347_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1625437_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1632457_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1632547_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1634257_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1634527_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1635247_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1635427_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1642357_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 2*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1642537_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1643257_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1643527_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1645237_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - 4*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1645327_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 5*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1652347_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 3*std::pow(Nc,2) - 3*std::pow(Nc,4) + std::pow(Nc,6)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1652437_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1653247_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 4*std::pow(Nc,2) - std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1653427_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 5*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1654237_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 5*std::pow(Nc,2) - 2*std::pow(Nc,4)))/std::pow(Nc,4))); 
SM->add_tree(cached_cross_term_md(T_1654327_up,T_1234567_up,one,4*((-1 + std::pow(Nc,2))*(1 + 6*std::pow(Nc,2) + std::pow(Nc,4)))/std::pow(Nc,4))); 
}

};};};};};};

return SM;
};


#if Zboson==1
Virtual_SME* vsme_2q5g2l(std::vector<int> indext,int ns,int nf,int nc, bool up_down_quark, int photonZW,int color, int tree_color,QCDorder lo_or_nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
#endif

#if Zboson==0
Virtual_SME_with_prefactor* vsme_2q5g2l(std::vector<int> indext,int ns,int nf,int nc, bool up_down_quark, int photonZW,int color, int tree_color,QCDorder lo_or_nlo){
	Virtual_SME_with_prefactor* VSM=new Virtual_SME_with_prefactor() ;
#endif



	vector<int> ind,ind_b, ind98, ind98_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);
	ind.push_back(indext[6]);
	ind.push_back(indext[7]);
	ind.push_back(indext[8]);

	ind98.push_back(indext[0]);
	ind98.push_back(indext[1]);
	ind98.push_back(indext[2]);
	ind98.push_back(indext[3]);
	ind98.push_back(indext[4]);
	ind98.push_back(indext[5]);
	ind98.push_back(indext[6]);
	ind98.push_back(indext[8]);
	ind98.push_back(indext[7]);

	ind_b.push_back(indext[6]);
	ind_b.push_back(indext[5]);
	ind_b.push_back(indext[4]);
	ind_b.push_back(indext[3]);
	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[7]);
	ind_b.push_back(indext[8]);

	ind98_b.push_back(indext[6]);
	ind98_b.push_back(indext[5]);
	ind98_b.push_back(indext[4]);
	ind98_b.push_back(indext[3]);
	ind98_b.push_back(indext[2]);
	ind98_b.push_back(indext[1]);
	ind98_b.push_back(indext[0]);
	ind98_b.push_back(indext[8]);
	ind98_b.push_back(indext[7]);


// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b,_ph_type98, _ph_type98_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

	_ph_type98.push_back(qp);
	_ph_type98.push_back(lp);
	_ph_type98_b.push_back(qm);
	_ph_type98_b.push_back(lp);



//report setup:
if(Zboson==0){
	cout<<"W+5j: using Virtual_SME with prefactor"<<endl;
}

if(BH::settings::BH_interface_settings::s_w5jet_number_of_threads==0) cout<<"threads tune: many"<<endl;
else cout<<"threads tune: "<<BH::settings::BH_interface_settings::s_w5jet_number_of_threads<<endl;


//only for Ws
#if Zboson==0
//--------------------- propagator correction
	prop_hel_fn* _prop_hel_fn= new prop_hel_fn(up_down_quark,photonZW,0,indext[7],indext[8],_ph_type_b);
	VSM->add_amplitude_prefactor(_prop_hel_fn);
//-------------------------------------------
#endif





int possible_threads;
switch(BH::settings::BH_interface_settings::s_w5jet_number_of_threads){
	case 0: possible_threads=32;break;
	case 3: possible_threads=3;break;
	case 4: possible_threads=4;break;
	case 6: possible_threads=6;break;
}

// for 6 threads and uneven workload
if(photonZW!=3) possible_threads*=4;

std::vector<Squared_ME*> SM;
for(int i=0;i<possible_threads;i++){
		SM.push_back(new Squared_ME(lo_or_nlo));
}

#if Zboson==0
std::vector<int> thread(32);

switch(BH::settings::BH_interface_settings::s_w5jet_number_of_threads){
	case 0: for(int i=0;i<32;i++) thread[i]=i;break;
	case 3: {
		thread[0]=0; 
		for(int i=1;i<6;i++) thread[i]=1;	
		for(int i=6;i<16;i++) thread[i]=0;	
		thread[15]=2;
		for(int i=16;i<22;i++) thread[i]=1;
		for(int i=22;i<32;i++) thread[i]=2;
		} break;	
	case 4: {
		thread[0]=0; 
		for(int i=1;i<6;i++) thread[i]=0;	
		for(int i=6;i<16;i++) thread[i]=1;	
		thread[15]=2;
		for(int i=16;i<22;i++) thread[i]=2;
		for(int i=22;i<32;i++) thread[i]=3;
		} break;	

	case 6: {
		thread[0]=0; 
		for(int i=1;i<6;i++) thread[i]=1;	
		for(int i=6;i<16;i++) thread[i]=2;	
		thread[15]=3;
		for(int i=16;i<22;i++) thread[i]=4;
		for(int i=22;i<32;i++) thread[i]=5;
		} break;	
}
#endif


#if Zboson==1
std::vector<int> thread(128);

switch(BH::settings::BH_interface_settings::s_w5jet_number_of_threads){
	case 0: for(int i=0;i<32;i++) thread[i]=i; break;
	case 3: {
		thread[0]=0; 
		for(int i=1;i<6;i++) thread[i]=1;	
		for(int i=6;i<16;i++) thread[i]=0;	
		thread[15]=2;
		for(int i=16;i<22;i++) thread[i]=1;
		for(int i=22;i<32;i++) thread[i]=2;
		} break;	
	case 4: {
		thread[0]=0; 
		for(int i=1;i<6;i++) thread[i]=0;	
		for(int i=6;i<16;i++) thread[i]=1;	
		thread[15]=2;
		for(int i=16;i<22;i++) thread[i]=2;
		for(int i=22;i<32;i++) thread[i]=3;
		} break;	

	case 6: {
		thread[0]=0; 
		for(int i=1;i<6;i++) thread[i]=1;	
		for(int i=6;i<16;i++) thread[i]=2;	
		thread[15]=3;
		for(int i=16;i<22;i++) thread[i]=4;
		for(int i=22;i<32;i++) thread[i]=5;
		} break;	
}
for(int i=32;i<128;i++) thread[i]=thread[i%32];
#endif

//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLESMNTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.
int pos(0);
#if Zboson==1	 
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);

	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);

	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);

	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo,SM[thread[pos++]]);


#endif
//----

	 
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;

	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;

	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;

	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo,SM[thread[pos]]);pos++;

#if Zboson==1	 
//----

	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);

	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);


	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);

	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,m,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,p,qbm,lm,lbp),ind98,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98,lo_or_nlo,SM[thread[pos++]]);


//---

	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);

	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);

	 
	 A_loop_2q_5g_2l_M2(process(qp,m,m,m,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);

	 A_loop_2q_5g_2l_M2(process(qp,p,p,p,m,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,p,m,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,p,m,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,p,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,p,m,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,p,m,m,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,p,m,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);
	 A_loop_2q_5g_2l_M2(process(qp,m,m,p,p,p,qbm,lm,lbp),ind98_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type98_b,lo_or_nlo,SM[thread[pos++]]);


#endif



for(int i=0;i<possible_threads;i++) VSM->add(SM[i]);
return VSM;
}

BH_Ampl_2q5g2l::BH_Ampl_2q5g2l(bool up_down_quark, int photonZW, int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q5g2l(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				up_down_quark,
				photonZW*settings::BH_interface_settings::s_photon_only,
				color,
				tree_color,
				lo_or_nlo),
			bhi,  // parent
			9,    // NbrExtParticles
			5,	  // NbrPowersOfAlphaS
			2,    // NbrPowersOfAlphaQED
			10,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
{}



}
