/* polylog_T.cpp */
/* adapted from: */
/*  David A. Kosower, November 19, 2007  */

/*  Implementation of real and complex polylogarithms.  The approach is
    based on the one described in 't Hooft and Veltman [NPB153:365 (1979)],
    and also documented in hep-ph/0502165v2 [Hameren, Vollinga and Weinzierl].

*/

//#include "c++standard.h"
#include <cassert>
#include "qd/dd_real.h"
#include "qd_suppl.h"
#include "str_constants.h"
#include "BH_debug.h"
using namespace qd;
namespace BH {
template<class T> inline T sq(const T& x) {return(x*x);}


//#define Li2TermLimit_HP 28  // 7 for double
#define Li2LesserTermThreshold 0.2
#define CLi2LesserTermThreshold 0.05
#define CLi2LesserTermLimit_HP 20
#define Cl2TermLimit_HP 29 // 14 for double
#define Cl2LesserTermLimit_HP 18 // for x < Pi/3

template <class T> int CLi2TermLimit();

template <> int CLi2TermLimit<R>(){return 8;};
template <> int CLi2TermLimit<RHP>(){return 19;};
template <> int CLi2TermLimit<RVHP>(){return 37;};
#if BH_USE_GMP
template <> int CLi2TermLimit<RGMP>(){ return RGMP::get_current_nbr_digits()/2 +3;}
#endif

template <class T> int CLi2LesserTermLimit();

template <> int CLi2LesserTermLimit<R>(){return 8;}; // probably less
template <> int CLi2LesserTermLimit<RHP>(){return 20;}; // probably less
template <> int CLi2LesserTermLimit<RVHP>(){return 20;};
#if BH_USE_GMP
template <> int CLi2LesserTermLimit<RGMP>(){ return RGMP::get_current_nbr_digits()/2 +3;}
#endif

template <class T> int Li2TermLimit();
template <> int  Li2TermLimit<R>(){ return 7;}
template <> int  Li2TermLimit<RHP>(){ return 15;}
template <> int  Li2TermLimit<RVHP>(){ return 29;}
#if BH_USE_GMP
template <> int  Li2TermLimit<RGMP>(){ return RGMP::get_current_nbr_digits()/2 +3;}
#endif

template <class T> int Li2LesserTermLimit();
template <> int  Li2LesserTermLimit<R>(){ return 4;}
template <> int  Li2LesserTermLimit<RHP>(){ return 9;}
template <> int  Li2LesserTermLimit<RVHP>(){ return 19;}
#if BH_USE_GMP
template <> int  Li2LesserTermLimit<RGMP>(){ return RGMP::get_current_nbr_digits()/3 ;}
#endif



const char* pi_str="3.14159265358979323846264338327950288419716939937510582097494459230781\
6406286208998628034825342117067982148086513282306647093844609550582231\
7253594081284811174502841027019385211055596446229489549303819644288109\
7566593344612847564823378678316527120190914564856692346034861045432664\
8213393607260249141273724587006606315588174881520920962829254091715364\
3678925903600113305305488204665213841469519415116094330572703657595919\
5309218611738193261179310511854807446237996274956735188575272489122793\
8183011949129833673362440656643086021394946395224737190702179860943702\
7705392171762931767523846748184676694051320005681271452635608277857713\
4275778960917363717872146844090122495343014654958537105079227968925892\
3542019956112129021960864034418159813629774771309960518707211349999998\
3729780499510597317328160963185950244594553469083026425223082533446850\
3526193118817101000313783875288658753320838142061717766914730359825349\
0428755468731159562863882353787593751957781857780532171226806613001927\
876611195909216420199";
const char* pisquaredover3_str="3.28986813369645287294483033329205037843789980241359687547111645874001\
4940806401747667257801239517410608008637924674381359257449374010015575\
8702058926617325536634666187355252101905020137442801095936231175897807\
2164655552383968151291175392647127341942019389780417186401610327295775\
6677692088890368119650290501366775262845531758785917612640894439581695\
4681821180416756579098556527780759527166687884090318241636199186908897\
5491759300176177881740222326942138632292368577596309724896718183668975\
1477485678816552057512642869202002715324196409744138080014765327120604\
8045689259260649132194343902855442631902511359972381743863087907048212\
7608814428426793095011603174463316798952486982844866967258097740193301\
1972452606821919347310562274334065382299756806871432321155335266613450\
5473788476833281779072455190801545589496254205040996756866460034331489\
6206057208699337695884334568671945639955876200169331215610755657718945\
5725186323732917658432131638771846483065161292356240376929955525196995\
512121876921210293537";
const char* pisquaredover6_str="1.64493406684822643647241516664602518921894990120679843773555822937000\
7470403200873833628900619758705304004318962337190679628724687005007787\
9351029463308662768317333093677626050952510068721400547968115587948903\
6082327776191984075645587696323563670971009694890208593200805163647887\
8338846044445184059825145250683387631422765879392958806320447219790847\
7340910590208378289549278263890379763583343942045159120818099593454448\
7745879650088088940870111163471069316146184288798154862448359091834487\
5738742839408276028756321434601001357662098204872069040007382663560302\
4022844629630324566097171951427721315951255679986190871931543953524106\
3804407214213396547505801587231658399476243491422433483629048870096650\
5986226303410959673655281137167032691149878403435716160577667633306725\
2736894238416640889536227595400772794748127102520498378433230017165744\
8103028604349668847942167284335972819977938100084665607805377828859472\
7862593161866458829216065819385923241532580646178120188464977762598497\
75606093846060514676";
template <class T> const T* get_Bernoulli();

#if BH_USE_GMP
template <> const RGMP* get_Bernoulli<RGMP>(){
	static int last_precision=0;
	static RGMP bern[250];
  if (last_precision < RGMP::get_current_precision()){
	  last_precision=RGMP::get_current_precision();
	  for (int i=0;i<250;i++){
		  bern[i]=RGMP(Bernoulli_mantisse[i])*RGMP(Bernoulli_exponent[i]);
	  }
  }

	return bern;
};
#endif

template <class T> const T* get_Bernoulli(){
	static bool isDone=false;
	static T bern[250];
  if (!isDone){
	  isDone=true;
	  for (int i=0;i<250;i++){
		  bern[i]=T(Bernoulli_mantisse[i])*T(Bernoulli_exponent[i]);
	  }
  }

	return bern;
};

template <class T> T ReLi2_tpl(const T& x_arg)
{

	unsigned int old_cw;
   fpu_fix_start(&old_cw);

   BH_DEBUG_MESSAGE2("ReLi2 called with argument: ",x_arg);

   T result;

   T added = T("0.0"), factor = T("1.0");
   // Improper rounding at 10^-20 or so seems to occur if the array is
   // outside the function
   // B_{2 index}
   const T B0 = T("1.0");
   const T B1over2 = T("-0.25");
   const T PI = T(pi_str);
   const T PiSquaredOver3 = T(pisquaredover3_str);
   const T PiSquaredOver6 = T(pisquaredover6_str);
   const T* Bernoulli2 = get_Bernoulli<T>();

//   assert (sizeof(Bernoulli2)/sizeof(Bernoulli2[0]) >= Li2TermLimit<T>());

	  T x=x_arg;
   // Map the argument into the range [0,1/2]
 if (x >= T(2))
   // Re[PolyLog[2, x]] ->  Pi^2/3 - Log[x]^2/2 - PolyLog[2, x^(-1)]}
    {added = PiSquaredOver3 - T("0.5")*sq(log(x));
     factor = T("-1.0");
     x = T("1.0")/x;}
 else if (x > T(1))
   // Re[PolyLog[2, x]] ->
   //           Pi^2/6 - Log[x-1]*Log[x] + Log[x]^2/2 + PolyLog[2, (x-1)/x]
    {T lnx = log(x);
     added = PiSquaredOver6 +(T("0.5")*lnx - log(x-T("1.0")))*lnx;
     factor = T("1.0");
     x = (x-T("1.0"))/x;}
 else if (x > T("0.5"))
   // PolyLog[2,x] -> -PolyLog[2,1-x] + Pi^2/6 - Log[x] Log[1-x]
    {added = PiSquaredOver6 - log(x)*log(1-x);
     factor = T("-1.0");
     x = T("1.0")-x;}
 else if (x > 0) {}
 else if (x >= -1)
   // PolyLog[2, y] ->  -Log[1 - y]^2/2 - PolyLog[2, y/(-1 + y)]
    {added = T("-0.5")*sq(log(T("1.0")-x));
     factor = T("-1.0");
     x = x/(x-T("1.0"));}
 else
   // PolyLog[2, y] ->  -Pi^2/6 + Log[1-y]^2/2 - Log[1-y]*Log[-y]
   //                   + PolyLog[2, 1/(1-y)]
    {T ln1x = log(T("1.0")-x);
     added = -PiSquaredOver6 + (T("0.5")*ln1x-log(-x))*ln1x;
     factor = T("1.0");
     x = T("1.0")/(T("1.0")-x);}
 // Compute PolyLog[2,x], x now in [0,1/2], by Bernoulli series
 T z = -log(T("1.0")-x);
 T li2 = (B0+B1over2*z)*z;
 T term = z, zsq = z*z;

 int limit = Li2TermLimit<T>();
 if (x < Li2LesserTermThreshold) limit = Li2LesserTermLimit<T>() ; //4 in double

 if (limit >= NbrBernoulliTerms){
	 std::cerr<< "Maximum precision reached in ReLi2_tpl. The result might not have the accuracy requested." << std::endl;
	 limit=NbrBernoulliTerms-1;
 };
 for (int j = 1;  j <= limit;  j += 1)
    {
	 term *= zsq/T((2*j+1)*2*j);
     li2 += Bernoulli2[j]*term;
}

result = factor*li2+added;
fpu_fix_end(&old_cw);

// return(factor*li2+added);
 return(result);
}



// Complex version of the above
template <class T> std::complex<T> Li2_tpl(const std::complex<T>& z_arg)
  /* Transformations may not be for optimal regions yet; following
     't Hooft & Veltman, put the argument inside the circle |z|=1,
     and make the real part < 0.5 */
{unsigned int old_cw;
 fpu_fix_start(&old_cw);
 typedef std::complex<T> CT;

 CT added(0,0);
 T factor(1);
 const T B0 = T("1.0");
 const T B1over2 = T("-0.25");
 const T PI = T(pi_str);
 const T PiSquaredOver3 = T(pisquaredover3_str);
 const T PiSquaredOver6 = T(pisquaredover6_str);
 const T* Bernoulli2 = get_Bernoulli<T>();

BH_DEBUG_MESSAGE2("Li2 called with argument: ",z_arg);
   CT z=z_arg;


 if (Re(z*conj(z)) > 1)
    {added = -PiSquaredOver6-sq(log(-z))/T(2);
     factor = T("-1");
     z = T("1.")/z;}
 if (Re(z) > T("0.5"))
    {added += factor*(PiSquaredOver6-log(z)*log(T("1.")-z));
     factor *= T("-1");
     z = T("1.")-z;}

 // Compute PolyLog[2,z], z now in unit disc with Re z < 1/2, by Bernoulli series
 CT w = -log(T("1.")-z);

 CT li2 = (B0+B1over2*w)*w;
 CT term = w, wsq = w*w;
 int limit = CLi2TermLimit<T>();
 if (Re(w*conj(w)) < CLi2LesserTermThreshold)
    limit = CLi2LesserTermLimit<T>();
 if (limit >= NbrBernoulliTerms){
 	 std::cerr<< "Maximum precision reached in Li2_tpl. The result might not have the accuracy requested." << std::endl;
 	 limit=NbrBernoulliTerms-1;
  };
 for (int j = 1;  j <= limit;  j += 1)
    {term *= wsq/T((2*j+1)*2*j);
    li2 += Bernoulli2[j]*term;}

 CT result = (factor*li2+added);
 fpu_fix_end(&old_cw);
 return(result);
}

// Clausen function, from hep-ph/0502165v2 [Hameren, Vollinga and Weinzierl]

template <class T> T Cl2(T x)
{unsigned int old_cw;
 fpu_fix_start(&old_cw);
 // Having these definitions (even the first!) outside the routine (therefore
 // outside the fpu_fix_start/fpu_fix_end pair) seems to cause errors at the
 // 10^-19 level...

 BH_DEBUG_MESSAGE2("Cl2 called with argument: ",x);

 const T PI_HP = T(pi_str);
 const T TwoPi = T(2)*PI_HP;
 const T TwoPiOver3 = T(2)*PI_HP/T(3);
 const T PiOver3 = PI_HP/T(3);

 const T B0_HP = T(1);
 const T B1over2_HP = T(-1)/T(4);
 const T Bernoulli2[]= {
 		T(Bernoulli_mantisse[0])*T(Bernoulli_exponent[0]),
 		T(Bernoulli_mantisse[1])*T(Bernoulli_exponent[1]),
 		T(Bernoulli_mantisse[2])*T(Bernoulli_exponent[2]),
 		T(Bernoulli_mantisse[3])*T(Bernoulli_exponent[3]),
 		T(Bernoulli_mantisse[4])*T(Bernoulli_exponent[4]),
 		T(Bernoulli_mantisse[5])*T(Bernoulli_exponent[5]),
 		T(Bernoulli_mantisse[6])*T(Bernoulli_exponent[6]),
 		T(Bernoulli_mantisse[7])*T(Bernoulli_exponent[7]),
 		T(Bernoulli_mantisse[8])*T(Bernoulli_exponent[8]),
 		T(Bernoulli_mantisse[9])*T(Bernoulli_exponent[9]),
 		T(Bernoulli_mantisse[10])*T(Bernoulli_exponent[10]),
 		T(Bernoulli_mantisse[11])*T(Bernoulli_exponent[11]),
 		T(Bernoulli_mantisse[12])*T(Bernoulli_exponent[12]),
 		T(Bernoulli_mantisse[13])*T(Bernoulli_exponent[13]),
 		T(Bernoulli_mantisse[14])*T(Bernoulli_exponent[14]),
 		T(Bernoulli_mantisse[15])*T(Bernoulli_exponent[15]),
 		T(Bernoulli_mantisse[16])*T(Bernoulli_exponent[16]),
 		T(Bernoulli_mantisse[17])*T(Bernoulli_exponent[17]),
 		T(Bernoulli_mantisse[18])*T(Bernoulli_exponent[18]),
 		T(Bernoulli_mantisse[19])*T(Bernoulli_exponent[19]),
 		T(Bernoulli_mantisse[20])*T(Bernoulli_exponent[20]),
 		T(Bernoulli_mantisse[21])*T(Bernoulli_exponent[21]),
 		T(Bernoulli_mantisse[22])*T(Bernoulli_exponent[22]),
 		T(Bernoulli_mantisse[23])*T(Bernoulli_exponent[23]),
 		T(Bernoulli_mantisse[24])*T(Bernoulli_exponent[24]),
 		T(Bernoulli_mantisse[25])*T(Bernoulli_exponent[25]),
 		T(Bernoulli_mantisse[26])*T(Bernoulli_exponent[26]),
 		T(Bernoulli_mantisse[27])*T(Bernoulli_exponent[27]),
 		T(Bernoulli_mantisse[28])*T(Bernoulli_exponent[28]),
 		T(Bernoulli_mantisse[29])*T(Bernoulli_exponent[29]),
 		T(Bernoulli_mantisse[30])*T(Bernoulli_exponent[30]),
 		T(Bernoulli_mantisse[31])*T(Bernoulli_exponent[31]),
 		T(Bernoulli_mantisse[32])*T(Bernoulli_exponent[32]),
 		T(Bernoulli_mantisse[33])*T(Bernoulli_exponent[33]),
 		T(Bernoulli_mantisse[34])*T(Bernoulli_exponent[34]),
 		T(Bernoulli_mantisse[35])*T(Bernoulli_exponent[35]),
 		T(Bernoulli_mantisse[36])*T(Bernoulli_exponent[36]),
 		T(Bernoulli_mantisse[37])*T(Bernoulli_exponent[37]),
 		T(Bernoulli_mantisse[38])*T(Bernoulli_exponent[38]),
 		T(Bernoulli_mantisse[39])*T(Bernoulli_exponent[39]),
 		T(Bernoulli_mantisse[40])*T(Bernoulli_exponent[40]),
 		T(Bernoulli_mantisse[41])*T(Bernoulli_exponent[41]),
 		T(Bernoulli_mantisse[42])*T(Bernoulli_exponent[42]),
 		T(Bernoulli_mantisse[43])*T(Bernoulli_exponent[43]),
 		T(Bernoulli_mantisse[44])*T(Bernoulli_exponent[44]),
 		T(Bernoulli_mantisse[45])*T(Bernoulli_exponent[45]),
 		T(Bernoulli_mantisse[46])*T(Bernoulli_exponent[46]),
 		T(Bernoulli_mantisse[47])*T(Bernoulli_exponent[47]),
 		T(Bernoulli_mantisse[48])*T(Bernoulli_exponent[48]),
 		T(Bernoulli_mantisse[49])*T(Bernoulli_exponent[49]),
 		T(Bernoulli_mantisse[50])*T(Bernoulli_exponent[50]),
 		T(Bernoulli_mantisse[51])*T(Bernoulli_exponent[51]),
 		T(Bernoulli_mantisse[52])*T(Bernoulli_exponent[52]),
 		T(Bernoulli_mantisse[53])*T(Bernoulli_exponent[53]),
 		T(Bernoulli_mantisse[54])*T(Bernoulli_exponent[54]),
 		T(Bernoulli_mantisse[55])*T(Bernoulli_exponent[55]),
 		T(Bernoulli_mantisse[56])*T(Bernoulli_exponent[56]),
 		T(Bernoulli_mantisse[57])*T(Bernoulli_exponent[57]),
 		T(Bernoulli_mantisse[58])*T(Bernoulli_exponent[58]),
 		T(Bernoulli_mantisse[59])*T(Bernoulli_exponent[59]),
 		T(Bernoulli_mantisse[60])*T(Bernoulli_exponent[60]),
 		T(Bernoulli_mantisse[61])*T(Bernoulli_exponent[61]),
 		T(Bernoulli_mantisse[62])*T(Bernoulli_exponent[62]),
 		T(Bernoulli_mantisse[63])*T(Bernoulli_exponent[63]),
 		T(Bernoulli_mantisse[64])*T(Bernoulli_exponent[64]),
 		T(Bernoulli_mantisse[65])*T(Bernoulli_exponent[65]),
 		T(Bernoulli_mantisse[66])*T(Bernoulli_exponent[66]),
 		T(Bernoulli_mantisse[67])*T(Bernoulli_exponent[67]),
 		T(Bernoulli_mantisse[68])*T(Bernoulli_exponent[68]),
 		T(Bernoulli_mantisse[69])*T(Bernoulli_exponent[69]),
 		T(Bernoulli_mantisse[70])*T(Bernoulli_exponent[70]),
 		T(Bernoulli_mantisse[71])*T(Bernoulli_exponent[71]),
 		T(Bernoulli_mantisse[72])*T(Bernoulli_exponent[72]),
 		T(Bernoulli_mantisse[73])*T(Bernoulli_exponent[73]),
 		T(Bernoulli_mantisse[74])*T(Bernoulli_exponent[74]),
 		T(Bernoulli_mantisse[75])*T(Bernoulli_exponent[75]),
 		T(Bernoulli_mantisse[76])*T(Bernoulli_exponent[76]),
 		T(Bernoulli_mantisse[77])*T(Bernoulli_exponent[77]),
 		T(Bernoulli_mantisse[78])*T(Bernoulli_exponent[78]),
 		T(Bernoulli_mantisse[79])*T(Bernoulli_exponent[79]),
 		T(Bernoulli_mantisse[80])*T(Bernoulli_exponent[80]),
 		T(Bernoulli_mantisse[81])*T(Bernoulli_exponent[81]),
 		T(Bernoulli_mantisse[82])*T(Bernoulli_exponent[82]),
 		T(Bernoulli_mantisse[83])*T(Bernoulli_exponent[83]),
 		T(Bernoulli_mantisse[84])*T(Bernoulli_exponent[84]),
 		T(Bernoulli_mantisse[85])*T(Bernoulli_exponent[85]),
 		T(Bernoulli_mantisse[86])*T(Bernoulli_exponent[86]),
 		T(Bernoulli_mantisse[87])*T(Bernoulli_exponent[87]),
 		T(Bernoulli_mantisse[88])*T(Bernoulli_exponent[88]),
 		T(Bernoulli_mantisse[89])*T(Bernoulli_exponent[89]),
 		T(Bernoulli_mantisse[90])*T(Bernoulli_exponent[90]),
 		T(Bernoulli_mantisse[91])*T(Bernoulli_exponent[91]),
 		T(Bernoulli_mantisse[92])*T(Bernoulli_exponent[92]),
 		T(Bernoulli_mantisse[93])*T(Bernoulli_exponent[93]),
 		T(Bernoulli_mantisse[94])*T(Bernoulli_exponent[94]),
 		T(Bernoulli_mantisse[95])*T(Bernoulli_exponent[95]),
 		T(Bernoulli_mantisse[96])*T(Bernoulli_exponent[96]),
 		T(Bernoulli_mantisse[97])*T(Bernoulli_exponent[97]),
 		T(Bernoulli_mantisse[98])*T(Bernoulli_exponent[98]),
 		T(Bernoulli_mantisse[99])*T(Bernoulli_exponent[99]),
 		T(Bernoulli_mantisse[100])*T(Bernoulli_exponent[100])
 };
					;


   assert (sizeof(Bernoulli2)/sizeof(Bernoulli2[0]) >= Cl2TermLimit_HP);

 T factor = T("1.0");
 T result;
 if (x < 0) {x = -x;  factor = T("-1.0");}

 while (x > TwoPi) x -= TwoPi;
 // shift to range [0,2 Pi/3] (larger upper limit prevents too many recursions)
 if (x > TwoPiOver3)
   {result = (T(2)*factor*(Cl2(T("0.5")*x)-Cl2(PI_HP-T("0.5")*x)));
    fpu_fix_end(&old_cw);
    return result;}
 // Compute by Bernoulli series
 T cl2 = x*(T("1.0")-log(x));
 T term = -x, mxsq = -x*x;
 int limit = Cl2TermLimit_HP;
 if (x < PiOver3) limit = Cl2LesserTermLimit_HP; //8 for double
 for (int j = 1;  j <= limit;  j += 1)
    {term *= mxsq/T((2*j+1)*2*j);
     cl2 += Bernoulli2*term/T(2*j);}

 result = (factor*cl2);
 fpu_fix_end(&old_cw);
 return result;
}
} /* BH */
