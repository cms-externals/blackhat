/*!\file mom_conf.cpp
\brief Implementation of the classe mom_conf
*/

#include "mom_conf.h"
#include "BH_error.h"
#include <cstdarg>
#include <stdio.h>
#include <string.h>
#include "mom_conf.hpp"


using std::string;
using std::vector;
using std::cout;using std::endl;


namespace BH {


#define Lw1(c1,c2) c1
#define Lw2(c1,c2) c2






void mathprint(momentum_configuration<RVHP>& momc)
{
    int n = (int) momc.n();
    cout << "{";
    for(int i = 1; i < n; i++)
    {    cout << "{" << real(momc.p(i).E()) << " + I " <<
imag(momc.p(i).E()) << ",";
        cout << real(momc.p(i).X()) << " + I " << imag(momc.p(i).X()) <<
",";
        cout << real(momc.p(i).Y()) << " + I " << imag(momc.p(i).Y()) <<
",";
        cout << real(momc.p(i).Z()) << " + I " << imag(momc.p(i).Z());
        cout << "},";
    }
    cout << "{" << real(momc.p(n).E()) << " + I " << imag(momc.p(n).E())
<< ",";
            cout << real(momc.p(n).X()) << " + I " <<
imag(momc.p(n).X()) << ",";
            cout << real(momc.p(n).Y()) << " + I " <<
imag(momc.p(n).Y()) << ",";
            cout << real(momc.p(n).Z()) << " + I " <<
imag(momc.p(n).Z());
            cout << "}";
    cout << "}" << endl;
    return;
}



// EXPLICIT INSTANCIATION

template class momentum_configuration<R>;
template class momentum_configuration<RHP>;
template class momentum_configuration<RVHP>;
template class sub_momentum_configuration<R>;
template class sub_momentum_configuration<RHP>;
template class sub_momentum_configuration<RVHP>;
template class mom_conf_reader<R>;
template class mom_conf_reader<RHP>;
template class mom_conf_reader<RVHP>;

template <> template <> momentum_configuration<RHP> momentum_configuration<R>::extend<RHP>(const vector<int>& indices);
template <> template <> momentum_configuration<RVHP> momentum_configuration<R>::extend<RVHP>(const vector<int>& indices);
template <> template <> momentum_configuration<RVHP> momentum_configuration<RHP>::extend<RVHP>(const vector<int>& indices);


template std::ostream& operator<<(std::ostream& s, const momentum_configuration<R>& mc);
template std::ostream& operator<<(std::ostream& s, const momentum_configuration<RHP>& mc);
template std::ostream& operator<<(std::ostream& s, const momentum_configuration<RVHP>& mc);

long int  momentum_configuration_base::mom_conf_next_ID=0;
//template<> long int  momentum_configuration_base<BH::RHP>::mom_conf_next_ID=0;
//template<> long int  momentum_configuration_base<BH::RVHP>::mom_conf_next_ID=0;

#if _WITH_CACHING
#if _CACHE_STATISTICS
std::map<std::string,momentum_configuration<R>::hash_stat> dummy;
std::map<std::string,momentum_configuration<RHP>::hash_stat> dummyHP;
std::map<std::string,momentum_configuration<RVHP>::hash_stat> dummyVHP;
template<> std::map<std::string,momentum_configuration<R>::hash_stat> momentum_configuration<R>::hash_stat_table=dummy;
template<>  std::map<std::string,momentum_configuration<RHP>::hash_stat> momentum_configuration<RHP>::hash_stat_table=dummyHP;
template<>  std::map<std::string,momentum_configuration<RVHP>::hash_stat> momentum_configuration<RVHP>::hash_stat_table=dummyVHP;
#endif
#endif

//template class momentum_configuration<R>* momentum_configuration<R>::default_momentum_configuration<R>;
//template class momentum_configuration<RHP>* momentum_configuration<RHP>::default_momentum_configuration<RHp>;

}
