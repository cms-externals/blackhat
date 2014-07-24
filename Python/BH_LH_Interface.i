%module BHLH

%include "std_string.i"
%include "std_vector.i"
%include "carrays.i"
%include "typemaps.i"
%include "std_vector.i"

namespace std {
   %template(vector_d) vector<double>;
}; 

%{
#include "Interface/BH_LH_interface.h"
namespace BH {
namespace LesHouches {
std::vector<double> EvalSubprocessForPython(int Label,std::vector<double>& Momenta,double mu,double alpha_s,double alpha_ew);
void SignContract(char* filename1,char* filename2);

}
}


 %}




namespace BH {
namespace LesHouches {

int Init(const char* filename);
void SignContract(char* filename1,char* filename2);
std::vector<double> EvalSubprocessForPython(int Label,std::vector<double>& Momenta,double mu,double alpha_s,double alpha_ew);
}
}


%pythoncode %{
def EvalSubprocess(i,momenta,mu,als,alew):
	return EvalSubprocessForPython(i,vector_d(momenta),mu,als,alew);
%}
