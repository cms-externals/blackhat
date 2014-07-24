#include <iostream>

using namespace std;

extern "C" {
extern void OLP_EvalSubProcess(int* Label,double* Momenta,double *mu,double *parameters,double *result);
extern void OLP_Start(const char* filename,int *status);
}

int main(){
	int success;
	OLP_Start("contract_file.lh",&success);
	if (success!=1) { return 1; }
	double result[4];
	double parameters[]={1.0,1.0};
	double momenta[]={
//   E      x      y      z        m
   500.0, 0.0,    0.0,  500.0,    0.0,
   500.0, 0.0,    0.0, -500.0,    0.0,
   500.0, 0.0, -500.0,    0.0,    0.0,
   500.0, 0.0,  500.0,    0.0,    0.0,
   };

  double mu=500; //
  int label=1;
  OLP_EvalSubProcess(&label,momenta,&mu,parameters,result);

  cout << "1/e^2 : " << result[0] << endl;
  cout << "1/e   : " << result[1] << endl;
  cout << "finite: " << result[2] << endl;

}



