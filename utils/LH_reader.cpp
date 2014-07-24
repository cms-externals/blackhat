/*
 * LH_reader.cpp
 *
 *  Created on: Jun 16, 2009
 *      Author: daniel
 */



#include <iostream>

#include "Interface/BH_LH_interface.h"

using namespace std;
using namespace BH;
using namespace LesHouches;

namespace BH {
namespace LesHouches {
void SignContract(char* filename1,char* filename2);
void PrintHelp(std::ostream&);
}
}

int main(int argc, char* argv[]) {
	   if (argc == 2  ){
		   cerr << argv[1] << endl;
		   string arg(argv[1]);
		   if ( (arg=="--help" || arg=="-h" || arg=="-H" || arg=="help" )) {
			   cerr << "Usage: " << argv[0] << " ORDERFILE CONTRACTFILE" << endl;
			   PrintHelp(cerr);
			   cerr << "\n";
		   return 1;
		   }
	   }
   if (argc != 3){
	   cerr << "Usage: " << argv[0] << " ORDERFILE CONTRACTFILE" << endl;
	   cerr << " or    " << argv[0] << " help    for a description of the options" << endl;
	   cerr << "\n";
	   return 1;
   }

//	cout << "argc = " << argc << endl;
//   for(int i = 0; i < argc; i++)
//      cout << "argv[" << i << "] = " << argv[i] << endl;

	SignContract(argv[1],argv[2]);

   return 0;



}
