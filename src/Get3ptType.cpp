#include "partitions_fwd.h"


// need this file because it is needed for the public version and I could
// not find a place that was appropriate

namespace BH {

//This function returns 1 if the three-point amplitude is mostly plus or -1 if the amplitude is mostely minus
int Get3ptType_new(corner_type corner)
{
	if(corner==a_type){
		return -1;
	}
	else if(corner==b_type){
		return 1;
	}

	return 0;
}




}
