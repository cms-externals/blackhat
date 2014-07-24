/*
 * BH_debug.h
 *
 *  Created on: Oct 13, 2009
 *      Author: daniel
 */

#ifndef BH_DEBUG_H_
#define BH_DEBUG_H_



#include "BH_utilities.h"

// the two functions are necessary to make sure __LINE__ is replaced before being concatenated with unique_ID_
#ifdef BH_DEBUG_ON
#define UNIQUE_ID(X)  unique_ID_ ## X
#define UNIQUE_STR_ID(X)  unique_ID_str ## X

#define BH_DEBUG_MESSAGEIMPL(X,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE2(UNIQUE_STR_ID(UNIQUE),X);}
#define BH_DEBUG_MESSAGE(X) BH_DEBUG_MESSAGEIMPL(X,__LINE__)
#define BH_DEBUG_MESSAGEIMPL2(X,Y,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE3(UNIQUE_STR_ID(UNIQUE),X,Y);}
#define BH_DEBUG_MESSAGE2(X,Y) BH_DEBUG_MESSAGEIMPL2(X,Y,__LINE__)
#define BH_DEBUG_MESSAGEIMPL3(X,Y,Z,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE4(UNIQUE_STR_ID(UNIQUE),X,Y,Z);}
#define BH_DEBUG_MESSAGE3(X,Y,Z) BH_DEBUG_MESSAGEIMPL3(X,Y,Z,__LINE__)
#define BH_DEBUG_MESSAGEIMPL4(X,Y,Z,A,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE5(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A);}
#define BH_DEBUG_MESSAGE4(X,Y,Z,A) BH_DEBUG_MESSAGEIMPL4(X,Y,Z,A,__LINE__)
#define BH_DEBUG_MESSAGEIMPL5(X,Y,Z,A,B,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE6(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B);}
#define BH_DEBUG_MESSAGE5(X,Y,Z,A,B) BH_DEBUG_MESSAGEIMPL5(X,Y,Z,A,B,__LINE__)
#define BH_DEBUG_MESSAGEIMPL6(X,Y,Z,A,B,C,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE7(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C);}
#define BH_DEBUG_MESSAGE6(X,Y,Z,A,B,C) BH_DEBUG_MESSAGEIMPL6(X,Y,Z,A,B,C,__LINE__)
#define BH_DEBUG_MESSAGEIMPL7(X,Y,Z,A,B,C,D,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE8(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D);}
#define BH_DEBUG_MESSAGE7(X,Y,Z,A,B,C,D) BH_DEBUG_MESSAGEIMPL7(X,Y,Z,A,B,C,D,__LINE__)
#define BH_DEBUG_MESSAGEIMPL8(X,Y,Z,A,B,C,D,E,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE9(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D,E);}
#define BH_DEBUG_MESSAGE8(X,Y,Z,A,B,C,D,E) BH_DEBUG_MESSAGEIMPL8(X,Y,Z,A,B,C,D,E,__LINE__)
#define BH_DEBUG_MESSAGEIMPL9(X,Y,Z,A,B,C,D,E,F,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE10(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D,E,F);}
#define BH_DEBUG_MESSAGE9(X,Y,Z,A,B,C,D,E,F) BH_DEBUG_MESSAGEIMPL9(X,Y,Z,A,B,C,D,E,F,__LINE__)
#define BH_DEBUG_MESSAGEIMPL10(X,Y,Z,A,B,C,D,E,F,G,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){_MESSAGE11(UNIQUE_STR_ID(UNIQUE),X,Y,Z,A,B,C,D,E,F,G);}
#define BH_DEBUG_MESSAGE10(X,Y,Z,A,B,C,D,E,F,G) BH_DEBUG_MESSAGEIMPL10(X,Y,Z,A,B,C,D,E,F,G,__LINE__)

#define BH_DEBUG_PRINTIMPL(X,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){ std::cout << UNIQUE_STR_ID(UNIQUE); _PRINT(X);}
#define BH_DEBUG_PRINT(X) BH_DEBUG_PRINTIMPL(X,__LINE__)
#define BH_DEBUG_IMPL(X,UNIQUE) static bool  UNIQUE_ID(UNIQUE) = need_debug(__FILE__,__FUNCTION__) ; std::string UNIQUE_STR_ID(UNIQUE)=get_info_str(__FILE__,__FUNCTION__,UNIQUE); if (UNIQUE_ID(UNIQUE)){ if (UNIQUE_STR_ID(UNIQUE).size()!=0){ _MESSAGE2("::begin of ",UNIQUE_STR_ID(UNIQUE)); }; X; if (UNIQUE_STR_ID(UNIQUE).size()!=0){ _MESSAGE2("::end of ",UNIQUE_STR_ID(UNIQUE)); }; }
#define BH_DEBUG(X) BH_DEBUG_IMPL(X,__LINE__)

#else
#define BH_DEBUG_MESSAGE(X)
#define BH_DEBUG_MESSAGE2(X,Y)
#define BH_DEBUG_MESSAGE3(X,Y,Z)
#define BH_DEBUG_MESSAGE4(X,Y,Z,A)
#define BH_DEBUG_MESSAGE5(X,Y,Z,A,B)
#define BH_DEBUG_MESSAGE6(X,Y,Z,A,B,C)
#define BH_DEBUG_MESSAGE7(X,Y,Z,A,B,C,D)
#define BH_DEBUG_MESSAGE8(X,Y,Z,A,B,C,D,E)
#define BH_DEBUG_MESSAGE9(X,Y,Z,A,B,C,D,E,F)
#define BH_DEBUG_MESSAGE10(X,Y,Z,A,B,C,D,E,F,G)

#define BH_DEBUG_PRINT(X)
#define BH_DEBUG(X)
#endif

namespace BH {

bool need_debug(const char* FileName,const char* FuncName);
std::string get_info_str(const char* FileName,const char* FuncName,int line);

} /* BH */

#endif /* BH_DEBUG_H_ */
