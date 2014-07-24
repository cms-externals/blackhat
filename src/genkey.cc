/* New key generation routines */

/* David A. Kosower, May 22, 2008 */

#include "mom_conf.h"
#include "key.h"
#include "BH_error.h"
#include <cstdarg>
#include "genkey.h"

using namespace std;

#define Lw1(c1,c2) c1
#define Lw2(c1,c2) c2
const int firstdigitofpair16[] = {
  Lw1('0','0'), Lw1('0','1'), 
  Lw1('0','2'), Lw1('0','3'), Lw1('0','4'), Lw1('0','5'), Lw1('0','6'), Lw1('0','7'), Lw1('0','8'), Lw1('0','9'), 
  Lw1('0','a'), Lw1('0','b'), Lw1('0','c'), Lw1('0','d'), Lw1('0','e'), Lw1('0','f'),
  Lw1('1','0'), Lw1('1','1'), 
  Lw1('1','2'), Lw1('1','3'), Lw1('1','4'), Lw1('1','5'), Lw1('1','6'), Lw1('1','7'), Lw1('1','8'), Lw1('1','9'), 
  Lw1('1','a'), Lw1('1','b'), Lw1('1','c'), Lw1('1','d'), Lw1('1','e'), Lw1('1','f'),
  Lw1('2','0'), Lw1('2','1'), 
  Lw1('2','2'), Lw1('2','3'), Lw1('2','4'), Lw1('2','5'), Lw1('2','6'), Lw1('2','7'), Lw1('2','8'), Lw1('2','9'), 
  Lw1('2','a'), Lw1('2','b'), Lw1('2','c'), Lw1('2','d'), Lw1('2','e'), Lw1('2','f'),
  Lw1('3','0'), Lw1('3','1'), 
  Lw1('3','2'), Lw1('3','3'), Lw1('3','4'), Lw1('3','5'), Lw1('3','6'), Lw1('3','7'), Lw1('3','8'), Lw1('3','9'), 
  Lw1('3','a'), Lw1('3','b'), Lw1('3','c'), Lw1('3','d'), Lw1('3','e'), Lw1('3','f'),
  Lw1('4','0'), Lw1('4','1'), 
  Lw1('4','2'), Lw1('4','3'), Lw1('4','4'), Lw1('4','5'), Lw1('4','6'), Lw1('4','7'), Lw1('4','8'), Lw1('4','9'), 
  Lw1('4','a'), Lw1('4','b'), Lw1('4','c'), Lw1('4','d'), Lw1('4','e'), Lw1('4','f'),
  Lw1('5','0'), Lw1('5','1'), 
  Lw1('5','2'), Lw1('5','3'), Lw1('5','4'), Lw1('5','5'), Lw1('5','6'), Lw1('5','7'), Lw1('5','8'), Lw1('5','9'), 
  Lw1('5','a'), Lw1('5','b'), Lw1('5','c'), Lw1('5','d'), Lw1('5','e'), Lw1('5','f'),
  Lw1('6','0'), Lw1('6','1'), 
  Lw1('6','2'), Lw1('6','3'), Lw1('6','4'), Lw1('6','5'), Lw1('6','6'), Lw1('6','7'), Lw1('6','8'), Lw1('6','9'), 
  Lw1('6','a'), Lw1('6','b'), Lw1('6','c'), Lw1('6','d'), Lw1('6','e'), Lw1('6','f'),
  Lw1('7','0'), Lw1('7','1'), 
  Lw1('7','2'), Lw1('7','3'), Lw1('7','4'), Lw1('7','5'), Lw1('7','6'), Lw1('7','7'), Lw1('7','8'), Lw1('7','9'), 
  Lw1('7','a'), Lw1('7','b'), Lw1('7','c'), Lw1('7','d'), Lw1('7','e'), Lw1('7','f'),
  Lw1('8','0'), Lw1('8','1'), 
  Lw1('8','2'), Lw1('8','3'), Lw1('8','4'), Lw1('8','5'), Lw1('8','6'), Lw1('8','7'), Lw1('8','8'), Lw1('8','9'), 
  Lw1('8','a'), Lw1('8','b'), Lw1('8','c'), Lw1('8','d'), Lw1('8','e'), Lw1('8','f'), 
  Lw1('9','0'), Lw1('9','1'), 
  Lw1('9','2'), Lw1('9','3'), Lw1('9','4'), Lw1('9','5'), Lw1('9','6'), Lw1('9','7'), Lw1('9','8'), Lw1('9','9'), 
  Lw1('9','a'), Lw1('9','b'), Lw1('9','c'), Lw1('9','d'), Lw1('9','e'), Lw1('9','f'), 
  Lw1('a','0'), Lw1('a','1'), 
  Lw1('a','2'), Lw1('a','3'), Lw1('a','4'), Lw1('a','5'), Lw1('a','6'), Lw1('a','7'), Lw1('a','8'), Lw1('a','9'), 
  Lw1('a','a'), Lw1('a','b'), Lw1('a','c'), Lw1('a','d'), Lw1('a','e'), Lw1('a','f'), 
  Lw1('b','0'), Lw1('b','1'), 
  Lw1('b','2'), Lw1('b','3'), Lw1('b','4'), Lw1('b','5'), Lw1('b','6'), Lw1('b','7'), Lw1('b','8'), Lw1('b','9'), 
  Lw1('b','a'), Lw1('b','b'), Lw1('b','c'), Lw1('b','d'), Lw1('b','e'), Lw1('b','f'), 
  Lw1('c','0'), Lw1('c','1'), 
  Lw1('c','2'), Lw1('c','3'), Lw1('c','4'), Lw1('c','5'), Lw1('c','6'), Lw1('c','7'), Lw1('c','8'), Lw1('c','9'), 
  Lw1('c','a'), Lw1('c','b'), Lw1('c','c'), Lw1('c','d'), Lw1('c','e'), Lw1('c','f'), 
  Lw1('d','0'), Lw1('d','1'), 
  Lw1('d','2'), Lw1('d','3'), Lw1('d','4'), Lw1('d','5'), Lw1('d','6'), Lw1('d','7'), Lw1('d','8'), Lw1('d','9'), 
  Lw1('d','a'), Lw1('d','b'), Lw1('d','c'), Lw1('d','d'), Lw1('d','e'), Lw1('d','f'), 
  Lw1('e','0'), Lw1('e','1'), 
  Lw1('e','2'), Lw1('e','3'), Lw1('e','4'), Lw1('e','5'), Lw1('e','6'), Lw1('e','7'), Lw1('e','8'), Lw1('e','9'), 
  Lw1('e','a'), Lw1('e','b'), Lw1('e','c'), Lw1('e','d'), Lw1('e','e'), Lw1('e','f'), 
  Lw1('f','0'), Lw1('f','1'), 
  Lw1('f','2'), Lw1('f','3'), Lw1('f','4'), Lw1('f','5'), Lw1('f','6'), Lw1('f','7'), Lw1('f','8'), Lw1('f','9'), 
  Lw1('f','a'), Lw1('f','b'), Lw1('f','c'), Lw1('f','d'), Lw1('f','e'), Lw1('f','f')
};

const int seconddigitofpair16[] = {
  Lw2('0','0'), Lw2('0','1'), 
  Lw2('0','2'), Lw2('0','3'), Lw2('0','4'), Lw2('0','5'), Lw2('0','6'), Lw2('0','7'), Lw2('0','8'), Lw2('0','9'), 
  Lw2('0','a'), Lw2('0','b'), Lw2('0','c'), Lw2('0','d'), Lw2('0','e'), Lw2('0','f'),
  Lw2('1','0'), Lw2('1','1'), 
  Lw2('1','2'), Lw2('1','3'), Lw2('1','4'), Lw2('1','5'), Lw2('1','6'), Lw2('1','7'), Lw2('1','8'), Lw2('1','9'), 
  Lw2('1','a'), Lw2('1','b'), Lw2('1','c'), Lw2('1','d'), Lw2('1','e'), Lw2('1','f'),
  Lw2('2','0'), Lw2('2','1'), 
  Lw2('2','2'), Lw2('2','3'), Lw2('2','4'), Lw2('2','5'), Lw2('2','6'), Lw2('2','7'), Lw2('2','8'), Lw2('2','9'), 
  Lw2('2','a'), Lw2('2','b'), Lw2('2','c'), Lw2('2','d'), Lw2('2','e'), Lw2('2','f'),
  Lw2('3','0'), Lw2('3','1'), 
  Lw2('3','2'), Lw2('3','3'), Lw2('3','4'), Lw2('3','5'), Lw2('3','6'), Lw2('3','7'), Lw2('3','8'), Lw2('3','9'), 
  Lw2('3','a'), Lw2('3','b'), Lw2('3','c'), Lw2('3','d'), Lw2('3','e'), Lw2('3','f'),
  Lw2('4','0'), Lw2('4','1'), 
  Lw2('4','2'), Lw2('4','3'), Lw2('4','4'), Lw2('4','5'), Lw2('4','6'), Lw2('4','7'), Lw2('4','8'), Lw2('4','9'), 
  Lw2('4','a'), Lw2('4','b'), Lw2('4','c'), Lw2('4','d'), Lw2('4','e'), Lw2('4','f'),
  Lw2('5','0'), Lw2('5','1'), 
  Lw2('5','2'), Lw2('5','3'), Lw2('5','4'), Lw2('5','5'), Lw2('5','6'), Lw2('5','7'), Lw2('5','8'), Lw2('5','9'), 
  Lw2('5','a'), Lw2('5','b'), Lw2('5','c'), Lw2('5','d'), Lw2('5','e'), Lw2('5','f'),
  Lw2('6','0'), Lw2('6','1'), 
  Lw2('6','2'), Lw2('6','3'), Lw2('6','4'), Lw2('6','5'), Lw2('6','6'), Lw2('6','7'), Lw2('6','8'), Lw2('6','9'), 
  Lw2('6','a'), Lw2('6','b'), Lw2('6','c'), Lw2('6','d'), Lw2('6','e'), Lw2('6','f'),
  Lw2('7','0'), Lw2('7','1'), 
  Lw2('7','2'), Lw2('7','3'), Lw2('7','4'), Lw2('7','5'), Lw2('7','6'), Lw2('7','7'), Lw2('7','8'), Lw2('7','9'), 
  Lw2('7','a'), Lw2('7','b'), Lw2('7','c'), Lw2('7','d'), Lw2('7','e'), Lw2('7','f'),
  Lw2('8','0'), Lw2('8','1'), 
  Lw2('8','2'), Lw2('8','3'), Lw2('8','4'), Lw2('8','5'), Lw2('8','6'), Lw2('8','7'), Lw2('8','8'), Lw2('8','9'), 
  Lw2('8','a'), Lw2('8','b'), Lw2('8','c'), Lw2('8','d'), Lw2('8','e'), Lw2('8','f'), 
  Lw2('9','0'), Lw2('9','1'), 
  Lw2('9','2'), Lw2('9','3'), Lw2('9','4'), Lw2('9','5'), Lw2('9','6'), Lw2('9','7'), Lw2('9','8'), Lw2('9','9'), 
  Lw2('9','a'), Lw2('9','b'), Lw2('9','c'), Lw2('9','d'), Lw2('9','e'), Lw2('9','f'), 
  Lw2('a','0'), Lw2('a','1'), 
  Lw2('a','2'), Lw2('a','3'), Lw2('a','4'), Lw2('a','5'), Lw2('a','6'), Lw2('a','7'), Lw2('a','8'), Lw2('a','9'), 
  Lw2('a','a'), Lw2('a','b'), Lw2('a','c'), Lw2('a','d'), Lw2('a','e'), Lw2('a','f'), 
  Lw2('b','0'), Lw2('b','1'), 
  Lw2('b','2'), Lw2('b','3'), Lw2('b','4'), Lw2('b','5'), Lw2('b','6'), Lw2('b','7'), Lw2('b','8'), Lw2('b','9'), 
  Lw2('b','a'), Lw2('b','b'), Lw2('b','c'), Lw2('b','d'), Lw2('b','e'), Lw2('b','f'), 
  Lw2('c','0'), Lw2('c','1'), 
  Lw2('c','2'), Lw2('c','3'), Lw2('c','4'), Lw2('c','5'), Lw2('c','6'), Lw2('c','7'), Lw2('c','8'), Lw2('c','9'), 
  Lw2('c','a'), Lw2('c','b'), Lw2('c','c'), Lw2('c','d'), Lw2('c','e'), Lw2('c','f'), 
  Lw2('d','0'), Lw2('d','1'), 
  Lw2('d','2'), Lw2('d','3'), Lw2('d','4'), Lw2('d','5'), Lw2('d','6'), Lw2('d','7'), Lw2('d','8'), Lw2('d','9'), 
  Lw2('d','a'), Lw2('d','b'), Lw2('d','c'), Lw2('d','d'), Lw2('d','e'), Lw2('d','f'), 
  Lw2('e','0'), Lw2('e','1'), 
  Lw2('e','2'), Lw2('e','3'), Lw2('e','4'), Lw2('e','5'), Lw2('e','6'), Lw2('e','7'), Lw2('e','8'), Lw2('e','9'), 
  Lw2('e','a'), Lw2('e','b'), Lw2('e','c'), Lw2('e','d'), Lw2('e','e'), Lw2('e','f'), 
  Lw2('f','0'), Lw2('f','1'), 
  Lw2('f','2'), Lw2('f','3'), Lw2('f','4'), Lw2('f','5'), Lw2('f','6'), Lw2('f','7'), Lw2('f','8'), Lw2('f','9'), 
  Lw2('f','a'), Lw2('f','b'), Lw2('f','c'), Lw2('f','d'), Lw2('f','e'), Lw2('f','f')
};

const char digit64[] = {
  '0', '1', 
  '2', '3', '4', '5', '6', '7', '8', '9', 
  'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 
  'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 
  'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 
  '$', '@'
};

namespace BH 
{
  
// Convert a 16-bit integer to a base-16 representation using the
// digits 0-9, lower-case a-f.
inline void ToBase16(char *digits,int value) {
  int ix;
  ix = (value >> 8)&0xFF;
  *digits++ = firstdigitofpair16[ix];
  *digits++ = seconddigitofpair16[ix];
  ix = value&0xFF;
  *digits++ = firstdigitofpair16[ix];
  *digits = seconddigitofpair16[ix];
}

// Convert a 16-bit integer to a base-64 representation using the
// digits 0-9, lower-case a-z, upper-case A-Z, and $ and @.
inline void ToBase64(char *digits,int value) {
  *digits++ = digit64[(value >> 12)&0x3F];
  *digits++ = digit64[(value >> 6)&0x3F];
  *digits = digit64[value&0x3F];
}


char separator = ':';

//! Key generation
string GenKey1(const string& tag,int i1){
  char res[_MAX_KEY_LENGTH];
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  //  sprintf(res,"%s:%d",tag,i1);
  ToBase64(out,i1);
  out += 3;  *out = 0;
  return res;
}

string GenKey1(const string& tag,int i1, int i2){
  char res[_MAX_KEY_LENGTH];
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  //  sprintf(res,"%s:%d",tag,i1);
  ToBase64(out,i1);  out += 3;
  ToBase64(out,i2);
  out += 3;  *out = 0;
  return res;
}

string GenKey1(const string& tag,int i1, int i2, int i3){
  char res[_MAX_KEY_LENGTH];
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  //  sprintf(res,"%s:%d",tag,i1);
  ToBase64(out,i1);  out += 3;
  ToBase64(out,i2);  out += 3;
  ToBase64(out,i3);
  out += 3;  *out = 0;
  return res;
}

string GenKey1(const string& tag,int i1, int i2, int i3,int i4){
  char res[_MAX_KEY_LENGTH];
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  //  sprintf(res,"%s:%d",tag,i1);
  ToBase64(out,i1);  out += 3;
  ToBase64(out,i2);  out += 3;
  ToBase64(out,i3);  out += 3;
  ToBase64(out,i4);
  out += 3;  *out = 0;
  return res;
}

string GenKey1(const string& tag,const vector<int>& v){
  if (v.size()*3 + tag.length()+1 >= _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  for (size_t j = 0;  j < v.size();  j++)
     {ToBase64(out,v[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey1(const string& tag, const vector<int>& v1, 
                 const vector<int>& v2) 
{ if ((v1.size()+v2.size())*3 + tag.length()+1 >= _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  for (size_t j = 0;  j < v1.size();  j++)
     {ToBase64(out,v1[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v2.size();  j++)
     {ToBase64(out,v2[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey1(const string& tag, const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3)
{ if ((v1.size()+v2.size()+v3.size())*3 + tag.length()+1 >= _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  for (size_t j = 0;  j < v1.size();  j++)
     {ToBase64(out,v1[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v2.size();  j++)
     {ToBase64(out,v2[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v3.size();  j++)
     {ToBase64(out,v3[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey1(const string& tag, const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3,
               const vector<int>& v4)
{ if ((v1.size()+v2.size()+v3.size()+v4.size())*3 + tag.length()+1 >= 
      _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  for (size_t j = 0;  j < v1.size();  j++)
     {ToBase64(out,v1[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v2.size();  j++)
     {ToBase64(out,v2[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v3.size();  j++)
     {ToBase64(out,v3[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v4.size();  j++)
     {ToBase64(out,v4[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey1(const string& tag, int i1, int i2,
               const vector<int>& v1)
{ if ((v1.size()+2)*3 + tag.length()+4 >= 
      _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  ToBase64(out,i1);
  out += 3;
  *out++ = separator;
  ToBase64(out,i2);
  out += 3;
  *out++ = separator;
  for (size_t j = 0;  j < v1.size();  j++)
     {ToBase64(out,v1[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey1(const string& tag, int i1, int i2,
               const vector<int>& v1, 
               const vector<int>& v2)
{ if ((v1.size()+v2.size()+2)*3 + tag.length()+5 >= 
      _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  ToBase64(out,i1);
  out += 3;
  *out++ = separator;
  ToBase64(out,i2);
  out += 3;
  *out++ = separator;
  for (size_t j = 0;  j < v1.size();  j++)
     {ToBase64(out,v1[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v2.size();  j++)
     {ToBase64(out,v2[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey1(const string& tag, int i1, int i2,
               const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3)
{ if ((v1.size()+v2.size()+v3.size()+2)*3 + tag.length()+6 >= 
      _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  ToBase64(out,i1);
  out += 3;
  *out++ = separator;
  ToBase64(out,i2);
  out += 3;
  *out++ = separator;
  for (size_t j = 0;  j < v1.size();  j++)
     {ToBase64(out,v1[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v2.size();  j++)
     {ToBase64(out,v2[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v3.size();  j++)
     {ToBase64(out,v3[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey1(const string& tag, int i1, int i2,
               const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3,
               const vector<int>& v4)
{ if ((v1.size()+v2.size()+v3.size()+v4.size()+2)*3 + tag.length()+7 >= 
      _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  ToBase64(out,i1);
  out += 3;
  *out++ = separator;
  ToBase64(out,i2);
  out += 3;
  *out++ = separator;
  for (size_t j = 0;  j < v1.size();  j++)
     {ToBase64(out,v1[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v2.size();  j++)
     {ToBase64(out,v2[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v3.size();  j++)
     {ToBase64(out,v3[j]);  out += 3;}
  *out++ = separator;
  for (size_t j = 0;  j < v4.size();  j++)
     {ToBase64(out,v4[j]);  out += 3;}
  *out = 0;
  return res;
}

string GenKey2(const string& tag,int i1){
  char res[_MAX_KEY_LENGTH];
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  //  sprintf(res,"%s:%d",tag,i1);
  ToBase16(out,i1);
  out += 4;  *out = 0;
  return res;
}

string GenKey2(const string& tag,int i1,int i2){
  char res[_MAX_KEY_LENGTH];
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  //  sprintf(res,"%s:%d",tag,i1);
  ToBase16(out,i1);
  out += 4;
  ToBase16(out,i2);
  out += 4;  *out = 0;
  return res;
}

string GenKey2(const string& tag,int i1,int i2,int i3){
  char res[_MAX_KEY_LENGTH];
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  //  sprintf(res,"%s:%d",tag,i1);
  ToBase16(out,i1);
  out += 4;
  ToBase16(out,i2);
  out += 4;
  ToBase16(out,i3);
  out += 4;  *out = 0;
  return res;
}
 
string GenKey2(const string& tag,const vector<int>& v){
  if (v.size()*3 + tag.length()+1 >= _MAX_KEY_LENGTH)
     {throw BHerror("Overflow in key generation.");}
  char res[_MAX_KEY_LENGTH];
  //	sprintf(res,"%s:%s",tag,VectorToString(v).c_str());
  //  while (*tag) *out++ = *tag++;
  tag.copy(res,sizeof(res));
  char *out = res+tag.length();
  *out++ = separator;
  for (size_t j = 0;  j < v.size();  j++)
     {ToBase16(out,v[j]);  out += 4;}
  *out = 0;
  return res;
}
}
