/*
 * key.h
 *
 *  Created on: 8 Mar 2010
 *      Author: daniel
 */

#ifndef KEY_H_
#define KEY_H_

#include <string>
#include <vector>

namespace BH {

const unsigned int _MAX_KEY_LENGTH=256;

#define Lw1(c1,c2) c1
#define Lw2(c1,c2) c2


void ToBase16(char *digits,int value) ;

void ToBase64(char *digits,int value);

std::string GenKey1(const char* tag,int i1);

std::string GenKey1(const char* tag,int i1, int i2);

std::string GenKey1(const char* tag,int i1, int i2, int i3);

std::string GenKey1(const char* tag,int i1, int i2, int i3,int i4);

std::string GenKey1(const char* tag,const std::vector<int>& v);

std::string GenKey2(const char* tag,int i1);

std::string GenKey2(const char* tag,int i1,int i2);

std::string GenKey2(const char* tag,int i1,int i2,int i3);
std::string GenKey2(const char* tag,const std::vector<int>& v);






// Key for antisymmetrized (i,j), i, j in [1..n]
int orderless_key2(int i, int j);

// Key for symmetrized (i,j), i, j in [1..n]
int SymmetrizedPairKey(int i, int j);


// works for three DIFFERENT i,j,k. Returns the 0-based index of the element (i1,i2,i3)
int orderless_key3(int i, int j, int k) ;

int key_spab(int i, int j , int k);



std::string GenKey(const char* tag,int i1);
std::string GenKey(const char* tag,int i1,int i2);
std::string GenKey(const char* tag,int i1,int i2,int i3);
std::string GenKey(const char* tag,int i1,int i2,int i3,int i4);
std::string GenKey(const char* tag,int i1,int i2,int i3, int i4,int i5);
std::string GenKey(const char* tag,int i1,int i2,int i3,int i4,int i5,int i6);


//! Key generation
std::string GenKey(const char* tag,const std::vector<int>& v);
std::string GenKey(const char* tag,const std::vector<int>& v1,const std::vector<int>& v2);
std::string GenKey(const char* tag,const std::vector<int>& v1,const std::vector<int>& v2,const std::vector<int>& v3);
std::string GenKey(const char* tag,const std::vector<int>& v1,const std::vector<int>& v2,const std::vector<int>& v3,const std::vector<int>& v4);
std::string GenKey(const char* tag,int i1,int i2,const std::vector<int>& v1);
std::string GenKey(const char* tag,int i1,int i2,int i3,const std::vector<int>& v1);
std::string GenKey(const char* tag,int i1,int i2,const std::vector<int>& v1,const std::vector<int>& v2);
std::string GenKey(const char* tag,int i1,int i2,const std::vector<int>& v1,const std::vector<int>& v2,const std::vector<int>& v3);
std::string GenKey(const char* tag,int i1,int i2,const std::vector<int>& v1,const std::vector<int>& v2,const std::vector<int>& v3,const std::vector<int>& v4);


}
#endif /* KEY_H_ */
