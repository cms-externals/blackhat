/* Key generation for caching */

/*  David A. Kosower, November 19, 2008

    Include file for routines implemented in genkey.cc, fast key generation

 */

#ifndef GenKeyDefined

using namespace std;

namespace BH {
string GenKey1(const string& tag,int i1);
string GenKey1(const string& tag,int i1, int i2);
string GenKey1(const string& tag,int i1, int i2, int i3);
string GenKey1(const string& tag,int i1, int i2, int i3,int i4);
string GenKey1(const string& tag,const vector<int>& v);
string GenKey1(const string& tag, const vector<int>& v1, 
               const vector<int>& v2);
string GenKey1(const string& tag, const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3);
string GenKey1(const string& tag, const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3,
               const vector<int>& v4);
string GenKey1(const string& tag, int i1, int i2, const vector<int>& v1);
string GenKey1(const string& tag, int i1, int i2, const vector<int>& v1, 
               const vector<int>& v2);
string GenKey1(const string& tag, int i1, int i2, const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3);
string GenKey1(const string& tag, int i1, int i2, const vector<int>& v1, 
               const vector<int>& v2,
               const vector<int>& v3,
               const vector<int>& v4);
string GenKey2(const string& tag,int i1);
string GenKey2(const string& tag,int i1,int i2);
string GenKey2(const string& tag,int i1,int i2,int i3);
string GenKey2(const string& tag,const vector<int>& v);
}

#define GenKeyDefined 1
#endif /* GenKeyDefined */
