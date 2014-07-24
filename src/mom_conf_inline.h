#ifndef MOM_CONF_INLINE_H_
#define MOM_CONF_INLINE_H_

//! inserts a new momentum
/**
 \param m complex momentum to be inserted
\return integer label of the inserted momentum
*/
int orderless_key2(int,int);

	template <class T> inline  const Cmom<T>& momentum_configuration<T>::p(size_t n) const {
	if (n<=momentum_configuration<T>::nbr){
		if ( n>_offset ) {

			return momentum_configuration<T>::ps[n-_offset-1];
		}
		else
		{
			return _parent->p(n);
		}
	}
	else{_WARNING5("Too large momentum index in sub_momentum_configuration::p: ",n," (max=",momentum_configuration<T>::nbr,")"  );throw BHerror("Mom_conf error"); }
}

template <class T> inline std::complex<T>  momentum_configuration<T>::m2(size_t n) const {
	if (n<=momentum_configuration<T>::nbr){
		if ( n>_offset ) {

			return momentum_configuration<T>::ms[n-_offset-1];
		}
		else
		{
			return _parent->m2(n);
		}
	}
	else{_WARNING5("Too large momentum index in sub_momentum_configuration::ms: ",n," (max=",momentum_configuration<T>::nbr,")"  );throw BHerror("Mom_conf error"); }
}




template <class T> inline int momentum_configuration<T>::insert(const lambdat<T>& l1,const lambda<T>& l2){
	return insert(Cmom<T>(l1,l2));
}
template <class T> inline int momentum_configuration<T>::insert(const momentum<std::complex<T> >& p){
	return insert(Cmom<T>(p));
}
template <class T> inline int momentum_configuration<T>::insert(const momentum<std::complex<T> >& p,momentum_type type){
	return insert(Cmom<T>(p,type));
}

template <class T> inline int momentum_configuration<T>::insert(const lambda<T>& l1,const lambdat<T>& l2){
	return insert(Cmom<T>(l2,l1));
}

//! spinor product \<i j\>
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
\return spinor product \<i j\> of the momenta with labels i and j
*/
// DAK 070925 Track sign, use orderless_key2
template <class T> inline std::complex<T> momentum_configuration<T>::spa(int i,int j){
#if _WITH_CACHING
	if (i==j) return std::complex<T>(0.,0.);
	int key=orderless_key2(i,j);
	std::complex<T> res;
	typename std::map<int,std::complex<T> >::iterator pos=Mspa.find(key);
	T sign = (i < j) ? 1 : -1;
#if _CACHE_STATISTICS
	find_and_stat<int,T>("Spa",key,Mspa, sign*p(i).L()*p(j).L() ,p(i).L()*p(j).L(),res );

#else
	if (pos != Mspa.end()){
       res = sign * pos->second ;
	}
	else {
 		res=p(i).L()*p(j).L();
		Mspa.insert(pair<int,std::complex<T> >(key,sign*res));
	}
	#endif
return res;
#else
	return p(i).L()*p(j).L();
#endif


}

//! spinor product [i j]
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
\return spinor product [i j] of the momenta with labels i and j
*/
// DAK 070925 Track sign
template <class T> inline std::complex<T> momentum_configuration<T>::spb(int i,int j){
#if _WITH_CACHING
	if (i==j) return std::complex<T>(0.,0.);
	int key=orderless_key2(i,j);
	std::complex<T> res;
	T sign = (i < j) ? 1 : -1;
#if _CACHE_STATISTICS==0
find_and_stat<int,T>("Spb",key,Mspb, sign*p(i).Lt()*p(j).Lt() ,p(i).Lt()*p(j).Lt(),res );
#else
	typename std::map<int,std::complex<T> >::iterator pos=Mspb.find(key);
	if (pos != Mspb.end()){
       res = sign * pos->second ;
	}
	else {
 		res=p(i).Lt()*p(j).Lt();
		Mspb.insert(pair<int,std::complex<T> >(key,sign*res));
		}

#endif

	return res;
#else
	return p(i).Lt()*p(j).Lt();
#endif
	}



//! spinor product \<i|j|k]
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
\return spinor product <i|j|k] of the momenta with labels i, j and k
*/
template <class T> std::complex<T> momentum_configuration<T>::spab(int i,int j,int k){
	if ((i==j)||(j==k)) return std::complex<T>(0.,0.);
#if _WITH_CACHING
	std::complex<T> res;
	std::string key=GenKey("spab",i,j,k);
#if _CACHE_STATISTICS
find_and_stat<std::string,T>("Spab",key,cache, L(i)*Sm(j)*Lt(k) ,L(i)*Sm(j)*Lt(k),res );
#else
	typename hash_map<std::string, std::std::complex<T>, hash<std::string> >::iterator pos=cache.find(key);
	if (pos != cache.end()) {
		res=pos->second;
	}
	else {
		res=L(i)*Sm(j)*Lt(k);
		cache[key]=res;
	}
#endif
	return res;
	#else
	return L(i)*Sm(j)*Lt(k);
#endif


}
//! spinor product \<i|p1+p2+...+pn|k] for a sum of vectors
/**
The slashed matrix inserted in the spinor product corresponds to the sum of the momenta with the labels from the vector v. The values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param v vector of  integer labels
 \param k integer label of the third momentum
\return spinor product <i|p1+...+pn|k] of the momenta with labels i,{p1,...,pn} and k
*/
template <class T> inline std::complex<T> momentum_configuration<T>::spab(int i,const std::vector<int>& v,int k){
	return this->spab(i,Sum(v),k);
}
//! spinor product \<i|p1+p2+...+pn|k] for a sum of the vectors represented by the plabels in v
/**
The slashed matrix inserted in the spinor product corresponds to the sum of the momenta with the labels from the plabel vector v. The values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param v vector of  plabels
 \param k integer label of the third momentum
\return spinor product <i|p1+...+pn|k] of the momenta with labels i,{p1,...,pn} and k
*/
template <class T> inline std::complex<T> momentum_configuration<T>::spab(int i,const std::vector<plabel>& v,int k){
	return this->spab(i,Sum(v),k);
}
//! spinor product [i|p1+p2+...+pn|k\> for a sum of vectors
/**
The slashed matrix inserted in the spinor product corresponds to the sum of the momenta with the labels from the vector v. The values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param v vector of  integer labels
 \param k integer label of the third momentum
\return spinor product [i|p1+...+pn|k\> of the momenta with labels i,{p1,...,pn} and k
*/
template <class T> inline std::complex<T> momentum_configuration<T>::spba(int i,const std::vector<int>& v,int k){
	return this->spba(i,Sum(v),k);
}
//! spinor product [i|p1+p2+...+pn|k\> for a sum of the vectors represented by the plabels in v
/**
The slashed matrix inserted in the spinor product corresponds to the sum of the momenta with the labels from the plabel vector v. The values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param v vector of  plabels
 \param k integer label of the third momentum
\return spinor product [i|p1+...+pn|k\> of the momenta with labels i,{p1,...,pn} and k
*/
template <class T> inline std::complex<T> momentum_configuration<T>::spba(int i,const std::vector<plabel>& v,int k){
	return this->spba(i,Sum(v),k);
}


//! spinor product \[i|j|k>
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
\return spinor product [i|j|k> of the momenta with labels i, j and k
*/

template <class T> inline std::complex<T> momentum_configuration<T>::spba(int i,int j,int k){
	if ((i==j)||(j==k)) return std::complex<T>(0.,0.);
	return momentum_configuration<T>::spab(k,j,i);
}


//! spinor product \<i|j|k|l\>
/**
 the values are cached
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the fourth momentum
\return spinor product \<i|j|k|l\> of the momenta with labels i, j, k and l
*/

template <class T> inline std::complex<T> momentum_configuration<T>::spaa(int i,int j,int k ,int l){
	if ((i==j)||(k==l)) return std::complex<T>(0.,0.);
#if _WITH_CACHING
	std::string key=GenKey("spaa",i,j,k,l);
	complex<T> res;
#if _CACHE_STATISTICS
	find_and_stat("Spaa4", key,this , L(i)*Sm(j)*Sm(k)*L(l),L(i)*Sm(j)*Sm(k)*L(l),  res );
#else
	if ( !get_value(key,res) ){
		res=L(i)*Sm(j)*Sm(k)*L(l);
		put_value(key,res);
	}
#endif    /* _CACHE_STATISTICS */
	return res;
#else
	return L(i)*Sm(j)*Sm(k)*L(l);
#endif
}

//! spinor product [i|j|k|l]
/**
 the values are cached
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the fourth momentum
\return spinor product [i|j|k|l] of the momenta with labels i, j, k and l
*/

template <class T> inline std::complex<T> momentum_configuration<T>::spbb(int i,int j,int k ,int l){
	if ((i==j)||(k==l)) return std::complex<T>(0.,0.);
#if _WITH_CACHING
	std::string key=GenKey("spbb",i,j,k,l);
	std::complex<T> res;
#if _CACHE_STATISTICS
	find_and_stat("Spbb4", key,this ,p(i).Lt()*p(j).Sm()*p(k).Sm()*p(l).Lt(),p(i).Lt()*p(j).Sm()*p(k).Sm()*p(l).Lt(),  res );
#else
	if ( !get_value(key,res) ){
		res= p(i).Lt()*p(j).Sm()*p(k).Sm()*p(l).Lt();
		put_value(key,res);
	}
#endif /* _CACHE_STATISTICS */
	return res;
#else
	return p(i).Lt()*p(j).Sm()*p(k).Sm()*p(l).Lt();

#endif
}

//! spinor product <i|j|k|l|m]
/**
 the values are cached
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the fourth momentum
 \param m integer label of the fifth momentum
\return spinor product <i|j|k|l|m] of the momenta with labels i, j, k, l and m
*/

template <class T> inline std::complex<T> momentum_configuration<T>::spab(int i,int j,int k ,int l,int m){
	if ((i==j)||(m==l)) return std::complex<T>(0.,0.);
#if _WITH_CACHING
	std::complex<T> res;
	std::string key=GenKey("spab",i,j,k,l,m);

#if _CACHE_STATISTICS
	find_and_stat("Spab5", key,this ,L(i)*Sm(j)*Sm(k)*Sm(l)*Lt(m),L(i)*Sm(j)*Sm(k)*Sm(l)*Lt(m),  res );
#else

	typename hash_map<std::string, std::complex<T>, hash<std::string> >::iterator pos=cache.find(key);
	if (pos != cache.end()) {
		res=pos->second;
	}
	else {
		res=L(i)*Sm(j)*Sm(k)*Sm(l)*Lt(m);
		cache[key]=res;
	}
#endif /* _CACHE_STATISTICS */
	return res;
#else
	return L(i)*Sm(j)*Sm(k)*Sm(l)*Lt(m);
#endif
}

//! spinor product [i|j|k|l|m>
/**
 the values are cached
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the fourth momentum
 \param m integer label of the fifth momentum
\return spinor product [i|j|k|l|m> of the momenta with labels i, j, k, l and m
*/

template <class T> inline std::complex<T> momentum_configuration<T>::spba(int i,int j,int k ,int l,int m){
	if ((i==j)||(m==l)) return std::complex<T>(0.,0.);
	return spab(m,l,k,j,i);
}

//! spinor product \<i|j|k|l|m|n\>
/**
 the values are cached
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the fourth momentum
 \param m integer label of the fifth momentum
 \param n integer label of the fifth momentum
\return spinor product \<i|j|k|l|m|n\> of the momenta with labels i, j, k, l, m and n
*/

template <class T> inline std::complex<T> momentum_configuration<T>::spaa(int i,int j,int k ,int l,int m,int n){
	if ((i==j)||(m==n)) return std::complex<T>(0.,0.);

#if _WITH_CACHING
	std::complex<T> res;
	std::string key=GenKey("spaa",i,j,k,l,m,n);

#if _CACHE_STATISTICS
	find_and_stat("Spaa6", key,this ,L(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*L(n),L(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*L(n),  res );
#else

	typename hash_map<std::string, std::complex<T>, hash<std::string> >::iterator pos=cache.find(key);
	if (pos != cache.end()) {
		res=pos->second;
	}
	else {
		res=L(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*L(n);
		cache[key]=res;
	}
#endif /* _CACHE_STATISTICS */
	return res;
#else
	return L(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*L(n);
#endif
//	if ((i==j)||(m==n)) return 0;
//	return p(i).L()*p(j).Sm()*p(k).Sm()*p(l).Sm()*p(m).Sm()*p(n).L();
}

//! spinor product [i|j|k|l|m|n]
/**
 the values are cached
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the fourth momentum
 \param m integer label of the fifth momentum
 \param n integer label of the fifth momentum
\return spinor product [i|j|k|l|m|n] of the momenta with labels i, j, k, l, m and n
*/

template <class T> inline std::complex<T> momentum_configuration<T>::spbb(int i,int j,int k ,int l,int m,int n){
	if ((i==j)||(m==l)) return std::complex<T>(0.,0.);
#if _WITH_CACHING
	std::complex<T> res;
	std::string key=GenKey("spbb",i,j,k,l,m,n);
#if _CACHE_STATISTICS
	find_and_stat("Spbb6", key,this ,Lt(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*Lt(n),Lt(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*Lt(n),  res );
#else
	typename hash_map<std::string, std::complex<T>, hash<std::string> >::iterator pos=cache.find(key);
	if (pos != cache.end()) {
		res=pos->second;
	}
	else {
		res=Lt(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*Lt(n);
		cache[key]=res;
	}
#endif /* _CACHE_STATISTICS */
	return res;
#else
	return Lt(i)*Sm(j)*Sm(k)*Sm(l)*Sm(m)*Lt(n);
#endif
}

//! invariant s(i,j)
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
\return invariant (pi+pj)^2 of the momenta with labels i, j
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(int i,int j){
// it makes only sense to push the sum of the four vectors into the mom_conf if it is cachedd and we have a chance to reuse it
#if _WITH_CACHING
	return m2(Sum(i,j));
#if _CACHE_STATISTICS
	find_and_stat<int,T>("S(i,j)",key,Mspa, sign*p(i).L()*p(j).L() ,p(i).L()*p(j).L(),res );
#endif
#else
momentum<std::complex<T> > p=mom(i)+mom(j); return p*p;
#endif
}
//! invariant s(i,j,k)
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
\return invariant (pi+pj+pk)^2 of the momenta with labels i, j and k
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(int i,int j,int k){
	// it makes only sense to push the sum of the four vectors into the mom_conf if it is cachedd and we have a chance to reuse it
	#if _WITH_CACHING
		return 	m2(Sum(i,j,k));
	#if _CACHE_STATISTICS
		find_and_stat<int,T>("S(i,j,k)",key,Mspa, sign*p(i).L()*p(j).L() ,p(i).L()*p(j).L(),res );
	#endif
	#else
	momentum<std::complex<T> > p=mom(i)+mom(j)+mom(k); return p*p;
	#endif
}
//! invariant s(i,j,k,l)
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the third momentum
\return invariant (pi+pj+pk+pl)^2 of the momenta with labels i, j, k and l
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(int i,int j,int k,int l){
	// it makes only sense to push the sum of the four vectors into the mom_conf if it is cachedd and we have a chance to reuse it
	#if _WITH_CACHING
		return m2(Sum(i,j,k,l));
	#if _CACHE_STATISTICS
		find_and_stat<int,T>("S(i,j,k,l)",key,Mspa, sign*p(i).L()*p(j).L() ,p(i).L()*p(j).L(),res );
	#endif
	#else
	momentum<std::complex<T> > p=mom(i)+mom(j)+mom(k)+mom(l); return p*p;
	#endif
}
//! invariant s(i,j,k,l,m)
/**
 the values are computed once and cached for later use.
 \param i integer label of the first momentum
 \param j integer label of the second momentum
 \param k integer label of the third momentum
 \param l integer label of the third momentum
\return invariant (pi+pj+pk+pl+pm)^2 of the momenta with labels i, j, k, l and m
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(int i,int j,int k,int l,int m){
	// it makes only sense to push the sum of the four vectors into the mom_conf if it is cachedd and we have a chance to reuse it
	#if _WITH_CACHING
	return m2(Sum(i,j,k,l,m));
	#if _CACHE_STATISTICS
		find_and_stat<int,T>("S(i,j,k,l)",key,Mspa, sign*p(i).L()*p(j).L() ,p(i).L()*p(j).L(),res );
	#endif
	#else
	momentum<std::complex<T> > p=mom(i)+mom(j)+mom(k)+mom(l)+mom(m); return p*p;
	#endif
}
//! invariant s for a set of momenta
/**
 the values are computed once and cached for later use.
 \param v vector of integer label of the momenta
\return invariant (p1+...+pn)^2 of the momenta labelled by the integers in the vector v.
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(const std::vector<int>& v){
#if _WITH_CACHING
	return m2(Sum(v));
#else
	momentum<std::complex<T> > sum(std::complex<T>(0.,0.),std::complex<T>(0.,0.),std::complex<T>(0.,0.),std::complex<T>(0.,0.));
		for (size_t i=0;i< v.size();i++){
			sum+=mom(v[i]);
		}
	return sum*sum;
#endif
}
//! invariant s for two sets of momenta
/**
 the values are computed once and cached for later use.
 \param v1 vector of integer label of the momenta
 \param v2 vector of integer label of the momenta
\return invariant (p1+...+pn)^2 of the momenta labelled by the integers in the vector v.
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(const std::vector<int>& v1,const std::vector<int>& v2){
	return m2(Sum(v1,v2));
}
//! invariant s for a set of momenta
/**
 the values are computed once and cached for later use.
 \param v1 vector of plabel of momenta
 \param v2 vector of plabel of momenta
\return invariant (p1+...+pn)^2 of the momenta labelled by the plabels in the vector v.
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(const std::vector<plabel>& v1,const std::vector<plabel>& v2){
	return m2(Sum(v1,v2));
}
//! invariant s for a set of momenta
/**
 the values are computed once and cached for later use.
 \param v vector of plabel of momenta
\return invariant (p1+...+pn)^2 of the momenta labelled by the plabels in the vector v.
*/
template <class T> inline std::complex<T> momentum_configuration<T>::s(const std::vector<plabel>& v){
	return m2(Sum(v));
}

//! scalar product p.q
/**
 \param i integer label of the first momentum
 \param j integer label of the second momentum
\return the scalar product pi.pj of the momenta with labels i, j
*/
template <class T> inline std::complex<T> momentum_configuration<T>::sp(int i,int j){
	return p(i).P()*p(j).P();
}

//! puts a std::complex value in the cache
/** Puts a value into the cache. It there is already a value, the old value is overwritten. The key can be generated wit GenKey \sa GenKey \sa put_label \sa get_value*/
template <class T> inline void momentum_configuration<T>::put_value(const std::string& key, std::complex<T> &value){
	cache[key]=value;
}

//! puts a label in the cache
/** Puts the new label into the label cache. It there is already a label with this key, the old label is overwritten. The key can be generated wit GenKey \sa GenKey \sa put_label \sa get_label*/
template <class T> inline void momentum_configuration<T>::put_label(const std::string& key,size_t &value){
	labelscache[key]=value;
}

#endif /*MOM_CONF_INLINE_H_*/
