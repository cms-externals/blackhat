/*
 * iterators.h
 *
 *  Created on: Sep 18, 2008
 *      Author: daniel
 */

#ifndef ITERATORS_H_
#define ITERATORS_H_



namespace BH {
namespace iterators {


template <class T,class Container> class cyclic_iterator;

template <class T,class Container> bool operator==(cyclic_iterator<T,Container> ci1,cyclic_iterator<T,Container> ci2);
template <class T, class Container> cyclic_iterator<T,Container> operator+(cyclic_iterator<T,Container> ci1,int offset);

template <class T,class Container> class cyclic_iterator : public std::iterator <std::bidirectional_iterator_tag,T>{
private:
	size_t m_length;
	size_t m_pos;    // one-based
	size_t m_rep;
	size_t m_maxrep;
	size_t m_base;
	bool m_end;
	const Container* m_cont;
	friend bool operator==<>(cyclic_iterator<T,Container> ci1,cyclic_iterator<T,Container> ci2);
	friend cyclic_iterator<T,Container> operator+<>(cyclic_iterator<T,Container> ci1,int offset);
public:
	cyclic_iterator(const Container& c,int n,size_t base=1): m_cont(&c), m_pos(1), m_rep(1),m_maxrep(n), m_end(false), m_length(c.size()), m_base(base-1) {};
	cyclic_iterator(const Container& c,int n,typename Container::const_iterator iter): m_cont(&c), m_pos(1), m_rep(1),m_maxrep(n), m_end(false), m_length(c.size()), m_base(iter-c.begin()) { if (m_base==m_length) m_end=true;};
	cyclic_iterator(const Container& c,int n,const cyclic_iterator<T,Container>& iter): m_cont(&c), m_pos(1), m_rep(1),m_maxrep(n), m_end(false), m_length(c.size()), m_base(iter.position()-1) { if (m_base==m_length) m_end=true;};
	cyclic_iterator(const Container& c): m_cont(&c), m_end(true) {};
	const T& operator*() const {return (*m_cont)[(m_pos+m_base-1)%m_length];};
	cyclic_iterator<T,Container> operator++();
	cyclic_iterator<T,Container> operator--();
	bool operator!=(const cyclic_iterator<T,Container>& other);
	size_t position() const {return (m_pos+m_base-1)%m_length+1;};
	void restart(){m_pos=1;m_rep=1;};
};

template <class T,class Container> bool cyclic_iterator<T,Container>::operator!=(const cyclic_iterator<T,Container>& other){
	if ( other.m_end == true &&  m_end==true) return false;
	if ( other.m_end != m_end) return true;
	if ( (other.m_pos+other.m_base)%other.m_length == (m_pos+m_base)%m_length ) return false;

return true;
}


template <class T,class Container> cyclic_iterator<T,Container> cyclic_iterator<T,Container>::operator++(){
	if ( m_pos < m_length ) {m_pos++; return *this;} else {
		if (m_rep<m_maxrep ){
			m_rep++; m_pos=1;
			return *this;
		}
		else {
			m_end=true;
			return *this;
		}
	}
}

template <class T,class Container> cyclic_iterator<T,Container> cyclic_iterator<T,Container>::operator--(){
	if ( m_pos > 1 ) {m_pos--; return *this;} else {
		if (m_rep<m_maxrep ){
			m_rep++; m_pos=m_end;
			return *this;
		}
		else {
			m_end=true;
			return *this;
		}
	}
}

template <class T,class Container> bool operator==(cyclic_iterator<T,Container> ci1,cyclic_iterator<T,Container> ci2){
	return ci1.m_end == ci2.m_end;
}
//could be done more efficiently, but don't bother now
template <class T, class Container> cyclic_iterator<T,Container> operator+(cyclic_iterator<T,Container> ci,int offset){
	cyclic_iterator<T,Container> res(ci);
	if (offset<0) { std::cerr << "Negative offset nt operator+(cyclyc_iterator,int)" << std::endl;return res;};
	for (int i=1;i<=offset;i++){
		++res;
	}
	return res;

}

}

}

using BH::iterators::cyclic_iterator ;

#endif /* ITERATORS_H_ */
