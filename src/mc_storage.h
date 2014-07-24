#ifndef _H_MC_STORAGE_H
#define _H_MC_STORAGE_H

#if BH_USE_OMP
#include <omp.h>
#include "BH_omp.h"

#endif


namespace BH {
namespace Tools {

struct FSArray_settings {
	static int s_initialMaxChunkNumber;
};

template <class Type,int ChunkSize> class FSArray
{
	typedef Type* StorageChunk;
	volatile long d_size;
	volatile long d_capacity;
	volatile int d_chunkNumber;
	volatile int d_maxChunkNumber;
	StorageChunk * volatile d_storage;
#if BH_USE_OMP
	HasOMPLock *d_lock_p;
#endif 
public:
	FSArray();
	~FSArray();
	volatile int size() const volatile {	return d_size;};
	const Type& operator[](int n) const volatile ; //zero based
	int push_back(Type value) volatile ;
	void reserve(int requestedSize) volatile ;
	void clear() volatile ;
private:
	void AddStorage() volatile;
};



} /* Tools */
} /* BH */




#endif /* _H_MC_STORAGE_H */

