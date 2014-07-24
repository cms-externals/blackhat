#ifndef _H_MC_STORAGE_HPP
#define _H_MC_STORAGE_HPP

#include "mc_storage.h"
#include <iostream>

#include <cassert>


namespace BH {
namespace Tools {


template <class Type,int ChunkSize> void FSArray<Type,ChunkSize>::clear() volatile  {
	d_size=0;
}

template <class Type,int ChunkSize> void FSArray<Type,ChunkSize>::reserve(int requestedSize) volatile {
#if BH_USE_OMP
	d_lock_p->get_lock();
#endif
	int nbrRequested=(requestedSize+ChunkSize-1)/ChunkSize;
	int nbrAdditional=nbrRequested-d_chunkNumber;
	for (int i=1;i<=nbrAdditional;i++) {
		AddStorage();
	}
#if BH_USE_OMP
#pragma omp flush (d_storage,d_maxChunkNumber)
	d_lock_p->release_lock();
#endif
}


template <class Type,int ChunkSize> void FSArray<Type,ChunkSize>::AddStorage() volatile {
#if BH_USE_OMP
	if (  d_lock_p->test_lock() ) {
		std::cerr << "Should not be here without a  lock!" << std::endl;
	}
#endif

	if (d_chunkNumber == d_maxChunkNumber ) {
		d_maxChunkNumber+=d_chunkNumber;
		Type** newStorage=new Type*[d_maxChunkNumber];
		for (int i=0;i<d_chunkNumber;i++) {
			newStorage[i]=d_storage[i];
		}
		//memcpy(newStorage,d_storage, d_chunkNumber*sizeof(Type*));
		delete[] d_storage;
		d_storage=newStorage;
#if BH_USE_OMP
#pragma omp flush (d_storage,d_maxChunkNumber)
#endif
	}
	d_storage[d_chunkNumber]=new Type[ChunkSize];
	d_capacity+=ChunkSize;
	d_chunkNumber++;
#if BH_USE_OMP
#pragma omp flush (d_capacity,d_storage,d_chunkNumber)
#endif
}


template <class Type,int ChunkSize> FSArray<Type,ChunkSize>::FSArray() :
		d_size(0),
		d_capacity(0),
		d_chunkNumber(0),
		d_maxChunkNumber(FSArray_settings::s_initialMaxChunkNumber) 	{
	d_storage=new Type*[FSArray_settings::s_initialMaxChunkNumber];
#if BH_USE_OMP
	d_lock_p=new HasOMPLock();
	d_lock_p->get_lock();
#endif
	AddStorage();
#if BH_USE_OMP
	d_lock_p->release_lock();
#endif
}

template <class Type,int ChunkSize> const Type& FSArray<Type,ChunkSize>::operator[](int n) const volatile {
	assert( n<d_size);
	int chunk=n/ChunkSize;
	int position=n%ChunkSize;
	return d_storage[chunk][position];
};


template <class Type,int ChunkSize> int FSArray<Type,ChunkSize>::push_back(Type value) volatile {
#if BH_USE_OMP
	d_lock_p->get_lock();
#endif
	int chunk=d_size/ChunkSize;
	int position=d_size%ChunkSize;
	if ( d_size == d_capacity) {
		AddStorage();
	}
	d_storage[chunk][position]=value;
	int currentSize=++d_size;  // this is needed because another thread could update d_size before the current thread returns
#if BH_USE_OMP
#pragma omp flush (d_size)
	d_lock_p->release_lock();
#endif
	return currentSize;
};

template <class Type,int ChunkSize> FSArray<Type,ChunkSize>::~FSArray() {
	for (int i=0;i<d_chunkNumber;++i) {
		delete[] d_storage[i];
	}
	delete[] d_storage;
#if BH_USE_OMP
	delete d_lock_p;
#endif
}

} /* Tools */
} /* BH */

#endif /* _H_MC_STORAGE_HPP */
