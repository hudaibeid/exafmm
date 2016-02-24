#ifndef HILBERT_WRAPPER
#define HILBERT_WRAPPER

#include <unordered_set>
#include "hilbert.h"
#include "logger.h"
#include "types.h"
#include "cell_mutex.h"
using namespace exafmm;
#define cast_uint32(V)               static_cast<uint32_t>(V)
#define COUNT_COLLISIONS 0	
#define DIM 3
#define ACCURACY 10000							  							// multiplier to maintain accuracy 		


//! Output one 64-bit Hilbert order from the 3D transposed key as per Skilling's paper
template<typename ArrayType> 
inline hilbert_t flattenTransposedKey(ArrayType const& arr, int order) {
		hilbert_t key = 0;
		int shifts = order-1;
		for(int i = shifts; i >= 0; --i) { 													// flatten transposed key		
			for(int j = 0; j < DIM; ++j) {
				if(arr[j]>>i & 1)																				// check x-bit
					key |= 1ull << ((DIM * i) + (DIM - j - 1));				   	// place x-bit in position 
			}
		}
		return key;
}

template<typename ArrayType>
inline ArrayType unflattenKey(hilbert_t const& key, int order) {	  
	  ArrayType arr{{0,0,0}};
		hilbert_t mask = 1ull<< (DIM * order - 1);		
		for(int i = 0; i < order; ++i) { 														// flatten transposed key							
			uint32_t inner_mask = (1ul << (order - i - 1));
			for(int j = 0; j < DIM; ++j) {
				if(key & mask)
					arr[j] |= inner_mask;
				mask >>= 1;
			}		
		}
		return arr;
}

inline bool isNeighbor(hilbert_t const& key1, hilbert_t const& key2, int order, int distance = 1) {
	auto&& arr1 = unflattenKey<hilbert_arr>(key1,order);
	auto&& arr2 = unflattenKey<hilbert_arr>(key2,order);
	for (int i = 0; i < DIM; ++i)
		if (std::abs(arr1[i] - arr2[i]) > distance)
			return false;
	//std::cout << key1 << " and " << key2 << " were detected as neigbors" << std::endl;
	return true;
}

//! Next power of 2 such that 2^p > val
template<typename T>
uint32_t nextPow2(T&& val) {
	return uint32_t(ceil(log(val)/log(2)));
}
//! Get Morton order for a specific posision
template<typename ArrayType>
inline uint64_t getMoronKey(ArrayType position, size_t order, size_t dim) {	
	uint64_t id = 0;
	for (size_t i = 0; i < order; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			id += (position[j] & 1) << (3 * i + j); 
			position[j] >>= 1;
		}
	}
	return id;	
}

//! Generate Hilbert order from given the particles structure 
//! This function also maps the floating point numbers to unsigned integers 
//! in order to fit in the space filling curve semantics 
//! and returns min/max orders
template<typename BodyType, typename BoundsType>
std::pair<hilbert_t,hilbert_t> generateHilbertKey(BodyType& bodies, BoundsType const& bounds, uint32_t& order) {	
	logger::startTimer("Hkey Generation");													// start Hilbert generation timer
	auto&& _min = min(bounds.Xmin);																	// get min of all dimensions
	auto&& _max = max(bounds.Xmax);																	// get max of all dimensions 	
	hilbert_t min_h = 0ull;																					// initialize min Hilbert order
	hilbert_t max_h = 0ull;																					// initialize max Hilbert order	
	order  = nextPow2(ACCURACY);			  					 							    // get needed bits to represent Hilbert order in each dimension 	
	real_t diameter = _max - _min;
	assert(order <= 21);																						// maximum key is 63 bits	
#if COUNT_COLLISIONS
  std::unordered_set<hilbert_t> myset;
	size_t collisions = 0;	
#endif 
#pragma omp parallel shared(max_h, min_h)
{
#pragma omp for reduction (max: max_h), reduction(min: min_h)
	for(size_t i = 0; i < bodies.size(); ++i) {											// loop over bodies
		auto&& body = bodies[i];
		hilbert_arr position 	  											// initialize shifted position 
		{{ cast_uint32((body.X[0] - _min)/diameter * ACCURACY),
		 	 cast_uint32((body.X[1] - _min)/diameter * ACCURACY), 
		   cast_uint32((body.X[2] - _min)/diameter * ACCURACY)}};  
#if HILBERT_CODE
		AxestoTranspose(position,order,DIM);														// get transposed Hilbert order
		body.ICELL = flattenTransposedKey(position,order);			      	// generate 1 flat Hilbert order (useful for efficient sorting)		
#else 
		body.ICELL = getMoronKey(position,order,DIM);			      				// generate 1 flat Morton order (useful for efficient sorting)	
#endif		
#if TEST_SFC		
		body.MORTON= getMoronKey(position,order,DIM);
#endif
#if COUNT_COLLISIONS
lock_mutex		
		auto&& iter = myset.find(body.ICELL);
		if(iter == myset.end())
			myset.insert(body.ICELL);
		else 
			collisions++;
unlock_mutex					
#endif 		
		if(min_h != 0ull) {
			if(body.ICELL > max_h)	  																// update min/max hilbert orders 
				max_h = body.ICELL;																				
			else if (body.ICELL < min_h)
				min_h = body.ICELL;
	  }
	  else {
	  	min_h = max_h = body.ICELL;
	  	assert(min_h != 0ull);
	  }
	}
}																																// end loop																														 
	logger::stopTimer("Hkey Generation");														// stop Hilbert generation timer
	return std::make_pair(min_h,max_h);														  // return min/max tuple
}

#endif
