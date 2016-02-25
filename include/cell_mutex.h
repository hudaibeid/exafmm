#ifndef CELL_MUTEX
#define CELL_MUTEX
#include "base_mpi.h"
#include "minicircle.h"
#include <mutex>
#include <iostream>
#include <queue>
using namespace exafmm;

#define NULLTAG 1
#define CELLTAG 2
#define CHILDCELLTAG 3
#define TIMEOUT 4
#define MULTIPOLETAG 5
#define PARTICLETAG 6
#define LEVELTAG 7
#define BODYTAG 8
#define FLUSHTAG 9
#define MAXTAG  15

#define LEVELSHIFT 5
#define REQUESTSHIFT 4
#define DIRECTIONSHIFT 1
#define GRAINSHIFT 16

#define LEVELMASK 0x1F
#define REQUESTMASK 0xF
#define DIRECTIONMASK 0x1
#define GRAINMASK 0xFFFF

#define IDMASK 0x2F

#define RAWPTR(vec) &vec[0]
#define NEWTRIGGER new int()
#define TRIGGER(task) *task = 1;
#define FINALIZE(task) *task = 2;
#define UNTRIGGER(task) *task = 0; 
#define GETTRIG(task) *task
#define DELTRIGGER(task) delete task;
#define TOGGLEDIRECTION(tag) tag ^= DIRECTIONMASK;

std::mutex cell_mutex;
#define lock_mutex cell_mutex.lock();
#define unlock_mutex cell_mutex.unlock();


typedef int* TaskTrigger;
typedef B_iter PointIterator;
typedef const real_t* CoordIterator;
typedef Miniball::Miniball <Miniball::XCoordAccessor<PointIterator, CoordIterator>> MB;


Cell __cell___; 
Body __body___;
Multipole __multipole___;
const int CELLWORD = sizeof(__cell___) / 4;
const int BODYWORD = sizeof(__body___) / 4;
const int MULTIPOLEWORD = sizeof(__multipole___) / 4;
const uint8_t SENDBIT = 1;
const uint8_t RECEIVEBIT = 0;

inline int encryptMessage(uint16_t grainSize, uint8_t requestType, uint8_t level, uint8_t direction) {
	int tag = int(grainSize);
	tag<<=REQUESTSHIFT; tag|=requestType;        
  tag<<=LEVELSHIFT; tag|=level;
  tag<<=DIRECTIONSHIFT; tag|= direction;
  return tag;  
}

inline void decryptMessage(int tag, uint16_t& grainSize, uint8_t& requestType, uint8_t& level, uint8_t& direction) {
	direction = tag & DIRECTIONMASK;
	tag >>= DIRECTIONSHIFT;
	level = tag & LEVELMASK;
	tag >>= LEVELSHIFT;
	requestType = tag & REQUESTMASK;	
	tag >>= REQUESTSHIFT;
	grainSize = tag & GRAINMASK;
}

inline uint16_t getGrainSize(int const& tag) {
	return ((tag >> (REQUESTSHIFT + DIRECTIONSHIFT + LEVELSHIFT)) & GRAINMASK);
}

inline uint8_t getMessageType(int const& tag) {
	return ((tag >> (DIRECTIONSHIFT + LEVELSHIFT)) & REQUESTMASK);
}

inline uint8_t getMessageLevel(int const& tag) {
	return ((tag >> DIRECTIONSHIFT) & LEVELMASK);
}

inline uint8_t getMessageDirection(int const& tag) {
	return (tag & DIRECTIONMASK);
}


class CellCounter {
	static size_t cellQueue;
public:
	inline static void incrementQueue() {
		lock_mutex;
		cellQueue++;
		unlock_mutex;
	}

	inline static void decrementQueue() {
		lock_mutex;
		cellQueue--;
		unlock_mutex;
	}

	inline static size_t getQueueSize() {
		return cellQueue;
	}
};
size_t CellCounter::cellQueue = 0;





#if SEND_MULTIPOLES

template <typename CellVec>
inline Multipoles getMultipoles(CellVec const& cells) {
	size_t const& size = cells.size();
	Multipoles multipoles(size);
	for (size_t i = 0; i < size; ++i)
	{
		multipoles[i].X 		= cells[i].X;
		multipoles[i].M 		= cells[i].M;
		multipoles[i].ICELL = cells[i].ICELL;
		multipoles[i].NCHILD = cells[i].NCHILD;
		multipoles[i].R = cells[i].R;
		multipoles[i].LEVEL = cells[i].LEVEL;
		multipoles[i].NBODY = cells[i].NBODY;
		multipoles[i].SCALE = cells[i].SCALE;

	}	
	return multipoles;
}

template <typename Iter>
inline Multipoles getMultipoles(Iter const& begin, Iter const& end) {
	size_t const& size = end - begin; 
	assert(size > 0);
	Multipoles multipoles(size);
	for (size_t i = 0; i < size; ++i) {	
		auto&& current =  begin + i;
		multipoles[i].X 		= current->X;
		multipoles[i].M 		= current->M;
		multipoles[i].ICELL = current->ICELL;
		multipoles[i].NCHILD = current->NCHILD;
		multipoles[i].R = current->R;
		multipoles[i].LEVEL = current->LEVEL;
		multipoles[i].NBODY = current->NBODY;
		multipoles[i].SCALE = current->SCALE;
	}	
	return multipoles;
}
#if DFS
template <typename VecType, typename Iter, typename GrainType>
inline void updateChildMultipoles(VecType& multipoles, Iter const& C0, Cell const& C, GrainType const& grainSize, size_t& index) {
	auto&& begin = C0 + C.ICHILD;
  auto&& end   = C0 + C.ICHILD + C.NCHILD;
  typedef typename VecType::value_type val_type;
	size_t const& size = C.NCHILD;
	for (size_t i = 0; i < size; ++i) {	
		auto&& current =  begin + i;
		val_type m;
		m.X 		= current->X;
		m.M 		= current->M;
		m.ICELL = current->ICELL;
		m.NCHILD = current->NCHILD;
		m.R = current->R;				
		m.LEVEL = current->LEVEL;
		m.NBODY = current->NBODY;	
		m.SCALE = current->SCALE;	
		multipoles.push_back(m);
		index++;
	}	
	for(auto&& cc = begin; cc<end; ++cc) 
		if(index < grainSize) updateChildMultipoles(multipoles,C0,*cc,grainSize,index);				
}
#else
template <typename VecType, typename Iter, typename GrainType>
inline void updateChildMultipoles(VecType& multipoles, Iter const& C0, Cell const& root, GrainType const& grainSize, size_t& index) {		
	std::queue<Cell> cell_queue;
	auto begin = C0 + root.ICHILD;  
	auto size = root.NCHILD;
	for (size_t i = 0; i < size; ++i) {							
	  auto&& current =  begin + i;	  
		cell_queue.push(*current);				
	}	
	typedef typename VecType::value_type val_type;
	while(cell_queue.size() > 0) {
		auto C = cell_queue.front();
		cell_queue.pop();
		begin = C0 + C.ICHILD;  	
		size = C.NCHILD;		
		val_type m;
		m.X 		= C.X;
		m.M 		= C.M;
		m.ICELL = C.ICELL;
		m.NCHILD = C.NCHILD;
		m.R = C.R;				
		m.LEVEL = C.LEVEL;
		m.NBODY = C.NBODY;				
		m.SCALE = C.SCALE;			
		m.ICHILD = index;
		m.IPARENT = C.IPARENT;
		multipoles.push_back(m);	
		index++;
		if(index < grainSize) {
			for (size_t i = 0; i < size; ++i) {							
			  auto&& current =  begin + i;			  
				cell_queue.push(*current);	
				cell_queue.back().IPARENT = m.ICHILD;	
			}	
		}
	}	
}
#endif

#endif

inline void setMinimumCircle(Cell& cell, real_t default_r) {
	if(cell.NBODY > 1){
		MB mb(3,cell.BODY,cell.BODY+cell.NBODY);
		cell.R = sqrt(mb.squared_radius());
		auto&& cen = mb.center();
		for (int i = 0; i < 3; ++i)
		{
			cell.X[i] = cen[i];
		}
	} else {
		cell.X = cell.BODY->X;
		cell.R = default_r*0.5;
	}
}

#endif