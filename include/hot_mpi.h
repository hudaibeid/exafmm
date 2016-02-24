#ifndef HOT_MPI
#define HOT_MPI

#include "logger.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include "cell_mutex.h"
#include "cell_dispatcher.h"
#include <queue>

using namespace std;

//!< TODO:: make the CELL type generic
namespace exafmm{
//!< class to handle MPI communication for the HOT
template <typename BodiesType, typename ReceiveCellsType>
class HOTMPI {	
public:
  typedef typename ReceiveCellsType::value_type CellType;
	typedef std::vector<std::array<real_t,3>>     GlobalBoundsVector;        //!< Type of Vector of Global Bounds
  typedef std::vector<hilbert_t>                GlobalHilbertBoundsVector; //!< Type of Global Hilbert Bounds Vector 
  typedef std::pair<hilbert_t,hilbert_t>        KeyPair;                   //!< Type of Key Pair
  typedef std::vector<KeyPair>                  KeyPairVector;             //!< Type of Key Pair Vector
  typedef std::vector<size_t>                   RankList;                  //!< Type of Rank list
  typedef std::vector<KeyPairVector>            BoundsBreadthMap;
  typedef std::vector<std::unordered_map<hilbert_t,BodiesType>> VectorBodyMap;
  typedef std::vector<std::unordered_map<hilbert_t,ReceiveCellsType>> VectorCellMap;
  typedef std::vector<std::unordered_map<hilbert_t,CellType>> VectorMap;
  typedef std::vector<VectorMap>                VectorVectorMap;
  typedef std::vector<VectorCellMap>            VectorVectorCellMap;
  typedef std::vector<VectorBodyMap>            VectorVectorBodyMap;
  typedef DispatcherWrapper<Cells,TreeType>     Dispatcher;
  typedef std::vector<int> VecInt;

  
  const int mpirank;                                       			           //!< Rank of MPI communicator
  const int mpisize;                                                   	   //!< Size of MPI communicator
  const int order;                                                         //!< Order of HOT tree
  const uint16_t grainSize;                                                //!< Granularity of communication for cells
private: 
  KeyPair globalBounds;																					  	     	 //!< Global Key Bounds
  KeyPair localBounds;																							     	 //!< Local Key Bounds
  VecInt sendBodyCount;                                                    //!< Send count
  VecInt sendBodyDispl;                                                    //!< Send displacement
  VecInt recvBodyCount;                                                    //!< Receive count
  VecInt recvBodyDispl;                                                    //!< Receive displacement
  GlobalBoundsVector allBoundsXmin;                                        //!< Array for local Xmin for all ranks
	GlobalBoundsVector allBoundsXmax;                                        //!< Array for local Xmax for all ranks
  GlobalHilbertBoundsVector allHilbertBoundsMin;                           //!< Vector of Hilbert bounds min of all ranks
  GlobalHilbertBoundsVector allHilbertBoundsMax;                           //!< Vector of Hilbert bounds max of all ranks
  VectorVectorMap cellsMap;                                                //!< A cache for the requested keys 
  VectorVectorCellMap childrenMap;
  VectorVectorBodyMap bodyMap;                        
  Dispatcher* dispatcher;
  BoundsBreadthMap boundsPerBreadth;                                       //!< Vector of Bounds per Level  
  int receivingTask;                                                       //!< Set current receiving thread
  
	
  inline void prepareBodyCounts(BodiesType const& bodies, VecInt& counts) {
		auto&& compare1 = [](int const& val,Body const& body){return (val < body.IRANK) ;};
		auto&& compare2 = [](Body const& body,int const& val){return (body.IRANK < val) ;};
		auto&& begin = bodies.begin();
		auto&& end 	 = bodies.end();
		for (int i = 0; i < mpisize; ++i) {																	   //!<bodies are already grouped based on rank (sorted ICELL)
			 auto&& lower = std::lower_bound(begin,end,i,compare2);
			 auto&& upper = std::upper_bound(begin,end,i,compare1);
			 counts[i] = upper - lower;
		}
	}

	inline void syncBodyCounts() {
		MPI_Alltoall(RAWPTR(sendBodyCount), 1, MPI_INT,                        //!<Communicate send count to get receive count
      RAWPTR(recvBodyCount), 1, MPI_INT, MPI_COMM_WORLD);
    sendBodyDispl[0] = recvBodyDispl[0] = 0;                               //!<Initialize send/receive displacements
    for (int irank=0; irank<mpisize-1; ++irank) {                          //!<Loop over ranks
      sendBodyDispl[irank+1] = sendBodyDispl[irank] + sendBodyCount[irank];//!< Set send displacement
      recvBodyDispl[irank+1] = recvBodyDispl[irank] + recvBodyCount[irank];//!<  Set receive displacement
    }                                                                      //!< End loop over ranks		
	}

	inline BodiesType syncBodies(BodiesType& bodies) {
    auto size = recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1];
		BodiesType recvBodies(size);
		assert( (sizeof(bodies[0]) & 3) == 0 );                     
    int word = sizeof(bodies[0]) / 4;                           
    for (int irank=0; irank<mpisize; ++irank) {                 
      sendBodyCount[irank] *= word;                             
      sendBodyDispl[irank] *= word;                             
      recvBodyCount[irank] *= word;                             
      recvBodyDispl[irank] *= word;                             
    }                                                           
    MPI_Alltoallv((int*)&bodies[0], RAWPTR(sendBodyCount),	    
    RAWPTR(sendBodyDispl), MPI_INT,  (int*)&recvBodies[0], 
    RAWPTR(recvBodyCount), RAWPTR(recvBodyDispl), MPI_INT, 
    MPI_COMM_WORLD);
    for (int irank=0; irank<mpisize; ++irank) {                 
      sendBodyCount[irank] /= word;                             
      sendBodyDispl[irank] /= word;                             
      recvBodyCount[irank] /= word;                             
      recvBodyDispl[irank] /= word;                             
    }                                                           
    return recvBodies;
	}

	//! merge the sorted bodies to form a big sorted list
	inline BodiesType mergeBodies(BodiesType& bodies) {
    auto dataSize = recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1];
    if(bodies.size() == sendBodyCount[mpirank] && dataSize == sendBodyCount[mpirank])   // Nothing to send to other ranks
      return bodies;                                            // return current value
    int irank = 0;
    while(recvBodyCount[irank] == 0 && irank < mpisize) 
      irank++;
    if(irank < mpisize) { 
      auto&& compare   = [&](Body body1, Body body2){return (body1.ICELL<body2.ICELL);};
      auto&& recvBegin = bodies.begin();   
      BodiesType sortedBodies(recvBegin + recvBodyDispl[irank], recvBegin + recvBodyDispl[irank] + recvBodyCount[irank]);
      irank++;
      for (; irank<mpisize; ++irank) {
      	size_t currentSize = recvBodyDispl[irank] + recvBodyCount[irank];
      	BodiesType output(currentSize);
      	std::merge(sortedBodies.begin(), sortedBodies.end(), recvBegin + recvBodyDispl[irank], recvBegin + currentSize, output.begin(),compare);
      	sortedBodies = std::move(output);
      }
      localBounds.first = sortedBodies.front().ICELL;
      localBounds.second = sortedBodies.back().ICELL;
      return sortedBodies; 
    } else {
      return bodies;
    }
	}


	inline BodiesType asyncBodies(BodiesType const& bodies){
    auto dataSize = recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1];
    if(bodies.size() == sendBodyCount[mpirank] && dataSize == sendBodyCount[mpirank])  
      return bodies;                                            
		BodiesType recvBodies(dataSize);
		assert( (sizeof(bodies[0]) & 3) == 0 );                     
    int word = sizeof(bodies[0]) / 4;                           
    for (int irank=0; irank<mpisize; ++irank) {                 
      sendBodyCount[irank] *= word;                             
      sendBodyDispl[irank] *= word;                             
      recvBodyCount[irank] *= word;                             
      recvBodyDispl[irank] *= word;                             
    }                                                           
    auto&& sendBuff = (int*)&bodies[0];													
    auto&& recvBuff = (int*)&recvBodies[0];											
    size_t sendSize = 0;
    size_t recvSize = 0;
    MPI_Request rreq[mpisize-1];                                
    MPI_Request sreq[mpisize-1];                                
    MPI_Status rstatus[mpisize-1];															
    MPI_Status sstatus[mpisize-1];															

		for (int irank=0; irank<mpisize; ++irank) {                 
	    if(irank != mpirank) {
        if(recvBodyCount[irank] > 0){   		
    		  MPI_Irecv(recvBuff + recvBodyDispl[irank], 
    		  recvBodyCount[irank],MPI_INT, irank, irank, MPI_COMM_WORLD, &rreq[recvSize]);
          recvSize++;
        }
			}
    }
    for (int irank=0; irank<mpisize; ++irank) {                 
	    if(irank != mpirank) {	    	
    		if(sendBodyCount[irank] > 0) {
          MPI_Isend(sendBuff + sendBodyDispl[irank], 
    		  sendBodyCount[irank],MPI_INT, irank, mpirank, MPI_COMM_WORLD, &sreq[sendSize]); 
          sendSize++;
        } 		    
			}
    }
    for (int irank=0; irank<mpisize; ++irank) {                
      sendBodyCount[irank] /= word;                            
      sendBodyDispl[irank] /= word;                            
      recvBodyCount[irank] /= word;                            
      recvBodyDispl[irank] /= word;                            
    }     
    auto&& localBuffer = bodies.begin() + sendBodyDispl[mpirank];
    std::copy(localBuffer, localBuffer + sendBodyCount[mpirank],recvBodies.begin() + recvBodyDispl[mpirank]);		
    MPI_Waitall(sendSize,&sreq[0], &sstatus[0]);
		MPI_Waitall(recvSize,&rreq[0], &rstatus[0]);
    return recvBodies;
	}



public:	
	HOTMPI(int rank, int size, size_t _order, uint16_t granularity)
		:mpirank(rank)
    ,mpisize(size),order(_order), grainSize(granularity), sendBodyCount(size)
		,sendBodyDispl(size),recvBodyCount(size),recvBodyDispl(size)
    ,allBoundsXmin(size),allBoundsXmax(size),allHilbertBoundsMin(size*(_order+1))
    ,allHilbertBoundsMax(size*(_order+1)),cellsMap(size,VectorMap(_order+1)), childrenMap(size,VectorCellMap(_order+1))
    ,bodyMap(size,VectorBodyMap(_order+1)) { }

  ~HOTMPI() { }	

	KeyPair allReduceHilbertBounds(KeyPair& local) {		
    logger::startTimer("Hilbert bounds");
		localBounds = local;
		hilbert_t min, max;
		MPI_Allreduce(&local.first,  &min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);// Reduce domain Xmin
    MPI_Allreduce(&local.second, &max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmax
    globalBounds = std::make_pair(min,max);
    logger::stopTimer("Hilbert bounds");
    return globalBounds;
	}

  template <typename RangeType>
  inline RangeType loadBalanceCounts(BodiesType const& bodies, RangeType range, int balancedCount) {
    const int halfCount = balancedCount >> 1;
    VecInt sendCounts(mpisize-1);
    VecInt recvCounts(mpisize-1);
    prepareBodyCounts(bodies,sendCounts);
    MPI_Allreduce(RAWPTR(sendCounts),RAWPTR(recvCounts),mpisize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    auto&& boundary = range[0];
    int sizeRange = mpisize;
    size_t current = 0;
    std::vector<int> erasable;
    erasable.reserve(mpisize);
    auto result = range;
    for(int i = 1; i <= sizeRange; ++i) {
      assert(current<mpisize);
      if(recvCounts[current] > balancedCount) {
        hilbert_t spiltSize = recvCounts[current]/hilbert_t(balancedCount);
        if(recvCounts[current]%balancedCount)
          spiltSize++;
        if(spiltSize > 1) {
          hilbert_t rangeSize = result[i] - boundary; 
          hilbert_t stepSize  = rangeSize/spiltSize;
          hilbert_t innerRange= boundary + stepSize;
          while((result[i] - innerRange) > hilbert_t(balancedCount)){
            result.insert(result.begin() + i, innerRange);
            innerRange+=stepSize; i++; sizeRange++;
            break;
          }
        }
      } else if(recvCounts[current] <= halfCount) {
        erasable.push_back(i);
      }
      boundary = result[i];
      current++;
    }
    size_t erasedCount = 0;
    size_t erasableSize = erasable.size();
    while(result.size() > mpisize + 1){
      if(erasedCount<erasableSize) {
        result.erase(result.begin()+erasable[erasedCount]-erasedCount);
        erasedCount++;
      } else { 
        result.erase(result.begin() + (result.size() - 2));
      }
    }
    return result;
  }

  //! Allgather bounds from all ranks
  void allgatherBounds(Bounds bounds) {
    float Xmin[3], Xmax[3];
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      Xmin[d] = bounds.Xmin[d];                                 //  Convert Xmin to float
      Xmax[d] = bounds.Xmax[d];                                 //  Convert Xmax to float
    }                                                           // End loop over dimensions
    MPI_Allgather(Xmin, 3, MPI_FLOAT, RAWPTR(allBoundsXmin), 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
    MPI_Allgather(Xmax, 3, MPI_FLOAT, RAWPTR(allBoundsXmax), 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
  }

  //! Allgather bounds from all ranks  
  void allgatherHilbertBounds(GlobalHilbertBoundsVector localMin, GlobalHilbertBoundsVector localMax) {
    const int depth = order + 1;
    MPI_Allgather(RAWPTR(localMin) , depth, MPI_UNSIGNED_LONG_LONG, RAWPTR(allHilbertBoundsMin), depth, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);// Gather all domain bounds
    MPI_Allgather(RAWPTR(localMax) , depth, MPI_UNSIGNED_LONG_LONG, RAWPTR(allHilbertBoundsMax), depth, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);// Gather all domain bounds
    KeyPairVector keyPairs(mpisize);    
    boundsPerBreadth.reserve(depth);
    for(int i=0; i<depth;++i) {
      KeyPairVector vec(mpisize);
      for(int rank=0; rank<mpisize;++rank) {
        vec[rank].first = allHilbertBoundsMin[rank * depth + i];
        vec[rank].second= allHilbertBoundsMax[rank * depth + i];
      }
      boundsPerBreadth.push_back(vec);
    }
  }

  template <typename HOT>
  void initDispatcher(Cells const& cells, HOT const& indexer) {
    dispatcher = new Dispatcher(mpirank, mpisize, cells, indexer);      
  }

  void startDispatcher(){
    dispatcher->callDispatcher();
  }

  void flushDispatcher() {
    MPI_Request request;    
    int tag = encryptMessage(1,FLUSHTAG,0,SENDBIT);                
    char null;   
    for(int i=0; i < mpisize; ++i)
      if(i!=mpirank)
        MPI_Isend(&null, 1, MPI_CHAR,i,tag,MPI_COMM_WORLD,&request);
  }
   
  void stopDispatcher() {    
    delete dispatcher;
  }
  void finalizeDispatcher(){
     dispatcher->finalizeDispatcher();
     delete dispatcher;
  }

  BodiesType getBodies(hilbert_t key, size_t level, int rank, int requestType, double& commtime) {
    assert(rank!=mpirank);
#if THREADED_REMOTE    
    lock_mutex;
#endif        
    commtime = 0.0;
    BodiesType recvData;
    auto&& bodyRankMap = bodyMap[rank][level];
    if(requestType == BODYTAG && bodyRankMap.find(key) ==  bodyRankMap.end()) { 
      MPI_Request request;
      int tag = encryptMessage(1,requestType,level,SENDBIT);
#if CALC_COM_COMP      
      logger::startTimer("Communication");
#endif    
      MPI_Isend(&key,1,MPI_UNSIGNED_LONG_LONG, rank, tag, MPI_COMM_WORLD,&request); 
      MPI_Status status;         
      int recvRank = rank;
      int receivedTag;      
      TOGGLEDIRECTION(tag)

      int ready = 0;
      while(!ready) {
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&ready,&status); 
        if(ready) {
          if(tag != status.MPI_TAG || rank != status.MPI_SOURCE) {           
            dispatcher->processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);          
            ready = 0;
          }
        }
      }    
      receivedTag = status.MPI_TAG;  
      int responseType = getMessageType(receivedTag);
      int recvCount = 0;
      assert(responseType == BODYTAG || responseType == NULLTAG);
      if(responseType == BODYTAG) { 
        MPI_Get_count(&status, MPI_INT, &recvCount);
        int bodyCount = recvCount/BODYWORD;
        recvData.resize(bodyCount);
        MPI_Recv(RAWPTR(recvData), recvCount, MPI_INT, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if CALC_COM_COMP      
        commtime = logger::stopTimer("Communication",0);
#endif    
        bodyRankMap.emplace(key,recvData);
      } else {
          std::cout << "bodies not found for cell "<<key << std::endl;
          char null;
          MPI_Recv(&null, 1, MPI_CHAR, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
#if CALC_COM_COMP      
        commtime = logger::stopTimer("Communication",0);
#endif    
      }
    }
    else if(requestType == BODYTAG) {
      recvData = bodyRankMap[key];
    }
#if THREADED_REMOTE    
    unlock_mutex;
#endif
    return recvData;
  }

#if DFS
  template <typename MapperType, typename CellsType>
  void appendDataToCellMapper(MapperType& map, CellsType const& cells, int& index, int id){    
    auto parent = cells[id];
    if(index < cells.size()) {    
      int level = parent.LEVEL;
      id = index;
      if(parent.NCHILD > 0) {        
        map[level].emplace(parent.ICELL,CellsType(cells.begin()+index,cells.begin()+index+parent.NCHILD));        
        index += parent.NCHILD;        
        for(size_t i = 0; i < parent.NCHILD; ++i) {
          appendDataToCellMapper(map, cells,index, id++);
        }
      }
    }    
  }
#else 
  template <typename MapperType, typename CellsType>
  void appendDataToCellMapper(MapperType& map, CellsType const& cells, int rootLevel, int nchild, uint64_t rootID){            
    map[rootLevel].emplace(rootID,CellsType(cells.begin(),cells.begin()+nchild));                        
    for(int i = nchild; i < cells.size(); ++i)
    {
      auto&& currentCell = cells[i];
      auto&& parentCell = cells[currentCell.IPARENT];      
      auto&& parentLevel = parentCell.LEVEL;
      auto&& parentID = parentCell.ICELL;
      if(map[parentLevel].find(parentID) == map[parentLevel].end()) {
        map[parentLevel][parentID] = CellsType(); 
        map[parentLevel][parentID].reserve(parentCell.NCHILD); 
      }
      map[parentLevel][parentID].emplace_back(currentCell);  
    }
  }

#endif
  ReceiveCellsType getCell(hilbert_t key, int nchild, size_t level, int rank, int requestType, double& commtime) {
    assert(rank!=mpirank);    
#if THREADED_REMOTE        
    lock_mutex;
#endif
    commtime = 0.0;
    ReceiveCellsType recvData;
    auto&& cellRankMap = cellsMap[rank][level];
    auto&& childRankMap = childrenMap[rank][level];    
    if((requestType == CHILDCELLTAG && childRankMap.find(key) == childRankMap.end()) ||
       (requestType == CELLTAG      &&  cellRankMap.find(key) ==  cellRankMap.end()) || 
        requestType ==LEVELTAG) { 
      MPI_Request request;      
      assert(requestType <= MAXTAG);      
      int tag = encryptMessage(grainSize,requestType,level,SENDBIT);
#if CALC_COM_COMP      
        logger::startTimer("Communication");
#endif 
      MPI_Isend(&key,1,MPI_UNSIGNED_LONG_LONG, rank, tag, MPI_COMM_WORLD,&request); 
      MPI_Status status;          
      int recvRank = rank;
      int receivedTag;       
      TOGGLEDIRECTION(tag)   
      int ready = 0;
      while(!ready) {
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&ready,&status); 
        if(ready) {
          if(tag != status.MPI_TAG || rank != status.MPI_SOURCE) {           
            dispatcher->processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);          
            ready = 0;
          }
        }
      }
      receivedTag = status.MPI_TAG;                
      int responseType = getMessageType(receivedTag);
      int recvCount = 0;
      assert(responseType != BODYTAG);
      if(responseType == CHILDCELLTAG || responseType == CELLTAG || responseType == LEVELTAG) { 
        MPI_Get_count(&status, MPI_INT, &recvCount);
#if SEND_MULTIPOLES        
        int cellCount = recvCount/MULTIPOLEWORD;
#else 
        int cellCount = recvCount/CELLWORD;
#endif         
        recvData.resize(cellCount);
        MPI_Recv(RAWPTR(recvData), recvCount, MPI_INT, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if CALC_COM_COMP      
        commtime = logger::stopTimer("Communication",0);
#endif          
        if(responseType == CELLTAG) cellRankMap.emplace(key, recvData[0]);
        else if(responseType == CHILDCELLTAG)  { 
          if(grainSize > 1) {
#if DFS            
            childRankMap.emplace(key,ReceiveCellsType(recvData.begin(), recvData.begin()+nchild));  
            int index = nchild;
            for (int i = 0; i < nchild; ++i) appendDataToCellMapper(childrenMap[rank],recvData,index,i);  
            recvData = childRankMap[key];  
#else
            appendDataToCellMapper(childrenMap[rank],recvData, level,nchild, key);        
            recvData = childRankMap[key];  
#endif            
          }
          else childRankMap.emplace(key,recvData);
        }
      } else if(responseType == NULLTAG) {
          char null;
          MPI_Recv(&null, 1, MPI_CHAR, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
#if CALC_COM_COMP      
          commtime = logger::stopTimer("Communication",0);
#endif                      
      }
    }
    else if(requestType == CELLTAG) {
      recvData.push_back(cellRankMap[key]);
    }
    else if(requestType == CHILDCELLTAG) {
      recvData = childRankMap[key];
    }
#if THREADED_REMOTE    
    unlock_mutex;
#endif    
    return recvData;
  }

  void clearCellCache(int rank) {
    cellsMap[rank].clear();
    bodyMap[rank].clear();  
    childrenMap[rank].clear();  
  }

  //! Global synchronization for particles in all partitions 
	BodiesType asyncGlobalCommunication(BodiesType& bodies, int printTime = 1, bool sort = false) {		
		if(mpisize == 1) return bodies;
    prepareBodyCounts(bodies,sendBodyCount);
		logger::startTimer("syncing counts");
		syncBodyCounts();
		logger::stopTimer("syncing counts", printTime);
		logger::startTimer("syncing bodies");
		auto&& recvBodies = asyncBodies(bodies);
    logger::stopTimer("syncing bodies", printTime);
    if(sort && recvBodies.size() > 0) {
      logger::startTimer("merging bodies");
  		recvBodies = mergeBodies(recvBodies);			
      logger::stopTimer("merging bodies", printTime);
    }
		return std::move(recvBodies);                                                 // move because RVO is not guaranteed here
	}

  //! Send bodies to next rank (round robin)
  void shiftBodies(Bodies & bodies) {
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    const int word = sizeof(bodies[0]) / 4;                     // Word size of body structure
    const int isend = (mpirank + 1          ) % mpisize;        // Send to next rank (wrap around)
    const int irecv = (mpirank - 1 + mpisize) % mpisize;        // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, receive request handles

    MPI_Isend(&oldSize, 1, MPI_INT, irecv, 0, MPI_COMM_WORLD, &sreq);// Send current number of bodies
    MPI_Irecv(&newSize, 1, MPI_INT, isend, 0, MPI_COMM_WORLD, &rreq);// Receive new number of bodies
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                         // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                         // Wait for receive to complete

    BodiesType recvBodies(newSize);                                 // Resize buffer to new number of bodies
    MPI_Isend((int*)&bodies[0], oldSize*word, MPI_INT, irecv,   // Send bodies to next rank
              1, MPI_COMM_WORLD, &sreq);
    MPI_Irecv((int*)&recvBodies[0], newSize*word, MPI_INT, isend,// Receive bodies from previous rank
              1, MPI_COMM_WORLD, &rreq);
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                         // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                         // Wait for receive to complete
    bodies = recvBodies;                                        // Copy bodies from buffer
  }



};
}

#endif