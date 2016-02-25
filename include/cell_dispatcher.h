#ifndef CELL_DISPATCHER
#define CELL_DISPATCHER
#include <chrono>
#include <thread>
#include "logger.h"
#include "cell_mutex.h"
#include "lightweight_vector.h"
#define INTERNAL_THREADING 0
template <typename CellVec, typename HOT>
class DispatcherWrapper  {

  class CellDispatcher {    
    typedef typename CellVec::const_iterator iter_type;
    int rank;
    int size;
    TaskTrigger trigger;
    CellVec const& cells;
    HOT const& hot; 
    uint64_t timeout;
    int terminated;
    iter_type begin;
    size_t hitCount; 
    size_t index;
    lightweight_vector<Multipole> sendingBuffer;     

  public:    
    CellDispatcher(int _rank, int _size, TaskTrigger _taskRunning, CellVec const& _cells, HOT const& _hot, uint64_t _timeout)
      :rank(_rank),size(_size),trigger(_taskRunning),cells(_cells),hot(_hot), timeout(_timeout), terminated(0), begin(_cells.begin()), hitCount(0),index(0){ }    
    
    void operator()() {     
      int ready; 
      MPI_Status status;      
      while(GETTRIG(trigger) && terminated < (size - 1)) { 
        ready = 0;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&ready,&status); 
        if(ready) {
          tryReceive(status.MPI_TAG, status.MPI_SOURCE);          
        }
      }
      logger::logFixed("hit count", hitCount, std::cout);      
    }

    inline bool tryReceive(int tag, int source) {
      if((tag & DIRECTIONMASK) == RECEIVEBIT)
        return false; 
      hitCount++;
      auto msgType = getMessageType(tag);
      auto nullTag = encryptMessage(1,NULLTAG,0,RECEIVEBIT);     
      hilbert_t recvBuff;
      MPI_Request request;
      char null; 
      if (msgType==CELLTAG || msgType==CHILDCELLTAG || msgType==BODYTAG) { 
        int level = getMessageLevel(tag);
        MPI_Recv(&recvBuff,1,MPI_UNSIGNED_LONG_LONG, source, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);  
        TOGGLEDIRECTION(tag)
        auto&& levelData = hot.at(level);
        auto&& end       = levelData.end();
        auto&& location  = levelData.find(recvBuff);
        if(location != end) {
          auto&& cell = cells[location->second];
          if(msgType == CELLTAG) {
#if SEND_MULTIPOLES
            Multipole multipole = cell;
            MPI_Isend(&multipole,MULTIPOLEWORD,MPI_INT,source,tag,MPI_COMM_WORLD,&request);      
#else 
            MPI_Isend(&cell,CELLWORD,MPI_INT,source,tag,MPI_COMM_WORLD,&request);
#endif                          
            return true;                               
          } else if (msgType == CHILDCELLTAG && cell.NCHILD > 0) {
            auto grainSize = getGrainSize(tag);
#if SEND_MULTIPOLES
            size_t sendingSize = 0;
            updateChildMultipoles(sendingBuffer,begin, cell,grainSize,sendingSize);            
            MPI_Isend(&sendingBuffer[index],  MULTIPOLEWORD*sendingSize,MPI_INT,source,tag,MPI_COMM_WORLD,&request);
            index += sendingSize;            
            // assuming MPI will deallocate sendingBuffer[0]
            //delete[] sendingBuffer;
#else
            MPI_Isend(&cells[cell.ICHILD],      CELLWORD*cell.NCHILD,MPI_INT,source,tag,MPI_COMM_WORLD,&request);            
#endif             
            return true;
          } else if(msgType == BODYTAG && cell.NBODY > 0) {
            MPI_Isend(&(*cell.BODY), BODYWORD*cell.NBODY,MPI_INT,source,tag,MPI_COMM_WORLD,&request);
          } else {
            MPI_Isend(&null,1,MPI_CHAR,source,nullTag,MPI_COMM_WORLD,&request);
          }
        } else {              
          MPI_Isend(&null,1,MPI_CHAR,source,nullTag,MPI_COMM_WORLD,&request);                          
        }
      }
      else if(msgType == LEVELTAG) {
        int level = (tag >> DIRECTIONSHIFT) & LEVELMASK;
        MPI_Recv(&recvBuff,1,MPI_UNSIGNED_LONG_LONG, source, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        TOGGLEDIRECTION(tag) 
        auto&& levelData = hot.at(level);
        auto&& dataSize = levelData.size();
#if SEND_MULTIPOLES
        Multipoles sendCells;
#else 
        Cells sendCells;
#endif  
        sendCells.reserve(dataSize);      
        assert(dataSize > 0);        
        for(auto&& pairData: levelData) 
          sendCells.push_back(cells[pairData.second]);            
#if SEND_MULTIPOLES
        MPI_Isend(RAWPTR(sendCells), MULTIPOLEWORD*dataSize,MPI_INT,source,tag,MPI_COMM_WORLD,&request);
#else 
        MPI_Isend(RAWPTR(sendCells),      CELLWORD*dataSize,MPI_INT,source,tag,MPI_COMM_WORLD,&request);
#endif        
        return true;                     
      }
      else if(msgType == FLUSHTAG) {
        MPI_Recv(&null,1,MPI_CHAR, source, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        terminated++;
      }
      return false;
    }

    int get_rank() const {
      return rank;
    }
  };

public:
  TaskTrigger taskTrigger;
  CellDispatcher dispatcher;
  bool taskCreated;
#if INTERNAL_THREADING  
  std::thread thread;
#endif
  bool taskKilled;
  DispatcherWrapper(int _rank, int _size, CellVec const& _cells, HOT const& _hot)
      :taskTrigger(NEWTRIGGER),dispatcher(_rank,_size,taskTrigger,_cells,_hot,TIMEOUT)
      ,taskCreated(false), taskKilled(false){ 
      }

  ~DispatcherWrapper() {
    stopDispatcher();
  }

  void stopDispatcher(){
    if(!taskKilled) { 
      MPI_Barrier(MPI_COMM_WORLD);
      DELTRIGGER(taskTrigger);      
#if INTERNAL_THREADING      
      thread.join();
#endif      
      taskKilled = true;
    }
  }

  void waitDispatcher(){
#if INTERNAL_THREADING      
    if(!taskKilled) {
      thread.join();
      taskKilled = true;
    }
#endif  
  }

  void processIncomingMessage(int tag, int source) {
    dispatcher.tryReceive(tag,source);
  }

  void callDispatcher() {
    taskCreated = true;
    TRIGGER(taskTrigger); 
    dispatcher();
    taskCreated = false;
  }

  void untriggerDispatcher() {
    UNTRIGGER(taskTrigger);
  }

  void startDispatcher() {
    if(!taskCreated) {  //remove when multiple threads supported 
      TRIGGER(taskTrigger); 
#if INTERNAL_THREADING      
      thread = std::thread(dispatcher);
#else       
      dispatcher();
#endif
      taskCreated = true;                         
    }
  }
  void finalizeDispatcher() {
    #if INTERNAL_THREADING
      UNTRIGGER(taskTrigger);
      try {
        thread.join();
        std::cout<<"joined successfully"<<std::endl;
      }
      catch(...)
      {

      }
      TRIGGER(taskTrigger); 
      dispatcher();
      MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }
};






#endif 