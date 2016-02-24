#ifndef HOT_PARTITIONER
#define HOT_PARTITIONER
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <set>

namespace exafmm{
class HOTPartition{
	typedef hilbert_t KeyType;
	typedef std::vector<KeyType> RangeType;
	typedef std::vector<size_t> WeightVecType;
	typedef typename RangeType::iterator 	RangeIteratorType;
	typedef std::pair<KeyType,KeyType> BoundsType; 
	typedef double WorkloadType;
	BoundsType localBounds;																												 //!< local partition boundaries 
	BoundsType globalBounds;																											 //!< global domain boundaries 
	RangeType range;																															 //!< vector of equidistant ranges as per key
	std::vector<size_t> workload;																									 //!< workload for nonuniform sampling 
	const int commrank;																														 //!< rank of MPI communicator  
	const int commsize;																														 //!< size of MPI communicator 
	RangeIteratorType rangeBegin;																						 			 //!< begin iterator of range
	RangeIteratorType rangeEnd;																		 				 				 //!< end iterator of range
	const int64_t hilbertDistance = 750000000000;                 								 //!< key acceptance criteria 
	const float sampleRate = 0.001f;																							 //!< sample rate for workload sampling 
	const float lamda  = 0.05f;

	void partitionSort(RangeType& keys, hilbert_t lbound, hilbert_t rbound, RangeType&rq, hilbert_t pivot=0, size_t depth = 0) {
		RangeType right;
		RangeType left;		
		for(auto&& key : keys) {
			if(key >= pivot)
				right.push_back(key);
			else left.push_back(key);
		}
		size_t sendSize[2];
		sendSize[0] = right.size();
		sendSize[1] = left.size();
		size_t size_d[2];
		MPI_Allreduce(sendSize,size_d, 2, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);// Reduce domain SUM
		float diff = 0.0f;
		bool recurse = true;		
		auto max = std::max(size_d[0],size_d[1]);
		auto min = std::min(size_d[0],size_d[1]);
		diff = (max - min)/float(max);		
		auto templ = lbound;
		auto tempr = rbound;
		while(diff > lamda) {
			if(size_d[0] > size_d[1]) {
				templ = pivot;
				pivot += (tempr - pivot)>>1;
			}
			else if(size_d[1] > size_d[0]) {
				tempr = pivot;
				pivot -= (pivot - templ)>>1;
			}
			right.clear();
			left.clear();			
			for(auto&& key : keys) {
				if(key >= pivot)
					right.push_back(key);
				else left.push_back(key);
			}
			sendSize[0] = right.size();
			sendSize[1] = left.size();
			MPI_Allreduce(sendSize,size_d, 2, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);// Reduce domain SUM
			max = std::max(size_d[0],size_d[1]);
		  min = std::min(size_d[0],size_d[1]);
			diff = (max - min)/float(max);	
		}
		if(recurse) {		
			rq.push_back(pivot);
			if(2<<depth < commsize) {
				partitionSort(left ,lbound,pivot,rq,lbound + ((pivot-lbound)>>1),depth+1);
				partitionSort(right,pivot,rbound,rq,pivot  + ((rbound-pivot)>>1),depth+1);	
			}
		}				
	}

	void weightedPartitionSort(Bodies& bodies, hilbert_t lbound, hilbert_t rbound, RangeType&rq, hilbert_t pivot = 0, size_t depth = 0) {
		Bodies right;
		Bodies left;		
		double leftWeight = 0.0;
		double rightWeight = 0.0;
		for(auto&& body : bodies) {
			if(body.ICELL >= pivot){
				right.push_back(body);
				rightWeight += body.WEIGHT;
			}				
			else {
			 	left.push_back(body);
			 	leftWeight += body.WEIGHT;
			}
		}
		double weightSum[2];
		weightSum[0] = rightWeight;
		weightSum[1] = leftWeight;
		double size_d[2];
		MPI_Allreduce(weightSum,size_d, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// Reduce domain SUM		
		bool recurse = true;		
		auto max = std::max(size_d[0],size_d[1]);
		auto min = std::min(size_d[0],size_d[1]);
		auto diff = (max - min)/(max);		
		auto templ = lbound;
		auto tempr = rbound;		
		while(diff > lamda) {
			if(size_d[0] > size_d[1]) {
				templ = pivot;
				auto&& offset = (tempr - pivot)>>1;
				pivot += offset;
			}
			else if(size_d[1] > size_d[0]) {
				tempr = pivot;
				auto&& offset = (pivot - templ)>>1;				
				pivot -= offset;
			}
			right.clear();
			left.clear();	
			leftWeight = 0;
			rightWeight = 0;		
			for(auto&& body : bodies) {
				if(body.ICELL >= pivot){
					right.push_back(body);
					rightWeight += body.WEIGHT;
				}				
				else {
				 	left.push_back(body);
				 	leftWeight += body.WEIGHT;
				}
			}
			weightSum[0] = rightWeight;
			weightSum[1] = leftWeight;
			MPI_Allreduce(weightSum,size_d, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// Reduce domain SUM
			max = std::max(size_d[0],size_d[1]);
		  min = std::min(size_d[0],size_d[1]);
			diff = (max - min)/(max);	
		}
		if(recurse) {		
			rq.push_back(pivot);
			if(2<<depth < commsize) {
				weightedPartitionSort(left ,lbound,pivot,rq,lbound + ((pivot-lbound)>>1),depth+1);
				weightedPartitionSort(right,pivot,rbound,rq,pivot  + ((rbound-pivot)>>1),depth+1);	
			}
		}				
	}

	//! Partition uniformly distributed particles using HOT 
	template <typename BodyType>
	void partitionUniform(BodyType& bodies, bool generateRange = true, bool sort = true) {		
		if(commsize <=1) return;																											// partitioning not needed in this case
		if(generateRange) {
			range.resize(commsize+1);
			rangeBegin = range.begin();
			rangeEnd 	 = range.end();
			auto&& min = globalBounds.first;																							// global min Hilbert
			auto&& max = globalBounds.second;																							// global max Hilbert
			KeyType val(min);																														  // initial range value
			KeyType step((max-min)/(commsize+1));																				  // step size
			std::generate(rangeBegin+1,rangeEnd-1,[&]{ return val+=step; });			  			// generate equi-distant range 
			range[0]				= min;
			range[commsize] = max;
		}
		if(sort) {
			logger::startTimer("Sort Hilbert");																  					// start HOT partitioning timer
			std::sort(bodies.begin(),bodies.end(),[&](Body a,Body b){return (a.ICELL<b.ICELL);});
			logger::stopTimer("Sort Hilbert");																						// start HOT partitioning timer
		}
		for(auto&& body: bodies) {																										// loop over bodies																								
			body.IRANK = getUniformPartition(body.ICELL);								  					    // get partition
			assert(body.IRANK >=0 && body.IRANK <commsize);															// make sure rank is within size
		}
		//assert(verifyParitioning(bodies));
	}


	bool verifyParitioning(Bodies const& bodies) const {
		std::vector<size_t> counts(commsize);
		size_t prev = 0;
		for (size_t i = 0; i < bodies.size(); ++i){
			auto&& body = bodies[i];
			if(size_t(body.IRANK) < prev)
				return false;
			prev = body.IRANK;
			counts[prev]++;
		}
		size_t sum_of_elements(0);
		std::for_each(counts.rbegin(), counts.rend(), [&](size_t n) { sum_of_elements += n; });
		return (sum_of_elements == bodies.size());		
	}

	void updateRangeParameters(RangeType const& rg) {
		range = rg;
		rangeBegin = range.begin();
	  rangeEnd = range.end();
	  globalBounds.first = range[0];
	  globalBounds.second = range[commsize];
	}

	//! Get the rank based on the key
	//! works only for particles partitioned uniformly 
	inline int getUniformPartition(KeyType const& key) {
	  auto&& low = std::lower_bound(rangeBegin, rangeEnd, key); 										// get iterator of first element that is >= input key		
		if(low > rangeBegin) {
			return (low - rangeBegin) - 1;
		}
		return 0;
	}

public:				
	HOTPartition(BoundsType local, BoundsType global, int rank, int size, uint32_t numBodies): 
		localBounds(local),globalBounds(global),workload(numBodies,0),
		commrank(rank),commsize(size) { }

	void partitionSort(Bodies& bodies) { 		
		if(commsize == 1) return;	  
	  auto lbound = globalBounds.first;
		auto rbound = globalBounds.second;
	  RangeType rq;		
		rq.reserve(commsize+1);
		rq.push_back(lbound); rq.push_back(rbound);
		auto startPivot = lbound + ((rbound - lbound)>>1);					
		RangeType keys(bodies.size());
		int index = 0;
		for(auto&& body : bodies) 
			keys[index++] = body.ICELL;				
		partitionSort(keys,lbound,rbound,rq,startPivot);				
		std::sort(rq.begin(),rq.end());
		updateRangeParameters(rq);		
		for(auto&& body: bodies) {																										// loop over bodies																								
			body.IRANK = getUniformPartition(body.ICELL);								  					    // get partition
			assert(body.IRANK >=0 && body.IRANK <commsize);															// make sure rank is within size
		}
	  logger::startTimer("Sort");
		std::sort(bodies.begin(),bodies.end(),[&](Body a,Body b){return (a.IRANK<b.IRANK);});
		logger::stopTimer("Sort");
	}

	void weightedPartitionSort(Bodies& bodies) {		
		auto lbound = globalBounds.first;
		auto rbound = globalBounds.second;
	  RangeType rq;		
		rq.reserve(commsize+1);
		rq.push_back(lbound); rq.push_back(rbound);
		auto startPivot = lbound + ((rbound - lbound)>>1);							
		weightedPartitionSort(bodies,lbound,rbound,rq,startPivot);
		std::sort(rq.begin(),rq.end());
		updateRangeParameters(rq);		
		for(auto&& body: bodies) {																									// loop over bodies																								
			body.IRANK = getUniformPartition(body.ICELL);							  					    // get partition
			assert(body.IRANK >=0 && body.IRANK <commsize);														// make sure rank is within size
		}
	  logger::startTimer("Sort");
		std::sort(bodies.begin(),bodies.end(),[&](Body a,Body b){return (a.IRANK<b.IRANK);});
		logger::stopTimer("Sort");		
	}

	template <class VecType>
	void traverse(C_iter C0, C_iter C, VecType& wq ){
		for(auto cc = C0+C->ICHILD; cc < C0+C->ICHILD+C->NCHILD;++cc) {
			cc->WEIGHT += C0->WEIGHT;
			if(cc->NCHILD == 0){
				for(int i =0; i < cc->NBODY; ++i)
					wq[i] = cc->BODY->WEIGHT;
			}
			traverse(C0,cc,wq);
		}
	}
	
	void balanceTraverse(C_iter C0, C_iter C, B_iter B0, WorkloadType* wbalanced, WorkloadType* wcomm){
		for(auto cc = C0+C->ICHILD; cc < C0+C->ICHILD+C->NCHILD;++cc) {						
			if(cc->NBODY <= 1024 && cc->WEIGHT > 0) {
				for (size_t j = 1; j < commsize + 1; ++j) {	    		
	    		if(wcomm[commrank] <= wbalanced[j]) {
	    			wcomm[commrank] += cc->WEIGHT;
	    			for(auto&& b = B0 + cc->IBODY; b < B0 + cc->IBODY + cc->NBODY; ++b) {
	    				b->IRANK = j-1;
	    			}
	    			break;
	    		}
				}	    			
		  }		
		  else balanceTraverse(C0, cc,B0,wbalanced,wcomm);
		}
	}
	
	// void balanceTraverse(C_iter C0, C_iter C, B_iter B0, WorkloadType* wbalanced, WorkloadType* wcomm){
	// 	for(auto cc = C0+C->ICHILD; cc < C0+C->ICHILD+C->NCHILD;++cc) {									
	// 	  if(cc->NCHILD > 0) balanceTraverse(C0, cc,B0, wbalanced,wcomm);
	// 	  if(cc->WEIGHT > 0) {
	// 		  for (size_t j = 1; j < commsize + 1; ++j) {	    		
	// 	    	if(wcomm[commrank] <= wbalanced[j]) {
	// 	    		wcomm[commrank] += cc->WEIGHT;
	// 	    		for(auto&& b = B0 + cc->IBODY; b < B0 + cc->IBODY + cc->NBODY; ++b) {
	// 	    			b->IRANK = j-1;
	// 	    		}
	// 	    		break;
	// 	    	}
	// 			}
	// 		}	
	// 	}
	// }



	template <typename CellType, typename BodyType>	
	void migrateWork(CellType&& cells, BodyType&& bodies, double commtime) {						
		WorkloadType wq[commsize];
		MPI_Allgather(&commtime, 1, MPI_DOUBLE, wq, 1, MPI_DOUBLE, MPI_COMM_WORLD);// Gather all domain bounds	
		if(bodies.size() > 0) {
	    WorkloadType wglobal(0);
	    WorkloadType wbalanced[commsize + 1];
	    for (size_t i = 0; i < commsize; ++i) {
	    	wglobal+=wq[i];
	    }
	    for (size_t i = 0; i < commsize + 1; ++i) {
	    	wbalanced[i] = (i * wglobal)/commsize;
	    }
	    WorkloadType wg[commsize + 1];
	    wg[0] = 0;
	    for (size_t i = 1; i < commsize + 1; ++i) 
	    	wg[i] = wq[i-1] + wg[i-1];
	    balanceTraverse(cells.begin(),cells.begin(), bodies.begin(),wbalanced,wg);

	    logger::startTimer("Sort");
		  std::sort(bodies.begin(),bodies.end(),[&](Body a,Body b){return (a.IRANK<b.IRANK);});
		  logger::stopTimer("Sort");
  	}
	}

	template <typename BodyType>	
	void migrateWork(BodyType& bodies) {		
		uint32_t bodyCount = bodies.size();  
		std::vector<WorkloadType> workPerLeaf(bodyCount, 0);		
		for(int i = 0;	 i < bodyCount; ++i) workPerLeaf[i] = bodies[i].WEIGHT;			
		WorkloadType wlocal(0);
		std::for_each(workPerLeaf.rbegin(), workPerLeaf.rend(), [&](WorkloadType load) { wlocal += load; });		
		WorkloadType wq[commsize];
		MPI_Allgather(&wlocal, 1, MPI_DOUBLE, wq, 1, MPI_DOUBLE, MPI_COMM_WORLD);// Gather all domain bounds	
		if(bodyCount > 0) {
	    WorkloadType wglobal(0);
	    WorkloadType wbalanced[commsize + 1];
	    for (size_t i = 0; i < commsize; ++i) {
	    	wglobal+=wq[i];
	    }
	    for (size_t i = 0; i < commsize + 1; ++i) {
	    	wbalanced[i] = (i * wglobal)/commsize;
	    }
	    WorkloadType wg[commsize + 1];
	    wg[0] = 0;
	    for (size_t i = 1; i < commsize + 1; ++i) 
	    	wg[i] = wq[i-1] + wg[i-1];
			int zeroWeightCount = 0;

	   	for (size_t i = 0; i < bodyCount; ++i) {
	    	for (size_t j = 1; j < commsize + 1; ++j) {	    		
	    		if(wg[commrank] <= wbalanced[j]) {
	    			if(workPerLeaf[i] > 0) {	    			
							wg[commrank] += workPerLeaf[i];
	    				for(int k = i - zeroWeightCount; k<=i; ++k) bodies[k].IRANK = j-1;	    	
	    				zeroWeightCount = 0;
	    			} else zeroWeightCount++;	    				    				    			
	    			break;
	    		}
	    	}	
	    }

	    logger::startTimer("Sort");
		  std::sort(bodies.begin(),bodies.end(),[&](Body a,Body b){return (a.IRANK<b.IRANK);});
		  logger::stopTimer("Sort");
  	}
	}

	template <typename CellType, typename BodyType>
	void loadBalance(CellType& cells, BodyType& bodies){		
		size_t bodyCount = bodies.size();  
    WorkloadType wlocal(0);
    RangeType rq(commsize + 1,0);
		WorkloadType wq[commsize];
		assert(bodyCount > 1);
		logger::startTimer("Sort");
		std::sort(bodies.begin(),bodies.end(),[&](Body a,Body b){return (a.ICELL<b.ICELL);});
		logger::stopTimer("Sort");		
		std::vector<WorkloadType> workPerLeaf(bodyCount, 0);		
		for(int i = 0; i < bodyCount; ++i)
			workPerLeaf[i] = bodies[i].WEIGHT;
		std::for_each(workPerLeaf.rbegin(), workPerLeaf.rend(), [&](WorkloadType load) { wlocal += load; });
    MPI_Allgather(&wlocal, 1, MPI_DOUBLE, wq, 1, MPI_DOUBLE, MPI_COMM_WORLD);// Gather all domain bounds
    if(bodyCount > 0) {
	    WorkloadType wglobal(0);
	    WorkloadType wbalanced[commsize + 1];
	    for (size_t i = 0; i < commsize; ++i) {
	    	wglobal+=wq[i];
	    }
	    for (size_t i = 0; i < commsize + 1; ++i) {
	    	wbalanced[i] = (i * wglobal)/commsize;
	    }
	    WorkloadType wg[commsize + 1];
	    wg[0] = 0;
	    for (size_t i = 1; i < commsize + 1; ++i) 
	    	wg[i] = wq[i-1] + wg[i-1];

	   	for (size_t i = 0; i < bodyCount; ++i) {
	    	for (size_t j = 1; j < commsize + 1; ++j) {
	    		if(wg[commrank] <= wbalanced[j]) {
	    			wg[commrank] += workPerLeaf[i];
	    			//bodies[i].IRANK = j-1;
	    			rq[j] = bodies[i].ICELL;
	    			break;
	    		}
	    	}	
	    }
  	}
  //	bool v = verifyParitioning(bodies);
  	rq[0] = 0;
  	rq[commsize] = globalBounds.second;
    std::vector<hilbert_t> rq_max(commsize + 1);
		MPI_Allreduce(RAWPTR(rq),RAWPTR(rq_max), commsize + 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmin
		updateRangeParameters(rq_max);	 
		if(bodyCount > 0) {
		  for(auto&& body: bodies) {																										// loop over bodies																								
				body.IRANK = getUniformPartition(body.ICELL);								  					    // get partition
				assert(body.IRANK >=0 && body.IRANK <commsize);															// make sure rank is within size
			}
	  }		
		//bodyCount = bodies.size();
		//logger::stopTimer("Loadbalance HTree");																				// start HOT partitioning timer			
	}

	template <typename BodyType>
	void partitionMasterSlave(BodyType& bodies) {	
		auto&& begin = bodies.begin();																								// bodies begin iterator 
		auto&& end   = bodies.end();																									// bodies end iterator 
		logger::startTimer("Sort Hilbert");																						// start sort timer		
		std::sort(begin,end,[&](Body a, Body b){return (a.ICELL<b.ICELL);});				  // sort based on Hilbert key
		logger::stopTimer("Sort Hilbert");																		        // stop sort timer
		logger::startTimer("HOT Partitioning");																				// start HOT partitioning timer
		auto&& numBodies  = bodies.size();																					  // number of bodies
		if(numBodies < commsize) {
			std::cout<< "insufficient problem size";
			MPI_Abort(MPI_COMM_WORLD,0);
		}							
		int partitionSize = numBodies/commsize;
		int remainderValue= numBodies%commsize;		
		int currentRank   = 0;
		int takenRemainders = 0;		
		for(int i=0; i < numBodies; ++i) {
			bodies[i].IRANK = currentRank;
			if((i-takenRemainders+1)%partitionSize == 0) {
				int j = i+1;
				if(remainderValue > 0) {
					remainderValue--;
					takenRemainders++;							
					if(j < numBodies) bodies[j].IRANK = currentRank;					
					currentRank++;
					i = j;
				}
				else if(j < numBodies) currentRank++;
			} 
		}				
		assert(currentRank == commsize - 1);
		logger::stopTimer("HOT Partitioning");																				// end HOT partitioning timer		
	}

	//! Partition non-uniformly distributed particles using HOT and workload sampling  
	template <typename BodyType>
	void partitionNonUniform(BodyType& bodies) {					
		auto&& begin = bodies.begin();																								// bodies begin iterator 
		auto&& end   = bodies.end();																									// bodies end iterator 
		logger::startTimer("Sort Hilbert");																						// start sort timer		
		std::sort(begin,end,[&](Body a, Body b){return (a.ICELL<b.ICELL);}); // sort based on Hilbert key
		logger::stopTimer("Sort Hilbert");																		        // stop sort timer
		logger::startTimer("HOT Partitioning");																				// start HOT partitioning timer
		auto&& numBodies  = bodies.size();																					  // number of bodies
		if(numBodies < commsize) {
			std::cout<< "insufficient problem size";
			MPI_Abort(MPI_COMM_WORLD,0);
		}						
		auto&& compareLower = [&,this](Body body1,Body body2){return (std::abs(body1.ICELL-body2.ICELL) > hilbertDistance);}; // lambda expression for lower bounds of neighbors
		auto&& compareUpper = [&,this](Body body1,Body body2){return (std::abs(body1.ICELL-body2.ICELL) < hilbertDistance);}; // lambda expression for upper bounds of neighbors 
		size_t numSamples = numBodies * sampleRate;																		// number of samples
		size_t stepSize = 1;																													// step size between samples
		if(numSamples > 1) stepSize = numBodies / numSamples;													// evaluate step size 
		int count, start;																														  // markers for interaction list 																	
		for(int i = 0; i < numBodies; i+=stepSize)																		// start workload sampling loop
		{
			auto&& body  = bodies[i];																											// take body alias
			auto&& pivot = begin + (i + 1);																								// pivot point
			count = 0; start = 0;
			auto&& leftBound = std::lower_bound(begin,pivot,body,compareLower);				  	// get first element that is <= Hilbert distance 			
			start = leftBound - begin;																								  	// mark the start of interaction list
			auto&& rightBound = std::lower_bound(pivot,  end,body,compareUpper);					// get first element that is >= Hilbert distance 
			if(rightBound == end && leftBound != pivot) count = pivot - leftBound;				// in case the there no right neighbors 
			else if(rightBound!=end) count = rightBound - leftBound;	  									// otherwise

			logger::startTimer("Workload Accum");																			 	  // start workload Accum 
			for (int j = start; j < start + count; ++j)		
				workload[j]++;
			logger::stopTimer("Workload Accum",0);														 			   	  // end workload Accum 
		}		
		logger::printTime("Workload Accum");																					// print time for workload accumulation 
		std::partial_sum(workload.begin(),workload.end(),workload.begin());						// build prefix sum for workload list
		double loadSize = workload[numBodies - 1] / double(commsize);		  						// calculate expected workload per partition
		size_t oldWorkload = 0;																												// previous work load value
		size_t currentRank = 0;																												// current rank value
		if(loadSize > 1.0) {
			for (size_t i = 0; i < numBodies; ++i) {
				bodies[i].IRANK = currentRank;
				if((workload[i] - oldWorkload) > loadSize && currentRank + 1 < commsize ) {															// if workload exceeds limits for current rank
					currentRank++;
					oldWorkload = workload[i];
				}
			}
		} else {			
			int partitionSize = numBodies/commsize;
			int remainderValue= numBodies%commsize;		
			int takenRemainders = 0;		
			for(int i=0; i < numBodies; ++i) {
				bodies[i].IRANK = currentRank;
				if((i-takenRemainders+1)%partitionSize == 0) {
					int j = i+1;
					if(remainderValue > 0) {
						remainderValue--;
						takenRemainders++;							
						if(j < numBodies) bodies[j].IRANK = currentRank;					
						currentRank++;
						i = j;
					}
					else if(j < numBodies) currentRank++;
				} 
			}					
		}	
		assert(verifyParitioning(bodies));
		assert(currentRank == commsize - 1);
		logger::stopTimer("HOT Partitioning");																				// end HOT partitioning timer		
	}

};
}
#endif 