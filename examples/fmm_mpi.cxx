#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "dataset.h"
#include "logger.h"
#include "up_down_pass.h"
#include "verify.h"
#include "hilbert_wrapper.h"
#include "hot_partitioner.h"
#include "traversal.h"
#include "build_tree_tbb.h"
#include "hot_mpi.h"

#define DIRECT 1


using namespace exafmm;
using namespace std;


#if SEND_MULTIPOLES
typedef HOTMPI<Bodies,Multipoles> BasicMPI;
#else 
typedef HOTMPI<Bodies,Cells> BasicMPI;
#endif
typedef Traversal<BasicMPI> TreeTraversal;


int main(int argc, char ** argv) {
  const real_t eps2 = 0.0;
  const real_t cycle = 2 * M_PI;
  Args args(argc, argv);
  Dataset data;
  BaseMPI baseMPI;
  Verify verify;
  BoundBox boundBox(args.nspawn);  
  int rank = baseMPI.mpirank;
  int size = baseMPI.mpisize;  
  args.verbose &= rank == 0;
  args.numBodies /= baseMPI.mpisize;  
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters"); 
  args.print(logger::stringLength, P, std::cout);  
  double commtime = 0.0;
  num_threads(args.threads);  
  Cells cells,jcells;
  Bodies jbodies;
  kernel::eps2 = 0.0;
#if EXAFMM_HELMHOLTZ
  kernel::wavek = complex_t(10.,1.) / real_t(2 * M_PI);
#endif
  kernel::setup();
  auto&& bodies = data.initBodies(args.numBodies, args.distribution,rank, size); 
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");    
    auto&& localBounds  = boundBox.getBounds(bodies);
    auto&& globalBounds = baseMPI.allreduceBounds(localBounds);
    uint32_t depth; 
    auto&& hilbertBounds = generateHilbertKey(bodies,globalBounds,depth);         
    BasicMPI* hotMPI = new BasicMPI(rank, size, depth,args.grain);    
    auto&& globalHilbertBounds = hotMPI->allReduceHilbertBounds(hilbertBounds);  
    HOTPartition partitioner(hilbertBounds,globalHilbertBounds,rank,size,args.numBodies);
    logger::stopTimer("Total FMM",0);    
    logger::startTimer("Partition");
    logger::startTimer("Load-balance");
    if(args.balance == 0 || t == 0) {
      logger::stopTimer("Load-balance",0);
      partitioner.partitionSort(bodies);
      bodies = hotMPI->asyncGlobalCommunication(bodies,1,false);
      logger::stopTimer("Partition");
    }
    if(args.balance != 0 && t > 0) {
      logger::stopTimer("Partition",0);      
      partitioner.migrateWork(bodies);
      bodies = hotMPI->asyncGlobalCommunication(bodies,1,false);
      logger::stopTimer("Load-balance");
    }
    logger::startTimer("Total FMM");  
    localBounds = boundBox.getBounds(bodies);
    BuildTree buildTree(args.ncrit,depth);
    cells = buildTree(bodies,localBounds);  
    localBounds = boundBox.getBounds(cells, localBounds);
    UpDownPass upDownPass(args.theta, args.useRmax, false);
    upDownPass.upwardPass(cells);
    TreeTraversal traversal(args.nspawn, 0, hotMPI);    
    traversal.traverse(cells, cells, cycle, true, args.mutual);
    traversal.dualTreeTraversalRemote(cells,rank,size,BuildTree::indexer,upDownPass);        

#if CALC_COM_COMP
    logger::printTime("Communication");
#endif
#if DIRECT
    Bodies buffer;
    if(args.repeat > 1) buffer.assign(bodies.begin(),bodies.end());
    logger::printTitle("MPI direct sum");
    jbodies = Bodies(bodies.begin(),bodies.end());
    const int numTargets = 100;
    data.sampleBodies(bodies, numTargets);
    Bodies bodies2(bodies.begin(),bodies.end());
    data.initTarget(bodies);
    logger::startTimer("Total Direct");
    for (int i=0; i<baseMPI.mpisize; i++) {
      if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << baseMPI.mpisize << std::endl;      
      hotMPI->shiftBodies(jbodies);      
      traversal.direct(bodies, jbodies, cycle);      
    }
    traversal.normalize(bodies);
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
    logger::stopTimer("Total Direct");
    logger::resetTimer("Total Direct");    
    double potDif = verify.getDifScalar(bodies, bodies2);
    double potNrm = verify.getNrmScalar(bodies);
    double accDif = verify.getDifVector(bodies, bodies2);
    double accNrm = verify.getNrmVector(bodies);
    double potDifGlob, potNrmGlob, accDifGlob, accNrmGlob;
    MPI_Reduce(&potDif, &potDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potNrm, &potNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accDif, &accDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accNrm, &accNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0) {
      logger::printTitle("FMM vs. direct");
      verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
      verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
    }
    if(args.repeat > 1) { 
      traversal.printTraversalData();
      bodies = buffer;
      data.initTarget(bodies);   
    } else {
      traversal.printTraversalData();
    } 

#else 
      traversal.printTraversalData();      
      if(args.repeat > 1) data.initTarget(bodies); 
#endif  
#if CALC_COM_COMP    
    commtime = logger::getTime("Communication");
#endif

#if WRITE_TIME              
    logger::writeTime(rank);
#endif
    logger::zeroTimer("Traverse Remote");
    logger::zeroTimer("Total FMM");
    logger::zeroTimer("Traverse");
#if CALC_COM_COMP    
    logger::zeroTimer("Communication");
#endif    
    logger::zeroTimer("Partition");
    logger::zeroTimer("Load-balance");
    buildTree.printTreeData(cells);
  }

  if(rank == 0) {
    std::ofstream configFile("config.dat");                 // Open list log file
    args.print(logger::stringLength, P, configFile);  
    logger::logFixed("Communication Size",size,configFile);
#if HILBERT_CODE
    configFile<<"using Hilbert partitioning"<<std::endl;
#else 
    configFile<<"using Morton partitioning"<<std::endl;
#endif    
#if WEIGH_COM
    configFile<<"using communication & workload for balancing"<<std::endl;
#else
    configFile<<"using workload only for balancing"<<std::endl;
#endif
#if CALC_COM_COMP   
    configFile<<"calculating communication time"<<std::endl;        
#endif    
#if USE_WEIGHT
    configFile<<"using weights for partitioning"<<std::endl;        
#else
    configFile<<"not using weights for partitioning"<<std::endl;        
#endif    
#if COUNT_KERNEL
    configFile<<"using counting kernel"<<std::endl;            
#else
    configFile<<"not using counting kernel"<<std::endl;        
#endif    
    configFile.close();
  }
}
