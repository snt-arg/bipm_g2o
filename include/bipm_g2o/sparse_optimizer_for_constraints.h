#ifndef G2O_GRAPH_OPTIMIZER_FOR_CONSTRAINTS_H_
#define G2O_GRAPH_OPTIMIZER_FOR_CONSTRAINTS_H_

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/g2o_core_api.h"
#include "termination_creteria.h"


namespace g2o {

//forwad declaration
template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeIneq;


template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeEq;

// Specialized optimizer inheriting from SparseOptimizer
template <typename Derived>
class G2O_CORE_API SparseOptimizerForConstraints : public SparseOptimizer {

 public: 
  SparseOptimizerForConstraints(); 
  virtual ~SparseOptimizerForConstraints(); // Virtual destructor


   template <int D, typename E, typename... VertexTypes>
   bool addEdgeIneq(BaseFixedSizedEdgeIneq<D, E, VertexTypes...>* e);
     /**
   * add the inequlaity edge to the optimizer and the set of inequality edges
   * @param e: the edge to be added
   * @returns false if somethings goes wrong
   */

   template <int D, typename E, typename... VertexTypes>
   bool addEdgeEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...>* e);
     /**
   * add the inequlaity edge to the optimizer and the set of inequality edges
   * @param e: the edge to be added
   * @returns false if somethings goes wrong
   */

  virtual void resetLagrangeMultiplierEq() = 0;

  
   virtual bool initializeOptimization(int level = 0) override;
   /**
   * Initializes the structures for optimizing the whole graph.
   * Before calling it be sure to invoke marginalized() and fixed() to the
   * vertices you want to include in the schur complement or to set as fixed
   * during the optimization.
   * @param level: is the level (in multilevel optimization)
   * @returns false if somethings goes wrong
   */

   virtual bool initializeOptimization(HyperGraph::EdgeSet& eset) override;
   /**
   * Initializes the structures for optimizing a portion of the graph specified
   * by a subset of vertices. Before calling it be sure to invoke marginalized()
   * and fixed() to the vertices you want to include in the schur complement or
   * to set as fixed during the optimization.
   * @param vset: the subgraph to be optimized.
   * @param level: is the level (in multilevel optimization)
   * @returns false if somethings goes wrong
   */
   virtual bool initializeOptimization(HyperGraph::VertexSet& vset,
                                      int level = 0) override;

   /**
   * Initializes the structures for optimizing the whole graph.
   * Before calling it be sure to invoke marginalized() and fixed() to the
   * vertices you want to include in the schur complement or to set as fixed
   * during the optimization.
   * @param level: is the level (in multilevel optimization)
   * @returns false if somethings goes wrong
   */


   virtual bool updateInitialization(HyperGraph::VertexSet& vset,
                                    HyperGraph::EdgeSet& eset) override;

   /** 
   * Propagates an initial guess from the vertex specified as origin.
   * It should be called after initializeOptimization(...), as it relies on the
   * _activeVertices/_edges structures. It constructs a set of trees starting
   * from the nodes in the graph which are fixed and eligible for preforming
   * optimization. The trees are constructed by utlizing a cost-function
   * specified.
   * @param costFunction: the cost function used
   * @patam maxDistance: the distance where to stop the search
   */


    /**
    * Verifies the termination condition for the BIPM algorithm.
    * The termination condition is based on the Newton decrement or the gradient norm.
    * @returns true if the termination condition is satisfied, false otherwise.
    */
    
   virtual double backtrackingAlgorithm([[maybe_unused]] const double* update) {return _alphaBacktracking[0];};

  

   //Termination
    
    TerminationCriteria& terminationCriteria(); 

    virtual void setConvergenceCriterion(int criterion);
    virtual bool verifyConvergence(double epsilon);
    virtual bool verifyEqFeasibility(EdgeContainer& edgeVec, double epsilon);
    virtual bool verifyIneqFeasibility(EdgeContainer& edgeVec, double epsilon);
    void setWarmStartLagrangeMultiplierFlag(bool flag);
    void setInnerIterationsMax(int maxInnerIterations);
    void setLagrangeMultiplierInitial(double lagrangeMultiplierInitial);
    void setEpsilonConstraint(double epsilon);
    void setEpsilonConvergence(double epsilon);
    void setAlgorithmEq(std::string algType);


    bool getWarmStartLagrangeMultiplierFlag();
    int getInnerIterationsMax();
    double getLagrangeMultiplierInitial();
    void setAlphaBacktrackingValue(double alpha);


 

   //virtual void defineAlphaBacktracking() {      }

   
  
  
 
 
 protected: 
   std::unordered_map<void*, std::function<void()>> updateMultipliersEqFunctionMap;
   TerminationCriteria _terminationCriteria;  // New termination criteria object
   EdgeContainer _activeEdgesIneq;       ///< sorted according to EdgeIDCompare
   std::set<OptimizableGraph::Edge*> _edgeIneqSet; // Set of inequality edges
   EdgeContainer _activeEdgesEq;       ///< sorted according to EdgeIDCompare
   std::set<OptimizableGraph::Edge*> _edgeEqSet; // Set of equality edges
   double _lagrange_multiplier_initial; // Initial value for the Lagrangian vertex
   int _num_inner_iterations_max; // Maximum number of inner iterations for the AL algorithm
   bool _warm_start_lagrange_multiplier_flag = false; // Flag to indicate if the Lagrangian multiplier should be warmed up
  // solver parameter
   std::vector<double> _alphaBacktracking;
   double _alpha = 0.9;

 
   void defineAlphaBacktracking(double alpha);
  };




}  // namespace g2o

#include "sparse_optimizer_for_constraints.hpp"


#endif  // G2O_GRAPH_OPTIMIZER_FOR_CONSTRAINTS_H_

 