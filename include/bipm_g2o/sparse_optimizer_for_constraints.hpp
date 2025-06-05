#include "g2o/core/optimization_algorithm.h"
#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/solver.h"
#include "g2o/stuff/logger.h"
#include "sparse_optimizer_for_constraints.h"
#include <Eigen/Core>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

#ifndef NDEBUG
namespace {
/**
 * tests whether there is a NaN in the array
 */
bool arrayHasNaN(const double *array, int size, int *nanIndex = 0) {
  for (int i = 0; i < size; ++i)
    if (std::isnan(array[i])) {
      if (nanIndex)
        *nanIndex = i;
      return true;
    }
  return false;
}
} // namespace
#endif

namespace g2o {
using namespace std;

// Default constructor
template <typename Derived>
SparseOptimizerForConstraints<Derived>::SparseOptimizerForConstraints() =
    default;

// Default destructor
template <typename Derived>
SparseOptimizerForConstraints<Derived>::~SparseOptimizerForConstraints() =
    default;

// Setter solver paramters
template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::initializeOptimization(
    HyperGraph::VertexSet &vset, int level) {
  if (edges().size() == 0) {
    G2O_WARN("Attempt to initialize an empty graph");
    return false;
  }
  preIteration(-1);
  bool workspaceAllocated = _jacobianWorkspace.allocate();
  (void)workspaceAllocated;
  assert(workspaceAllocated &&
         "Error while allocating memory for the Jacobians");
  clearIndexMapping();
  _activeVertices.clear();
  _activeVertices.reserve(vset.size());
  _activeEdges.clear();
  set<Edge *> auxEdgeSet; // temporary structure to avoid duplicates
  for (HyperGraph::VertexSet::iterator it = vset.begin(); it != vset.end();
       ++it) {
    OptimizableGraph::Vertex *v = (OptimizableGraph::Vertex *)*it;
    const OptimizableGraph::EdgeSet &vEdges = v->edges();
    // count if there are edges in that level. If not remove from the pool
    int levelEdges = 0;
    for (OptimizableGraph::EdgeSet::const_iterator it = vEdges.begin();
         it != vEdges.end(); ++it) {
      OptimizableGraph::Edge *e =
          reinterpret_cast<OptimizableGraph::Edge *>(*it);
      if (level < 0 || e->level() == level) {
        bool allVerticesOK = true;
        for (vector<HyperGraph::Vertex *>::const_iterator vit =
                 e->vertices().begin();
             vit != e->vertices().end(); ++vit) {
          if (vset.find(*vit) == vset.end()) {
            allVerticesOK = false;
            break;
          }
        }
        if (allVerticesOK && !e->allVerticesFixed()) {
          auxEdgeSet.insert(e);
          levelEdges++;
        }
      }
    }
    if (levelEdges) {
      _activeVertices.push_back(v);

      // test for NANs in the current estimate if we are debugging
#ifndef NDEBUG
      int estimateDim = v->estimateDimension();
      if (estimateDim > 0) {
        VectorX estimateData(estimateDim);
        if (v->getEstimateData(estimateData.data()) == true) {
          int k;
          bool hasNan = arrayHasNaN(estimateData.data(), estimateDim, &k);
          if (hasNan)
            G2O_WARN("Vertex {} contains a nan entry at index {}", v->id(), k);
        }
      }
#endif
    }
  }

  _activeEdges.reserve(auxEdgeSet.size());
  _activeEdgesIneq.reserve(_edgeIneqSet.size());
  for (set<Edge *>::iterator it = auxEdgeSet.begin(); it != auxEdgeSet.end();
       ++it) {
    _activeEdges.push_back(*it);
    //-----------------------------------------
    // add the active edge (it) to the active inequality edges if it is the
    // inequality set
    if (_edgeIneqSet.find(*it) != _edgeIneqSet.end())
      // Element found, add it to _activeEdgesIneq
      _activeEdgesIneq.push_back(
          *it); // Dereference the iterator to get the element
    // add the active edge (it) to the active equality edges if it is the
    // equality set
    if (_edgeEqSet.find(*it) != _edgeEqSet.end())
      // Element found, add it to _activeEdgesIneq
      _activeEdgesEq.push_back(
          *it); // Dereference the iterator to get the element
    //-----------------------------------------
  }
  sortVectorContainers();
  bool indexMappingStatus = buildIndexMapping(_activeVertices);
  postIteration(-1);
  return indexMappingStatus;
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::initializeOptimization(int level) {
  HyperGraph::VertexSet vset;
  for (VertexIDMap::iterator it = vertices().begin(); it != vertices().end();
       ++it)
    vset.insert(it->second);
  return initializeOptimization(vset, level);
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::updateInitialization(
    HyperGraph::VertexSet &vset, HyperGraph::EdgeSet &eset) {
  std::vector<HyperGraph::Vertex *> newVertices;
  newVertices.reserve(vset.size());
  _activeVertices.reserve(_activeVertices.size() + vset.size());
  _activeEdgesIneq.reserve(_edgeIneqSet.size());

  _activeEdges.reserve(_activeEdges.size() + eset.size());
  for (HyperGraph::EdgeSet::iterator it = eset.begin(); it != eset.end();
       ++it) {
    OptimizableGraph::Edge *e = static_cast<OptimizableGraph::Edge *>(*it);
    if (!e->allVerticesFixed()) {
      _activeEdges.push_back(e);
      //-----------------------------------------
      // add the inequality edges to the active edges if they are in the
      // inequality set
      if (_edgeIneqSet.find(e) != _edgeIneqSet.end())
        // Element found, add it to _activeEdgesIneq
        _activeEdgesIneq.push_back(
            e); // Dereference the iterator to get the element
      //-----------------------------------------
    }
  }
  // update the index mapping
  size_t next = _ivMap.size();
  for (HyperGraph::VertexSet::iterator it = vset.begin(); it != vset.end();
       ++it) {
    OptimizableGraph::Vertex *v = static_cast<OptimizableGraph::Vertex *>(*it);
    if (!v->fixed()) {
      if (!v->marginalized()) {
        v->setHessianIndex(next);
        _ivMap.push_back(v);
        newVertices.push_back(v);
        _activeVertices.push_back(v);
        next++;
      } else // not supported right now
        abort();
    } else {
      v->setHessianIndex(-1);
    }
  }

  if (newVertices.size() != vset.size()) {
    G2O_ERROR("something went wrong, size mismatch {} != {}", vset.size(),
              newVertices.size());
  }
  return _algorithm->updateStructure(newVertices, eset);
}



template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::initializeOptimization(
    HyperGraph::EdgeSet &eset) {
  preIteration(-1);
  bool workspaceAllocated = _jacobianWorkspace.allocate();
  (void)workspaceAllocated;
  assert(workspaceAllocated &&
         "Error while allocating memory for the Jacobians");
  clearIndexMapping();
  _activeVertices.clear();
  _activeEdges.clear();
  _activeEdges.reserve(eset.size());
  //-----------------------------------------
  _activeEdgesIneq.reserve(_edgeIneqSet.size());
  //-----------------------------------------
  set<Vertex *> auxVertexSet; // temporary structure to avoid duplicates
  for (HyperGraph::EdgeSet::iterator it = eset.begin(); it != eset.end();
       ++it) {
    OptimizableGraph::Edge *e = (OptimizableGraph::Edge *)(*it);
    if (e->numUndefinedVertices())
      continue;
    for (vector<HyperGraph::Vertex *>::const_iterator vit =
             e->vertices().begin();
         vit != e->vertices().end(); ++vit) {
      auxVertexSet.insert(static_cast<OptimizableGraph::Vertex *>(*vit));
    }
    _activeEdges.push_back(reinterpret_cast<OptimizableGraph::Edge *>(*it));
    //-----------------------------------------
    // add the inequality edges to the active edges if they are in the
    // inequality set
    if (_edgeIneqSet.find(static_cast<OptimizableGraph::Edge *>(*it)) !=
        _edgeIneqSet.end())
      // Element found, add it to _activeEdgesIneq
      _activeEdgesIneq.push_back(static_cast<OptimizableGraph::Edge *>(
          *it)); // Dereference the iterator to get the element

    //-----------------------------------------
  }

  _activeVertices.reserve(auxVertexSet.size());
  for (set<Vertex *>::iterator it = auxVertexSet.begin();
       it != auxVertexSet.end(); ++it)
    _activeVertices.push_back(*it);

  sortVectorContainers();
  bool indexMappingStatus = buildIndexMapping(_activeVertices);
  postIteration(-1);
  return indexMappingStatus;
}

// Implementation of `addEdgeIneq`
template <typename Derived>
template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerForConstraints<Derived>::addEdgeIneq(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e) {
  return static_cast<Derived *>(this)->addEdgeIneqImpl(e);
}

// Implementation of `addEdgeEq`
template <typename Derived>
template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerForConstraints<Derived>::addEdgeEq(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e) {
  return static_cast<Derived *>(this)->addEdgeEqImpl(e);
}

// Termination functions
template <typename Derived>
TerminationCriteria &
SparseOptimizerForConstraints<Derived>::terminationCriteria() {
  return _terminationCriteria;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setConvergenceCriterion(
    int criterion) {
  /*
  the termination criterion is based on the Newton decrement, the gradient norm,
  or the update norm where GradientNorm =0
  */
  _terminationCriteria.setConvergenceCriterion(criterion);
};

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::verifyConvergence(
    double epsilon) {
  OptimizationAlgorithmWithHessian *algorithm =
      static_cast<OptimizationAlgorithmWithHessian *>(_algorithm);

  // Get the gradient vector b and the vector size
  const double *b = algorithm->solver().b();
  const double *update = algorithm->solver().x();
  size_t xSize = algorithm->solver().vectorSize();

  // Map the b vector to an Eigen vector
  Eigen::Map<const Eigen::VectorXd> bVec(b, xSize);
  Eigen::Map<const Eigen::VectorXd> updateVec(update, xSize);

  if (terminationCriteria().verifyConvergence(bVec, updateVec,
                                                     epsilon)) {
    return true;
  }
  return false;
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::verifyIneqFeasibility(
    EdgeContainer &edgeVec, double epsilon) {
  return terminationCriteria().verifyIneqFeasibility(edgeVec, epsilon);
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::verifyEqFeasibility(
    EdgeContainer &edgeVec, double epsilon) {
  return terminationCriteria().verifyEqFeasibility(edgeVec, epsilon);
}

template<typename Derived>
void SparseOptimizerForConstraints<Derived>::setWarmStartLagrangeMultiplierFlag(bool flag) {
      _warm_start_lagrange_multiplier_flag = flag;
}

template<typename Derived>
void SparseOptimizerForConstraints<Derived>::setInnerIterationsMax(int maxInnerIterations) {
  _num_inner_iterations_max = maxInnerIterations;
}

template<typename Derived>
void SparseOptimizerForConstraints<Derived>::setLagrangeMultiplierInitial(double lagrangeMultiplierInitial) {
  _lagrange_multiplier_initial = lagrangeMultiplierInitial;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setEpsilonConstraint(
    double epsilon_feas) {
  _terminationCriteria.epsilon_constraint = epsilon_feas;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setEpsilonConvergence(double epsilon) {
  _terminationCriteria.epsilon_convergence = epsilon;
}

template<typename Derived>
bool SparseOptimizerForConstraints<Derived>::getWarmStartLagrangeMultiplierFlag() {
    return  _warm_start_lagrange_multiplier_flag ;
}

template<typename Derived>
int SparseOptimizerForConstraints<Derived>::getInnerIterationsMax() {
  return _num_inner_iterations_max;
}


template<typename Derived>
double SparseOptimizerForConstraints<Derived>::getLagrangeMultiplierInitial() {
  return _lagrange_multiplier_initial;
}


} // namespace g2o
