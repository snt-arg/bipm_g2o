#ifndef G2O_SPARSE_OPTIMIZER_AL_HPP
#define G2O_SPARSE_OPTIMIZER_AL_HPP

#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"
#include "optimization_algorithm.h"
#include "solver.h"
#include "sparse_optimizer_al.h"
#include <Eigen/Core>
#include <iostream>
#include <string>
#include <vector>

namespace g2o {
using namespace std;

// Default constructor
SparseOptimizerAL::SparseOptimizerAL() {
  _lagrange_multiplier_initial = 0;
  _num_inner_iterations_max= 10;
}

// Default destructor
SparseOptimizerAL::~SparseOptimizerAL() = default;

// Setter solver paramters
void SparseOptimizerAL::setRhoInitial(double rhoInitial) {
  _rho_initial = rhoInitial;
}
void SparseOptimizerAL::setRhoUpdateFactor(double rhoUpdateFactor) {
  _rho_update_factor = rhoUpdateFactor;
}
void SparseOptimizerAL::setRhoMax(double rhoMax) { _rho_max = rhoMax; }

// getter for the solver parameters
double SparseOptimizerAL::getRhoInitial() { return _rho_initial; }
double SparseOptimizerAL::getRhoMax() { return _rho_max; }
double SparseOptimizerAL::getRhoUpdateFactor() { return _rho_update_factor; }


// Setter solver paramters
void SparseOptimizerAL::resetLagrangeMultiplierEq() {
  for (auto &vertex : _vEqLagrangeMultipliers) {
    vertex->setToOrigin();
  }

  // update the multiplier for Equality constraints
  for (auto &edge : _activeEdgesEq) {
    executeEqMultiplierUpdate(edge);
  }
}

// Call the stored function later
void SparseOptimizerAL::executeIneqMultiplierUpdate(void *edgePtr) {
  auto it = updateMultipliersIneqFunctionMap.find(edgePtr);
  if (it != updateMultipliersIneqFunctionMap.end()) {
    it->second(); // Calls updateMultipliers on the correct edge
  } else {
    std::cerr << "[Error] Edge not found in function map" << std::endl;
  }
}

void SparseOptimizerAL::executeEqMultiplierUpdate(void *edgePtr) {
  auto it = updateMultipliersEqFunctionMap.find(edgePtr);
  if (it != updateMultipliersEqFunctionMap.end()) {
    it->second(); // Calls updateMultipliers on the correct edge
  } else {
    std::cerr << "[Error] Edge not found in function map" << std::endl;
  }
}

void SparseOptimizerAL::setAlphaBacktracking(
    std::vector<double> alphaBacktracking) {
  _alphaBacktracking = alphaBacktracking;
}

int SparseOptimizerAL::optimize(int iterations, bool online) {
  if (_ivMap.size() == 0) {
    G2O_WARN("0 vertices to optimize, maybe forgot to call "
             "initializeOptimization()");
    return -1;
  }

  int cjIterations = 0;
  double cumTime = 0;
  bool innerLoopStop = false;
  bool outerLoopStop = false;
  bool ok = true;

  if (_reset_lagrange_multipliers) {
    resetLagrangeMultiplierEq();
    for (auto &edge : _activeEdgesIneq) {
      executeIneqMultiplierUpdate(edge);
    }
  }

  ok = _algorithm->init(online);
  if (!ok) {
    G2O_ERROR("Error while initializing");
    return -1;
  }

  _batchStatistics.clear();
  if (_computeBatchStatistics)
    _batchStatistics.resize(iterations);

  OptimizationAlgorithm::SolverResult result = OptimizationAlgorithm::OK;

  while (!outerLoopStop) {
    int i = 0;
    innerLoopStop = false;
    // BIPM - Inner loop
    while (!innerLoopStop) {
      preIteration(cjIterations);

      if (_computeBatchStatistics) {
        G2OBatchStatistics &cstat = _batchStatistics[cjIterations];
        G2OBatchStatistics::setGlobalStats(&cstat);
        cstat.iteration = cjIterations;
        cstat.numEdges = _activeEdges.size();
        cstat.numVertices = _activeVertices.size();
      }

      double ts = get_monotonic_time();
      result = _algorithm->solve(cjIterations, online);
      ok = (result == OptimizationAlgorithm::OK);

      bool errorComputed = false;
      if (_computeBatchStatistics) {
        computeActiveErrors();
        errorComputed = true;
        _batchStatistics[cjIterations].chi2 = activeRobustChi2();
        _batchStatistics[cjIterations].timeIteration =
            get_monotonic_time() - ts;
      }

      if (verbose()) {
        double dts = get_monotonic_time() - ts;
        cumTime += dts;
        if (!errorComputed)
          computeActiveErrors();
        cerr << "iteration= " << cjIterations
             << "\t chi2= " << FIXED(activeRobustChi2()) << "\t time= " << dts
             << "\t cumTime= " << cumTime
             << "\t edges= " << _activeEdges.size();
        _algorithm->printVerbose(cerr);
        cerr << endl;
      }
      ++cjIterations;
      postIteration(cjIterations);

      // termination criteria
      i++;
      innerLoopStop = i >= _num_inner_iterations_max ||
                      cjIterations >= iterations || terminate() ||
                      verifyConvergence(-1 * 10) || !ok;

      if (innerLoopStop) {
        // update the multiplier for Equality constraints
        for (auto &edge : _activeEdgesEq) {
          executeEqMultiplierUpdate(edge);
        }
      }
    }

    if (result == OptimizationAlgorithm::Fail) {
      return 0;
    }

    if (cjIterations >= iterations ||
        terminate()) { // check if the maximum number of iterations is reached
      break;
    }
    // check if the termination condition is satisfied
    outerLoopStop = verifyConvergence(-1) &&
                    verifyEqFeasibility(_activeEdgesEq, -1.0) &&
                    verifyIneqFeasibility(_activeEdgesIneq, -1.0);

    // update the multiplier for inequality constraints
    for (auto &edge : _activeEdgesIneq) {
      executeIneqMultiplierUpdate(edge);
    }

  } // End of outer loop

  return cjIterations;
}

} // namespace g2o

#endif // G2O_SPARSE_OPTIMIZER_AL_HPP