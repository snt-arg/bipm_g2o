#ifndef G2O_SPARSE_OPTIMIZER_BIPM_HPP
#define G2O_SPARSE_OPTIMIZER_BIPM_HPP

#include "sparse_optimizer_bipm.h"
#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"
#include "optimization_algorithm.h"
#include "solver.h"
#include <Eigen/Core>
#include <iostream>
#include <vector>

namespace g2o {
using namespace std;

// Default constructor
SparseOptimizerBIPM::SparseOptimizerBIPM() {
   _lagrange_multiplier_initial = 0;
   _num_inner_iterations_max = 100;
}

// Default destructor
SparseOptimizerBIPM::~SparseOptimizerBIPM() = default;

// Setter solver paramters
void SparseOptimizerBIPM::setKappa(double kappa) {_kappa = kappa;}
void SparseOptimizerBIPM::setKappaInitial(double kappaInitial) {_kappa_initial = kappaInitial;}
void SparseOptimizerBIPM::setKappaUpdateFactor(double kappaUpdateFactor) {_kappa_update_factor = kappaUpdateFactor;} 
void SparseOptimizerBIPM::setKappaFinal(double kappaFinal) {_kappa_final = kappaFinal;}

 


// getter for the solver parameters
double SparseOptimizerBIPM::getKappa(){ return _kappa;}
double SparseOptimizerBIPM::getKappaInitial() { return _kappa_initial; }
double SparseOptimizerBIPM::getKappaFinal() { return _kappa_final; }
double SparseOptimizerBIPM::getKappaUpdateFactor() { return _kappa_update_factor; }

void SparseOptimizerBIPM::resetLagrangeMultiplierEq() {
  for (auto &vertex : _vEqLagrangeMultipliers) {
    vertex->setToOrigin();
  }
}


void SparseOptimizerBIPM::setAlphaBacktracking(
    std::vector<double> alphaBacktracking) {
  _alphaBacktracking = alphaBacktracking;
}

double SparseOptimizerBIPM::backtrackingAlgorithm(const double *update) {
  // Initialize variables
  size_t c = 0;
  int offset = 0;
  bool _inside_interior_flag = true;
  bool debugFlag = false;

  if (debugFlag) {
    // cout the vertex estimate before backtracking
    std::cout << " before backtracking = ";
    std::cout << "  xy: "
              << static_cast<VertexLagrangeMultiplier<2> *>(this->vertex(0))
                     ->estimate()
                     .transpose();
    std::cout << "  z: "
              << static_cast<VertexLagrangeMultiplier<1> *>(this->vertex(1))
                     ->estimate()
                     .transpose()
              << std::endl;
  }

  // Loop over all inequality edges
  for (auto &edge : _activeEdgesIneq) {
    // cast the edge from OptimizableGraph to BaseFixedSizedEdgeIneq
    while (c < _alphaBacktracking.size() - 1) {
      _inside_interior_flag = true;
      // 1- Scale the update and apply it to the vertices
      for (size_t i = 0; i < edge->vertices().size(); ++i) {
        // Get the vertex
        OptimizableGraph::Vertex *vertex =
            static_cast<OptimizableGraph::Vertex *>(edge->vertex(i));
        offset = vertex->colInHessian(); // Assumes hessianIndex() gives the
                                         // position of the vertex in _x

        // Scale the update for this step size
        Eigen::Map<const Eigen::VectorXd> updateVec(update + offset,
                                                    vertex->dimension());
        Eigen::VectorXd scaledUpdate = _alphaBacktracking[c] * updateVec;

        // Apply the scalled update
        vertex->oplus(scaledUpdate.data());
      }

      // 2- Check if the inequality constraint is satisfied
      edge->computeError();
      auto error = edge->errorData();

      // 3- Reverse the update and apply it to the vertices
      for (size_t i = 0; i < edge->vertices().size(); ++i) {
        // Get the vertex
        OptimizableGraph::Vertex *vertex =
            static_cast<OptimizableGraph::Vertex *>(edge->vertex(i));
        offset = vertex->colInHessian(); // Assumes hessianIndex() gives the
                                         // position of the vertex in _x

        // Scale the update for this step size to reverse
        Eigen::Map<const Eigen::VectorXd> updateVec(update + offset,
                                                    vertex->dimension());
        Eigen::VectorXd reverseUpdate = -_alphaBacktracking[c] * updateVec;

        // Apply the scalled update
        vertex->oplus(reverseUpdate.data());

        if (debugFlag) {
          std::cout << " After  local  update reversed  scaled update  "
                    << endl;
          if (offset == 0) {
            std::cout << "  xy: "
                      << static_cast<VertexLagrangeMultiplier<2> *>(
                             this->vertex(0))
                             ->estimate()
                             .transpose()
                      << std::endl;
          } else {
            std::cout << "  z: "
                      << static_cast<VertexLagrangeMultiplier<1> *>(
                             this->vertex(1))
                             ->estimate()
                             .transpose()
                      << std::endl;
          }
          std::cout << "-------------------------------------------------------"
                       "----------"
                    << endl;
        }
      }

      // 4- Check if the inequality constraint is satisfied
      for (int i = 0; i < edge->dimension(); ++i) {
        if (error[i] > 0) {
          _inside_interior_flag = false;
          if (debugFlag) {
            std::cout << "-----------------------------------------------------"
                         "------------"
                      << endl;
            std::cout << " The inequality constraint is not satisfied at c = "
                      << c << endl;
          }
          break;
        }
      }

      // 5- Break if the inequality constraint is satisfied
      if (_inside_interior_flag)
        break;

      // 6- Increment the step size index if the inequality constraint is not
      // satisfied
      ++c;
    }
  }

  return _alphaBacktracking[c]; // Return the index of _alphaBacktracking for
                                // the corresponding successful step size
}

void SparseOptimizerBIPM::update(const double *update) {

  OptimizationAlgorithmWithHessian *algorithm =
      static_cast<OptimizationAlgorithmWithHessian *>(_algorithm);
  size_t xSize = algorithm->solver().vectorSize();

  Eigen::Map<const Eigen::VectorXd> updateVec(update, xSize);

  // 1- find the scaling index for the scaling factors vector _alphaBacktracking
  double scaleFactor = backtrackingAlgorithm(update);
  Eigen::VectorXd scaledUpdateVec;
  const double *requiredUpdate = nullptr;
  // 2- Calculate the required update step
  if (scaleFactor == 1) { // No scaling required
    requiredUpdate = update;
  } else {
    scaledUpdateVec = scaleFactor * updateVec;
    requiredUpdate = scaledUpdateVec.data();
  }

  // 3- Update the decision variables
  SparseOptimizer::update(requiredUpdate);
}

int SparseOptimizerBIPM::optimize(int iterations, bool online) {
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
  _kappa = _kappa_initial; // reset _kappa to the initial value
  if (~_warm_start_lagrange_multiplier_flag){
  resetLagrangeMultiplierEq();
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

  if (_activeEdgesIneq.empty()) {
    // set mu to 0 if there are no inequality edges
    this->setKappa(_kappa_final);
  }

  if (_kappa_update_factor <= 1 ) {
    std::cerr << "[WARNING] Value of_nu = " << _kappa_update_factor
              << "; It should be greater than 1." << std::endl;
  }
  while (!outerLoopStop) { 
     
    int i = 1;
    innerLoopStop = false;
    // BIPM - Inner loop
    while (!innerLoopStop) {
      i++;
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
      innerLoopStop = i >= _num_inner_iterations_max ||
                      cjIterations >= iterations || terminate() ||
                      verifyConvergence(-1 * 10) || !ok;  
    }

    if (result == OptimizationAlgorithm::Fail) {
      return 0;
    }

    if (cjIterations >= iterations ||
        terminate()) { // check if the maximum number of iterations is reached
      break;
    }
    // check if the termination condition is satisfied
    outerLoopStop =
        verifyConvergence(-1) &&
        verifyEqFeasibility(_activeEdgesEq, -1.0) &&
        verifyIneqFeasibility(_activeEdgesIneq, -1.0); // defualt * 1
    outerLoopStop = _kappa >= _kappa_final && (outerLoopStop);

    // update _kappa
    _kappa = (_kappa < _kappa_final) ? (_kappa * _kappa_update_factor) : _kappa;

  } // End of inner loop

  return cjIterations;
}

} // namespace g2o


#endif // G2O_SPARSE_OPTIMIZER_BIPM_HPP