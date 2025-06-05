#ifndef G2O_GRAPH_OPTIMIZER_BIPM_H_
#define G2O_GRAPH_OPTIMIZER_BIPM_H_

#include "g2o/core/g2o_core_api.h"
#include "sparse_optimizer_for_constraints.h"
#include "vertex_lagrange_multiplier.h"
#include <iostream>
#include <memory>
#include <optimizable_graph.h>
#include <vector>

namespace g2o {

// forwad declaration
template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeIneq;

template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeEq;

// Specialized optimizer inheriting from SparseOptimizer
class G2O_CORE_API SparseOptimizerBIPM
    : public SparseOptimizerForConstraints<SparseOptimizerBIPM> {
  /**
   * @class SparseOptimizerBIPM
   * @brief Implements the Barrier Interior Point Method (BIPM) for efficiently
   * handling constrained optimization problems.
   *
   * This class extends `SparseOptimizer` and introduces functionality for
   * optimization problems with inequality constraints using the Barrier
   * Interior Point Method. It includes various parameters and algorithms to
   * ensure efficient handling of such problems.
   *
   * ### Key Features:
   * - Implements a barrier method where the inequality constraints are
   * incorporated using a logarithmic penalty.
   * - Supports backtracking line search for scaling the update steps to
   * maintain interior feasibility.
   * - Allows customization of optimization parameters to balance between outer
   * and inner iteration performance.
   *
   * ### Parameters:
   * - `_alphaBacktracking`: A vector of scaling factors for the update step
   * during backtracking.
   *   - The update step is computed as:
   *     `x+ = x + _alphaBacktracking[i] * update_step`.
   *   - Index 0 corresponds to a full step (scaling factor of 1), while higher
   * indices correspond to smaller steps.
   *   - `_kappa`: The approximation parameter for the logarithmic barrier
   * function.
   *   - The optimization objective is adjusted as:
   *     `J = f(x) - 1/_kappa * log(g(x))`, where `g(x)` represents the
   * inequality constraints.
   * - `_kappaUpdateFactor`: The update factor for `_kappa` in the outer loop.
   * Typically between 10 and 100 (default is 20).
   *   - A larger `_kappaUpdateFactor` reduces the number of outer iterations
   * but increases the inner iterations.
   * - `_kappa`: Default value for the barrier function parameter `_kappa`
   * (default: 1).
   * - `_epsilon`: Stopping criterion for the outer loop (default: 1e-3).
   *   - Optimization stops when ` _kappa >= _kappaFinal`.
   * - `_epsilon_feas`: Stopping criterion for the inner loop (default: 0.01).
   *
   * ### Default Values:
   * - `_kappa = 1`: Default barrier function parameter.
   * - `_kappaUpdateFactor = 8`: Default updating factor for `_mu`.
   * - `_epsilon = 1e-3`: Outer loop stopping criterion.
   * - `_epsilon_feas = 0.01`: Inner loop stopping criterion.
   * - `_alphaBacktracking`: `{1, 0.9, 0.8, 0.6, 0.5, 0.4, 0.31, 0.19, 0.08,
   * 0.02, 1e-2, 1e-3, 1e-4, 1e-5, 1e-7, 1e-9, 1e-10, 1e-12, 1e-14, 1e-16,
   * 1e-18, 1e-20}`.
   *
   * ### Algorithms:
   * - **Backtracking Algorithm**: Implements a backtracking line search using
   * the `_alphaBacktracking` vector to ensure the decision variable remains
   * within the interior of the feasible region.
   * - **Barrier Method**: Incorporates a logarithmic barrier for inequality
   * constraints with `_mu` as the approximation parameter.
   *
   * @note Proper initialization of optimization structures is required before
   * invoking methods like `optimize()` or `initializeOptimization()`.
   */

public:
  SparseOptimizerBIPM();  // Default constructor
  ~SparseOptimizerBIPM(); // Virtual destructor

  // Function to construct the quadratic form implementation
  template <int D, typename E, typename... VertexTypes>
  void constructQuadraticFormIneq(
      BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge);
  // Function to construct the quadratic form implementation for Augmented
  // Lagrangian Eqaulity Algorithm
  template <int D, typename E, typename... VertexTypes>
  void
  constructQuadraticFormEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...> &edge);

  template <int D, typename E, typename... VertexTypes>
  bool addEdgeIneqImpl(BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e);
  /**
   * add the inequlaity edge to the optimizer and the set of inequality edges
   * @param e: the edge to be added
   * @returns false if somethings goes wrong
   */

  template <int D, typename E, typename... VertexTypes>
  bool addEdgeEqImpl(BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e);
  /**
   * add the inequlaity edge to the optimizer and the set of equality edges
   * @param e: the edge to be added
   * @returns false if somethings goes wrong
   */

  virtual void update(const double *update) override;

  virtual int optimize(int iterations, bool online = false) override;
  /**
   * computes the blocks of the inverse of the specified pattern.
   * the pattern is given via pairs <row, col> of the blocks in the hessian
   * @param blockIndices: the pattern
   * @param spinv: the sparse block matrix with the result
   * @returns false if the operation is not supported by the solver
   */

  virtual double backtrackingAlgorithm(const double *update) override;
  /**
  backtracking algorithm for scaling the update step such that the decision
  variable remains within the interior of the feasible region
  @param update: the update step
  @returns the the scaling factor for a successful step size
  */

  // virtual void defineAlphaBacktracking() {      }

  // solver parameter

public:
  // Setter for the mu value
  void setKappa(double kappa);
  void setKappaInitial(double kappaInitial);
  void setKappaFinal(double kappaFinal);
  void setKappaUpdateFactor(double kappaUpdateFactor);

  // getter for the solver parameters
  double getKappa();
  double getKappaInitial();
  double getKappaFinal();
  double getKappaUpdateFactor();


  void setAlphaBacktracking(std::vector<double> alphaBacktracking);
  void resetLagrangeMultiplierEq() override;

protected:
  std::vector<std::shared_ptr<g2o::OptimizableGraph::Vertex>>
      _vEqLagrangeMultipliers; // Set of Lagrangian vertices
  std::unordered_map<void *, std::function<void()>>
      updateMultipliersIneqFunctionMap;
  // solver parameter
public:
  double _kappa; // Default value for the barrier function parameter _kappa
  double _kappa_initial =
      1; // Default value for the barrier function parameter _kappa0
  double _kappa_update_factor =
      8; // Default value for the barrier function parameter _kappaUpdate (\nu).
         // Typically between 5 and 100 (default is 20).
  double _kappa_final =
      1500; // Maximum value for the barrier function parameter
  bool _reset_lagrange_multipliers = false;
};

template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerBIPM::addEdgeIneqImpl(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e) {
  bool eresult = OptimizableGraph::addEdge(e);
  if (!eresult) {
    std::cerr << "[Error] adding Ineq edge to the optimizer" << std::endl;
    return false;
  }

  // add teh edge to the set of inequality set
  _edgeIneqSet.insert(e);

  // set the ConstructQuadraticFormImpl function

  auto lambda = [this](BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge) {
    this->constructQuadraticFormIneq(edge);
  };

  e->setConstructQuadraticFormImpl(lambda);
  return true;
}

template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerBIPM::addEdgeEqImpl(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e) {
  // create new vertex for the lagrangian
  auto vLagrangian = std::make_shared<VertexLagrangeMultiplier<D>>();
  // auto* nu = new VertexLagrangeMultiplier<D>();
  vLagrangian->setId(e->getVertexLagrangeMultiplierId());
  vLagrangian->setInitialValue(_lagrange_multiplier_initial);

  // Add the Nu vertex to the optimizer
  this->addVertex(vLagrangian.get());

  // assign the vertex to the edge
  e->setVertex(e->vertices().size() - 1, vLagrangian.get());

  bool eresult = OptimizableGraph::addEdge(e);
  if (!eresult) {
    std::cerr << "[Error] adding Eq edge to the optimizer" << std::endl;
    return false;
  }

  vLagrangian->setFixed(false);
  e->initializeInformationMatrix();

  // add teh edge to the set of inequality set
  _edgeEqSet.insert(e);
  _vEqLagrangeMultipliers.push_back(vLagrangian);

  return true;
}

// Implementation of the quadratic form construction
template <int D, typename E, typename... VertexTypes>
void SparseOptimizerBIPM::constructQuadraticFormIneq(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge) {
  // Inverse of the error (if needed)
  auto error_inv = edge.error().array().inverse();
  // Weighted error with the barrier function parameter _mu
  double _mu = 1.0 / _kappa;
  auto weightedError = error_inv * _mu;
  // Diagonal matrix (omega) calculation
  auto omega = weightedError.array() * error_inv.array();
  auto omega_matrix =
      omega.matrix().asDiagonal(); // Convert into a diagonal matrix
  // Pass the result to the edge method for constructing the quadratic form
  static const std::size_t _nr_of_vertices = sizeof...(VertexTypes);
  edge.constructQuadraticFormNs(omega_matrix, weightedError,
                                std::make_index_sequence<_nr_of_vertices>());
}

// In case we use AL alg for Eq with BIPM
template <int D, typename E, typename... VertexTypes>
void SparseOptimizerBIPM::constructQuadraticFormEq(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> &edge) {

  auto error = edge.error();
  auto gamma = edge.lagrangeMultiplier();
  auto rho = edge.rho();
  Eigen::Matrix<double, 2 * D, 1> weightedError =
      Eigen::Matrix<double, 2 * D, 1>::Zero();
  Eigen::DiagonalMatrix<double, 2 * D> omega_matrix =
      Eigen::DiagonalMatrix<double, 2 * D>(
          Eigen::Matrix<double, 2 * D, 1>::Zero());
  for (int i = 0; i < D; ++i) {
    omega_matrix.diagonal()[i] = rho[i];
    weightedError[i] =
        -(omega_matrix.diagonal()[i] * error[i] + 0.5 * gamma[i]);
  }

  static const std::size_t _nr_of_vertices =
      sizeof...(VertexTypes) + 1; // add the lagrange multiplier vertex
  edge.constructQuadraticFormNs(omega_matrix, weightedError,
                                std::make_index_sequence<_nr_of_vertices>());
}

} // namespace g2o

#include "sparse_optimizer_bipm.hpp"

#endif // G2O_GRAPH_OPTIMIZER_BIPM_H_
