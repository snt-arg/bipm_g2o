#ifndef G2O_GRAPH_OPTIMIZER_AL_H_
#define G2O_GRAPH_OPTIMIZER_AL_H_

#include "sparse_optimizer_for_constraints.h"
#include "vertex_lagrange_multiplier.h"

#include "g2o/core/g2o_core_api.h"
#include <Eigen/src/Core/DiagonalMatrix.h>
#include <Eigen/src/Core/Matrix.h>
#include <iostream>
#include <memory>
#include <optimizable_graph.h>
#include <unordered_map>

#include <type_traits> //

namespace g2o {

// forwad declaration
template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeIneq;

// Specialized optimizer inheriting from SparseOptimizer
class G2O_CORE_API SparseOptimizerAL
    : public SparseOptimizerForConstraints<SparseOptimizerAL> {
  /**
   * @class SparseOptimizerAL
   * @brief Implements the Augmented Lagrangian (AL) method for efficiently
   * handling constrained optimization problems.
   *
   * This class extends `SparseOptimizerForConstraints` and introduces
   * functionality for optimization problems with equanu lity and inequality
   * constraints using the Augmented Lagrangian method.
   *
   * ### Key Features:
   * - Implements an Augmented Lagrangian approach where constraints are
   * incorporated by modifying the objective function.
   * - Supports dynamic updates of Lagrange multipliers to improve convergence.
   * - Allows customization of penalty parameters for efficient constraint
   * handling.
   *
   * ### Parameters:
   * - `_rho`: Penalty parameter for the augmented Lagrangian method
   *   - Determines the weight of constraint violations in the objective
   * function.
   * - `_rho_update_factor`: Multiplicative factor for updating `_rho`
   *   - Controls how quickly the penalty increases to enforce constraints.
   * - `_epsilon`: Stopping criterion for optimization (default: 1e-3).
   * - `_epsilon_feas`: Stopping criterion for constraint satisfaction (default:
   * 0.01).
   *
   * ### Default Values:
   * - `_rho = 1.0`: Initial penalty parameter.
   * - `_rho_update_factor = 10.0`: Default factor for penalty updates.
   * - `_epsilon = 1e-3`: Optimization stopping criterion.
   * - `_epsilon_feas = 0.01`: Constraint satisfaction stopping criterion.
   *
   * ### Algorithms:
   * - **Multiplier Update Algorithm**: Updates the Lagrange multipliers
   * iteratively to improve constraint satisfaction.
   * - **Penalty Method**: Incorporates squared constraint violations into the
   * objective function to guide optimization.
   *
   * @note Proper initialization of optimization structures is required before
   * invoking methods like `optimize()` or `initializeOptimization()`.
   */

public:
  SparseOptimizerAL();  // Default constructor
  ~SparseOptimizerAL(); // Virtual destructor

  // Function to construct the quadratic form implementation
  template <int D, typename E, typename... VertexTypes>
  void constructQuadraticFormInEq(
      BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge);
  // Function to construct the quadratic form implementation
  template <int D, typename E, typename... VertexTypes>
  void
  constructQuadraticFormEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...> &edge);
  // Note that both the Eq and InEq cosntraints use the same edge type
  // The difference is in the implementation of the constructQuadraticForm
  // function

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

  virtual int optimize(int iterations, bool online = false) override;
  /**
   * computes the blocks of the inverse of the specified pattern.
   * the pattern is given via pairs <row, col> of the blocks in the hessian
   * @param blockIndices: the pattern
   * @param spinv: the sparse block matrix with the result
   * @returns false if the operation is not supported by the solver
   */

  /*
  virtual double backtrackingAlgorithm (const double* update);
  / / * *
  backtracking algorithm for scaling the update step such that the decision
  variable remains within the interior of the feasible region
  @param update: the update step
  @returns the the scaling factor for a successful step size
  */

  void executeIneqMultiplierUpdate(void *edgePtr);
  void executeEqMultiplierUpdate(void* edgePtr);


  // virtual void defineAlphaBacktracking() {      }
  template <int D, typename E, typename... VertexTypes, typename EdgeType>
  void updateMultipliers(EdgeType &edge);
  ;

  template <int D, typename E, typename... VertexTypes>
  void updateMultipliersInEq(OptimizableGraph::Edge *edge);

  // solver parameter

public:
  // Setter for the mu value
   void setRho(double rho);
   void setRhoInitial(double rhoInitial);
   void setRhoMax(double rhoMax);
   void setRhoUpdateFactor(double rhoUpdateFactor);
   
   //getter for the solver parameters
   double getRho();
   double getRhoInitial();
   double getRhoMax();
   double getRhoUpdateFactor();




  // Setter function
  void setAlphaBacktracking(std::vector<double> alphaBacktracking);
  void resetLagrangeMultiplierEq() override;


protected:
  std::unordered_map<void *, std::function<void()>>
      updateMultipliersIneqFunctionMap;


  double  _rho_min = .00001;  // Minimum value for the penalty parameter
  double _rho_initial = 10; // Initial value for the penalty parameter
  double _rho_max = 500; // Maximum value for the penalty parameter
  double _rho_update_factor = 5; // Update factor for the penalty parameter
 
  bool _useSlackVariables = true; // Using slack variables for the inequality
                                  // constraints. it does not appear directly
  bool _reset_lagrange_multipliers = true;

  std::vector<std::shared_ptr<g2o::OptimizableGraph::Vertex>>
      _vEqLagrangeMultipliers; // Set of Lagrangian vertices

  // _rho_bar is saved in the variable called slack variable
};

template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerAL::addEdgeIneqImpl(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e) {
  bool eresult = OptimizableGraph::addEdge(e);
  if (!eresult) {
    std::cerr << "[Error] Failed adding Ineq edge to the optimizer"
              << std::endl;
    return false;
  }

  // use the slack variable member to save /rho_bar
  Eigen::Matrix<double, D, 1> rho_bar =
      Eigen::Matrix<double, D, 1>::Constant(_rho_initial);
  e->setSlackVariable(rho_bar);
  e->setRho(rho_bar);
  // save the constraint variation in the information matrix

  e->setConstraintViolationPrev(e->error().array().min(0.0));
  // add teh edge to the set of inequality set

  // set the ConstructQuadraticFormImpl function
  auto lambda = [this](BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge) {
    this->constructQuadraticFormInEq(edge);
  };
  e->setConstructQuadraticFormImpl(lambda);

  // Store the function in a map instead
  updateMultipliersIneqFunctionMap[e] = [this, e]() {
    this->updateMultipliers<D, E, VertexTypes...>(*e);
  };

  _edgeIneqSet.insert(e);
  return true;
}

template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerAL::addEdgeEqImpl(
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

 
    vLagrangian->setFixed(true);

    Eigen::Matrix<double, D, 1> rho_bar =
        Eigen::Matrix<double, D, 1>::Constant(_rho_initial);

    // In equality constraints, the information matrix 2D is used to store the
    // constraint violation and the rho_bar we initialize the the constraint
    // violation to error and rho_bar to as initial value Assign first D
    // diagonal entries from `error.segment(0, D)
    e->setConstraintViolationPrev(e->error().segment(0, D).cwiseAbs());
    e->setConstraintViolationPrev(std::move(rho_bar), D);

    e->setRho(rho_bar);

    // set the ConstructQuadraticFormImpl function
    auto lambda = [this](BaseFixedSizedEdgeEq<D, E, VertexTypes...> &edge) {
      this->constructQuadraticFormEq(edge);
    };
    e->setConstructQuadraticFormImpl(lambda);

    // Store the function in a map instead
    updateMultipliersEqFunctionMap[e] = [this, e]() {
      this->updateMultipliers<D, E, VertexTypes...>(*e);
    };
  

  // add teh edge to the set of inequality set
  _edgeEqSet.insert(e);
  _vEqLagrangeMultipliers.push_back(vLagrangian);

  return true;
}

template <int D, typename E, typename... VertexTypes, typename EdgeType>
void SparseOptimizerAL::updateMultipliers(EdgeType &edge) {
  // Ensure EdgeType is either BaseFixedSizedEdgeEq or BaseFixedSizedEdgeIneq
  static_assert(
      std::is_base_of<BaseFixedSizedEdgeEq<D, E, VertexTypes...>,
                      EdgeType>::value ||
          std::is_base_of<BaseFixedSizedEdgeIneq<D, E, VertexTypes...>,
                          EdgeType>::value,
      "updateMultipliers can only accept BaseFixedSizedEdgeEq or "
      "BaseFixedSizedEdgeIneq");

  constexpr bool is_eq =
      std::is_base_of<BaseFixedSizedEdgeEq<D, E, VertexTypes...>,
                      EdgeType>::value;
  constexpr bool is_ineq =
      std::is_base_of<BaseFixedSizedEdgeIneq<D, E, VertexTypes...>,
                      EdgeType>::value;
  if (_reset_lagrange_multipliers) {
    if constexpr (is_ineq) {
      auto multiplier = edge.lagrangeMultiplier();
      multiplier.setConstant(_lagrange_multiplier_initial);
      edge.setLagrangeMultiplier(multiplier);

      Eigen::Matrix<double, D, 1> rho_bar =
          Eigen::Matrix<double, D, 1>::Constant(_rho_initial);
      edge.setSlackVariable(rho_bar);
      edge.setRho(rho_bar);
      // save the constraint variation in the information matrix

      edge.setConstraintViolationPrev(edge.error().array().min(0.0));
      // add teh edge to the set of inequality set

      return;
    }

    if constexpr (is_eq) {
      auto multiplier = edge.lagrangeMultiplier();
      multiplier.setConstant(_init_lagrange_multiplier);
      edge.setLagrangeMultiplier(multiplier);

      Eigen::Matrix<double, D, 1> rho_bar =
          Eigen::Matrix<double, D, 1>::Constant(_rho_initial);

      edge.setRho(rho_bar);
      // In equality constraints, the information matrix 2D is used to store the
      // constraint violation and the rho_bar we initialize the the constraint
      // violation to error and rho_bar to as initial value Assign first D
      // diagonal entries from `error.segment(0, D)
      edge.setConstraintViolationPrev(edge.error().segment(0, D).cwiseAbs());
      edge.setConstraintViolationPrev(std::move(rho_bar), D); 
      return;
    }
  }

  auto error = edge.error();
  auto rho = edge.rho();
  auto multiplier = edge.lagrangeMultiplier();
  Eigen::Matrix<double, D, 1> rho_bar;
  auto constraints_violation_prev = edge.information().diagonal().segment(0, D);
  auto constraints_violation = error.segment(0, D);

  double rho_min;
  double rho_max;
  double rho_upate_factor;
  if constexpr (is_ineq) {
    rho_min = _rho_min;
    rho_max = _rho_max;
    rho_upate_factor = _rho_update_factor;
  } else if constexpr (is_eq) {
    rho_min = _rho_min;
    rho_max = _rho_max;
    rho_upate_factor = _rho_update_factor;
  }

  for (int i = 0; i < D; ++i) {
    [[maybe_unused]] double x_d = 0, x_i = 0;

    switch (2) {
    case 1:

      if constexpr (is_ineq) {
        constraints_violation = constraints_violation.array().min(0.0);
      } else if constexpr (is_eq) {
        constraints_violation = constraints_violation.cwiseAbs();
      }

      if constexpr (is_ineq) {
        rho_bar = edge.slackVariable();
      } else if constexpr (is_eq) {
        rho_bar = constraints_violation_prev.segment(D, D);
      }

      if (constraints_violation_prev[i] > constraints_violation[i])
        x_d = 1 - constraints_violation[i] / constraints_violation_prev[i];

      if (constraints_violation[i] > constraints_violation_prev[i])
        x_i = 1 - constraints_violation_prev[i] / constraints_violation[i];
      rho[i] = rho_bar[i] + x_d * (rho_max - rho_bar[i]);
      rho_bar[i] = rho_bar[i] + x_d * (rho_max - rho_bar[i]);

      if constexpr (is_ineq) {
        edge.setConstraintViolationPrev(constraints_violation);
        edge.setSlackVariable(rho_bar);
      } else if constexpr (is_eq) {
        edge.setConstraintViolationPrev(constraints_violation);
        edge.setConstraintViolationPrev(std::move(rho_bar), D);
      }
      break;
    case 2:
      rho[i] = rho_upate_factor * rho[i];
      if (rho[i] < rho_min)
        rho[i] = rho_min;
      if (rho[i] > rho_max)
        rho[i] = rho_max;
      break;
    }

    // update mulitplier
    if constexpr (is_ineq) {
      multiplier[i] = std::max(0.0, multiplier[i] + 2 * rho[i] * error[i]);
    } else if constexpr (is_eq) {
      multiplier[i] = multiplier[i] + 2 * rho[i] * error[i];
    }
  }

  edge.setLagrangeMultiplier(multiplier);
  edge.setRho(rho);
}

// Implementation of the quadratic form construction
template <int D, typename E, typename... VertexTypes>
void SparseOptimizerAL::constructQuadraticFormInEq(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge) {
  auto error = edge.error();
  auto lambda = edge.lagrangeMultiplier();
  auto rho = edge.rho();
  Eigen::Matrix<double, D, 1> weightedError;
  Eigen::DiagonalMatrix<double, D> omega_matrix =
      Eigen::DiagonalMatrix<double, D>(Eigen::Matrix<double, D, 1>::Zero());
  for (int i = 0; i < D; ++i) {
    int active = 0;
    if (!_useSlackVariables && error[i] > 0)
      active = 1;

    if (_useSlackVariables && error[i] > -lambda[i] / (2 * rho[i]))
      active = 1;

    omega_matrix.diagonal()[i] = active * rho[i];

    weightedError[i] =
        (omega_matrix.diagonal()[i] * error[i] + 0.5 * lambda[i]);
    if (_useSlackVariables)
      weightedError[i] = -active * weightedError[i];
  }
  static const std::size_t _nr_of_vertices = sizeof...(VertexTypes);
  edge.constructQuadraticFormNs(omega_matrix, weightedError,
                                std::make_index_sequence<_nr_of_vertices>());
}

// Implementation of the quadratic form construction
template <int D, typename E, typename... VertexTypes>
void SparseOptimizerAL::constructQuadraticFormEq(
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

#include"sparse_optimizer_al.hpp"

#endif // G2O_GRAPH_OPTIMIZER_BIPM_H_
