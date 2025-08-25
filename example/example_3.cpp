/*
Copyright (c) 2023, University of Luxembourg
All rights reserved.

Redistributions and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*/

/*
     minimize     0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 + 2*x2 + 1||^2
    subject to    x1 + x1^3 + x2 + x2^2  = 0
                   
*/

#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>

#include <Eigen/Core>
#include <iostream>

#include "bipm_g2o/linear_solver_eigen_qr.h"

#ifdef USE_BIPM
#include "bipm_g2o/sparse_optimizer_bipm.h" // Barrier interior point method
#else
#include "bipm_g2o/sparse_optimizer_al.h" // using the Augmented Lagrangian
#endif

#include "include/example_3_vertices_edges.h" // for defining the edges

int main(int argc, char **argv) {
  std::cout << "OP: min  0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 + 2*x2 + 1||^2 s.t.    "
               " x1 + x1^3 + x2 + x2^2  = 0"
            << std::endl;
  std::cout
      << "Usage: " << argv[0]
      << "solverType:<0=GN,1=LV,2=DL> iterations<int> a<int> b<int> c<int>"
      << std::endl;
  std::cout << "x_start<int> y_start<int> z_start<int> " << std::endl;
  std::cout << "./example -0.2 -0.2 0 150" << std::endl;

#ifdef USE_BIPM
  std::cout << " The Barriar Interior Point Solver is Used\n";
#else
  std::cout << "The Augumented Lagrangian Solver is Used\n";
#endif

  int argCount = 1;
  double xStart = (argc > argCount) ? std::atof(argv[argCount]) : -1.2 ;

  argCount++;
  double yStart = (argc > argCount) ? std::atof(argv[argCount]) : -1;
  argCount++;
  int solverType = (argc > argCount) ? std::atoi(argv[argCount]) : 0;

  argCount++;
  int numberOfIterations = (argc > argCount) ? std::atoi(argv[argCount]) : 150;

// Initialize optimizer
#ifdef USE_BIPM
  g2o::SparseOptimizerBIPM optimizer;
  // Add the solver settings, kappa is the the barrier parameter [ -1/kappa
  // log(-g(x))]
  optimizer.setKappaInitial(
      1.0); // initial value for the barrier function parameter
  optimizer.setKappaUpdateFactor(
      50.0); // update factor for the barrier function parameter
  optimizer.setKappaFinal(
      1.0e5); // maximum value for the barrier function parameter
  optimizer.setLagrangeMultiplierInitial(
      0.0); // initial value for the Lagrangian vertex of the Equality
            // constraints
  optimizer.setInnerIterationsMax(100); // maximum number of inner iterations
#else
  g2o::SparseOptimizerAL optimizer;
  optimizer.setRhoInitial(
      11.0); // initial value for the barrier function parameter
  optimizer.setRhoUpdateFactor(
      20.0);                // update factor for the barrier function parameter
  optimizer.setRhoMax(5e5); // maximum value for the barrier function parameter
  optimizer.setLagrangeMultiplierInitial(
      0.0); // initial value for the Lagrangian vertex of the Equality
            // constraints
  optimizer.setInnerIterationsMax(10); // maximum number of inner iterations
#endif

  optimizer.setVerbose(true);

  std::unique_ptr<g2o::BlockSolverX> blockSolver;
  std::unique_ptr<g2o::LinearSolver<g2o::BlockSolverX::PoseMatrixType>>
      linearSolver;
  std::unique_ptr<g2o::OptimizationAlgorithm> algorithm;

   linearSolver = std::make_unique<
      g2o::LinearSolverEigenQR<g2o::BlockSolverX::PoseMatrixType>>();

  blockSolver = std::make_unique<g2o::BlockSolverX>(std::move(linearSolver));

  // Set up the optimization algorithm
  switch (solverType) {
  case 1:
    algorithm = std::make_unique<g2o::OptimizationAlgorithmLevenberg>(
        std::move(blockSolver));
    break;
  case 2:
    algorithm = std::make_unique<g2o::OptimizationAlgorithmDogleg>(
        std::move(blockSolver));
    break;
  default:
    algorithm = std::make_unique<g2o::OptimizationAlgorithmGaussNewton>(
        std::move(blockSolver));
    break;
  }

  optimizer.setAlgorithm(algorithm.get());

  // Add the vertex to the optimizer
  auto xy = std::make_shared<VertexXY>();
  xy->setId(0);
  xy->setEstimate(Eigen::Vector2d(xStart, yStart)); // initial estimate
  optimizer.addVertex(xy.get());

  // Add the cost edgeXY to the optimizer: 0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 + 2*x2 + 1||^2
  auto edgeXY = std::make_shared<EdgeXY>();
  edgeXY->setVertex(0, xy.get());
  edgeXY->setMeasurement(Eigen::Vector2d(0, 0)); //
  edgeXY->setInformation(Eigen::Matrix2d::Identity() * 0.5);
  optimizer.addEdge(edgeXY.get());


  auto *edgeEq = new EdgeEq(); // add x1 + x1^3 + x2 + x2^2 = 0
  edgeEq->setVertexLagrangeMultiplierId(10);
  edgeEq->setVertex(0, xy.get());
  optimizer.addEdgeEq(edgeEq);

  // Optimize
  optimizer.initializeOptimization();
  class terminationCriterionType;

  optimizer.setEpsilonConvergence(1e-3); // Stopping criterion for update norm
  optimizer.setEpsilonConstraint(1e-3);

  optimizer.optimize(numberOfIterations);

  // Output the results
  std::cout << "Optimized x: " << xy->estimate().transpose() << std::endl;

  // std::cout << "Optimized nu: " << nu->estimate() << std::endl;

  return 0;
}
