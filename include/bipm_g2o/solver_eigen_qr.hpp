#ifndef G2O_LINEAR_SOLVER_EIGEN_QR_HPP
#define G2O_LINEAR_SOLVER_EIGEN_QR_HPP

#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_dogleg.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "linear_solver_eigen_qr.h"


using namespace std;
namespace g2o {

namespace {
template <int p, int l, bool blockorder>
std::unique_ptr<BlockSolverBase> AllocateSolverSparseQR() {
  G2O_DEBUG("Using EigenQR poseDim {} landMarkDim {} blockordering {}", p, l, blockorder);
  auto linearSolver = std::make_unique<
      LinearSolverEigenQR<typename BlockSolverPL<p, l>::PoseMatrixType>>();
  return std::make_unique<BlockSolverPL<p, l>>(std::move(linearSolver));
}
}  // namespace

static OptimizationAlgorithm* createSolverSparseQR(const std::string& fullSolverName) {
  static const std::map<std::string,
                        std::function<std::unique_ptr<BlockSolverBase>()>>
      solver_factories{
          {"var", &AllocateSolverSparseQR<-1, -1, true>},
          {"fix3_2", &AllocateSolverSparseQR<3, 2, true>},
          {"fix6_3", &AllocateSolverSparseQR<6, 3, true>},
          {"fix7_3", &AllocateSolverSparseQR<7, 3, true>},
          {"fix3_2_scalar", &AllocateSolverSparseQR<3, 2, false>},
          {"fix6_3_scalar", &AllocateSolverSparseQR<6, 3, false>},
          {"fix7_3_scalar", &AllocateSolverSparseQR<7, 3, false>},
      };

  string solverName = fullSolverName.substr(3);
  auto solverf = solver_factories.find(solverName);
  if (solverf == solver_factories.end()) return nullptr;

  string methodName = fullSolverName.substr(0, 2);

  if (methodName == "gn") {
    return new OptimizationAlgorithmGaussNewton(solverf->second());
  } else if (methodName == "lm") {
    return new OptimizationAlgorithmLevenberg(solverf->second());
  } else if (methodName == "dl") {
    return new OptimizationAlgorithmDogleg(solverf->second());
  }

  return nullptr;
}

class EigenSparseQR_SolverCreator : public AbstractOptimizationAlgorithmCreator {
 public:
  explicit EigenSparseQR_SolverCreator(const OptimizationAlgorithmProperty& p)
      : AbstractOptimizationAlgorithmCreator(p) {}
  virtual OptimizationAlgorithm* construct() {
    return createSolverSparseQR(property().name);
  }
};

G2O_REGISTER_OPTIMIZATION_LIBRARY(eigen_sparse_qr);

// Register solvers (similar to the dense ones)
G2O_REGISTER_OPTIMIZATION_ALGORITHM(gn_var_sparse_qr, new EigenSparseQR_SolverCreator(
    OptimizationAlgorithmProperty("gn_var_sparse_qr", "Gauss-Newton: Sparse QR solver (variable blocksize)", "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(lm_var_sparse_qr, new EigenSparseQR_SolverCreator(
    OptimizationAlgorithmProperty("lm_var_sparse_qr", "Levenberg: Sparse QR solver (variable blocksize)", "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));
// Add other registrations as needed


}  // namespace g2o


#endif //#define G2O_LINEAR_SOLVER_EIGEN_QR_HPP
