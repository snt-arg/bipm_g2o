#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <functional>
#include "optimizable_graph.h"
 
namespace g2o {


class TerminationCriteria {
public:
    enum class ConvergenceCriterion {
        UpdateNorm = 0,      
    };


    // Public termination thresholds (accessible in optimizer)
    double epsilon_convergence = 1e-2;    // Stopping criterion for update norm
    double epsilon_constraint = 1e-2;  // Stopping criterion for individual constraint check

private:
    ConvergenceCriterion _convergenceCriterion = ConvergenceCriterion::UpdateNorm; 

    using ConvergenceFunc = std::function<bool(const Eigen::VectorXd&, const Eigen::VectorXd&, double)>;
    using ConstraintFeasibilityFunc = std::function<bool(double)>;

    std::unordered_map<ConvergenceCriterion, ConvergenceFunc> ConvergenceStrategies;

public:
    TerminationCriteria() {
        // Initialize solution acceptance strategies
        ConvergenceStrategies[ConvergenceCriterion::UpdateNorm] = 
            [this](const Eigen::VectorXd&, const Eigen::VectorXd& updateVec, double epsilon) {
                if (epsilon <0) epsilon =  std::abs(epsilon) * epsilon_convergence;
                return updateVec.norm() < epsilon;
            }; 


    }  

    // Public methods to evaluate termination conditions
    bool verifyConvergence(const Eigen::VectorXd& bVec, const Eigen::VectorXd& updateVec, double epsilon) {
        auto it = ConvergenceStrategies.find(_convergenceCriterion);
        return it != ConvergenceStrategies.end() && it->second(bVec, updateVec,epsilon);
    }

    bool verifyEqFeasibility(OptimizableGraph::EdgeContainer& edgeVec, double epsilon) {
        if (epsilon <0) epsilon = std::abs(epsilon) * epsilon_constraint;

        for (auto& edge : edgeVec) {
            edge->computeError();             
            const double* error = edge->errorData();   
            int edgeDimension = edge->dimension() /2 ; // divide by 2 to account for equality constraints excluding the lagrange multiplier
            for (int i = 0; i < edgeDimension; ++i) {
                if (std::abs(error[i]) > epsilon) {
                    return false;  // If any error exceeds epsilon, return false
                }
            }
        }
        
        return true; // Return true only if all edges satisfy the condition
    }

    
    bool verifyIneqFeasibility(OptimizableGraph::EdgeContainer& edgeVec, double epsilon) {
        if (epsilon <0) epsilon = std::abs(epsilon) * epsilon_constraint;

        for (auto& edge : edgeVec) {
            edge->computeError();             
            const double* error = edge->errorData();   
            int edgeDimension = edge->dimension();
            for (int i = 0; i < edgeDimension; ++i) {
                if (std::max(0.0, error[i]) > epsilon) {
                    return false;  // If any error exceeds epsilon, return false
                }
            }
        }
        
        return true; // Return true only if all edges satisfy the condition
    }

 

    // Setters to configure termination criteria
    void setConvergenceCriteria(ConvergenceCriterion type) {
        _convergenceCriterion = type;
    }

    void setConvergenceCriterion(int criterion) {
  /*
  the termination criterion is based on the Newton decrement, the gradient norm, or the update norm
  where GradientNorm =0
  
  */
    _convergenceCriterion = static_cast<ConvergenceCriterion>(criterion);
    
 if (criterion == 0) {
      std::cout
          << "The stopping criterion of the update convergence is Gradient Norm."
          << std::endl;
    } else {
      throw std::invalid_argument(
          "Invalid stopping criterion. Valid options are 0:Gradient Norm, 1:Newton Decrement, or 2:Update Norm.");
    }

 };
 
 

};

}  // namespace g2o
