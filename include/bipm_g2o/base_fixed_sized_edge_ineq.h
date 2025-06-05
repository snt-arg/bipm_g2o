#ifndef G2O_BASE_FIXED_SIZED_EDGE_INEQ_H
#define G2O_BASE_FIXED_SIZED_EDGE_INEQ_H

#include "g2o/core/base_fixed_sized_edge.h"
#include <functional>  // Required for std::function

namespace g2o {

// Declare the template class
template <int D, typename E, typename... VertexTypes>
class G2O_CORE_API BaseFixedSizedEdgeIneq : public BaseFixedSizedEdge<D, E, VertexTypes...> {
public: 

    // Type alias for the function pointer using std::function
    using VoidEdgeFuncType = std::function<void(BaseFixedSizedEdgeIneq<D, E, VertexTypes...>&)>;
    using InformationType = typename BaseFixedSizedEdge<D, E, VertexTypes...>::InformationType;

    BaseFixedSizedEdgeIneq();

    // Virtual method to compute inequality; must be implemented by derived classes
    virtual void computeIneq() = 0;

    // Override methods
    void computeError() override;
    void constructQuadraticForm() override;
    bool read(std::istream& is) override;
    bool write(std::ostream& os) const override;


    // Setter for the function pointers
    void setConstructQuadraticFormImpl(const VoidEdgeFuncType& func);
    void setFindLagrangeMuliplierUpdate(const VoidEdgeFuncType& func);
    void setFindLagrangeMuliplierMaxUpdateSize(const VoidEdgeFuncType& func);
    void setFindSlackVariableUpdate(const VoidEdgeFuncType& func);
    void setFindSlackVariableMaxUpdateSize(const VoidEdgeFuncType& func);


    // Setter
    void setLagrangeMultiplier(const Eigen::Matrix<double, D, 1>& lagrangeMultiplier);
    void setSlackVariable(const Eigen::Matrix<double, D, 1>& slackVariable); 
    void setRho(const Eigen::Matrix<double, D, 1>& rho);

 
    // Getter functions
    Eigen::Matrix<double, D, 1> lagrangeMultiplier() const;
    Eigen::Matrix<double, D, 1> slackVariable() const;
    Eigen::Matrix<double, D, 1> rho() const; 
    
    void setConstraintViolationPrev(Eigen::Matrix<double, D, 1>&& input, size_t start = 0);  // in AL method we use the information matrix to save the constraint violation
    
    protected:
    void setInformation(const InformationType& information) ;
    // Member to store the function pointer (or lambda)
    typename BaseFixedSizedEdge<D, E, VertexTypes...>::ErrorVector _ineq;  // Vector to store the inequality values

    VoidEdgeFuncType _constructQuadraticFormImpl = nullptr;

    Eigen::Matrix<double, D, 1> _lagrangeMultiplier; 
    Eigen::Matrix<double, D, 1> _lagrangeMultiplierUpdate; 
    double _lagrangeMultiplierMaxUpdateSize; 
    VoidEdgeFuncType findLagrangeMuliplierUpdate = nullptr;
    VoidEdgeFuncType findLagrangeMuliplierMaxUpdateSize = nullptr;

     
    // In case the slack variable is used
    Eigen::Matrix<double, D, 1> _slackVariable;            
    Eigen::Matrix<double, D, 1> _slackVariableUpdate;
    double _slackVariableMaxUpdateSize;
    VoidEdgeFuncType findSlackVariableUpdate = nullptr;
    VoidEdgeFuncType findSlackVariableMaxUpdateSize = nullptr;

    // In case Rho is used
    Eigen::Matrix<double, D, 1> _rho;
    
 
 
    

 
     

};

// Forward declare the customQuadraticForm function template
template <int D, typename E, typename... VertexTypes>
void customQuadraticForm(BaseFixedSizedEdgeIneq<D, E, VertexTypes...>& edge);

}  // namespace g2o

namespace g2o {

template <int D, typename E, typename... VertexTypes>
BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::BaseFixedSizedEdgeIneq() 
    : BaseFixedSizedEdge<D, E, VertexTypes...>(), _constructQuadraticFormImpl(nullptr) {
        this->_information.setZero();
        // Constructor body (if needed)
}



// Override the constructQuadraticForm method
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::constructQuadraticForm() {
    // Ensure that a function has been assigned before calling it
    if (_constructQuadraticFormImpl) {
        _constructQuadraticFormImpl(*this);  // Call the assigned function (or lambda)
    } else {
        // If no function is assigned, call the base class implementation (if needed)
        BaseFixedSizedEdge<D, E, VertexTypes...>::constructQuadraticForm();
    }
}


template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::computeError() {
    this->computeIneq();
    this->_error = _ineq;
}


template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::write(std::ostream& os) const {
    for (int i = 0; i < D; ++i) {
        os << this->_error[i] << " ";
    }
    return os.good();
}

template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::read(std::istream& is) {
    for (int i = 0; i < D; ++i) {
        is >> this->_error[i];
    }
    return is.good();
}



// Setter functions pointers 
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setInformation(const InformationType& information) {
    this->_information = information;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setConstraintViolationPrev(Eigen::Matrix<double, D, 1>&& input, size_t start) {
      this->_information.diagonal().segment(start, D) = input; 
}  
  



template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setConstructQuadraticFormImpl(const VoidEdgeFuncType& func) {
    _constructQuadraticFormImpl = func;
}
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setFindLagrangeMuliplierUpdate(const VoidEdgeFuncType& func) {
    findLagrangeMuliplierUpdate = func;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setFindLagrangeMuliplierMaxUpdateSize(const VoidEdgeFuncType& func) {
    findLagrangeMuliplierMaxUpdateSize = func;
}

template <int D, typename E, typename... VertexTypes> 
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setFindSlackVariableUpdate(const VoidEdgeFuncType& func) {
    findSlackVariableUpdate = func;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setFindSlackVariableMaxUpdateSize(const VoidEdgeFuncType& func) {
    findSlackVariableMaxUpdateSize = func;
}


// Setter functions
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setLagrangeMultiplier(const Eigen::Matrix<double, D, 1>& lagrangeMultiplier) {
    _lagrangeMultiplier = lagrangeMultiplier;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setSlackVariable(const Eigen::Matrix<double, D, 1>& slackVariable) {
    _slackVariable = slackVariable;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setRho(const Eigen::Matrix<double, D, 1>& rho) {
    _rho = rho;
}


//Getter functions

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::lagrangeMultiplier() const {
    return _lagrangeMultiplier;
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::slackVariable() const {
    return _slackVariable;
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::rho() const {
    return _rho;
}



}  // namespace g2o
#endif  // G2O_BASE_FIXED_SIZED_EDGE_INEQ_H