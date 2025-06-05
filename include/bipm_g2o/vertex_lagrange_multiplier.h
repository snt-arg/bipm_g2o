#ifndef VERTEX_LAGRANGE_MULTIPLIER_H
#define VERTEX_LAGRANGE_MULTIPLIER_H

#include <g2o/core/base_vertex.h>
#include <Eigen/Core>
#include <iostream>

namespace g2o {

template <int D>
class VertexLagrangeMultiplier : public BaseVertex<D, Eigen::Matrix<double, D, 1>> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VertexLagrangeMultiplier();

    // Set the vertex state to the origin
    virtual void setToOriginImpl() override;

    // Apply an increment (delta) to the state
    virtual void oplusImpl(const double* update) override;

    // Write the vertex state to an output stream
    virtual bool write(std::ostream& os) const override;

    // Read the vertex state from an input stream
    virtual bool read(std::istream& is) override;
    

    void setInitialValue(double initialValue);      
    double getInitialValue() const ;
    
    
    private:
   
    double _initialValue = 0;
};

} // namespace g2o




namespace g2o {

template <int D> 
void VertexLagrangeMultiplier<D>::setInitialValue(double initialValue) {
        _initialValue = initialValue;
        setToOriginImpl();
    }
 

template <int D> 
double VertexLagrangeMultiplier<D>::getInitialValue() const { 
        return _initialValue;
    }

  
template <int D>
VertexLagrangeMultiplier<D>::VertexLagrangeMultiplier() {
    setToOriginImpl();    
}   

template <int D>
void VertexLagrangeMultiplier<D>::setToOriginImpl() {
    this->_estimate.setConstant(_initialValue);
}

template <int D>
void VertexLagrangeMultiplier<D>::oplusImpl(const double* update) {
    Eigen::Map<const Eigen::Matrix<double, D, 1>> delta(update);
    this->_estimate += delta;
}

template <int D>
bool VertexLagrangeMultiplier<D>::write(std::ostream& os) const {
    for (int i = 0; i < D; ++i) {
        os << this->_estimate[i] << " ";
    }
    return os.good();
}

template <int D>
bool VertexLagrangeMultiplier<D>::read(std::istream& is) {
    for (int i = 0; i < D; ++i) {
        is >> this->_estimate[i];
    }
    return is.good();
}



} // namespace g2o













#endif // VERTEX_LAGRANGE_MULTIPLIER_H
