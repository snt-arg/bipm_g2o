# BIPM_g2o: Barrier Interior Point Method for Constrained Factor Graph Optimization

[![Paper](https://img.shields.io/badge/Paper-IEEE%20RA--L-blue)](https://ieeexplore.ieee.org/document/11192685)

This repository extends the `g2o` framework with support for **inequality constraints** using the **Barrier Interior Point Method (BIPM)**. It introduces a new class of factor graph nodes that implement logarithmic barrier functions, enabling efficient and robust constrained optimization for robotics control problems such as Model Predictive Control (MPC).

---

## 🧠 Project Summary

This project provides:

- A g2o extension supporting both equality and inequality constraints  
- A BIPM-based optimization backend using barrier functions  
- A demonstration on adaptive cruise control (ACC) using constrained factor graphs  
- Optional fallback to the Augmented Lagrangian (AL) method  

---

## 📚 Reference

> *Barrier Method for Inequality Constrained Factor Graph Optimization With Application to Model Predictive Control*  
> Published in IEEE Robotics and Automation Letters (RA-L), 2025.

📄 Paper: https://ieeexplore.ieee.org/document/11192685

---

## 📥 Clone the Repository

```bash
git clone https://github.com/snt-arg/bipm_g2o.git
cd bipm_g2o
```

---

## 🚀 Quickstart with Docker

### 1. Build the Docker image

```bash
docker build -t bipm_g2o -f docker/Dockerfile .
```

### 2. Run the container

```bash
docker run -it --rm bipm_g2o
```

### 3. Inside the container

```bash
cd /workspace
cmake .. -DUSE_BIPM=ON   # or OFF to use AL
make
./example
```

---

## 🔧 Manual Installation (Without Docker)

### 1. Install Dependencies (Ubuntu 22.04)

```bash
sudo apt update
sudo apt install -y \
    git cmake build-essential libeigen3-dev libspdlog-dev \
    libsuitesparse-dev qtdeclarative5-dev qt5-qmake \
    libqglviewer-dev-qt5 libmetis-dev
```

---

### 2. Clone and Build g2o

```bash
cd /opt
git clone https://github.com/RainerKuemmerle/g2o
cd g2o
git checkout master
```

---

### 3. Apply Required Patch

Edit:

```
g2o/core/sparse_optimizer.h
```

Replace:

```cpp
void update(const double* update);
```

with:

```cpp
virtual void update(const double* update);
```

---

### 4. Build and Install g2o

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DG2O_USE_OPENMP=ON \
      -DCMAKE_CXX_FLAGS="-DEIGEN_USE_THREADS -DEIGEN_USE_OPENMP" ..
make -j$(nproc)
sudo make install
```

---

### 5. Build This Project

```bash
cmake .. -DUSE_BIPM=ON   # or OFF to use AL
make
./example
```

---

## ⚙️ Configuration Flags

- `USE_BIPM=ON` → Enables Barrier Interior Point Method (default)  
- `USE_BIPM=OFF` → Uses Augmented Lagrangian method  

---

## 📁 Project Structure

```
.
├── docker/
│   ├── Dockerfile
│   └── startup.sh
├── include/bipm_g2o/
├── examples/
├── CMakeLists.txt
├── LICENSE
└── README.md
```

---

## 🧪 Application: Adaptive Cruise Control

The example demonstrates a constrained MPC formulation for adaptive cruise control using factor graphs:

- Dynamic system constraints  
- Safety distance constraints  
- Force and velocity limits  

The BIPM solver shows:

- Faster convergence  
- Improved numerical stability  
- Better scalability for long prediction horizons  

---

## 📌 Citation

If you use this work, please cite:

```bibtex
@article{abdelkarim2025barrier,
  title={Barrier Method for Inequality Constrained Factor Graph Optimization With Application to Model Predictive Control},
  author={Abdelkarim, Anas and Görges, Daniel and Voos, Holger},
  journal={IEEE Robotics and Automation Letters},
  year={2025},
  publisher={IEEE}
}
```

---

## 📜 License

This project is open-source and distributed under the GPL-3.0 License.

---

## 📫 Acknowledgments

This work builds upon the `g2o` framework and is part of ongoing research on unifying perception and control through constrained factor graph optimization.