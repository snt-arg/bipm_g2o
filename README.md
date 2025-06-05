# BIPM-g2o: Barrier Interior Point Method for Constrained Factor Graph Optimization

This repository extends the `g2o` framework with support for **inequality constraints** using the **Barrier Interior Point Method (BIPM)**. It introduces a new class of factor graph nodes that implement logarithmic barrier functions, enabling efficient and robust constrained optimization for robotics control problems such as Model Predictive Control (MPC).

## 🧠 Project Summary

This project provides:

- A g2o extension supporting both equality and inequality constraints.
- A BIPM-based optimization backend using barrier functions for constraint handling.
- A demonstration on adaptive cruise control (ACC) using constrained factor graphs.
- Optional fallback to the Augmented Lagrangian (AL) method.

### Reference

> *From Inference to Optimal Control: Barrier-Based Factor Graphs for Model Predictive Control*  
> Submitted to IEEE RAL.

This project introduces BIPM as a robust and efficient alternative to Augmented Lagrangian (AL) for constrained optimization in factor graphs.

---

## 📥 Clone the Repository

```bash
git clone https://github.com/snt-arg/bipm_g2o.git
cd bipm_g2o
```

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

## 🔧 Manual Installation Without Docker

To set up the environment without Docker, follow these steps:

1. **Install Dependencies (Ubuntu 22.04)**

   ```bash
   sudo apt update
   sudo apt install -y \
       git cmake build-essential libeigen3-dev libspdlog-dev \
       libsuitesparse-dev qtdeclarative5-dev qt5-qmake \
       libqglviewer-dev-qt5 libmetis-dev
   ```

2. **Clone and Build g2o with Patch**

   ```bash
   cd /opt
   git clone https://github.com/RainerKuemmerle/g2o
   cd g2o
   git checkout master
   ```

3. **Apply Required Code Patch**

   Edit the file `g2o/core/sparse_optimizer.h`:

   Find the line:
   ```cpp
     void update(const double* update);
   ```

   Replace it with:
   ```cpp
     virtual void update(const double* update);
   ```

   This change enables overriding the `update` method required by the BIPM extension.

4. **Build and Install g2o**

   ```bash
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DG2O_USE_OPENMP=ON -DCMAKE_CXX_FLAGS="-DEIGEN_USE_THREADS -DEIGEN_USE_OPENMP" ..
   make -j$(nproc)
   sudo make install
   ```

5. **Build Your Application**

   Assuming you're in your app's build directory:

   ```bash
   cmake .. -DUSE_BIPM=ON   # or OFF to use AL
   make
   ./example
   ```

---

## ⚙️ Configuration Flags

- `USE_BIPM=ON`: Enables the Barrier Interior Point Method (default)
- `USE_BIPM=OFF`: Enables Augmented Lagrangian method instead

---

## 📁 Project Structure

```
.
├── docker/
│   ├── Dockerfile         # Reproducible build
│   └── startup.sh         # Entry point inside container
├── src/                   # Source files
├── examples/              # Demonstrations (e.g., ACC)
├── CMakeLists.txt         # Build configuration
└── README.md
```

---

## 🧪 Application: Adaptive Cruise Control

The example application solves a multi-objective adaptive cruise control (ACC) task using factor graph optimization with:

- Dynamic constraints
- Safety distance constraints
- Force and velocity limits

The BIPM solver shows faster convergence and better scalability compared to AL, especially for long prediction horizons.

---

## 📜 License

This project is released under an open-source license. See `LICENSE` for details.

---

## 📫 Acknowledgments

This project is based on the work *“Barrier-Based Factor Graphs for Model Predictive Control,”* submitted to IEEE RAL. It builds on the g2o framework and unifies robotic perception and control via constrained factor graphs.
