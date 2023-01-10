## [Snap Buckling in Overhand Knots](https://asmedigitalcollection.asme.org/appliedmechanics/article-abstract/doi/10.1115/1.4056478/1154422/Snap-Buckling-in-Overhand-Knots)

This work studies the snap buckling process when tightening an over hand knot. Uses [Discrete Elastic Rod (DER)](http://www.cs.columbia.edu/cg/pdfs/143-rods.pdf) framework and incorporates contact and friction. Comparison between simulations and experiments is seen below in Figure 1.

<p align="center">
<img src= "comparison.png" alt>
<br>
<em> Figure 1. Comparison between simulations and experiments for different overhand knots. </em>
</p>


## How to Use

### Dependencies
Install the following C++ dependencies:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  - Eigen is used for various linear algebra operations.
  - IMC is built with Eigen version 3.4.0 which can be downloaded [here](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). After downloading the source code, install through cmake as follows.
    ```bash
    cd eigen-3.4.0 && mkdir build && cd build
    cmake ..
    sudo make install
    ```
- [Pardiso Solver Project](https://www.pardiso-project.org/)
  - Necessary for access to Pardiso, which is used as a sparse matrix solver.
  - In linux, Put the .so file under /usr/lib. In Windows, check the link attached in the pardiso solve project title.

- [OpenGL / GLUT](https://www.opengl.org/)
  - OpenGL / GLUT is used for rendering the knot through a simple graphic.
  - Simply install through apt package manager:
      ```bash
    sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
    ```
- Lapack (*usually preinstalled on your computer*)

***
### Compiling
After completing all the necessary above steps, clone the source repository of IMC and then build the project through cmake.
```bash
mkdir build && cd build
cmake ..
make -j4
```

***

### Setting Parameters

All simulation parameters are set through a parameter file ```option.txt```. A template file ```template_option.txt``` is provided that can be used to construct ```option.txt```.

```bash
cp template_option.txt option.txt   # create option.txt
```
Specifiable parameters are as follows (we use SI units):
- ```RodLength``` - Contour length of the rod.
- ```numVertices``` - Number of nodes on the rod.
- ```rodRadius``` - Cross-sectional radius of the rod.
- ```density``` - Mass per unit volume.
- ```youngM``` - Young's modulus.
- ```Poisson``` - Poisson ratio.
- ```tol``` and ```stol``` - Small numbers used in solving the linear system. Fraction of a percent, e.g. 1.0e-3, is often a good choice.
- ```maxIter``` - Maximum number of iterations allowed before the solver quits. 
- ```gVector``` - 3x1 vector specifying acceleration due to gravity.
- ```viscosity``` - Viscosity for applying damping forces.
- ```render (0 or 1) ```- Flag indicating whether OpenGL visualization should be rendered.
- ```saveData (0 or 1)``` - Flag indicating whether pull forces and rod end positions should be reocrded.
- ```dataResolution``` - Rate of data recording in seconds. Applies to both ```saveData``` and ```recordNodes```.
- ```statictime``` Initial wait period duration.
- ```pullTime``` - Duration to pull for (*starts after ```waitTime``` is done*).
- ```pullSpeed``` - Speed at which to pull and/or loosen each end.
- ```deltaTime``` - Time step size.
- ```filename``` - File name for the initial knot configuration. Should be a txt file located in ```knot_configurations``` directory. Note that overhand knot configurations for ```n1, n2, n3, n4``` are provided with a discretization of 301 nodes.

***
### Running the Simulation
Once parameters are set to your liking, the simulation can be ran from the terminal by running the provided script:
```bash
./run.sh
```
If this doesn't work, execute ```chmod +x run.sh``` prior to running.

***

### Citation
If our work has helped your research, please cite the following paper.
```

@article{tong2021snap,
  title={Snap buckling in overhand knots},
  author={Tong, Dezhong and Choi, Andrew and Joo, Jungseock and Borum, Andy and Jawed, Mohammad Khalid},
  journal={Journal of Applied Mechanics},
  pages={1--23},
  year={2021}
}

```




Compile and build:
------------------

Instructions for Ubuntu:
(1) To run this code you need Eigen, OpenGL, Pardiso and Lapack. Lapack is usually preinstalled on your computer. 
Eigen can be found at http://eigen.tuxfamily.org/index.php?title=Main_Page


(3) Open a terminal, "cd" to this folder and run the command "make" (without the quotes).

(4) Compile Command:
e.g.  g++ -fopenmp -I /usr/local/include/eigen3/  main.cpp world.cpp elasticRod.cpp elasticStretchingForce.cpp elasticBendingForce.cpp elasticTwistingForce.cpp externalGravityForce.cpp inertialForce.cpp dampingForce.cpp timeStepper.cpp setInput.cpp collision.cpp -lGL -lglut -lGLU -lpardiso600-GNU720-X86-64 -llapack -lgfortran -fopenmp -lpthread -lm -Ofast -o simDER

(5) To start the simulation, run the command "export OMP_NUM_THREADS=1;./simDER option.txt" (without the quotes). More on option.txt later.