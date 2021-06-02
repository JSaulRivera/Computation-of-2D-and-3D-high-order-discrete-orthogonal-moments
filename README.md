# Computation-of-2D-and-3D-high-order-discrete-orthogonal-moments
 The content of this repository is related to the project that helps eliminate numerical instability and high-order orthogonal moment calculation error while preserving it's orthogonality, which is one of the most important properties of these types of functions.
### Kernel computation of discrete moments

The high-order orthogonal polynomial base causes a numerical instability that makes it difficult to calculate the moments. In this project, we analyze the calculation of the following discrete orthogonal polynomials,
  
* Tchebycheff polynomials
* Krawtchouk polynomials
* Hahn polynomials
* Meixner polynomials
* Charlier polynomials
  
The source codes of the project are described in the chapter 3: "Computation of 2D and 3D high-order discrete orthogonal moments" of the book "Recent Progress in Moments and Moment Invariants" GCSR Volume 7 Edited by Dr. George A. Papakostas. 

If you find these codes useful as a support tool.

You can access the information about the book through the following link:
https://sciencegatepub.com/sgp-books/gcsr/gcsr_vol7/

or you can directly download the chapter using the following link:
https://sciencegatepub.com/?smd_process_download=1&download_id=3179

### Please cite the following documents as:

*   José S. Rivera-Lopez, César Camacho-Bello, and Lucia Gutiérrez-Lazcano, Chapter 3. “Computation of 2D and 3D High-order Discrete Orthogonal Moments”. Recent Progress in Image     Moments and Moment Invariants, GCSR Volume 7 (2021), 53-74,  DOI: 10.15579/gcsr.vol7.ch3.

        Bibtex:
        @article{Rivera2021Recent,
        title={Recent Progress in Image Moments and Moment Invariants},
        author={José S. Rivera-Lopez and César Camacho-Bello and Lucia Gutiérrez-Lazcano},
        journal={Science Gate},
        volume={7},
        pages={53--74},
        year={2021},
        URL={https://sciencegatepub.com/sgp-books/gcsr/gcsr_vol7/},
        ISBN={2241-9063}
        }

*    C. Camacho-Bello and J. S. Rivera-Lopez, “Some computational aspects of tchebichef moments for higher orders,”Pattern Recognition Letters, vol. 112, pp. 332–339, 2018.
   
         Bibtex:
         @article{camacho2018some,
           title={Some computational aspects of Tchebichef moments for higher orders},
           author={Camacho-Bello, C{\'e}sar and Rivera-Lopez, Jos{\'e} S},
           journal={Pattern Recognition Letters},
           volume={112},
           pages={332--339},
           year={2018},
           publisher={Elsevier}
          }

# Getting Started

1. Download this repository.
2. To obtain discrete ortogonal moments of image test  run the file "Discrete_orthogonal_moments.m" using Matlab.
3. By default the value assigned for the order of calculation of the moment is the maximum possible. Change the value of "n" in line 66 from "Discrete_orthogonal_moments.m", to obtain different orders of discrete orthogonal moments form image test.
4. By default set discrete orthogonal polynomials is Tchebichef polynomial. For obtain different set polynomial modify line 69 from "Discrete_orthogonal_moments.m" corresponding to the calculation of the polynomial base.
