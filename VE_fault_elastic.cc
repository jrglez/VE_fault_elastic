/* VE_fault_elastic
 * Copyright (c) 2016, Juan Rodriguez-Gonzalez, University of Maryland
 * https://github.com/jrglez/VE_fault_elastic
 * jrglez@umd.edu

 * VE_fault_elastic was developed by
 * Juan Rodríguez-González and Laurent G. J. Montesi.
 * All work was funded by the National Science Foundation (#1419826).

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

 * ---------------------------------------------------------------------
 *
 * This file has been created using the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

/* TO DO:
 * - Obtain the integration constant A from Turcottes expression.
 *   + Obtain the analytical solution for viscoelastic using Laplace transform
 * - Decide on time step (adaptive or fixed) and modify comments accordingly
 * - Repalce:
 *   + velocity (v or V) --> displacement (u or U)
 *   + stress --> viscous strain (v_strain)
 *   + remove old_solution?
 */

/*! \mainpage
 *
 * \warning
 * This code is under development and has not been properly benchmarked.
 * Feel free to use it, but be warned that it is distributed WITHOUT ANY WARRANTY.
 * Please, contact Juan Rodriguez-Gonzalez (jrglez@umd.edu) if you make any
 * significant modifications to this code.
 *
 * \section intro Introduction
 * This code simulates the deformation of a viscoelastic medium in presence
 * of a discrete fault under the anti-plane shear approximation.
 * Under this assumption there is no displacement in the modeled plane
 * \f$(u_{x}=u_{y}=0)\f$ and there is no spatial variation in the direction
 * perpendicular to the modeled plane \f$(\partial z=0)\f$. We solve the equilibrium
 * equation for an elastic medium in with an isotropic constitutive relationship
 * in which the elastic strain is relaxed by viscous creep. Under this conditions,
 * the right-hand side depends on the displacement values from the previous time step.
 *
 *
 * \subsection equations Viscoelastic approximation
 *
 * We solve the equation of equilibrium for stress in a continuum medium:
 *
 * \f[\partial_{i} \sigma = f^g_j \f]
 *
 *
 * where \f$i\f$ and \f$j\f$ run from 1 to 3 and designate spatial coordinate,
 * \f$\sigma\f$ is the stress tensor and \f$\boldsymbol{f^g}\f$ is the total
 * external force, in this case, only the gravitational force
 * \f$f^g_j=f^g_j=\rho g \left(1-\alpha T \right) \delta_{iy}\f$.
 * We assume a linear maxwel model in which the total strain \f$ \varepsilon \f$
 * is the sum of the elastic strain \f$ \varepsilon^e \f$ and the viscous strain
 * \f$ \varepsilon^v \f$:
 *
 * \f[ \varepsilon_{ij} =
 * \varepsilon^e_{ij} + \varepsilon^v_{ij} + \varepsilon^0_{ij} \f]
 *
 * where \f$ \varepsilon^0 \f$ is the initial strain. We will use a semi-elastic
 * approach \cite zienkiewicz_cormeau_74 \cite yamasaki_houseman_12 for which
 * we assume a constitutive relation between stress and elastic strain:
 *
 * \f[ \sigma_{ij} = C_{ijkl} \varepsilon_{kl} + \sigma^0_{ij}\f]
 *
 * where \f$ C \f$ is the elastic stiffness tensor and \f$ \sigma^0 \f$ is the
 * initial stress. This can be interpreted as the stress being relaxed by viscous creep:
 *
 * \f[ \sigma_{ij} = C_{ijkl} \left( \varepsilon_{kl} -
 * \varepsilon^v_{kl} -\varepsilon^0_{kl} \right) + \sigma^0_{ij} \f]
 *
 * In an isotropic medium:
 *
 * \f[ \sigma_{ij} = \lambda \left( \epsilon_{kk} - \epsilon^0_{kk} \right) \delta_{ij} +
 * 2\mu \left( \varepsilon_{ij} - \varepsilon^v_{ij} - \varepsilon^0_{ij} \right)
 * + \sigma^0_{ij} \f]
 *
 * where \f$ \lambda \f$ is the first Lamé parameter and \f$\mu\f$ is the shear modulus. The
 * strain can be expressed in terms of the displacement \f$ u_i \f$:
 *
 * \f[ \varepsilon_{ij} =
 * \frac{1}{2} \left( \partial_j u_i + \partial_i u_j \right) \f]
 *
 * The viscous strain rate \f$ \dot \varepsilon^v \f$ depends on the stress tensor:
 *
 * \f[ \dot \varepsilon^v_{ij} = \beta_{ij}(\sigma) \f]
 *
 * where \f$ \beta \f$ are given functions. A relationship to link the viscous strain
 * rate and the elastic strain is given in the next section.
 *
 * We can now use this in the momentum equation and we have:
 *
 * \f[ \partial_i \left( \lambda \varepsilon_{kk} \delta_{ij} +
 * 2 \mu \varepsilon_{ij} \right) = f_j + f^v_j \f]
 *
 * where
 *
 * \f[ f_j =
 * f^g_j +
 * \partial_i \left( \lambda \varepsilon^0_{kk}\delta_{ij} + 2 \mu \varepsilon^0_{ij} -
 * \sigma^0_{ij} \right) \f]
 *
 * This equation is very similar to the elasticity equation but with an internal force
 * due to the viscous coupling:
 *
 * \f[ f^v_j = \partial_i \left( 2 \mu \varepsilon^v_{ij} \right) \f]
 *
 * \subsection num_approach Numerical approach
 * To use a semi-elastic approach we need to manipulate the constitutive equation
 * to find a relation between viscous and elastic strains. To achieve this, the first step
 * is to express the viscous strain in terms of the viscous strain rate. We replace the
 * strain rate by the forward derivative:
 *
 * \f[ \dot \varepsilon^v(t_n) = \frac{\varepsilon^v(t_{n+1}) - \varepsilon^v(t_n)}{\Delta t_n} \f]
 *
 * where the subscripts \f$ n \f$ and \f$ n+1 \f$ indicate the current and future time
 * steps and \f$ \Delta t_{n+1} = t_{n+1} - t_n \f$ is the time step increment. Therefore,
 * we can write the viscous strain at a given time step as:
 *
 * \f[ \varepsilon^v(t_{n+1}) = \varepsilon^v(t_{n}) + \Delta t_{n+1} \dot \varepsilon(t_{n}) \f]
 *
 * Following \cite yamasaki_houseman_12, we use the following relation between viscous
 * strain rate and elastic strain:
 *
 * \f[ \dot \varepsilon^v_{ij} = \frac{\mu}{\eta}
 * \left( \varepsilon_{ij} -
 * \frac {1}{3} \epsilon_{kk}\delta_{ij} \right) \f]
 *
 * where \f$ \eta \f$ is the viscosity. And at time step n+1 we will need to solve the
 * equation:
 *
 * \f[ \partial_i \left( \lambda \varepsilon_{kk}(t_{n+1}) \delta_{ij} +
 * 2 \mu \varepsilon_{ij}(t_{n+1}) \right) = f^g_j(t_{n+1}) + f^{v}_j(t_{n+1}) \f]
 *
 * where the internal elastic force now only depends on the strain and past viscous strain
 * rate:
 *
 * \f[ f^{v}_j(t_{n+1}) = \partial_i \left( 2\mu \varepsilon^v_{ij}(t_{n+1}) \right) =
 * \partial_i \left\{ 2\mu \left[\varepsilon^v_{ij}(t_n) +
 * \frac{\mu \Delta t_{n+1}}{\eta} \left(\varepsilon_{ij}(t_n) -
 * \frac{1}{3} \varepsilon_{kk}(t_n) \delta_{ij} \right) \right] \right\} \f]
 *
 *
 * Therefore, we will follow the next scheme to solve the problem:
 *   -# Set initial conditions \f$ (n = 0) \f$ for the stress \f$\sigma(t_0) = \sigma^0\f$,
 *   displacement \f$ \boldsymbol{u}(t_0) \f$ and viscous strain
 *   \f$ \varepsilon^v(t_0) \f$.
 *   -# Compute the viscous strain \f$ \varepsilon^v(t_{n+1}) \f$, which depends on
 *   \f$ \varepsilon(t_n) \f$ and \f$ \varepsilon^v(t_n) \f$.
 *   -# Solve the momentum equation to obtain the displacement
 *   \f$ \varepsilon (t_{n+1}) \f$
 *   -# Compute \f$ \sigma (t_{n+1}) \f$.
 *   -# Update \f$ n \leftarrow n+1 \f$ and repeat steps 2 to 5.
 *
 *
 *
 * \subsubsection anti_plane Anti-plane shear approximation
 * The equations presented until now are valid for viscoelastic problems in any
 * dimensions. However, solving it can be computationally demanding. For certain
 * problems it is enough to solve the equation under certain approximations.
 * In this case we are going to model the displacement in the \f$x-y\f$ plane
 * produced by a fault that runs parallel to the \f$y-z\f$ plane and is sliding
 * parallel to the \f$z\f$ direction.
 * \image html aps.png "Figure 1: Schematic representation of a fault and the modeled plane"
 * \image latex aps.eps "Figure 1: Schematic representation of a fault and the modeled plane"
 *
 * In this case, we can considered that displacement only occurs in the \f$z\f$
 * direction \f$ (u_{x}=u_{y}=0) \f$. Moreover, under this approximation no
 * magnitude can vary along the \f$ z \f$ direction \f$(\frac{\partial}{\partial z}=0)\f$.
 * Therefore, under this approximation
 * \f$ \varepsilon_{xx} = \varepsilon_{yy} = \varepsilon_{zz} = \varepsilon_{xy} =
 * \varepsilon_{yx}=0\f$, * \f$ \varepsilon^v_{xx} = \varepsilon^v_{yy} = \varepsilon^v_{zz}
 * = \varepsilon^v_{xy} = \varepsilon^v_{yx} = 0 \f$ and \f$ \sigma_{xx} = \sigma_{yy} =
 * \sigma_{zz} = \sigma_{xy} = \sigma_{yx} = 0 \f$. In addition, we neglect the gravitational
 * force. Therefore, the governing equations are reduced to:
 *
 * \f[ \partial_i \left( \mu \partial_i u_z \right) = f^v_z \f]
 *
 * where
 *
 * \f[ f^v_z(t_{n+1}) = \partial_i \left( 2\mu \varepsilon^v_{iz}(t_{n+1}) \right) =
 * \partial_i \left[ 2 \mu \left(\varepsilon^v_{iz}(t_n) +
 * \frac{\mu \Delta t_{n+1}}{2 \eta} \partial_i u_z(t_n) \right) \right] \f]
 *
 *
 * \section model Modeling
 * \subsection setup Setup
 * Here we will solve the visco-elastic equation in a two-dimensional domain
 * of dimensions \f$a\f$ by \f$b\f$. We will solve for the displacement in the
 * direction perpendicular to the plane \f$(u_z)\f$ caused by a fault also
 * perpendicular to the modeled plane.
 *
 * There are four possible boundary conditions:
 *
 * - \f$\Gamma_0\f$ - The locked part above the fault, with imposed zero displacement, \f$u_z=0\f$
 * (homogeneous Dirichlet boundary conditions).
 *
 * - \f$\Gamma_1\f$ - Regions of imposed displacement, \f$u_z=U\f$ (non-homogeneous Dirichlet
 * boundary conditions).
 *
 * - \f$\Gamma_2\f$ - Zero tangential stress boundaries, \f$ n_i \sigma_{iz} = 0 \f$
 * (homogeneous Neumann boundary conditions).
 *
 * - \f$\Gamma_3\f$ - Imposed tangential stress \f$ n_i \sigma_{iz} = S\f$
 * (non-homogeneous Neumann boundary conditions).
 *
 *
 * \subsection weak_form The Weak Form
 *
 * In order to solve the equation using the Finite Element Method we need to derive
 * the weak form. We will begin with the partial differential equation that governs
 * the dicplacement in a visco-elastic medium under the anti-plane shear approximation:
 *
 * \f[ \partial_i \left( \mu \partial_i u_z\right) = f^v_z \f]
 *
 * First we multiply the equation from the left by a test function \f$v\f$,
 * and then we integrate over the whole domain:
 *
 * \f[ \int_{\Omega} v \cdot \partial_i \left( \mu \partial_i u_z \right) =
 * \int_{\Omega} v \cdot \partial_i \left( 2 \mu \varepsilon^v_{iz} \right) \f]
 *
 * Integrating by parts both sides of the equation and using the Gauss theorem
 *
 * \f[ \int_{\Omega} \partial_i v \cdot \mu \partial_i u_z =
 * \int_{\Omega} \partial_i v \cdot 2 \mu \varepsilon^v_{iz}  +
 * \int_{\partial \Omega} v \cdot \mu n_i \partial_i u_z -
 * \int_{\partial \Omega} v \cdot n_i 2 \mu \varepsilon^v_{iz} \f]
 *
 * Considering that:
 *
 * \f[ \sigma_{iz} =
 * 2 \mu \left( \varepsilon_{iz} - \varepsilon^v_{iz} \right) \f]
 *
 * The boundary term can be rewriten as:
 *
 * \f[ \int_{\Omega} \partial_i v \cdot \mu \partial_i u_z =
 * \int_{\Omega} \partial_i v \cdot 2 \mu \varepsilon^v_{iz}  +
 * \int_{\partial \Omega} v \cdot n_i \sigma_{iz} \f]
 *
 * or:
 *
 * \f[ \left( \partial_i v, \mu \partial_i u_z \right)_{\Omega} =
 * \left( \partial_i v, 2 \mu \varepsilon^v_{iz} \right)_{\Omega} +
 * \left(v, n_i \sigma_{iz} \right)_{\partial \Omega} \f]
 *
 * where the operator \f$(a, b)\f$ is \f$ \int a \cdot b \f$.
 *
 * In those boundaries where the
 * displacement is imposed (Dirichlet boundary conditions), this is \f$\Gamma_0\f$ and
 * \f$\Gamma_1\f$,the test function \f$v = 0\f$. The contribution of those boundaries in
 * which the imposed tangential stress is zero (i.e., \f$\Gamma_2\f$) will not contribute to
 * the right hand side either. On the other hand,those boundaries in which we impose a
 * non-zero tangential stress (non-homogeneous Neumann boundary conditions), this is
 * \f$\Gamma_3\f$, the test functions \f$v \neq 0\f$, and that part of the boundary integral
 * has to be included in the weak form.
 *
 * \f[ \left( \partial_i v, \mu \partial_i u_z \right)_{\Omega} =
 * \left( \partial_i v, 2 \mu \varepsilon^v_{iz} \right)_{\Omega} +
 * \left(v, S \right)_{\Gamma^3} \f]
 *
 * where \f$S\f$ is the imposed tangential stress at the boundares \f$ \Gamma_3 \f$.
 *
 * We can now  replace the exact solution \f$ u_z \f$ for an approximate solution in terms
 * of the shape functions, \f$ \sum_{k} U_k v_k(x_i) \f$. Here \f$U_k\f$ are the
 * expansion coefficients we need to determine to find the approximate solution to the equation.
 * Then, the problem is reduced to solving
 *
 * \f[AU=B\f]
 *
 * where the matrix A and the vector B are defined as:
 *
 * \f[A_{kl} = (\partial_i v_k, \mu \partial_i v_l)_{\Omega}\f]
 * \f[B_k = (\partial_i v_k, 2 \mu \varepsilon^v_{iz})_\Omega
 * + (v_k, S)_{\Gamma_3} \f]
 *
 * where \f$ k \f$ and \f$ l \f$ run through all the degrees of freedom.
 *
 *
 * \subsection adaptive_refinement Adaptive Mesh Refinement in Time Dependent Problems
 * Every time step we need to use the solution and the viscous strain from the previous
 * time step. The old solution is stored in a
 * <a href="https://www.dealii.org/8.4.0/doxygen/deal.II/classVector.html">Vector</a> with as
 * many elements as the degrees of freedom the problem has, and is updated every time step.
 * The old viscous strain is stored at each quadrature and can be accessed through a user
 * pointer that each cell holds.
 *
 * Every time the mesh is refined, the number of cells and degrees of freedom changes,
 * therefore, as part of the mesh refinement we need to transfer the solution vector and the
 * old viscous strain from the old mesh to the new one. For the solution vectors we just need to use
 * <a href="https://www.dealii.org/8.4.0/doxygen/deal.II/classSolutionTransfer.html">SolutionTransfer</a>
 * (see for example
 * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_26.html#codeHeatEquationrefine_meshcode">Step-26</a>).
 * Transferring the old viscous strin is slightly more complicated, we first need to transfer the data
 * stored in the quadrature points to a finite element field that is defined everywhere so that we can
 * later transfer it to the new mesh (using
 * <a href="https://www.dealii.org/8.4.0/doxygen/deal.II/classSolutionTransfer.html">SolutionTransfer</a>
 * and then interpolate it to the new quadrature points. We need a discontinuous field that
 * matches the values in the quadrature points (we will use a Discontinuous Galerking finite element
 * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/classFE__DGQ.html">FE_DGQ</a>).
 */


/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth and Ralf Hartmann, University of Heidelberg, 2000
 */



#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

#include <typeinfo>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <math.h>

#define PI 3.14159265

namespace vsf
{
  using namespace dealii;

  /**
   * \struct PointHistory
   * We first define a structure where we will store the old data that we will need
   * to compute the right hand side of the equation. In this case it is enough if we save
   * the viscous strain from the previous time step.
   */
  template <int dim>
  struct PointHistory
  {
      Tensor<1,2>       new_v_strain;
  };
  /**
   * We declare and define some function classes that represent the body force, the
   * viscosity and the boundary values. For simplicity we have chosen to define all relevant parameters
   * in a base class, which later will serve as the base class for the rest of the classes.
   * These parameters are the width and height of the model (FunctionBase::width and
   * FunctionBase::height), lock depth \f$d\f$ (FunctionBase::locked_depth),
   * the imposed displacement at the boundary \f$U\f$ (FunctionBase::boundary_displacement) and the imposed
   * stress at the boundary \f$S\f$ (FunctionBase::boundary_stress).
   *
   * \note These parameters have to be hardcoded here because the solution and boundary
   * conditions classes are derived from
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/classFunction.html">Function</a>
   * and are used directly by other deal.II class
   * functions. Adding any extra input parameters in a class that inherits from
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/classFunction.html">Function</a>
   * would require modifying the source code. Another option would be to input these parameter
   * via an input file (see deal.II tutorial
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_33.html">step-33</a>).
   *
   * All the parameters used in the code are set here to keep things in order.
   */

  template <int dim>
  class FunctionBase
  {
  protected:
    static const double width; /**< Model's width. Will not be used here, but we want to keep input parameters together.*/
    static const double height; /**< Model's height.*/
    static const double locked_depth; /**< Depth of the locked region.*/
    static const double boundary_displacement; /**< Displacement imposed in some boundaries.*/
    static const double boundary_stress; /**< Stress imposed in some boundaries.*/
    static const double shear_modulus; /**< Reference shear modulus.*/
    static const double viscosity; /**< Reference viscosity.*/
  };

  template <int dim>
  const double FunctionBase<dim>::width = 2.0e5;

  template <int dim>
  const double FunctionBase<dim>::height = 1.0e5;

  template <int dim>
  const double FunctionBase<dim>::locked_depth = 0.25e5;

  template <int dim>
  const double FunctionBase<dim>::boundary_displacement = 10.0;

  template <int dim>
  const double FunctionBase<dim>::boundary_stress = 1.0e5;

  template <int dim>
  const double FunctionBase<dim>::shear_modulus = 1.0e11;

  template <int dim>
  const double FunctionBase<dim>::viscosity = 1.0e20;



  /** This class is used to compute the analytical solution of the surface displacement caused
   * by a fault in a completely elastic medium. we use Turcotte and Spence viscouse model\cite turcotte_spence_74.
   * and applyt the correspondence principle to obtain the viscoelastic solution.
   */
  template <int dim>
  class Solution_Turcotte : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Solution_Turcotte () : Function<dim>() {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
    virtual void get_time(const double &t);
  private:
    double time; /**< Time (s) at which the displacement is calculated.*/
  };

//  template <int dim>
//  double Solution_Turcotte<dim>::time = 0.0;

  /** Solution_Turcotte::value extracts the value of the analytical solution in one point ("p").
   */
  template <int dim>
  double Solution_Turcotte<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
  {
    /** This solution is only valid at the surface, therefore, if the point at which it
     * will be evaluated is not in the surface, an exception will be thrown.
     */
    Assert (std::fabs(p[1] - 0)<1e-12, ExcNotImplemented ());

    const double c_1 = PI/(2*this->height);
    const double c_2 = sin (c_1 * this->locked_depth);
    const double A = this->boundary_stress / (c_1 * this->shear_modulus);
    const  double u_x = A * log((sinh(c_1*p[0]) + sqrt(pow(sinh(c_1*p[0]),2) + pow(c_2,2))) / c_2);
    return u_x *(1 + time * this->shear_modulus / this->viscosity);
  }

  /** Solution_Turcotte::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Solution_Turcotte<dim>::value_list (const std::vector<Point<dim> >   &points,
                                           std::vector<double>              &values,
                                           const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Solution_Turcotte::value(points[i]);
  }

  /** Solution_Turcotte::get_time takes the given time and use it to populate the variable "time".
    */
  template <int dim>
  void Solution_Turcotte<dim>::get_time (const double &t)
  {
    time = t;
  }



  /** This class is used to compute the analytical solution of the surface displacement caused
   * by a fault in a completely elastic medium using \cite savage_burford_73.
   */
  template <int dim>
  class Solution_Savage : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Solution_Savage () : Function<dim>() {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Solution_Savage::value extracts the value of the analytical solution in one point ("p").
   */
  template <int dim>
  double Solution_Savage<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
  {
    /** This solution is only valid at the surface, therefore, if the point at which it
     * will be evaluated is not in the surface, an exception will be thrown.
     */
    Assert (std::fabs(p[1] - 0)<1e-12, ExcNotImplemented ());

    return (this->boundary_displacement/PI)*atan(p[0]/(this->locked_depth));
  }

  /** Solution_Savage::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Solution_Savage<dim>::value_list (const std::vector<Point<dim> >   &points,
                                  std::vector<double>              &values,
                                  const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Solution_Savage::value(points[i]);
  }



  /**
  * The template class Shear_Modulus gives the value of the shear modulus in the requested point/s
  */
  template <int dim>
  class Shear_Modulus : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Shear_Modulus () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Shear_Modulus::value extracts the value in one point ("p"). In the current implementation
   * we consider a uniform shear modulus. To use a spatially variable shear modulus modify this
   * method accordingly.
   */
  template <int dim>
  double Shear_Modulus<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
     return this->shear_modulus;
  }

  /** Shear_Modulus::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Shear_Modulus<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Shear_Modulus::value(points[i]);
  }



  /**
  * The template class Viscosity gives the value of the viscosity in the requested point/s
  */
  template <int dim>
  class Viscosity : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Viscosity () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Viscosity::value extracts the value in one point ("p"). In the current implementation
   * we consider a uniform shear modulus. To use a spatially variable shear modulus modify this
   * method accordingly.
   */
  template <int dim>
  double Viscosity<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
     return this->viscosity;
  }

  /** Viscosity::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Viscosity<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Viscosity::value(points[i]);
  }



   /**
   * The template class BodyForce gives the value of body force (on the right hand side of the
   * elastic equation) in the requested point/s
   */
   template <int dim>
   class BodyForce : public Function<dim>,
   protected FunctionBase<dim>
   {
   public:
     BodyForce () : Function<dim>() {}

     virtual double value (const Point<dim>   &p,
                           const unsigned int  component = 0) const;
     virtual void value_list (const std::vector<Point<dim> >      &points,
                              std::vector<double>                 &values,
                              const unsigned int                  component = 0) const;
   };

   /** RightHandSide::value extracts the value in one point ("p"). In the current implementation
   * we consider a zero body force. To use a spatially variable shear modulus modify this
   * method accordingly.
   */
   template <int dim>
   double BodyForce<dim>::value (const Point<dim>   &p,
                                     const unsigned int) const
   {
      return 0.0;
   }

   /** BodyForce::value_list allows to extract values for several points ("points") at once and
    * outputs them in a vector ("values").
    */
   template <int dim>
   void BodyForce<dim>::value_list (const std::vector<Point<dim> >   &points,
                                std::vector<double>              &values,
                                const unsigned int               component) const
   {
     const unsigned int n_points = points.size();
     Assert (values.size() == n_points,
             ExcDimensionMismatch (values.size(), n_points));
     Assert (component == 0,
             ExcIndexRange (component, 0, 1));

     for (unsigned int i=0; i<n_points; ++i)
       values[i] = BodyForce::value(points[i]);
   }



  /**
  * The template class Dirichlet_BC computes the values of the non-homogeneous Dirichlet BC
  * in each point.
  */
  template <int dim>
  class Dirichlet_BC : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
      Dirichlet_BC () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Dirichlet_BC::value extracts the value in one point ("p"). To use a spatially variable BC
   * modify this method accordingly.
   */
  template <int dim>
  double Dirichlet_BC<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
    /** The displacement imposed at the boundary is only half the total displacement because we are only
     * modeling half the domain.
     */
    return this->boundary_displacement/2;
  }

  /** Dirichlet_BC::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Dirichlet_BC<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Dirichlet_BC::value(points[i]);
  }



  /**
   * The template class Neumann_BC computes the values of the non-homogeneous BC in each point.
   */
  template <int dim>
  class Neumann_BC : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Neumann_BC () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Neumann_BC::value extracts the value in one point ("p"). To use a spatially variable BC
   * modify this method accordingly.
   */
  template <int dim>
  double Neumann_BC<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
    return this->boundary_stress;
  }

  /** Neumann_BC::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Neumann_BC<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Neumann_BC::value(points[i]);
  }



  /**
   * The main class of the code, contains all the functions that do the job of setting up and
   * solving the model.
   *
   * \note Since the model parameters have been set in FunctionBase, this new class
   * will also inherit from it.
   */
  template <int dim>
  class ApShear : protected FunctionBase<dim>
  {
  public:
    enum RefinementMode
    {
      global_refinement, adaptive_refinement
    };

    enum Model
    {
      turcotte, savage
    };

    ApShear (const FiniteElement<dim>  &fe, /**< The code is written so different elements can be used and are defined as input parameters.*/
             const FiniteElement<dim>  &history_fe, /**< FE to store history data*/
             const unsigned int        degree, /**< To be consistent, we define the degree of the finite elements as an input parameter. */
             const RefinementMode      refinement_mode,/**< The mesh can be refined either globally (global_refinement) or adaptively (adaptive_refinement), which one to use is defined as an input parameters */
             const Model               model /**< The code can be used to simmulate the deformation caused by a fault using two models.*/
                          );

    ~ApShear ();

    void run ();

  private:
    void do_first_time_step (); /**< Create the initial mesh and perform preliminary refinements and solve for the initial time step*/
    void do_time_step (); /**< Solve for successive time steps (includes mesh refinement every several time steps).*/

    void make_grid_turcotte (); /**<  Create initial mesh by refining globally several times. Prepare boundaries for Turcotte and Spence Model */
    void make_grid_savage (); /**<  Create initial mesh by refining globally several times. Prepare boundaries for Savage and Burford Model */
    void refine_grid (const unsigned int min_grid_level,
                      const unsigned int max_grid_level); /**<  Refines the mesh adaptively. */
    void qpoints_to_DG (); /**< The information stored in the quadrature points is moved to a Discontinuous Galerkin space*/
    void DG_to_qpoints (); /**< The information stored in a Discontinious Galerking space is moved to the quadrature points*/

    void setup_system (); /**<  Setup DoFs, renumbering and constraints. */
    void setup_quadrature_point_history (); /**< Initialize the data stored at quadrature points. */
    void update_quadrature_point_history (); /**< Update the data stored at quadrature points and sore it in a DG space. */
    void assemble_system (); /**<  Assemble stiffness matrix and RHS. */

    void solve (); /**< The solver an refinement steps are defined here.*/
    void get_time_step (); /**< Obtain the the time step increase*/

    void compare_solutions_turcotte (); /**< The numerical solution is compared against Turcotte and Spence Model*/
    void compare_solutions_savage (); /**< The numerical solution is compared against Savage and Burford Model*/
    void output_results (); /* Results are output in vtk format.*/

    const unsigned int                            degree;

    Triangulation<dim>                            triangulation;
    DoFHandler<dim>                               dof_handler,
                                                  history_dof_handler;
    const QGauss<dim>                             quadrature_formula;

    SmartPointer<const FiniteElement<dim> >       fe,
                                                  history_fe;

    ConstraintMatrix                              constraints;

    SparsityPattern                               sparsity_pattern;
    SparseMatrix<double>                          system_matrix;
    
    Vector<double>                                system_rhs,
                                                  solution;

    std::vector<PointHistory<dim> >               quadrature_point_history;
    std::vector<Vector<double> >                     history_field ;

    const RefinementMode                          refinement_mode;

    const Model                                   model;

    const unsigned int                            n_initial_global_refinement,
                                                  n_pre_refinment;
    unsigned int                                  step;

    double                                        time_step,
                                                  total_time;

    std::vector<double>                           surface_q_points,
                                                  surface_analytical_solution,
                                                  surface_numerical_solution;
    double                                        error;
  };


  template <int dim>
  ApShear<dim>::ApShear (const FiniteElement<dim> &fe,
                         const FiniteElement<dim> &history_fe,
                         const unsigned int       degree,
                         const RefinementMode     refinement_mode,
                         const                    Model model) :
    degree (degree),
    dof_handler (triangulation),
    history_dof_handler (triangulation),
    quadrature_formula(degree+1),
    fe (&fe),
    history_fe (&history_fe),
    history_field (2, Vector<double> ()),
    refinement_mode (refinement_mode),
    model (model),
    n_initial_global_refinement (2),
    n_pre_refinment ((refinement_mode==global_refinement)?2:5),
    step (0),
    time_step (1e-6),
    total_time (0),
    error (0)
  {}


  template <int dim>
  ApShear<dim>::~ApShear ()
  {
    dof_handler.clear ();
  }


  /**
   * The problem is solved in several steps:
   *   - Initial time step (do_first_time_step()):
   *     -# Create an initial uniform grid and set boundary conditions
   *        (make_grid_turcotte() and make_grid_savage()).
   *     -# Initialize the history values stored at quadrature points with zeros
   *        (setup_quadrature_point_history()).
   *     -# Setup DoFs and constraints (setup_system()).
   *     -# Assemble the stiffness matrix and RHS (assemble_system()).
   *     -# Solve (solve()).
   *     -# Refine the mesh adaptively (refine_grid()) and repeat 2-5.
   *
   *   - Successive time steps (do_time_step()):
   *     -# Calculate the time step and compute the total displacement (get_time_step_andisplacement()).
   *     -# Update the values stored in the quadrature points using the new solution
   *        (update_quadrature_point_history()).
   *     -# Assemble the stiffness matrix and RHS (setup_system()).
   *     -# Solve (solve()).
   *     -# Every several time steps refine the grid.
   *       + a) Transfer the quadrature point history to a Discontinuous Galerkin field
   *            (qpoints_to_DG()).
   *       + b) Refine the mesh and interpolate the old solution and viscous strain in the new mesh
   *            (refine_grid()).
   *       + c) Transfer the data stored in the DG field to the new quadrature points
   *          (DG_to_qpoints()).
   */
  template <int dim>
  void ApShear<dim>::run ()
  {
    do_first_time_step ();

    for (step=1; step<1; ++step)
      {
        do_time_step ();
     }
    output_results ();
  }


  /**
   * In the first time step we will create a uniform grid and set the boundary conditions depending
   * on the model that we want to solve. Then we will solve the problem and refine the mesh adaptively
   * several times
   */
  template <int dim>
  void ApShear<dim>::do_first_time_step ()
  {
    std::cout << "Step 0, t = 0" << std::endl;

    std::cout << "Refining initial mesh...";
    for (unsigned int refinement_step=0; refinement_step<=n_pre_refinment; ++refinement_step)
      {
        if (refinement_step==0)
          switch (model)
          {
          case turcotte:
            make_grid_turcotte ();
            break;
          case savage:
            make_grid_savage ();
            break;
          default:
            Assert (false, ExcNotImplemented());
          }
        else
          refine_grid (n_initial_global_refinement,
                       n_initial_global_refinement + n_pre_refinment);

        assemble_system ();
        solve ();
      }
    std::cout << triangulation.n_active_cells() << " cells." << std::endl;
    output_results ();
  }


  /**
   * In successive time steps we will update the viscous strain (using the
   * solution and the viscous strain from the last time step) and we will solve
   * again. The mesh will be refined adaptively every several time steps.
   *
   */
  template <int dim>
  void ApShear<dim>::do_time_step ()
  {
    get_time_step ();
    total_time += time_step;
    std::cout << "Step " << step << ", t = " << total_time << std::endl;
    update_quadrature_point_history ();
    assemble_system ();
    solve ();
    output_results ();

    if (step % 5 == 0)
      {
        std::cout << "Refining grid...";
        qpoints_to_DG ();
        refine_grid (n_initial_global_refinement,
                             n_initial_global_refinement + n_pre_refinment);
        DG_to_qpoints ();
        std::cout << triangulation.n_active_cells() << " cells." << std::endl;
      }
 }


  /**
   * The initial mesh will consist of regular rectangular elements. This mesh will be refined
   * in successive steps using refine_grid(). First a mesh is created with
   * elements of the highest possible quality (as close to square as possible). Then this mesh
   * is refined globally a given number of times.
   *
   * \note If the fault is very small it may be ignored (if it's size is smaller than the
   * element size). If this happens, include more global refinement steps.
   */
  template <int dim>
  void ApShear<dim>::make_grid_turcotte ()
  {
    /* The domain can be rectangle, but we want the cells go be as square as possible.*/
    Point<dim> corner_low_left(0.0,-this->height);
    Point<dim> corner_up_right(this->width,0.0);

    std::vector<unsigned int> repetitions(dim,1);
    repetitions[0] = std::max(1.0,round(this->width/this->height));
    GridGenerator::subdivided_hyper_rectangle (triangulation,
                                               repetitions,
                                               corner_low_left,
                                               corner_up_right);

    triangulation.refine_global (n_initial_global_refinement);

    /**
     * Once the mesh has been created and refined, we assign an indicator to each part of the
     * boundary, depending on what kind of boundary it is. There are four possible boundary
     * conditions:
     *
     * - 0) The locked part above the fault, with imposed zero displacement (homogeneous
     * Dirichlet boundary conditions).
     *
     * - 1) Regions of imposed displacement (non-homogeneous Dirichlet
     * boundary conditions).
     *
     * - 2) Free boundaries, with zero imposed stress (homogeneous Neumann boundary
     * conditions).
     *
     * - 3) Imposed stress (non-homogeneous boundary conditions.
     */
    typename Triangulation<dim>::cell_iterator cell = triangulation.begin (),
                                 endc = triangulation.end();
      for (; cell!=endc; ++cell)
        for (unsigned int face_number=0;
             face_number<GeometryInfo<dim>::faces_per_cell;
             ++face_number)
          /** For Turcotte and Spence's model \cite turcotte_spence_74, the
           * locked part (0) will be on the left boundary, from the top to a certain
           * depth \f$d\f$. From that depth until the bottom of the model the
           * boundary will be free (2). The top and bottom boundaries will also
           * be free (2). Finally, the leftmost boundary will have an imposed
           * stress \f$S\f$ (3).
           */
          if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Locked part
              &&
              (cell->face(face_number)->center()(1) > -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (0);
          else if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Fault
              &&
              (cell->face(face_number)->center()(1) <= -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (2);
          else if (std::fabs(cell->face(face_number)->center()(1) - (0)) < 1e-12) // Top
            cell->face(face_number)->set_boundary_indicator (2);
          else if  (std::fabs(cell->face(face_number)->center()(1) - (-this->height)) < 1e-12) // Bottom
            cell->face(face_number)->set_boundary_indicator (2);
          else if (std::fabs(cell->face(face_number)->center()(0) - (this->width)) < 1e-12) // Left
            cell->face(face_number)->set_boundary_indicator (3);

      setup_system ();
      setup_quadrature_point_history ();
  }


  /**
   * The initial mesh will consist of regular rectangular elements. This mesh will be refined
   * in successive steps using refine_grid(). First a mesh is created with
   * elements of the highest possible quality (as close to square as possible). Then this mesh
   * is refined globally a given number of times.
   *
   * \note If the fault is very small it can be ignored (if it's size is smaller than the
   * element size). If this happens, include more global refinement steps.
   */
  template <int dim>
  void ApShear<dim>::make_grid_savage ()
  {
    /* The domain can be rectangle, but we want the cells go be as square as possible. */
    Point<dim> corner_low_left(0.0,-this->height);
    Point<dim> corner_up_right(this->width,0.0);

    std::vector<unsigned int> repetitions(dim,1);
    repetitions[0] = std::max(1.0,round(this->width/this->height));
    GridGenerator::subdivided_hyper_rectangle (triangulation,
                                               repetitions,
                                               corner_low_left,
                                               corner_up_right);

    triangulation.refine_global (n_initial_global_refinement);

    /**
     * Once the mesh has been created and refined, we assign an indicator to each part of the
     * boundary, depending on what kind of boundary it is. There are four possible boundary
     * conditions:
     *
     * - 0) The locked part above the fault, with imposed zero displacement (homogeneous
     * Dirichlet boundary conditions).
     *
     * - 1) Regions of imposed displacement (non-homogeneous Dirichlet
     * boundary conditions).
     *
     * - 2) Free boundaries, with zero imposed stress (homogeneous Neumann boundary
     * conditions).
     *
     * - 3) Imposed stress (non-homogeneous boundary conditions.
     */
    typename Triangulation<dim>::cell_iterator cell = triangulation.begin (),
                                 endc = triangulation.end();
      for (; cell!=endc; ++cell)
        for (unsigned int face_number=0;
             face_number<GeometryInfo<dim>::faces_per_cell;
             ++face_number)

          /** For Savage and Burford model \cite savage_burford_73, the
           * locked part (0) will be on the left boundary, from the top to a certain
           * depth \f$d\f$. From that depth until the bottom of the model we impose
           * a displacement \f$UD\f$ (1). The bottom boundary will have the same imposed
           * displacement (1). Finally, the top and left boundaries will be free (2).
           */
          if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Locked part
              &&
              (cell->face(face_number)->center()(1) > -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (0);
          else if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Fault
              &&
              (cell->face(face_number)->center()(1) <= -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (1);
          else if (std::fabs(cell->face(face_number)->center()(1) - (0)) < 1e-12) // Top
            cell->face(face_number)->set_boundary_indicator (2);
          else if  (std::fabs(cell->face(face_number)->center()(1) - (-this->height)) < 1e-12) // Bottom
            cell->face(face_number)->set_boundary_indicator (1);
          else if (std::fabs(cell->face(face_number)->center()(0) - (this->width)) < 1e-12) // Left
            cell->face(face_number)->set_boundary_indicator (2);

      setup_system ();
      setup_quadrature_point_history ();
  }


  template <int dim>
  void ApShear<dim>::refine_grid (const unsigned int min_grid_level,
                                  const unsigned int max_grid_level)
  {
    switch (refinement_mode)
      {
      case global_refinement:
      {
        triangulation.refine_global (1);
        break;
      }

      /**
       * The first step is to identify which cells need to be refined, for which we use we use a
       * Kelly Error Estimator. We will refine 60% of the cells that have the highest error and coarsen
       * the rest).
       */
      case adaptive_refinement:
      {
        Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

        KellyErrorEstimator<dim>::estimate (dof_handler,
                                            QGauss<dim-1>(3),
                                            typename FunctionMap<dim>::type(),
                                            solution,
                                            estimated_error_per_cell);

        GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                         estimated_error_per_cell,
                                                         0.6, 0.4);
        /**
         * We do not want to make any cell too small, because that would mean (according to
         * the Courant–Friedrichs–Lewy condition) that we have to make the time steps very small.
         * Therefore, we will un-flag those cells that have already been refined a certain number of
         * times so they are not refined any further. A similar thing happens with cells that
         * are too big, therefore, we will also remove the flag of the cells that have reached a
         * certain maximum size so they are not coarsened any more.
         */

        if (triangulation.n_levels() > max_grid_level)
          for (typename Triangulation<dim>::active_cell_iterator
               cell = triangulation.begin_active(max_grid_level);
               cell != triangulation.end(); ++cell)
            cell->clear_refine_flag ();
        for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active(min_grid_level);
             cell != triangulation.end_active(min_grid_level); ++cell)
          cell->clear_coarsen_flag ();

        /**
         * Here we also need to transfer the #solution and the #history_field
         * (where the old viscous strain is stored) to the new mesh. First we have to prepare the
         * vectors that will be transferred to the new grid (the old grid will disappear
         * once we have done the refinement so we have to transfer the data at the same
         * time as the refinement).
         */

        SolutionTransfer<dim> solution_trans(dof_handler);
        Vector<double> unrefined_solution (dof_handler.n_dofs());
        unrefined_solution = solution;

        SolutionTransfer<dim> history_trans(history_dof_handler);
        std::vector<Vector<double> > unrefined_history_field (2, Vector<double>(history_dof_handler.n_dofs()));
        for (unsigned int i=0; i<2; ++i)
          unrefined_history_field [i] = history_field[i];

        triangulation.prepare_coarsening_and_refinement();
        solution_trans.prepare_for_coarsening_and_refinement(unrefined_solution);
        history_trans.prepare_for_coarsening_and_refinement(unrefined_history_field);

        /**
         * Once everything is ready we can actually refine and coarsen the mesh.
         */
        triangulation.execute_coarsening_and_refinement ();

        /**
         * Once the new mesh is ready, we need to resent the system and the data at the quadrature
         * points.
         */
        setup_system ();
        setup_quadrature_point_history ();

        /**
         * Finally, we can interpolate the #solution and the
         * #history_field in the new mesh.
         */
        Vector<double> refined_solution (dof_handler.n_dofs());
        solution_trans.interpolate(unrefined_solution, refined_solution);
        solution = refined_solution;

        std::vector<Vector<double> > refined_history_field (2, Vector<double>(history_dof_handler.n_dofs()));
        history_trans.interpolate(unrefined_history_field, refined_history_field);
        for (unsigned int i=0; i<2; ++i)
          history_field[i] = refined_history_field[i];

        break;
      }

      default:
      {
        Assert (false, ExcNotImplemented());
      }
      }
  }


  /**
   * As part of the refinement process, we need to store the old viscous strain in a FE field defined
   * everywhere. Once that has been done, we will be able to transfer the old viscous strain (or any
   * data stored in the quadrature point using PointHistory) from the old mesh to the refined
   * mesh. For more details see
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_18.html#Refinementduringtimesteps">Step-18</a>
   */
  template <int dim>
  void ApShear<dim>::qpoints_to_DG ()
  {
    std::vector<Vector<double> >
                 local_history_values_at_qpoints (2, Vector<double>(quadrature_formula.size())),
                 local_history_fe_values (2, Vector<double>(history_fe->dofs_per_cell));
    FullMatrix<double> qpoint_to_dof_matrix (history_fe->dofs_per_cell,
                                             quadrature_formula.size());
    FETools::compute_projection_from_quadrature_points_matrix
              (*history_fe,
               quadrature_formula, quadrature_formula,
               qpoint_to_dof_matrix);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();
    for (; cell!=endc; ++cell, ++dg_cell)
      {
        PointHistory<dim> *local_quadrature_points_history
               = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
        Assert (local_quadrature_points_history >=
                    &quadrature_point_history.front(),
                    ExcInternalError());
        Assert (local_quadrature_points_history <
                    &quadrature_point_history.back(),
                    ExcInternalError());
        for (unsigned int i=0; i<2; i++)
          {
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              local_history_values_at_qpoints[i](q)
                       = local_quadrature_points_history[q].new_v_strain[i];
            qpoint_to_dof_matrix.vmult (local_history_fe_values[i],
                                        local_history_values_at_qpoints[i]);
            dg_cell->set_dof_values (local_history_fe_values[i],
                                     history_field[i]);
          }
      }
  }


  /**
   * After the mesh has been refined, we need to transfer the viscous strain back to the
   * quadrature points. For more details see
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_18.html#Refinementduringtimesteps">Step-18</a>
   */
  template <int dim>
  void ApShear<dim>::DG_to_qpoints ()
  {
    std::vector<Vector<double> >
                   local_history_values_at_qpoints (2, Vector<double>(quadrature_formula.size())),
                   local_history_fe_values (2, Vector<double>(history_fe->dofs_per_cell));

    FullMatrix<double> dof_to_qpoint_matrix (quadrature_formula.size(),
                                             history_fe->dofs_per_cell);
    FETools::compute_interpolation_to_quadrature_points_matrix
              (*history_fe,
               quadrature_formula,
               dof_to_qpoint_matrix);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();
    for (; cell != endc; ++cell, ++dg_cell)
    {
     PointHistory<dim> *local_quadrature_points_history
             = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
      Assert (local_quadrature_points_history >=
                  &quadrature_point_history.front(),
                  ExcInternalError());
      Assert (local_quadrature_points_history <
                  &quadrature_point_history.back(),
                  ExcInternalError());
      for (unsigned int i=0; i<2; i++)
        {
          dg_cell->get_dof_values (history_field[i],
                                   local_history_fe_values[i]);
          dof_to_qpoint_matrix.vmult (local_history_values_at_qpoints[i],
                                      local_history_fe_values[i]);
          for (unsigned int q=0; q<quadrature_formula.size(); ++q)
            local_quadrature_points_history[q].new_v_strain[i]
                       = local_history_values_at_qpoints[i](q);
        }
    }
  }


  /**
   * Here we set up DoFs, renumber them for efficiency and set required constraints.
   */
  template <int dim>
  void ApShear<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (*fe);
    history_dof_handler.distribute_dofs (*history_fe);

    /**
     * Renumber DoFs to increase efficiency of the solver. In this case, it is not really necessary
     * for the solver and preconditioners that are being used
     */
    DoFRenumbering::Cuthill_McKee (dof_handler);

    /**
     * Set constraints for hanging nodes (due to adaptive refinement) and for homogeneous and
     * non-homogeneous BC (indicators 0 and 1 respectively).
     */
    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              Dirichlet_BC<dim>(),
                                              constraints);
    constraints.close ();


    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    c_sparsity,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);

    sparsity_pattern.copy_from(c_sparsity);

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    for (unsigned int i=0; i<2; i++)
      history_field[i].reinit(history_dof_handler.n_dofs());

    system_rhs.reinit (dof_handler.n_dofs());
  }


  template <int dim>
  void ApShear<dim>::setup_quadrature_point_history ()
  {

    triangulation.clear_user_data();

    /* Only necesary if implementing parallel code (the number of quadrature points
     * might change between processors
    {
      std::vector<PointHistory<dim> > tmp;
      tmp.swap (quadrature_point_history);
    }
    */
    quadrature_point_history.resize (triangulation.n_active_cells() *
                                     quadrature_formula.size());
    unsigned int history_index = 0;
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell->set_user_pointer (&quadrature_point_history[history_index]);
        history_index += quadrature_formula.size();
      }
    Assert (history_index == quadrature_point_history.size(),
            ExcInternalError());
  }


  /*
   * The data from previous time steps, which is stored in the quadrature points, is updated each
   * time step. In this case, the only data stored is the old viscous strain.
   *
   * \note
   * Because we are using the anti-plane shear approximation, the only non-zero components of
   * the strain are \f$ \varepsilon_{xz} = \varepsilon_{zx} \f$ and
   * \f$ \varepsilon_{yz} = \varepsilon_{zy} \f$. Therefore, we only need to store two variables,
   * which will be stored in a vector with two elements. In general, we will need the 9 elements
   * of the tensor, which will be stored in a 3 by 3 matrix (see VE_fault_viscous).
   */
  template <int dim>
  void ApShear<dim>::update_quadrature_point_history ()
  {
    FEValues<dim>                 fe_values (*fe, quadrature_formula,
                                             update_quadrature_points | update_values |
                                             update_gradients);

    const Shear_Modulus<dim>      shear_modulus;
    std::vector<double>           sm_values (quadrature_formula.size());
    const Viscosity<dim>          viscosity;
    std::vector<double>           visc_values (quadrature_formula.size());

    std::vector<Tensor<1,dim> >   new_solution_gradients (quadrature_formula.size());

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      PointHistory<dim> *local_quadrature_points_history
        = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());

      fe_values.reinit (cell);

      shear_modulus.value_list(fe_values.get_quadrature_points(),sm_values);
      viscosity.value_list(fe_values.get_quadrature_points(),visc_values);

      fe_values.get_function_gradients (solution,
                                        new_solution_gradients);

      for (unsigned int q=0; q<quadrature_formula.size(); ++q)
        {
          const double A = time_step * sm_values[q] / visc_values[q];

          Tensor<1,2> new_v_strain (2);

          Tensor<1,2> old_v_strain = local_quadrature_points_history[q].new_v_strain;

          //The values of the strain are only valid under the anti-plane shear approximation
          Tensor<1,2> old_strain (2);
          old_strain[0] = 0.5 * new_solution_gradients[q][0];
          old_strain[1] = 0.5 * new_solution_gradients[q][1];

          for (unsigned int i=0; i < 2; ++i)
            {
              new_v_strain[i] = old_v_strain[i] + A * old_strain [i];
            }
          local_quadrature_points_history[q].new_v_strain = new_v_strain;
        }
    }
  }


  /**
   * The stiffness matrix and RHS are assembled here.
   */
  template <int dim>
  void ApShear<dim>::assemble_system ()
  {
    system_rhs = 0;
    system_matrix = 0;

    QGauss<dim-1> face_quadrature_formula(degree+1);

    const unsigned int n_q_points    = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    FEValues<dim>  fe_values (*fe, quadrature_formula,
                              update_values   | update_gradients |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula,
                                      update_values         | update_quadrature_points  |
                                      update_normal_vectors | update_JxW_values);

    const Shear_Modulus<dim>     shear_modulus;
    std::vector<double>          sm_values (n_q_points);
    const Viscosity<dim>         viscosity;
    std::vector<double>          visc_values (n_q_points);
    const BodyForce<dim>         body_force;
    std::vector<double>          bf_values (n_q_points);
    const Neumann_BC<dim>        neumann_bc;
    std::vector<double>          nbc_values (n_face_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    /**
     * To assemble the stiffness matrix and RHS for the given problem we assemble
     * the local matrices for each element and the transfer it to the global
     * matrix.
     */
    unsigned int cell_n =0;
    for (; cell!=endc; ++cell)
      {
        ++cell_n;
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);

        PointHistory<dim> *local_quadrature_points_history
                          = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());

        shear_modulus.value_list (fe_values.get_quadrature_points(),sm_values);
        viscosity.value_list (fe_values.get_quadrature_points(),visc_values);
        body_force.value_list (fe_values.get_quadrature_points(),bf_values);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            /*
             * To compute the contribution of the internal elastic force in the right hand
             * side of the weak form, we need the viscous strain
             * from the previous time step (stored at each quadrature point).
             */

            Tensor<1,2> new_v_strain = local_quadrature_points_history[q].new_v_strain;

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += (sm_values [q] *
                                       fe_values.shape_grad(i,q) *
                                       fe_values.shape_grad(j,q)*
                                       fe_values.JxW(q));

                cell_rhs(i) += (fe_values.shape_grad(i,q) *
                                2.0 * sm_values[q] * new_v_strain *
                                fe_values.JxW(q)) -
                               (fe_values.shape_value(i,q) *
                                bf_values [q] *
                                fe_values.JxW(q));
              }
          }

        for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
          /**
           * Apply homogeneous Neumann boundary conditions. We loop over all the faces in
           * each cell. We check if that face is in a boundary and if it is a boundary of type 3
           * (non-homogeneous Neumann boundary conditions), to those faces we add the corresponing
           * contribution.
           */
          if (cell->face(face_number)->at_boundary()
              &&
              (cell->face(face_number)->boundary_indicator() == 3))
            {
              fe_face_values.reinit (cell, face_number);
              neumann_bc.value_list (fe_face_values.get_quadrature_points(), nbc_values);

              for (unsigned int q=0; q<n_face_q_points; ++q)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    cell_rhs(i) += (fe_face_values.shape_value(i,q) *
                                    nbc_values[q] *
                                    fe_face_values.JxW(q));
            }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (cell_matrix,
                                                cell_rhs,
                                                local_dof_indices,
                                                system_matrix,
                                                system_rhs);
      }
  }


  /**
   * The strategy to invert the stiffness matrix and solve the system is set here.
   * Because the problem is simple, we just use a Conjugate Gradient solver with the identity as
   * a preconditioner.
   */
  template <int dim>
  void ApShear<dim>::solve ()
  {
    SolverControl           solver_control (2000, 1e-12);
    SolverCG<>              solver (solver_control);

    solver.solve (system_matrix, solution, system_rhs,
                  PreconditionIdentity());

    constraints.distribute (solution);
  }


  /**
   * We need to calculate the size of each time step. We choose a size such that, the
   * algorithm is stable \cite cormeau_75 \cite yamasaki_houseman_12:
   *
   * \f[ \Delta t = \frac{\eta}{24 \mu} \f]
   *
   * We use the minimum value of \f$ \Delta t \f$ in the whole domain. Therefore,
   * we need to compute the time increment at every quadrature point.
   *
   * \note
   * In the current implementation, both the shear modulus and the viscosity are
   * uniform, and so will be the time step increment. However, we will still calculate
   * the time step increment at every point, so it is easier to implement non-uniform
   * parameters in the future.
   */
  template <int dim>
  void ApShear<dim>::get_time_step ()
  {
    const unsigned int   n_q_points
      = quadrature_formula.size();
    FEValues<dim> fe_values (*fe, quadrature_formula,
                             update_quadrature_points);


    time_step = 1e24; // We start with a very high value of the time step

    const Shear_Modulus<dim>     shear_modulus;
    std::vector<double>          sm_values (n_q_points);
    const Viscosity<dim>         viscosity;
    std::vector<double>          visc_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        shear_modulus.value_list (fe_values.get_quadrature_points(),sm_values);
        viscosity.value_list (fe_values.get_quadrature_points(),visc_values);

        for (unsigned int q=0; q<n_q_points; ++q)
          time_step = std::min (time_step,
                                visc_values[q] / sm_values[q] / 24.0);
      }
  }


  /**
   * This function will take care of computing and storing the difference between the numerical
   * solution and the analytical solution from Turcotte and Spence Model \cite turcotte_spence_74.
   * At the end we will output the results in several graphs. The comparison will stored in three
   * different variables:
   * - analytical_solution: Analytical solution computed at each quadrature point of the final
   * mesh.
   * - surface_numerical_solution: Numerical solution in each quadrature point at each refinement cycle.
   * - integrated_diff: Integrated difference between both solutions using the L1 norm.
   */
  template <int dim>
  void ApShear<dim>::compare_solutions_turcotte ()
  {

    QGauss<dim-1>   face_quadrature_formula(degree+1);
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
    const unsigned int n_face_q_points    = face_quadrature_formula.size();

    Solution_Turcotte<dim> exact_solution;

    FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula,
                                         update_values         | update_quadrature_points  |
                                         update_normal_vectors | update_JxW_values);

    exact_solution.get_time (total_time);

    surface_q_points.clear ();
    surface_analytical_solution.clear ();
    surface_numerical_solution.clear ();
    /**
     * First we need to run a loop over every cell and every quadrature point to determine if
     * they are located at the surface, if so, we will compute the analytical solution and extract
     * the numerical solution to store them in analytical_solution and numerical_solution.
     * Once we have both solutions at the surface, we will compute the integrated
     * difference between them and we will store it in integrated_diff.
     */
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int face_number=0; face_number<n_faces; ++face_number)
        if (std::fabs(cell->face(face_number)->center()(1)) < 1e-12)
          {
            fe_face_values.reinit (cell, face_number);
            for (unsigned int q=0; q<n_face_q_points; ++q)
              {
                surface_q_points.push_back(fe_face_values.quadrature_point(q)[0]);
                surface_analytical_solution.push_back(exact_solution.value(fe_face_values.quadrature_point(q)));
                surface_numerical_solution.push_back(VectorTools::point_value (dof_handler, solution,
                                                                    fe_face_values.quadrature_point(q)));
                error += (std::fabs(surface_analytical_solution.back()-surface_numerical_solution.back())*
                                                    fe_face_values.JxW(q));

              }
          }
  }


  /**
   * This function will take care of computing and storing the difference between the numerical
   * solution and the analytical solution from Turcotte and Spence Model \cite turcotte_spence_74.
   * At the end we will output the results in several graphs. The comparison will stored in three
   * different variables:
   * - analytical_solution: Analytical solution computed at each quadrature point of the final
   * mesh.
   * - surface_numerical_solution: Numerical solution in each quadrature point at each refinement cycle.
   * - integrated_diff: Integrated difference between both solutions using the L1 norm.
   */
  template <int dim>
  void ApShear<dim>::compare_solutions_savage ()
  {

    QGauss<dim-1>   face_quadrature_formula(degree+1);
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
    const unsigned int n_face_q_points    = face_quadrature_formula.size();

    const Solution_Savage<dim> exact_solution;

    FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula,
                                         update_values         | update_quadrature_points  |
                                         update_normal_vectors | update_JxW_values);
    /**
     * First we need to run a loop over every cell and every quadrature point to determine if
     * they are located at the surface, if so, we will compute the analytical solution and extract
     * the numerical solution to store them in analytical_solution and numerical_solution.
     * Once we have both solutions at the surface, we will compute the integrated
     * difference between them and we will store it in integrated_diff.
     */
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int face_number=0; face_number<n_faces; ++face_number)
        if (std::fabs(cell->face(face_number)->center()(1)) < 1e-12)
          {
            fe_face_values.reinit (cell, face_number);
            for (unsigned int q=0; q<n_face_q_points; ++q)
              {
                surface_q_points.push_back(fe_face_values.quadrature_point(q)[0]);
                surface_analytical_solution.push_back(exact_solution.value(fe_face_values.quadrature_point(q)));
                surface_numerical_solution.push_back(VectorTools::point_value (dof_handler, solution,
                                                                    fe_face_values.quadrature_point(q)));
                error += (std::fabs(surface_analytical_solution.back()-surface_numerical_solution.back())*
                                                    fe_face_values.JxW(q));

              }
          }
  }


  /**
   * The results will be output in this step. Three different outputs are generated.
   * - The final mesh in eps format.
   * - The final displacement in the 2D domain in vtk format
   * - Information regarding the conparison between the numerical and the analytical solutions.
   */
  template <int dim>
  void ApShear<dim>::output_results ()
  {
    std::ostringstream filename;
    // The output file name depends on the refinement mode and the element type.
    switch (model)
      {
      case turcotte:
        filename << "solution-turcotte";
        break;
      case savage:
        filename << "solution-savage";
        break;
      default:
        Assert (false, ExcNotImplemented());
      }

    switch (refinement_mode)
      {
      case global_refinement:
        filename << "-global";
        break;
      case adaptive_refinement:
        filename << "-adaptive";
        break;
      default:
        Assert (false, ExcNotImplemented());
      }

    switch (fe->degree)
      {
      case 1:
        filename << "-q1";
        break;
      case 2:
        filename << "-q2";
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    filename << "(step" << step << ")";

//    // Output final mesh
//    std::ostringstream eps_filename;
//    eps_filename << filename.str() << ".eps";
//    std::ofstream output_mesh (eps_filename.str().c_str());
//
//    GridOut grid_out;
//    grid_out.write_eps (triangulation, output_mesh);

    // Output numerical solution
    std::ostringstream vtk_filename;
    vtk_filename << filename.str() << ".vtk";
    std::ofstream output_num_sol (vtk_filename.str().c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "displacement");

    data_out.build_patches (fe->degree);
    data_out.write_vtk (output_num_sol);

    // Output point history
    std::ostringstream vtk_filename_history;
    vtk_filename_history << filename.str() << "_history.vtk";
    std::ofstream output_history (vtk_filename_history.str().c_str());

    DataOut<dim> data_out_history;
    data_out_history.attach_dof_handler (history_dof_handler);
    data_out_history.add_data_vector (history_field[0], "history_field_xz");
    data_out_history.add_data_vector (history_field[1], "history_field_yz");

    data_out_history.build_patches (history_fe->degree);
    data_out_history.write_vtk (output_history);


    // Output comparison
    compare_solutions_turcotte ();

    std::ostringstream txt_filename_error;
    txt_filename_error << filename.str() << "-error" << ".txt";
    std::ofstream out_compare_error(txt_filename_error.str().c_str());

    std::ostringstream txt_filename_analytical;
    txt_filename_analytical << filename.str() << "-compare" << ".txt";
    std::ofstream out_compare (txt_filename_analytical.str().c_str());

    const int n_q_points = surface_q_points.size ();
    for (unsigned int q = 0; q<n_q_points; ++q)
      {
        out_compare << surface_q_points [q] << "\t";
        out_compare << surface_analytical_solution[q] << "\t";
        out_compare << surface_numerical_solution[q] << std::endl;
      }
    out_compare_error << error << std::endl;
  }
}





int main ()
{
  const unsigned int dim = 2;
  const unsigned int degree = 1;

  try
    {
      using namespace dealii;
      using namespace vsf;

      deallog.depth_console (0);
      /**
       * \note The equations are only valid for 2D geometries (anti-plane shear approximation). If the
       * user tries to create a 3D geometry, an exception is thrown and the program stops.
       */
      Assert (dim == 2, ExcNotImplemented ());

      /**
       * When an ApShear is created the user must specify the kind of elements and
       * refinement (adaptive_refinement or global_refinement) and Model (turcotte or savage)
       * are going oto be used.
       */
      {
        std::cout << "Solving with Q1 elements, adaptive refinement for Turcotte and Spence Model" << std::endl
                  << "==========================================================================" << std::endl
                  << std::endl;

        FE_Q<dim> fe(degree);
        FE_DGQ<dim> history_fe(degree);
        ApShear<dim>
        antiplane_shear_problem_2d (fe,
                                    history_fe,
                                    degree,
                                    ApShear<dim>::adaptive_refinement,
                                    ApShear<dim>::turcotte);

        antiplane_shear_problem_2d.run ();

        std::cout << std::endl;
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
