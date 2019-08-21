/*!
 * @file diffusion_problem.hpp
 * @brief Contains implementation of the main object.
 * @author Konrad Simon
 * @date August 2019
 */

#ifndef INCLUDE_DIFFUSION_PROBLEM_HPP_
#define INCLUDE_DIFFUSION_PROBLEM_HPP_

// Deal.ii
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

// STL
#include <cmath>
#include <fstream>
#include <iostream>

// My Headers
#include "matrix_coeff.hpp"
#include "right_hand_side.hpp"
#include "neumann_bc.hpp"
#include "dirichlet_bc.hpp"

/*!
 * @namespace DiffusionProblem
 * @brief Contains implementation of the main object
 * and all functions to solve a
 * Dirichlet-Neumann problem on a unit square.
 */
namespace DiffusionProblem
{
using namespace dealii;

/*!
 * @class DiffusionProblem
 * @brief Main class to solve
 * Dirichlet-Neumann problem on a unit square.
 */
template <int dim>
class DiffusionProblem
{
public:
	DiffusionProblem (unsigned int n_refine);
	void run ();

private:
	void make_grid ();
	void setup_system ();
	void assemble_system ();
	void solve_iterative ();
	void output_results () const;

	Triangulation<dim>   			triangulation;
	FE_Q<dim>            			fe;
	DoFHandler<dim>      			dof_handler;

	AffineConstraints<double> 		constraints;

	SparsityPattern      			sparsity_pattern;
	SparseMatrix<double> 			system_matrix;

	/*!
	 * Current solution. Needed forping.
	 */
	Vector<double>       			solution;

	/*!
	 * Contains all parts of the right-hand side needed to
	 * solve the linear system.
	 */
	Vector<double>       			system_rhs;

	/*!
	 * Number of initial grid refinements.
	 */
	unsigned int n_refine;
};


/*!
 * Default constructor.
 */
template <int dim>
DiffusionProblem<dim>::DiffusionProblem (unsigned int n_refine) :
  fe (1),
  dof_handler (triangulation),
  n_refine (n_refine)
{}


/*!
 * @brief Set up the grid with a certain number of refinements.
 *
 * Generate a triangulation of \f$[0,1]^{\rm{dim}}\f$ with edges/faces
 * numbered form \f$1,\dots,2\rm{dim}\f$.
 */
template <int dim>
void DiffusionProblem<dim>::make_grid ()
{
	GridGenerator::hyper_cube (triangulation, 0, 1, /* colorize */ true);

	triangulation.refine_global (n_refine);

	std::cout << "Number of active cells: "
			<< triangulation.n_active_cells()
			<< std::endl;
}


/*!
 * @brief Setup sparsity pattern and system matrix.
 *
 * Compute sparsity pattern and reserve memory for the sparse system matrix
 * and a number of right-hand side vectors. Also build a constraint object
 * to take care of Dirichlet boundary conditions.
 */
template <int dim>
void DiffusionProblem<dim>::setup_system ()
{
	dof_handler.distribute_dofs (fe);

	std::cout << "Number of active cells: " << triangulation.n_active_cells()
			<< std::endl
			<< "Number of degrees of freedom: " << dof_handler.n_dofs()
			<< std::endl
			<< std::endl;


	constraints.clear();
	DoFTools::make_hanging_node_constraints(dof_handler, constraints);

	/*
	 * Set up Dirichlet boundary conditions.
	 */
	const Coefficients::DirichletBC<dim> dirichlet_bc;
	for (unsigned int i = 0; i<dim; ++i)
	{
		VectorTools::interpolate_boundary_values(dof_handler,
													/*boundary id*/ 2*i, // only even boundary id
													dirichlet_bc,
													constraints);
	}

	constraints.close();

	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler,
									dsp,
									constraints,
									/*keep_constrained_dofs =*/ true); // forping this is essential to be true

	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit (sparsity_pattern);

	solution.reinit (dof_handler.n_dofs());
	system_rhs.reinit (dof_handler.n_dofs());
}


/*!
 * @brief Assemble the system matrix and the static right hand side.
 *
 * Assembly routine to build the time-independent (static)part.
 * Neumann boundary conditions will be put on edges/faces
 * with odd number. Constraints are not applied here yet.
 */
template <int dim>
void DiffusionProblem<dim>::assemble_system ()
{
	QGauss<dim>  quadrature_formula(fe.degree + 1);
	QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

	FEValues<dim> 	fe_values (fe, quadrature_formula,
								update_values    |  update_gradients |
								update_quadrature_points  |  update_JxW_values);

	FEFaceValues<dim> 	fe_face_values(fe,
										face_quadrature_formula,
										update_values | update_quadrature_points |
										update_normal_vectors |
										update_JxW_values);

	const unsigned int   	dofs_per_cell = fe.dofs_per_cell;
	const unsigned int   	n_q_points    = quadrature_formula.size();
	const unsigned int 	n_face_q_points = face_quadrature_formula.size();

	FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>       cell_rhs (dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

	/*
	 * Matrix coefficient and vector to store the values.
	 */
	const Coefficients::MatrixCoeff<dim> 		matrix_coeff;
	std::vector<Tensor<2,dim>> 	matrix_coeff_values(n_q_points);

	/*
	 * Right hand side and vector to store the values.
	 */
	const Coefficients::RightHandSide<dim> 	right_hand_side;
	std::vector<double>      	rhs_values(n_q_points);

	/*
	 * Neumann BCs and vector to store the values.
	 */
	const Coefficients::NeumannBC<dim> 	neumann_bc;
	std::vector<double>  	neumann_values(n_face_q_points);

	/*
	 * Integration over cells.
	 */
	for (const auto &cell: dof_handler.active_cell_iterators())
	{
		cell_matrix = 0;
		cell_rhs = 0;

		fe_values.reinit (cell);

		// Now actually fill with values.
		matrix_coeff.value_list(fe_values.get_quadrature_points (),
						  	  	  matrix_coeff_values);
		right_hand_side.value_list(fe_values.get_quadrature_points(),
									   rhs_values);

		for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
		{
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{

					cell_matrix(i,j) += fe_values.shape_grad(i,q_index) *
										 matrix_coeff_values[q_index] *
										 fe_values.shape_grad(j,q_index) *
										 fe_values.JxW(q_index);
				} // end ++j

				cell_rhs(i) += fe_values.shape_value(i,q_index) *
								   rhs_values[q_index] *
								   fe_values.JxW(q_index);
			} // end ++i
		} // end ++q_index

		/*
		 * Boundary integral for Neumann values for odd boundary_id.
		 */
		for (unsigned int face_number = 0;
			 face_number < GeometryInfo<dim>::faces_per_cell;
			 ++face_number)
		{
			if (cell->face(face_number)->at_boundary() &&
					(
						(cell->face(face_number)->boundary_id() == 1) ||
						(cell->face(face_number)->boundary_id() == 3) ||
						(cell->face(face_number)->boundary_id() == 5)
					)
				)
			{
				fe_face_values.reinit(cell, face_number);

				/*
				 * Fill in values at this particular face.
				 */
				neumann_bc.value_list(fe_face_values.get_quadrature_points(),
										   neumann_values);

				for (unsigned int q_face_point = 0; q_face_point < n_face_q_points; ++q_face_point)
				{
					for (unsigned int i = 0; i < dofs_per_cell; ++i)
					{
						cell_rhs(i) += neumann_values[q_face_point] // g(x_q)
										* fe_face_values.shape_value(i, q_face_point) // phi_i(x_q)
										* fe_face_values.JxW(q_face_point); // dS
					} // end ++i
				} // end ++q_face_point
			} // end if
		} // end ++face_number


		// get global indices
		cell->get_dof_indices (local_dof_indices);
		/*
		 * Now add the cell matrix and rhs to the right spots
		 * in the global matrix and global rhs. Constraints will
		 * be taken care of later.
		 */
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
		{
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
			{
				system_matrix.add(local_dof_indices[i],
							local_dof_indices[j],
							cell_matrix(i, j));
			}
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	} // end ++cell
}


/*!
 * @brief Iterative solver.
 *
 * CG-based solver with SSOR-preconditioning.
 */
template <int dim>
void DiffusionProblem<dim>::solve_iterative ()
{
	SolverControl           solver_control (1000, 1e-12);
	SolverCG<>              solver (solver_control);

	PreconditionSSOR<> preconditioner;
	preconditioner.initialize(system_matrix, 1.2);

	solver.solve (system_matrix,
				solution,
				system_rhs,
				preconditioner);

	constraints.distribute (solution);

	std::cout << "   " << solver_control.last_step()
			<< " CG iterations needed to obtain convergence."
			<< std::endl;
}


/*!
 * @brief Write results to disk.
 *
 * Write results to disk in vtu-format.
 */
template <int dim>
void DiffusionProblem<dim>::output_results () const
{
	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, "solution");
	data_out.build_patches ();

	std::string filename = (dim == 2 ?
								"solution-std_2d" :
								"solution-std_3d" );
	filename += "_refinements-" + Utilities::int_to_string(n_refine, 1)
				+ ".vtu";

	std::ofstream output (filename.c_str());
	data_out.write_vtu (output);
}


/*!
 * @brief Run function of the object.
 *
 * Run the computation after object is built. Implements theping loop.
 */
template <int dim>
void DiffusionProblem<dim>::run ()
{
	std::cout << std::endl
				<< "===========================================" << std::endl;
	std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;

	make_grid ();

	setup_system ();

	assemble_system ();

	// Now solve
	constraints.condense(system_matrix, system_rhs);
	solve_iterative ();

	output_results ();

	std::cout << std::endl
			<< "===========================================" << std::endl;
}

} // end namespace DiffusionProblem


#endif /* INCLUDE_DIFFUSION_PROBLEM_HPP_ */
