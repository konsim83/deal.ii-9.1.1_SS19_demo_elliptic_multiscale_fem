/*!
 * @file  DiffusionProblem.cc
 * @brief Contains main function.
 * @author Konrad Simon
 * @date August 2019
 */

// My Headers
#include "diffusion_problem.hpp"
#include "diffusion_problem_ms.hpp"

using namespace dealii;


int main ()

{
    deallog.depth_console (0);
    dealii::MultithreadInfo::set_thread_limit ();

    const unsigned int n_refine = 3,
    		n_refine_local = 4;

    const bool compute_2d = true,
			compute_3d = false;

    if (compute_2d)
    {
    	DiffusionProblem::DiffusionProblem<2> diffusion_problem_2d_coarse (n_refine);
    	diffusion_problem_2d_coarse.run ();

    	DiffusionProblem::DiffusionProblem<2> diffusion_problem_2d_fine (n_refine + n_refine_local);
		diffusion_problem_2d_fine.run ();

		DiffusionProblem::DiffusionProblemMultiscale<2> diffusion_ms_problem_2d (n_refine, n_refine_local);
		diffusion_ms_problem_2d.run ();
    }

    if (compute_3d)
    {
    	DiffusionProblem::DiffusionProblem<3> diffusion_problem_3d_coarse (n_refine);
		diffusion_problem_3d_coarse.run ();

		DiffusionProblem::DiffusionProblem<3> diffusion_problem_3d_fine (n_refine + n_refine_local);
		diffusion_problem_3d_fine.run ();

		DiffusionProblem::DiffusionProblemMultiscale<3> diffusion_ms_problem_3d (n_refine, n_refine_local);
		diffusion_ms_problem_3d.run ();
    }

  return 0;
}
