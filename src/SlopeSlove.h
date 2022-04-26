//
// Created by huacheng on 4/21/22.
//

#ifndef SLOPEINIT_SLOPESLOVE_H
#define SLOPEINIT_SLOPESLOVE_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <fstream>
#include <iostream>

// import petsc warppers (linear algebra)

namespace LA
{
    using namespace dealii::LinearAlgebraPETSc;
#define USE_PETSC_LA
}

#include "SlopeGeometry.h"
#include "eos_cpp.h"

namespace SlopeEquationData
{
    using namespace dealii;
    template<int dim>
    class BoundaryForce : public Function<dim>
    {
    public:
        BoundaryForce();
        virtual double value(const Point<dim> &p,unsigned int component = 0) const override;
        virtual void vector_value(const Point<dim> &p,Vector<double> &values) const override;
    };

    template<int dim>
    BoundaryForce<dim>::BoundaryForce():Function<dim>(dim){}

    template<int dim>
    double BoundaryForce<dim>::value(const Point<dim> &,unsigned int) const
    {
        return 0.;
    }

    template<int dim>
    void BoundaryForce<dim>::vector_value(const Point<dim> &p,Vector<double> &values) const
    {
        for (unsigned int c = 0; c < this->n_components; ++c)
            values(c) = BoundaryForce<dim>::value(p, c);
    }


    template<int dim>
    class BoundaryValues : public Function<dim>
    {
    public:
        BoundaryValues();
        virtual double value(const Point<dim> &p,const unsigned int component = 0) const override;
    };

    template<int dim>
    BoundaryValues<dim>::BoundaryValues():Function<dim>(dim){}

    template<int dim>
    double BoundaryValues<dim>::value(const Point<dim> &,const unsigned int) const
    {
        return 0.;
    }


    template<int dim>
    class BodyForce : public Function<dim>
    {
    private:
        SlopeInfo * slope;
    public:
        BodyForce() = delete;
        explicit BodyForce(SlopeInfo * _slope);
        virtual double value(const Point<dim> &p,const unsigned int component = 0) const override;
        virtual void vector_value(const Point<dim> &p,Vector<double> &values) const override;
        virtual void vector_value_list(const std::vector<Point<dim>> &points, std::vector<Vector<double>> & value_list) const override;
    };

    template<int dim>
    BodyForce<dim>::BodyForce(SlopeInfo * _slope):Function<dim>(dim),slope(_slope){}

    template<int dim>
    double BodyForce<dim>::value(const Point<dim> &p, const unsigned int component) const
    {
        double pdepth = EvaluateDepth(slope,p[0],p[1],p[2]);
        int pindex = (int)(fabs(pdepth)/slope->dh);
        if(pindex >= slope->nstep) pindex = slope->nstep-1;
        if(component == (dim-1))
            return slope->Grav[pindex] * slope->Den[pindex];
        else
            return 0.;
    }

    template<int dim>
    void BodyForce<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
    {
        double pdepth = EvaluateDepth(slope,p[0],p[1],p[2]);
        int pindex = (int)(fabs(pdepth)/slope->dh);
        if(pindex >= slope->nstep) pindex = slope->nstep-1;

        for (unsigned int c = 0; c < this->n_components; ++c)
        {
            if((dim-1)==c)
                values(c) = slope->Grav[pindex] * slope->Den[pindex];
            else
                values(c) = 0.;
        }
    }

    template<int dim>
    void BodyForce<dim>::vector_value_list(const std::vector<Point<dim>> &points, std::vector<Vector<double>> &  value_list) const
    {
        const unsigned int n_points = points.size();

        AssertDimension(value_list.size(), n_points);

        for (unsigned int p = 0; p < n_points; ++p)
            this->vector_value(points[p], value_list[p]);
    }

    template <int dim>
    class ConstitutiveLaw
    {
    public:
        ConstitutiveLaw(SlopeInfo * _slope);
        void GetStrainStressTensor(const Point<dim> &p,SymmetricTensor<4, dim> & c);
        void GetStrainStressTensorList(const std::vector<Point<dim>> &pts,std::vector<SymmetricTensor<4, dim>> &clist);
    private:
        SymmetricTensor<4,dim> c_kappa;
        SymmetricTensor<4,dim> c_mu;
        SlopeInfo * slope;
    };

/*    template<int dim>
    ConstitutiveLaw<dim>::ConstitutiveLaw(double E, double nu):
            kappa(E/(3.*(1.0-2.0*nu))),mu(E/(2.0*(1.0+nu))),
            c_kappa(kappa* outer_product(unit_symmetric_tensor<dim>(),unit_symmetric_tensor<dim>())),
            c_mu(2.0*mu*(identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(),unit_symmetric_tensor<dim>())/3.0))
    {}*/

    template<int dim>
    ConstitutiveLaw<dim>::ConstitutiveLaw(SlopeInfo * _slope):
    slope(_slope),
    c_kappa(outer_product(unit_symmetric_tensor<dim>(),unit_symmetric_tensor<dim>())),
    c_mu((identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(),unit_symmetric_tensor<dim>())/3.0))
    {}

    template<int dim>
    void ConstitutiveLaw<dim>::GetStrainStressTensor(const Point<dim> &p, SymmetricTensor<4, dim> &c)
    {
        double pdepth = EvaluateDepth(slope,p[0],p[1],p[2]);
        int pindex = (int)(fabs(pdepth)/slope->dh);
        if(pindex >= slope->nstep) pindex = slope->nstep-1;

        double tmp_kappa = slope->Den[pindex]*slope->Cs[pindex]*slope->Cs[pindex];
        double tmp_mu = tmp_kappa*3.0*(1.-2.*slope->Nu[pindex])/(2.0*(1.0+slope->Nu[pindex]));
        c = tmp_kappa*c_kappa + 2.0*tmp_mu*c_mu;
    }

    template<int dim>
    void ConstitutiveLaw<dim>::GetStrainStressTensorList(const std::vector<Point<dim>> &pts,
                                                         std::vector<SymmetricTensor<4, dim>> &clist)
    {
        size_t n_points = pts.size();
        for(size_t k=0;k<n_points;++k)
        {
            this->GetStrainStressTensor(pts[k],clist[k]);
        }
    }

    template<int dim>
    class DomainTransform
    {
    public:
        DomainTransform() = delete;
        explicit DomainTransform(double (*pts)[dim])
        {
            for(size_t k=0;k<(1<<dim);++k)
            {
                for(size_t j=0;j<dim;++j)
                    this->corner_point[k][j] = pts[k][j];
            }
        }
        Point<dim> operator()(const Point<dim> &p) const
        {
            assert(3==dim);
            double tmp_p[dim],result_p[dim];
            for(size_t k=0;k<dim;++k) tmp_p[k] = p[k];
            TransformInterpolate(result_p,tmp_p,this->corner_point);
            return Point<dim>(result_p[0],result_p[1],result_p[2]);
        }
    private:
        double corner_point[1<<dim][dim];
    };

    template<int dim>
    uintptr_t GetMaterial(Point<dim> p, SlopeInfo * _s)
    {
        double &px = p[0];
        double &py = p[1];
        double &pz = p[2];

        double r = _s->equation[0]*px + _s->equation[1]*py + _s->equation[2]*pz + _s->equation[3];
        assert(r < 1.0e4);
        int k = 0;
        while(k<_s->nlayers)
        {
            if(_s->sum_depth[k] + r > 0) break;
            ++k;
        }
        return _s->mdata[k];
    }

    void InitSlopeProfile(SlopeInfo * _s);
    void InitSlopeEos(SlopeInfo * pSlopeInfo);

}

namespace OutHelper
{
    using namespace dealii;
    template <int dim>
    class StrainPostprocessor : public DataPostprocessorTensor<dim>
    {
    private:
        SlopeEquationData::ConstitutiveLaw<dim> constitutive_law;
    public:
        explicit StrainPostprocessor(SlopeInfo * sinfo):
        DataPostprocessorTensor<dim>("stress", update_gradients|update_quadrature_points),
        constitutive_law(sinfo){}

        virtual void evaluate_vector_field(
                const DataPostprocessorInputs::Vector<dim> &input_data,
                std::vector<Vector<double>> &               computed_quantities) const
        {
            AssertDimension(input_data.solution_gradients.size(),
                            computed_quantities.size());

            for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
            {
                AssertDimension(computed_quantities[p].size(),
                                (Tensor<2, dim>::n_independent_components));

                SymmetricTensor<4,dim> c_tensor;
                SymmetricTensor<2,dim> t_stress,t_strain;
                constitutive_law.GetStrainStressTensor(input_data.evaluation_points[p],c_tensor);

                for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int e = 0; e < dim; ++e)
                    {
                        t_strain[d][e] = (input_data.solution_gradients[p][d][e] +
                                          input_data.solution_gradients[p][e][d])/2;
                    }

                t_stress = c_tensor*t_strain;

                for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int e = 0; e < dim; ++e)
                    {
                        computed_quantities[p][Tensor<2, dim>::component_to_unrolled_index(TableIndices<2>(d, e))]
                        = t_stress[d][e];
                    }



            }
        }
    };
}


namespace Slope
{
    using namespace dealii;

    template<int dim>
    class SlopeProblem
    {
    public:
        SlopeProblem();
        void run();
        void set_info(SlopeInfo * sinfo);


    private:
        void setup_system();
        void assemble_system();
        void solve_system();
        void make_grid_trapezoid(int refine_time = 0);
        void refine_grid();
        void output_results(const unsigned int cycle) const;
        MPI_Comm mpi_communicator;
        parallel::distributed::Triangulation<dim> triangulation;
        DoFHandler<dim> dof_handler;
        const MappingQ<dim> mapping;
        FESystem<dim> fe;
        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;
        AffineConstraints<double> constraints;
        LA::MPI::SparseMatrix system_matrix;
        LA::MPI::Vector       locally_relevant_solution;
        LA::MPI::Vector       system_rhs;
        ConditionalOStream pcout;
        TimerOutput computing_timer;
        SlopeInfo * slope_info;
    };

    template<int dim>
    SlopeProblem<dim>::SlopeProblem()
            :mpi_communicator(MPI_COMM_WORLD),
             triangulation(mpi_communicator,typename Triangulation<dim>::MeshSmoothing(
                     Triangulation<dim>::smoothing_on_refinement
                     |Triangulation<dim>::smoothing_on_coarsening)),
             dof_handler(triangulation),mapping(2),fe(FE_Q<dim>(1), dim),
             pcout(std::cout,(Utilities::MPI::this_mpi_process(mpi_communicator)==0)),
             computing_timer(mpi_communicator,pcout,TimerOutput::summary,TimerOutput::wall_times)
    {}

    template<int dim>
    void SlopeProblem<dim>::set_info(SlopeInfo *sinfo)
    {
        this->slope_info = sinfo;
    }

    template<int dim>
    void SlopeProblem<dim>::setup_system()
    {
        TimerOutput::Scope t(computing_timer, "setup");

        dof_handler.distribute_dofs(fe);
        // information about local_owned and ghost cell
        locally_owned_dofs = dof_handler.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler,locally_relevant_dofs);
        locally_relevant_solution.reinit(locally_owned_dofs,locally_relevant_dofs,mpi_communicator);
        system_rhs.reinit(locally_owned_dofs,mpi_communicator);

        // hanging node and boundary value constraints
        constraints.clear();
        constraints.reinit(locally_relevant_dofs);
        DoFTools::make_hanging_node_constraints(dof_handler,constraints);

        // boundary_id:0 the bottom, fixed boundary
        VectorTools::interpolate_boundary_values(dof_handler,0,
                                                 SlopeEquationData::BoundaryValues<dim>(),constraints);

        // boundary_id:2 the left/right/front/back side face, robbin boundary
        const FEValuesExtractors::Scalar x_displacement(0);
        const FEValuesExtractors::Scalar y_displacement(1);
        const FEValuesExtractors::Scalar z_displacement(2);

        /*
        VectorTools::interpolate_boundary_values(dof_handler, 2,
                                                 SlopeEquationData::BoundaryValues<dim>(),constraints,
                                                 (fe.component_mask(x_displacement)|fe.component_mask(z_displacement)));
        */
        std::set<types::boundary_id> robbin_boundary_ids;
        robbin_boundary_ids.insert(2);
        VectorTools::compute_no_normal_flux_constraints(dof_handler,0,robbin_boundary_ids,constraints,mapping);
        constraints.close();

        // init the system matrix
        DynamicSparsityPattern dsp(locally_relevant_dofs);
        DoFTools::make_sparsity_pattern(dof_handler,dsp,constraints, false);
        SparsityTools::distribute_sparsity_pattern(dsp,
                                                   dof_handler.locally_owned_dofs(),
                                                   mpi_communicator,
                                                   locally_relevant_dofs);
        system_matrix.reinit(locally_owned_dofs,locally_owned_dofs,dsp,mpi_communicator);
    }

    template<int dim>
    void SlopeProblem<dim>::assemble_system()
    {
        TimerOutput::Scope t(computing_timer, "assembly");

        QGauss<dim> quadrature_formula(fe.degree + 1);
        QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

        FEValues<dim> fe_values(mapping,fe,quadrature_formula,
                                update_values | update_gradients
                                | update_quadrature_points|update_JxW_values);
        FEFaceValues<dim> fe_values_face(mapping,fe,face_quadrature_formula,
                                         update_values | update_gradients |
                                         update_quadrature_points | update_normal_vectors |
                                         update_JxW_values);

        const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
        const unsigned int n_q_points = quadrature_formula.size();
        const unsigned int n_face_q_points = face_quadrature_formula.size();

        // f and k on local cell
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        FullMatrix<double> t_cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs(dofs_per_cell);
        // information on cell
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        std::vector<double> lambda_values(n_q_points);
        std::vector<double> mu_values(n_q_points);
        std::vector<Vector<double>> rhs_values(n_q_points,Vector<double>(dim));
        std::vector<Vector<double>> boundary_force_values(n_face_q_points,Vector<double>(dim));

//        typename DoFHandler<dim>::active_cell_iterators
//        cell = dof_handler.begin_active(),endc = dof_handler.end();


        SlopeEquationData::ConstitutiveLaw<dim> constitutive_law(this->slope_info);
        std::vector<SymmetricTensor<4, dim>> c_tensor(n_q_points, identity_tensor<dim>());

        SlopeEquationData::BodyForce<dim> gravity_field(this->slope_info);

        for(const auto &cell: dof_handler.active_cell_iterators())
        {
            if(!cell->is_locally_owned()) continue;
            cell_matrix = 0.;
            cell_rhs = 0.;
            fe_values.reinit(cell);
            t_cell_matrix = 0.;
            // mass matrix

            /*
             * use tensor operation instead of shape_grad
             */

            constitutive_law.GetStrainStressTensorList(fe_values.get_quadrature_points(),c_tensor);
            const FEValuesExtractors::Vector displacement(0);
            for (const unsigned int q_point: fe_values.quadrature_point_indices())
                for (const unsigned int i: fe_values.dof_indices())
                {
                    const SymmetricTensor<2,dim> stress_phi_i = c_tensor[q_point]*fe_values[displacement].symmetric_gradient(i,q_point);
                    for(const unsigned int j:fe_values.dof_indices())
                    {
                        cell_matrix(i,j) += stress_phi_i * fe_values[displacement].symmetric_gradient(j,q_point)
                                            * fe_values.JxW(q_point);
                    }
                }




            // the body force term
            gravity_field.vector_value_list(fe_values.get_quadrature_points(),rhs_values);
            for(const unsigned int i: fe_values.dof_indices())
            {
                const unsigned int component_i = fe.system_to_component_index(i).first;
                for (const unsigned int q_point: fe_values.quadrature_point_indices()) {
                    cell_rhs(i) += fe_values.shape_value(i, q_point)
                                   * rhs_values[q_point][component_i] * fe_values.JxW(q_point);
                }
            }


            for (const auto &face : cell->face_iterators())
            {
                if(face->at_boundary() == false) continue;
                if(face->boundary_id() == 1)
                {
                    // process neumann boundary condition
                    fe_values_face.reinit(cell, face);
                    const SlopeEquationData::BoundaryForce<dim> f_boundary_force;
                    f_boundary_force.vector_value_list(fe_values_face.get_quadrature_points(),boundary_force_values);

                    for (unsigned int q_point = 0; q_point < n_face_q_points;++q_point)
                    {
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            const unsigned int component_i = fe.system_to_component_index(i).first;
                            cell_rhs(i) += fe_values.shape_value(i,q_point)
                                           * boundary_force_values[q_point][component_i]
                                           * fe_values_face.JxW(q_point);
                        }
                    }
                } else if(face->boundary_id() == 0)
                {
                    // no flux boundary has been implemented in affine_constrain
                    // flowing is the no shear stress condition
                    fe_values_face.reinit(cell, face);
                    t_cell_matrix = 0.;
                    for(unsigned int q_point=0;q_point<n_face_q_points;++q_point)
                    {
                        const Tensor<1,dim> &q_norm = fe_values_face.normal_vector(q_point);
                        const Tensor<2,dim> &q_norm2= outer_product(q_norm,q_norm);
                        for(const unsigned int j: fe_values.dof_indices())
                        {
                            const SymmetricTensor<2,dim> norm_stress = c_tensor[q_point]*fe_values_face[displacement].symmetric_gradient(j,q_point);
                            double shear_boundary = norm_stress * symmetrize(q_norm2);
                            for(const unsigned int i: fe_values.dof_indices())
                            {
                                double shape_w = fe_values_face[displacement].value(i,q_point)*q_norm;
                                cell_matrix(i,j) += shape_w * shear_boundary * fe_values_face.JxW(q_point);
                                t_cell_matrix(i,j) += shape_w * shear_boundary * fe_values_face.JxW(q_point);
                            }
                        }
                    }
                }
            }

            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices,
                                                   system_matrix, system_rhs);
        }

        system_matrix.compress(VectorOperation::add);
        system_rhs.compress(VectorOperation::add);
    }

    template<int dim>
    void SlopeProblem<dim>::solve_system()
    {
        TimerOutput::Scope t(computing_timer, "solve");

        LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,mpi_communicator);

        SolverControl solver_control(dof_handler.n_dofs(),1e-12);
        LA::SolverCG solver_cg(solver_control,mpi_communicator);
        LA::MPI::PreconditionAMG preconditioner;
        LA::MPI::PreconditionAMG::AdditionalData data;

        data.symmetric_operator = true;
        // other control on amg setting
        preconditioner.initialize(system_matrix,data);
        solver_cg.solve(system_matrix,
                        completely_distributed_solution,
                        system_rhs,preconditioner);

        constraints.distribute(completely_distributed_solution);
        locally_relevant_solution = completely_distributed_solution;
    }

    template<int dim>
    void SlopeProblem<dim>::refine_grid()
    {
        TimerOutput::Scope t(computing_timer, "refine");

        Vector<float> estimate_error_per_cell(triangulation.n_active_cells());
        KellyErrorEstimator<dim>::estimate(
                dof_handler,
                QGauss<dim-1>(fe.degree+1),
                {},
                locally_relevant_solution,
                estimate_error_per_cell);
        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(triangulation,estimate_error_per_cell,0.3,0.03);
        triangulation.execute_coarsening_and_refinement();
    }

    template<int dim>
    void SlopeProblem<dim>::make_grid_trapezoid(int refine_time)
    {
        std::vector< unsigned int > repetitions(dim, 16);
        const Point<dim> p1 = (dim == 3 ? Point<dim>(-1.0, -1.0, -1.0) : Point<dim>(-1.0, -1.0));
        const Point<dim> p2 = (dim == 3 ? Point<dim>( 1.0,  1.0,  1.0) : Point<dim>( 1.0,  1.0));
        GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,p1,p2);

        for(const auto &cell:triangulation.active_cell_iterators())
        {
            for(const auto &face:cell->face_iterators())
            {
                if(!face->at_boundary()) continue;
                if(std::fabs(face->center()[2]-1.0)<1e-10)
                    face->set_boundary_id(1);
                else if(std::fabs(face->center()[2]+1.0)<1e-10)
                    face->set_boundary_id(0);
                else
                    face->set_boundary_id(2);
            }
        }

        double t_corner[1<<dim][dim];
        GetTransformCorner(this->slope_info,t_corner);
        const SlopeEquationData::DomainTransform<dim> local_transformer(t_corner);

        GridTools::transform(local_transformer, triangulation);
//        GridTools::scale(1.0e5,triangulation);

        if(refine_time > 0)
            triangulation.refine_global(refine_time);
    }

    template<int dim>
    void SlopeProblem<dim>::output_results(const unsigned int cycle) const
    {
        DataOut<dim> data_out;
        OutHelper::StrainPostprocessor<dim> strain_calculator;
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
                data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(locally_relevant_solution,
                                 std::vector<std::string>(dim, "displacement"),
                                 DataOut<dim>::type_dof_data,
                                 data_component_interpretation);

        data_out.add_data_vector(locally_relevant_solution,
                                 strain_calculator);
        Vector<float> subdomain(triangulation.n_active_cells());
        for(unsigned int i=0;i<subdomain.size();++i)
        {
            subdomain(i) = triangulation.locally_owned_subdomain();
        }

        data_out.add_data_vector(subdomain,"subdomain");
        data_out.build_patches();
        data_out.write_vtu_with_pvtu_record("./","solution",cycle,mpi_communicator,2,4);
    }

    template<int dim>
    void SlopeProblem<dim>::run()
    {
        pcout << "Running with PETSc" << " on "
              << Utilities::MPI::n_mpi_processes(mpi_communicator)
              << " MPI rank(s)..." << std::endl;

        const unsigned int n_cycle = 2;
        for(unsigned int cycle=0;cycle<n_cycle;++cycle)
        {
            pcout<< "CYCLE " << cycle << ":" << std::endl;
            if(cycle == 0)
            {
                make_grid_trapezoid();
            } else
                refine_grid();

            setup_system();
            pcout << "   Number of active cells:       "
                  << triangulation.n_global_active_cells() << std::endl
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl;

            assemble_system();
            solve_system();

            if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
            {
                TimerOutput::Scope t(computing_timer, "output");
                output_results(cycle);
            }
            computing_timer.print_summary();
            computing_timer.reset();
            pcout << std::endl;
        }
    }
}



#endif //SLOPEINIT_SLOPESLOVE_H
