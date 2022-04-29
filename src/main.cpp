#include <iostream>
#include "InputParser.h"
#include "SlopeGeometry.h"
#include "eos_cpp.h"
#include "SlopeSlove.h"

int main(int argc,char * argv[])
{
    SlopeInfo pSlopeInfo;
    LoadSlopeInfo(&pSlopeInfo,"../SALEc.inp");
    SlopeEquationData::InitSlopeEos(&pSlopeInfo);
    SlopeEquationData::InitSlopeProfile(&pSlopeInfo);

   /*
    double Xc[8][3];
    GetTransformCorner(&pSlopeInfo,Xc);
    printf("%d\n",pSlopeInfo.noffset);
    */

    try {
        using namespace dealii;
        using namespace Slope;
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,argv,1);

        SlopeProblem<3> slope_problem;
        slope_problem.set_info(&pSlopeInfo);
        slope_problem.run();

       /* if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
        {
            for(int k=0;k<3;++k)
            {
                for(int j=0;j<3;++j)
                    std::cout << Tensor<2, 3>::component_to_unrolled_index(TableIndices<2>(k, j)) << std::endl;
            }
        }*/
    }
    catch(std::exception &exc)
    {
        std::cerr << std::endl
                  << std::endl
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
        std::cerr << std::endl
                  << std::endl
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
