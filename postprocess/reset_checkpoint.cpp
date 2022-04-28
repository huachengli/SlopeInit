//
// Created by huacheng on 4/27/22.
//

#include "vtkMathUtilities.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkTesting.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"

int main(int argc, char* argv[])
{
    vtkNew<vtkXMLPUnstructuredGridReader> reader;
    if(argc < 2)
        reader->SetFileName("solution_00.pvtu");
    else
        reader->SetFileName(argv[1]);
    reader->Update();

    vtkUnstructuredGrid * unstructuredGrid = vtkUnstructuredGrid::SafeDownCast(reader->GetOutputAsDataSet());

    std::cout << "number of points:" << unstructuredGrid->GetNumberOfPoints()
              << "\n" << "number of cells:"<< unstructuredGrid->GetNumberOfCells() << std::endl;

    double test_x[3] = {0.,0.,-2.0e4};
    double weight[24] = {0.},pcoord[3] = {0.,0.,-2.0e4};
    int subId = -1;
    int targetId = unstructuredGrid->FindCell(test_x, nullptr,-1,1.0e-6,subId,pcoord,weight);

    vtkCell * targetCell = unstructuredGrid->GetCell(targetId);

    size_t nc_points = targetCell->GetNumberOfPoints();
    for(size_t k=0;k<nc_points;k++)
    {
        double * tmp = targetCell->Points->GetPoint((vtkIdType)k);
        std::cout << k << " " << targetCell->PointIds->GetId((vtkIdType)k)
                  <<":  "<< tmp[0] << "," << tmp[1] << "," << tmp[2] << " = " << weight[k] << std::endl;

    }

    std::cout << "target_id:" << targetId <<std::endl;

    vtkPointData * vtuPointData = unstructuredGrid->GetPointData();
    int n_array = vtuPointData->GetNumberOfArrays();
    vtkAbstractArray * stress_array_abstract = vtuPointData->GetAbstractArray("stress");
    auto * stress_array = vtkArrayDownCast<vtkFloatArray>(stress_array_abstract);

    vtkNew<vtkFloatArray> target_stress_array;
    target_stress_array->SetNumberOfComponents(stress_array->GetNumberOfComponents());
    target_stress_array->SetNumberOfTuples(targetCell->GetNumberOfPoints());
    target_stress_array->SetName("cell_stress");
    stress_array->GetTuples(targetCell->GetPointIds(),target_stress_array);

    target_stress_array->PrintSelf(std::cout,vtkIndent(3));
    size_t n_component = target_stress_array->GetNumberOfComponents();
    for(size_t k=0;k<nc_points;k++)
    {
        std::cout << k << " " << targetCell->PointIds->GetId((vtkIdType)k) << " " << std::endl;
        for(size_t j=0;j<n_component/3;++j)
        {
            std::cout << "  " << target_stress_array->GetValue((vtkIdType)k*n_component + j*3 + 0)
                      << ", " << target_stress_array->GetValue((vtkIdType)k*n_component + j*3 + 1)
                      << ", " << target_stress_array->GetValue((vtkIdType)k*n_component + j*3 + 2) << std::endl;
        }
        std::cout << "----------\n" ;

    }

    return EXIT_SUCCESS;
}
