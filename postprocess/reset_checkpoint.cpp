//
// Created by huacheng on 4/27/22.
//

#include "InputParser.h"
#include "CheckPoint.h"
#include "memory"

int main(int argc, char* argv[])
{
    char VtuFile[200];
    char ChkPrefix[200];
    char OutPrefix[200];
    char Default[] = "0";
    int MaxRank;
    InputFile * ifp = OpenInputFile2("../postprocess.inp");
    GetValueS(ifp,"process.pvtu",VtuFile,Default);
    GetValueS(ifp,"process.chkprefix",ChkPrefix,Default);
    GetValueS(ifp,"process.outprefix",OutPrefix,Default);
    MaxRank = GetValueI(ifp,"process.nrank",Default);
    CloseInputFile(ifp);

    std::cout << VtuFile << std::endl;
    vtkNew<vtkXMLPUnstructuredGridReader> reader;
    reader->SetFileName(VtuFile);
    reader->Update();
    vtkUnstructuredGrid * slope_vtu_data = vtkUnstructuredGrid::SafeDownCast(reader->GetOutputAsDataSet());

    std::cout << "number of points:" << slope_vtu_data->GetNumberOfPoints()
              << "\n" << "number of cells:" << slope_vtu_data->GetNumberOfCells() << std::endl;


    std::allocator<Mesh> mesh_allocator;
    Mesh * mesh_per_rank = mesh_allocator.allocate(MaxRank);
    for(int k_rank=0;k_rank<MaxRank;++k_rank)
    {
        Mesh * local_mesh = mesh_per_rank + k_rank;
        mesh_allocator.construct(local_mesh,0,k_rank);
        local_mesh->SetPrefix(ChkPrefix);
        local_mesh->LoadChk();
        local_mesh->MergeVtu(slope_vtu_data);
        mesh_allocator.destroy(local_mesh);
    }
    mesh_allocator.deallocate(mesh_per_rank,MaxRank);



    double test_x[3] = {0.,0.,-2.0e4};
    double weight[24] = {0.},pcoord[3] = {0.,0.,-2.0e4};
    int subId = -1;
    int targetId = slope_vtu_data->FindCell(test_x, nullptr, -1, 1.0e-6, subId, pcoord, weight);

    vtkCell * targetCell = slope_vtu_data->GetCell(targetId);

    size_t nc_points = targetCell->GetNumberOfPoints();
    for(size_t k=0;k<nc_points;k++)
    {
        double * tmp = targetCell->Points->GetPoint((vtkIdType)k);
        std::cout << k << " " << targetCell->PointIds->GetId((vtkIdType)k)
                  <<":  "<< tmp[0] << "," << tmp[1] << "," << tmp[2] << " = " << weight[k] << std::endl;

    }

    std::cout << "target_id:" << targetId <<std::endl;

    vtkPointData * vtuPointData = slope_vtu_data->GetPointData();
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
