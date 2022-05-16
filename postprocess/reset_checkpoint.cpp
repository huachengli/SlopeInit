//
// Created by huacheng on 4/27/22.
//

#include "InputParser.h"
#include "CheckPoint.h"
#include <memory>
#include <future>
#include <thread>
#include <ctime>


int main(int argc, char* argv[])
{
    char VtuFile[200];
    char ChkPrefix[200];
    char OutPrefix[200];
    char Default[] = "0";
    int NumOfChk;
    int num_thread = 4;


    InputFile * ifp = NULL;
    if(argc == 2)
    {
        ifp = OpenInputFile2(argv[1]);
        char TmpVtuFile[200];
        GetValueS(ifp,"process.pvtu",TmpVtuFile,Default);
        sprintf(VtuFile,"%s/%s",argv[1],TmpVtuFile);
    }
    else
    {
        ifp = OpenInputFile2("../postprocess.inp");
        GetValueS(ifp,"process.pvtu",VtuFile,Default);
    }


    GetValueS(ifp,"process.chkprefix",ChkPrefix,Default);
    GetValueS(ifp,"process.outprefix",OutPrefix,Default);
    NumOfChk = GetValueI(ifp, "process.nchk", Default);
    num_thread = GetValueI(ifp, "process.nrank", "4");
    CloseInputFile(ifp);

    time_t current_time = time(0);
    char ascii_cur_time[50];
    ctime_r(&current_time,ascii_cur_time);
    *(ascii_cur_time+ strlen(ascii_cur_time) -1) = '\0';
    std::cout << "[" << ascii_cur_time << "]:process" << VtuFile << std::endl;
    vtkNew<vtkXMLPUnstructuredGridReader> reader;
    reader->SetFileName(VtuFile);
    reader->Update();
    vtkUnstructuredGrid * slope_vtu_data = vtkUnstructuredGrid::SafeDownCast(reader->GetOutputAsDataSet());

    std::allocator<Mesh> mesh_allocator;
    Mesh * mesh_per_rank = mesh_allocator.allocate(NumOfChk);


    std::vector<std::future<int>> thread_list;
    thread_list.reserve(num_thread);
    for(int k_thread=0;k_thread<num_thread;++k_thread)
    {
        thread_list.push_back(std::async(std::launch::async,[&](int max_mesh_offset,int user_thread_id,int max_thread_id){
            for(int mesh_offset=0;mesh_offset<NumOfChk;++mesh_offset)
            {
                if(user_thread_id != (mesh_offset%max_thread_id)) continue;
                Mesh * local_mesh = mesh_per_rank + mesh_offset;
                mesh_allocator.construct(local_mesh,0,mesh_offset);
                local_mesh->SetPrefix(ChkPrefix);
                local_mesh->LoadChk();
                local_mesh->MergeVtu(slope_vtu_data);
                local_mesh->ExportChk(OutPrefix);
                mesh_allocator.destroy(local_mesh);
            }
            return user_thread_id;
        },NumOfChk,k_thread,num_thread));
    }


    for(auto &thread_handle:thread_list)
    {
        std::cout << "join thread:" << thread_handle.get() << std::endl;
    }

   /* for(int k_rank=0; k_rank < NumOfChk; ++k_rank)
    {

        Mesh * local_mesh = mesh_per_rank + k_rank;
        mesh_allocator.construct(local_mesh,0,k_rank);
        local_mesh->SetPrefix(ChkPrefix);
        local_mesh->LoadChk();
        local_mesh->MergeVtu(slope_vtu_data);
        local_mesh->ExportChk(OutPrefix);
        mesh_allocator.destroy(local_mesh);

    }*/
    mesh_allocator.deallocate(mesh_per_rank, NumOfChk);
    return EXIT_SUCCESS;
}
