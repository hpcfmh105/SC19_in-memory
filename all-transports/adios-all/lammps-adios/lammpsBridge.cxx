#include "lammpsBridge.h"
#include "lammpsDataAdaptor.h"
#include <vtkSmartPointer.h>
#include <ConfigurableAnalysis.h>


namespace lammpsBridge
{
  static vtkSmartPointer<senseiLammps::lammpsDataAdaptor>     GlobalDataAdaptor;
  static vtkSmartPointer<sensei::ConfigurableAnalysis> GlobalAnalysisAdaptor;

void Initialize(MPI_Comm world,int step, int nlocal, double **x,
               double xsublo, double xsubhi,
               double ysublo, double ysubhi,
               double zsublo, double zsubhi,
               const std::string& config_file)
{ 
  
  timer::Initialize();

  GlobalDataAdaptor = vtkSmartPointer<senseiLammps::lammpsDataAdaptor>::New();
  GlobalDataAdaptor->SetCommunicator(world);
  GlobalDataAdaptor->SetDataTimeStep(-1);

  lammpsBridge::GlobalDataAdaptor->Initialize(step, nlocal, x, xsublo, xsubhi,
                                       xsublo, xsubhi, xsublo, xsubhi);

  GlobalAnalysisAdaptor = vtkSmartPointer<sensei::ConfigurableAnalysis>::New();
  GlobalAnalysisAdaptor->Initialize(config_file);
}

void update_bridge(int step, /*double time,*/ double *buffer)
{
  //lammpsBridge::GlobalDataAdaptor->SetDataTime(time);
  lammpsBridge::GlobalDataAdaptor->SetDataTimeStep(step);
  //lammpsBridge::GlobalDataAdaptor->AddArray("buffer", buffer);
  lammpsBridge::GlobalAnalysisAdaptor->Execute(lammpsBridge::GlobalDataAdaptor.GetPointer());
  lammpsBridge::GlobalDataAdaptor->ReleaseData();
}
//-----------------------------------------------------------------------------
void Finalize()
{
  lammpsBridge::GlobalAnalysisAdaptor->Finalize();

  lammpsBridge::GlobalAnalysisAdaptor = nullptr;
  lammpsBridge::GlobalDataAdaptor = nullptr;

  timer::Finalize();  
}

}	// namespace lammpsBridge

