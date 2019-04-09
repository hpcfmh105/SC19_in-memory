#pragma once

#include <DataAdaptor.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

namespace senseiLammps
{

class lammpsDataAdaptor : public sensei::DataAdaptor
{
public:
  static lammpsDataAdaptor* New();
  senseiTypeMacro(lammpsDataAdaptor *a, sensei::DataAdaptor);

  void Initialize( int nstep, int nlocal, double **x, 
                   double xsublo, double xsubhi, 
                   double ysublo, double ysubhi, 
                   double zsublo, double zsubhi);

  void GetBounds ( double &xsublo, double &xsubhi, 
                   double &ysublo, double &ysubhi, 
                   double &zsublo, double &zsubhi);

  void GetN ( int &nlocal );

  void GetPointers ( double **&x);

  void GetAtoms ( vtkDoubleArray *&buffers );

  /*void GetTypes ( vtkIntArray *&types );

  void GetIDs ( vtkIntArray *&ids );*/

// SENSEI API
  int GetNumberOfMeshes(unsigned int &numMeshes) override;

  int GetMeshName(unsigned int id, std::string &meshName) override;

  int GetMesh(const std::string &meshName, bool structureOnly,
    vtkDataObject *&mesh) override;

  int AddArray(vtkDataObject* mesh, const std::string &meshName,
    int association, const std::string &arrayName) override;

  int GetNumberOfArrays(const std::string &meshName, int association,
    unsigned int &numberOfArrays) override;

  int GetArrayName(const std::string &meshName, int association,
    unsigned int index, std::string &arrayName) override;

  int ReleaseData() override;

protected:
  lammpsDataAdaptor();
  ~lammpsDataAdaptor();

private:
  lammpsDataAdaptor(const lammpsDataAdaptor&); // not implemented.
  void operator=(const lammpsDataAdaptor&); // not implemented.

  struct DInternals;
  DInternals* Internals;
};

}
