#ifndef VTKMeshWriter_h
#define VTKMeshWriter_h

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkHexagonalPrism.h>
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkPentagonalPrism.h>
#include <vtkPixel.h>
#include <vtkPolyLine.h>
#include <vtkPolyVertex.h>
#include <vtkPolygon.h>
#include <vtkPyramid.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkTriangleStrip.h>
#include <vtkVertex.h>
#include <vtkVoxel.h>
#include <vtkWedge.h>

#include <string>
#include <memory>

namespace HEM{

class VTKMeshWriter
{
public:


public:
    VTKMeshWriter()
    {
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }

    template<class Mesh>
    void set_mesh(Mesh & mesh)
    {
      m_ugrid->Initialize(); //把网格清空
      set_points(mesh);
      set_cells(mesh);
    }


    template<class Mesh>
    void set_points(Mesh & mesh)
    {
      using Node = typename Mesh::Node;
      auto NN = mesh.number_of_nodes();
      auto points = vtkSmartPointer<vtkPoints>::New();
      points->Allocate(NN);
      auto GD = mesh.geo_dimension();

      if(GD == 3)
      {
        auto func = [&](const Node & pp)->void
        {
          auto coordinate = pp.coordinate();
          points->InsertNextPoint(coordinate[0], coordinate[1], coordinate[2]);
        };
        mesh.for_each_node(func);
      }
      else if(GD == 2)
      {
        auto func = [&](const Node & pp)->void
        {
          auto coordinate = pp.coordinate();
          points->InsertNextPoint(coordinate[0], coordinate[1], 0.0);
        };
        mesh.for_each_entity(func);
      }
      m_ugrid->SetPoints(points);
    }

    template<class Mesh>
    void set_cells(Mesh & mesh)
    {
      using Node = typename Mesh::Node;
      using Cell = typename Mesh::Cell;
      auto NC = mesh.number_of_cells();
      auto cells = vtkSmartPointer<vtkCellArray>::New();

      Node * c2n[128];
      auto & nidx = *mesh.get_node_indices();
      auto func = [&](const Cell & cell)->void
      {
        int N = cell.cell_to_node(c2n);
        cells->InsertNextCell(N); 
        for(int i = 0; i < N; i++)
          cells->InsertCellPoint(nidx[c2n->index()]);
      };
      mesh.for_each_entity(func);
      m_ugrid->SetCells(7, cells);
    }

    template<typename T>
    void set_point_data(std::vector<T> & data, int ncomponents, const std::string name)
    {
        int n = data.size()/ncomponents;
        auto vtkdata = vtkSmartPointer<vtkAOSDataArrayTemplate<T>>::New();
        vtkdata->SetNumberOfComponents(ncomponents);
        vtkdata->SetNumberOfTuples(n);
        vtkdata->SetName(name.c_str());
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < ncomponents; j ++)
                vtkdata->SetComponent(i, j, data[i*ncomponents + j]);
        }
        m_ugrid->GetPointData()->AddArray(vtkdata);
    }

    template<typename T>
    void set_cell_data(std::vector<T> & data, int ncomponents, const std::string name)
    {
        int n = data.size()/ncomponents;
        auto vtkdata = vtkSmartPointer<vtkAOSDataArrayTemplate<T>>::New();
        vtkdata->SetNumberOfComponents(ncomponents);
        vtkdata->SetNumberOfTuples(n);
        vtkdata->SetName(name.c_str());
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < ncomponents; j ++)
                vtkdata->SetComponent(i, j, data[i*ncomponents + j]);
        }
        m_ugrid->GetCellData()->AddArray(vtkdata);
    }

    void write(const std::string & fname)
    {
        m_writer->SetFileName(fname.c_str());
        m_writer->SetInputData(m_ugrid);
        m_writer->Write();
    }

private:
  vtkSmartPointer<vtkUnstructuredGrid> m_ugrid;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> m_writer;
};

} // end of namespace Mesh

#endif // end of VTKMeshWriter_h
