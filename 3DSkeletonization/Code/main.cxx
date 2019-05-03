#include <vtkPolyData.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyhedron.h>
#include <vtkPolyLine.h>
#include <vtkFeatureEdges.h>
#include <vtkDiskSource.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <sstream>
#include <map>

//#define DEBUG

void WriteMeshToVTK(vtkSmartPointer<vtkUnstructuredGrid> polyData,
			std::string filename);

vtkSmartPointer<vtkUnstructuredGrid> ReadMeshFromVTK(std::string filename);

int main ( int argc, char *argv[] )
{
  if(argc != 2)
  {
    std::cout << "Usage: " << argv[0] << "  Filename(.vtk)"
	      << " <Number of iterations>" << std::endl;
    return EXIT_FAILURE;
  }

  vtkSmartPointer<vtkPolyData> dualMesh = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> dualMeshPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> dualMeshCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> dualMeshLines = vtkSmartPointer<vtkCellArray>::New();
  dualMesh->SetPoints(dualMeshPoints);
  dualMesh->SetPolys(dualMeshCells);
  dualMesh->SetLines(dualMeshLines);
    
  std::string inputFilename = argv[1];
  vtkSmartPointer<vtkUnstructuredGrid> mesh = ReadMeshFromVTK(inputFilename);
  
  
  mesh->BuildLinks();
  mesh->ComputeBounds();

  vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
  
  unsigned int numberOfPoints = points->GetNumberOfPoints();
  for (unsigned int pointCounter = 0; pointCounter < numberOfPoints;
       ++pointCounter)
  {
    double* point = points->GetPoint(pointCounter);
    std::cout << pointCounter << std::endl;
    vtkSmartPointer<vtkIdList> adjacentCells = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(pointCounter, adjacentCells);
    if (adjacentCells->GetNumberOfIds() > 1)
    {
      // This loop extracts the neighbourhood of the current point. 
      std::vector<vtkIdType> neighbouringVertices;
      for (unsigned int i = 0;i < adjacentCells->GetNumberOfIds();++i)
      {
	vtkSmartPointer<vtkIdList> vertices = vtkSmartPointer<vtkIdList>::New();
	mesh->GetCellPoints(adjacentCells->GetId(i), vertices);
        for (unsigned int k = 0;k < vertices->GetNumberOfIds();++k)
	{
	  if (vertices->GetId(k) != pointCounter)
	    neighbouringVertices.push_back(vertices->GetId(k));
	}
      }
    }
  }

  /**
     Create the cell data to store the index to the dual vertex of the cell. 
     in the primal triangular cells. 
  */
  vtkSmartPointer<vtkIdTypeArray> cellData =
    vtkSmartPointer<vtkIdTypeArray>::New();
  mesh->GetCellData()->SetScalars(cellData);

  /*
   * Contains information related to the cells. 
   * - Index of the primal vertex at the origin of the dual cell.
   * - Value attached to the primal vertex (curvature, indication value, etc).
   *
   */
  vtkSmartPointer<vtkFloatArray> arrayData =
    vtkSmartPointer<vtkFloatArray>::New();
  arrayData->SetNumberOfComponents(5);
  arrayData->SetName ("ArrayData");
  mesh->GetCellData()->AddArray(arrayData);
  
  /**
   * Traverse all the 2-dimensional cells in the polyData in order to create 
   * the dual vertices. 
   * This loop also detects the primal cells (faces) that are only connected to 
   * one other face. 
   * 
   */
  vtkSmartPointer<vtkCellArray> primalCells = mesh->GetCells();
  vtkSmartPointer<vtkIdList> idListPoints = vtkSmartPointer<vtkIdList>::New();
  vtkIdType cellCounter = 0;
  for (cellCounter = 0, primalCells->InitTraversal();
       primalCells->GetNextCell(idListPoints) ;
       ++cellCounter)
  {
    std::cout << "Les ids des points de la cellule " << cellCounter << " sont ";
    for (unsigned int i = 0;i < idListPoints->GetNumberOfIds();++i)
    {
      std::cout << " " << idListPoints->GetId(i);
    }
    std::cout << std::endl;
    /* 
       Adding data to the primal and dual meshes. We have to store the coordinates 
       and the index of the barycenter in each cell of the primal mesh. 
    */
    
    double* cellCounterData = new double(cellCounter);
   
    vtkIdType position = mesh->GetCellData()->GetScalars()
      ->InsertNextTuple(cellCounterData);
    delete cellCounterData;

    double* cellData = new double[5];
    cellData[0] = 1;
    cellData[1] = 1;
    cellData[2] = -1.0;
    cellData[3] = -1.0;
    cellData[4] = 0.0;
    vtkIdType idNewCellData = arrayData->InsertNextTuple(cellData);
    delete[] cellData;    
  }
  
  //Write the dual dual PolyData.
  std::stringstream ssDualDual;
  ssDualDual << "ProcessedMesh"; 
  WriteMeshToVTK(mesh, ssDualDual.str().c_str());
 
  return EXIT_SUCCESS;
}


void WriteMeshToVTK(vtkSmartPointer<vtkUnstructuredGrid> polyData,
			std::string filename)
{
  //Write the mesh.
  vtkSmartPointer<vtkUnstructuredGridWriter> writer =  
    vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  std::stringstream ssFilename;
  ssFilename << filename << ".vtk";
  writer->SetFileName(ssFilename.str().c_str());  
  writer->SetInputData(polyData);
  writer->SetFileTypeToASCII();
  writer->Update();
}

vtkSmartPointer<vtkUnstructuredGrid> ReadMeshFromVTK(std::string filename)
{
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
    vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName ( filename.c_str() );
  reader->Update();
  
  //Obtain the polydata.
  vtkSmartPointer<vtkUnstructuredGrid> polyData = reader->GetOutput();
  return polyData;
}



 
