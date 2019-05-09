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
#include <vector>
#include <vtkProperty.h>

#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

//#define DEBUG

void WriteMeshToVTK(vtkSmartPointer<vtkUnstructuredGrid> polyData, std::string filename);

vtkSmartPointer<vtkUnstructuredGrid> ReadMeshFromVTK(std::string filename);

void findDualPoints(int count, vtkIdType cellCounter, vtkSmartPointer<vtkUnstructuredGrid> &mesh, vtkSmartPointer<vtkPoints> &dualMeshPoints,
	 vtkSmartPointer<vtkIdList> &idListPoints);

void setDuaLine(int count,vtkIdType cellCounter, vtkSmartPointer<vtkIdList> idListPoints, vtkSmartPointer<vtkUnstructuredGrid> &mesh, 
	 vtkSmartPointer<vtkUnstructuredGrid> &dualMesh);
	 
void WriteMeshToPolyVTK(vtkSmartPointer<vtkPolyData> polyData, std::string filename);

void getEdgeCells (vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId, vtkIdList * cellIdsNeighborsFromEdge 
	 /*vtkSmartPointer<vtkCellArray> & polyMeshCells*/);

int main ( int argc, char *argv[] )
{
  if(argc != 2)
  {
    std::cout << "Usage: " << argv[0] << "  Filename(.vtk)"
	      << " <Number of iterations>" << std::endl;
    return EXIT_FAILURE;
  }

  ///////// Color /////////

  vtkSmartPointer<vtkUnsignedCharArray> colors =
  vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName ("Colors");
  
  /*colors->InsertNextTupleValue(red);
  colors->InsertNextTupleValue(green);
  colors->InsertNextTupleValue(blue);
  */

  unsigned char red[3] = {255, 0, 0};
  unsigned char green[3] = {0, 255, 0};
  unsigned char blue[3] = {0, 0, 255};

  //polydata->GetPointData()->SetScalars(colors);
  ////////////////////
  
  vtkSmartPointer<vtkUnstructuredGrid> dualMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> dualMeshPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> dualMeshCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> dualMeshLines = vtkSmartPointer<vtkCellArray>::New();
  dualMesh->SetPoints(dualMeshPoints);
  
  //dualMesh->SetPolys(dualMeshCells);
  //dualMesh->SetLines(dualMeshLines);
  //dualMesh->SetCells(dualMeshLines);
  
  //dualMesh->GetPointData()->SetScalars(colors);

  std::string inputFilename = argv[1];
  vtkSmartPointer<vtkUnstructuredGrid> mesh = ReadMeshFromVTK(inputFilename);
  
  mesh->BuildLinks();
  mesh->ComputeBounds();

  vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
  
  
  unsigned int numberOfPoints = points->GetNumberOfPoints();
  //std:cout << "nombre de point avant ajout des sommet dual" << numberOfPoints <<endl;
  int timeOut = 5;
  
  for (unsigned int pointCounter = 0; pointCounter < numberOfPoints;
       ++pointCounter)
  {
		double* point = points->GetPoint(pointCounter);
		//std::cout << pointCounter << std::endl;
		vtkSmartPointer<vtkIdList> adjacentCells = vtkSmartPointer<vtkIdList>::New();
		mesh->GetPointCells(pointCounter, adjacentCells);

       	
		if (adjacentCells->GetNumberOfIds() > 1)
        {
		  std::vector<vtkIdType> neighbouringVertices;
          // This loop extracts the neighbourhood of the current point.
          for (unsigned int i = 0; i < adjacentCells->GetNumberOfIds(); ++i)
          {
            vtkSmartPointer<vtkIdList> vertices = vtkSmartPointer<vtkIdList>::New();
            //if (pointer_neightbors_vertices != nullptr)
			//	pointer_neightbors_verices =
            mesh->GetCellPoints(adjacentCells->GetId(i), vertices);
			//vectorNeighborsVertice.push_back(vertices);
			
            for (unsigned int k = 0; k < vertices->GetNumberOfIds(); ++k)
            {
              if (vertices->GetId(k) != pointCounter)
                neighbouringVertices.push_back(vertices->GetId(k));
                
            }
           }
		
		  
		 std::sort(neighbouringVertices.begin(), neighbouringVertices.end(), std::greater<int>());
		 auto last  = std::unique(neighbouringVertices.begin(), neighbouringVertices.end());
		 neighbouringVertices.erase(last, neighbouringVertices.end());
		 //std::cout << "Les sommets voisins du " << pointCounter << " sont ";
		 for (unsigned int i = 0;i < neighbouringVertices.size();++i)
		 {
			//std::cout << " " << neighbouringVertices[i];
								
			//mesh->IsEdge(pointCounter,vertices->GetId(k));
			//====== test ========
			if(pointCounter==timeOut)
			{
				vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();
				vtkIdType p1 = pointCounter;
				vtkIdType p2 = neighbouringVertices[i];
				line->InsertNextId(p1);
				line->InsertNextId(p2);
				dualMesh->InsertNextCell(VTK_LINE, line);
			}
		 }
		 
		 dualMesh->SetPoints(points);
		 //std::cout << std::endl;
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
  
  int count = 0;
  
  for (cellCounter = 0, primalCells->InitTraversal();
       primalCells->GetNextCell(idListPoints) ;
       ++cellCounter)
  {
		//std::cout << "Les ids des points de la cellule " << cellCounter << " sont ";
		std::vector<vtkIdType> CellPointsId;
		
		//======= set dual point  ========
		
		findDualPoints(count,cellCounter,mesh,dualMeshPoints,idListPoints);
		
		//======= sel dual line ======
		
		setDuaLine(count, cellCounter, idListPoints, mesh, dualMesh);
		
		
		count ++; 
		
		//dualMeshCells->InsertNextCell(... 
		
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

  std::stringstream ssDual2;
  ssDual2 << "ProcessedDualMesh"; 
  WriteMeshToVTK(dualMesh, ssDual2.str().c_str());
  
  return EXIT_SUCCESS;
}

void findDualPoints(int count, vtkIdType cellCounter, vtkSmartPointer<vtkUnstructuredGrid> &mesh, vtkSmartPointer<vtkPoints> &dualMeshPoints,
	 vtkSmartPointer<vtkIdList> &idListPoints)
{
	 	double p[3];
		double x=0;
		double y=0;
		double z=0;
		
		for (unsigned int i = 0; i < idListPoints->GetNumberOfIds(); ++i)
		{
			//std::cout << " " << idListPoints->GetId(i);
			vtkIdType v = idListPoints->GetId(i);
			mesh->GetPoint(v,p);
			x += p[0];
			y += p[1];
			z += p[2];
		}
		x /= (double)idListPoints->GetNumberOfIds();
		y /= (double)idListPoints->GetNumberOfIds();
		z /= (double)idListPoints->GetNumberOfIds();
		
		if(count < 10) //affiche les 10 premier points duals
		{
			std::cout << "point " << cellCounter <<" : ";		
			std::cout << x << " ";
			std::cout << y << " ";
			std::cout << z;
			std::cout << std::endl;
		}
		
		dualMeshPoints->InsertNextPoint(x,y,z);

}

void setDuaLine(int count,vtkIdType cellCounter, vtkSmartPointer<vtkIdList> idListPoints, vtkSmartPointer<vtkUnstructuredGrid> &mesh, 
	 vtkSmartPointer<vtkUnstructuredGrid> &dualMesh)
{

	//les cellules voisines d'une cellule primale correspond aux points voisins d'un voisin Dual
	
	//vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
	
	//On stocke les identifiant des cellules voisines de chaque cellule (correspondant à l'indice de nos sommets duals)
	//mesh->GetCellNeighbors(cellCounter,idListPoints, neighborCellIds);
	
	//std::cout << neighborCellIds->GetNumberOfIds() << " ";
	vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();
	vtkIdType p1 = cellCounter;
	line->InsertNextId(p1);
	
	//On parcourt chaque voisin
	for(int i = 0 ; i <idListPoints->GetNumberOfIds();i++)
	{
	
		vtkIdType p2 = idListPoints->GetId(i);

		line->InsertNextId(p2);
		
		//dualMesh->InsertNextCell(VTK_LINE, line);
		
		line->DeleteId(p2);
	}

	
}

void getEdgeCells (vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId, vtkIdList * cellIdsNeighborsFromEdge /*vtkSmartPointer<vtkCellArray> & polyMeshCells*/) {
    /**
        * Get the neightbors cells at an edge
        */

    /** -- Notes Doc
        *    void vtkPolyData::GetCellPoints    (    vtkIdType     cellId,
        *        vtkIdList *     ptIds
        *    )
        *
        *    void vtkPolyData::GetCellEdgeNeighbors    (    vtkIdType     cellId,
        *        vtkIdType     p1,
        *        vtkIdType     p2,
        *        vtkIdList *     cellIds
        *    )
        *    Notes Doc -- */

    // All points ids from a cell
    vtkSmartPointer<vtkIdList> ptIdsFromCell = vtkSmartPointer<vtkIdList>::New();
    _mesh->GetCellPoints(cellId, ptIdsFromCell);

    // two !first! points ids forms a segment of current cell
    vtkSmartPointer<vtkIdList> edgeByTwoPtIds = vtkSmartPointer<vtkIdList>::New();
    edgeByTwoPtIds->InsertNextId(ptIdsFromCell->GetId(0));
    edgeByTwoPtIds->InsertNextId(ptIdsFromCell->GetId(1));

    //vtkSmartPointer<vtkIdList> cellIdsNeighborsFromEdge = vtkSmartPointer<vtkIdList>::New();
    _mesh->GetCellNeighbors(cellId,  edgeByTwoPtIds,  cellIdsNeighborsFromEdge);
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

void WriteMeshToPolyVTK(vtkSmartPointer<vtkPolyData> polyData,
			std::string filename)
{
  //Write the mesh.
  vtkSmartPointer<vtkPolyDataWriter> writer =  
    vtkSmartPointer<vtkPolyDataWriter>::New();
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



 
