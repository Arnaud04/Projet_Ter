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

void getEdgeCells (vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId, vtkCellArray * cells);

bool compareCellsByFaces(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId1, vtkIdType cellId2);

void WriteMeshToPolyVTK(vtkSmartPointer<vtkPolyData> polyData, std::string filename);



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

void findDualPoints(int count, vtkIdType cellCounter, vtkSmartPointer<vtkUnstructuredGrid> &mesh, vtkSmartPointer<vtkPoints> &dualMeshPoints,
	 vtkSmartPointer<vtkIdList> &idListPoints)
{
	 	double p[3];
		double x=0;
		double y=0;
		double z=0;

		for (unsigned int i = 0; i < idListPoints->GetNumberOfIds(); ++i)
		{
			vtkIdType v = idListPoints->GetId(i);
			mesh->GetPoint(v,p);
			x += p[0];
			y += p[1];
			z += p[2];
		}
		x /= (double)idListPoints->GetNumberOfIds();
		y /= (double)idListPoints->GetNumberOfIds();
		z /= (double)idListPoints->GetNumberOfIds();

		if(count < 10) //print the 10th first duals points
		{
			std::cout << "point " << cellCounter <<" : ";
			std::cout << x << " ";
			std::cout << y << " ";
			std::cout << z;
			std::cout << std::endl;
		}

		dualMeshPoints->InsertNextPoint(x,y,z);

}

void setDuaLine(vtkCellArray cells, vtkSmartPointer<vtkUnstructuredGrid> &mesh,
     vtkSmartPointer<vtkUnstructuredGrid> &dualMesh)
{


    vtkSmartPointer<vtkIdList> idListPoints1 = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> idListPoints2 = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();

    //vtkIdType cellCounter1 = 0;
    //vtkIdType cellCounter2 = 0;

    //vtkSmartPointer <vtkPoint> pointsDual1 = mesh -> GetPoints();

    std::cout<<cells->GetNextCell(idListPoints1)<<" ";

    for (int cellCounter1=0;cellCounter1<8;cellCounter1++/cellCounter1 = 0 ; cells->GetNextCell(idListPoints1) ; cellCounter1 ++/)
    {
        //vtkIdType pointsDual1 = idListPoints1->GetId(i);
        //mesh->GetPoint(v,p);
        std::cout << "cell " << cellCounter1 <<" neigbors of : ";
        for (int cellCounter2=0; cellCounter2<8; cellCounter2++/cellCounter2 = 0 ; cells->GetNextCell(idListPoints2) ; cellCounter2 ++/)
        {
            if(cellCounter2 != cellCounter1)
                if(compareCellsByFaces(mesh, cellCounter1, cellCounter2));
                {

                    std::cout << cellCounter2 << " ";
                    //line->InsertNextId(/mesh->GetId(/cellCounter1/)/);
                    //line->InsertNextId(/mesh->GetId(/cellCounter2/)*/);
                    //dualMesh->InsertNextCell(VTK_LINE, line);
                }
        }
        std::cout<< endl;
    }
	/*for (cellCounter2 = 0, primalCells->InitTraversal();
    	primalCells->GetNextCell(idListPoints2) ;
        ++cellCounter2)
    {
        if(cellCounter2 != cellCounter)
            if(compareCellsByFaces(mesh, cellCounter, cellCounter2));
            {
                line->InsertNextId(cellCounter);
                line->InsertNextId(cellCounter2);
                dualMesh->InsertNextCell(VTK_LINE, line);
            }
    }*/

    /*
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
    */

}

void getEdgeCells (vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId, vtkCellArray * cells) {
	/**
		* Get the neightbors cells at an edge
		*/

  // All points ids from a cell
  vtkSmartPointer<vtkIdList> ptIdsFromCell = vtkSmartPointer<vtkIdList>::New();
  _mesh->GetCellPoints(cellId, ptIdsFromCell);


  // two points ids !defined in tough! forms a segment of current cell
  vtkSmartPointer<vtkIdList> edgeByTwoPtIds = vtkSmartPointer<vtkIdList>::New();
  edgeByTwoPtIds->InsertNextId(5);
  edgeByTwoPtIds->InsertNextId(6);

  vtkSmartPointer<vtkIdList> cellIdsNeighborsFromEdge = vtkSmartPointer<vtkIdList>::New();
  _mesh->GetCellNeighbors(cellId,  edgeByTwoPtIds,  cellIdsNeighborsFromEdge);

  // Stock all the points from Current Cell in vtkCellArray
	cells->InsertNextCell(_mesh->GetCell(cellId));
  for (int counterCurrentCell = 0; counterCurrentCell < cellIdsNeighborsFromEdge->GetNumberOfIds(); counterCurrentCell++) {
		vtkIdType idCurrentCell = cellIdsNeighborsFromEdge->GetId(counterCurrentCell);
		cells->InsertNextCell(_mesh->GetCell(idCurrentCell));
	}

}

bool compareCellsByFaces(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId1, vtkIdType cellId2) {
	/**
		*	compare two cell by their	faces, determine if they are adjacent.
		*/

	// Recovere all points ids from each cell
	vtkSmartPointer<vtkIdList> ptIdsFromCell1 = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> ptIdsFromCell2 = vtkSmartPointer<vtkIdList>::New();
	_mesh->GetCellPoints(cellId1, ptIdsFromCell1);
	_mesh->GetCellPoints(cellId2, ptIdsFromCell2);

	// Course all points ids from each cell and add the points ids common in list
	vtkSmartPointer<vtkIdList> resultPointsIdsCommonCompareFromCells = vtkSmartPointer<vtkIdList>::New();
	for (int counterPointIdFromCell1 = 0; counterPointIdFromCell1 < ptIdsFromCell1->GetNumberOfIds(); counterPointIdFromCell1++) {
		for (int counterPointIdFromCell2 = 0; counterPointIdFromCell2 < ptIdsFromCell2->GetNumberOfIds(); counterPointIdFromCell2++) {
			if (ptIdsFromCell1->GetId(counterPointIdFromCell1) == ptIdsFromCell2->GetId(counterPointIdFromCell2)) {
				resultPointsIdsCommonCompareFromCells->InsertUniqueId(ptIdsFromCell1->GetId(counterPointIdFromCell1));
			}
		}
	}

	// true if face common
	if (resultPointsIdsCommonCompareFromCells->GetNumberOfIds() >= 3)
		return true;

	return false;
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



int main ( int argc, char *argv[] )
{
  if(argc != 2)
  {
    std::cout << "Usage: " << argv[0] << "  Filename(.vtk)"
	      << " <Number of iterations>" << std::endl;
    return EXIT_FAILURE;
  }


  vtkSmartPointer<vtkUnstructuredGrid> dualMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> dualMeshPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> dualMeshCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> dualMeshLines = vtkSmartPointer<vtkCellArray>::New();
  dualMesh->SetPoints(dualMeshPoints);

  //dualMesh->SetPolys(dualMeshCells);
  //dualMesh->SetLines(dualMeshLines);
  //dualMesh->SetCells(dualMeshLines);

  std::string inputFilename = argv[1];
  vtkSmartPointer<vtkUnstructuredGrid> mesh = ReadMeshFromVTK(inputFilename);

  mesh->BuildLinks();
  mesh->ComputeBounds();

  vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

  unsigned int numberOfPoints = points->GetNumberOfPoints();

  int numVertex = 5;
  int numVertex2 = 0;

  for (unsigned int pointCounter = 0; pointCounter < numberOfPoints;
       ++pointCounter)
  {
		double* point = points->GetPoint(pointCounter);

		vtkSmartPointer<vtkIdList> adjacentCells = vtkSmartPointer<vtkIdList>::New();
		mesh->GetPointCells(pointCounter, adjacentCells);


		if (adjacentCells->GetNumberOfIds() > 1)
		{
			std::vector<vtkIdType> neighbouringVertices;
			// This loop extracts the neighbourhood of the current point.
			for (unsigned int i = 0; i < adjacentCells->GetNumberOfIds(); ++i)
			{
				vtkSmartPointer<vtkIdList> vertices = vtkSmartPointer<vtkIdList>::New();

				mesh->GetCellPoints(adjacentCells->GetId(i), vertices);

				for (unsigned int k = 0; k < vertices->GetNumberOfIds(); ++k)
				{
					if (vertices->GetId(k) != pointCounter)
						neighbouringVertices.push_back(vertices->GetId(k));

				}
			}


			std::sort(neighbouringVertices.begin(), neighbouringVertices.end(), std::greater<int>());
			auto last  = std::unique(neighbouringVertices.begin(), neighbouringVertices.end());
			neighbouringVertices.erase(last, neighbouringVertices.end());

			for (unsigned int i = 0;i < neighbouringVertices.size();++i)
			{

				//====== test : Etraire les arêtes voisines d'un sommet i ========
				if(pointCounter == numVertex)
				{
					vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();
					vtkIdType p1 = pointCounter;
					vtkIdType p2 = neighbouringVertices[i];
					line->InsertNextId(p1);
					line->InsertNextId(p2);
					dualMesh->InsertNextCell(VTK_LINE, line);

					//on extrait une arête et on regarde les cellules adjacatentes au deux points de celle ci
					//if(numVertex2 < 1)
					//{

						//vtkIdType pointId,vtkSmartPointer<vtkIdList> &cells,vtkSmartPointer<vtkUnstructuredGrid> &mesh
						//mes


						//mesh->GetPointCells(p2,cells);
						//dualMesh->InsertNextCell(VTK_QUAD, cells);
						//mesh->GetPointCells(p1,cells);
						//dualMesh->InsertNextCell(VTK_TRIANGLE, cells);



						//numVertex2 ++;
					//}

					//dualMesh->Set

				}

			}
			//dualMesh->SetPoints(points);

		}

  }

	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	getEdgeCells(mesh, 0, cells);

	// --------- Set DualMesh --------------
	// TODO : Test
	// Meth 1 :
	//dualMesh->SetCells(VTK_TETRA, cells);

	// Meth 2 : => Functional
	dualMesh->SetPoints(points);
	dualMesh->SetCells(VTK_TETRA, cells);

	// Meth 3 :
	/*for (int counterCurrentCell = 0; counterCurrentCell < cells->GetNumberOfCells(); counterCurrentCell++) {
		vtkSmartPointer<vtkIdList> ptsIdsFromCurrentCell = cells->;
		dualMesh->InsertNextCell(VTK_POLYGON, ptsIdsFromCurrentCell);
	}*/
	// -------------------------------------
	std::cout << "compare cells 0 and 1 by faces" << compareCellsByFaces(mesh, 0, 3) << endl;

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
		std::vector<vtkIdType> CellPoints;

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

  std::stringstream ssDualDual2;
  ssDualDual2 << "ProcessedDualMesh";
  WriteMeshToVTK(dualMesh, ssDualDual2.str().c_str());

  return EXIT_SUCCESS;
}
