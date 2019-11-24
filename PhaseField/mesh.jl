# Functions related to mesh generation
#
# Vinh Phu Nguyen
# (nvinhphu@gmail.com)
#

module FeMesh

import moduleMaterialPoint
using PyPlot

# Data structure to store mesh
type Mesh
  dim        :: Int
  lengthX    :: Float64
  lengthY    :: Float64
  elemCountX :: Int
  elemCountY :: Int
  elemCount  :: Int
  nodeCount  :: Int
  deltaX     :: Float64
  deltaY     :: Float64

  nodes      :: Array{Float64, 2}
  elements   :: Array{Int, 2}

  elem2MP    :: Dict{Int64,Array{Int64,1}}

  function Mesh(dim::Int, lx::Float64, ly::Float64, elemCountX::Int, elemCountY::Int)
    elemCount  = elemCountX * elemCountY
    nodeCount  = (elemCountX+1)*(elemCountY+1)
    nodes      = zeros(2,nodeCount)
    elements   = zeros(4,elemCount)

    dx         = lx / Float64(elemCountX)
    dy         = ly / Float64(elemCountY)

    for j=1:elemCountY+1
      for i=1:elemCountX+1
        index          = (elemCountX+1)*(j-1)+i
        nodes[1,index] = (i-1)*dx 
        nodes[2,index] = (j-1)*dy 
      end
    end

    index = 0
    
    for j=1:elemCountY
      for i=1:elemCountX
        index += 1
        elements[1,index] = (elemCountX+1)*(j-1)+i
        elements[2,index] = (elemCountX+1)*(j-1)+i+1
        elements[3,index] = (elemCountX+1)*(j  )+i+1
        elements[4,index] = (elemCountX+1)*(j  )+i
      end
    end

    new(dim, lx, ly, elemCountX, elemCountY, elemCount, nodeCount, dx, dy, nodes, elements)
  end
end


# Write mesh node and connectivity data to screen
function write_to_file(mesh::Mesh)
  for i = 1:mesh.nodeCount
    @printf("%f, %f\n", mesh.nodes[1,i], mesh.nodes[2,i])
  end

  for i = 1:mesh.elemCount
    @printf("%d %d %d %d\n", mesh.elements[1,i], mesh.elements[2,i],
            mesh.elements[3,i],mesh.elements[4,i])
  end
end

function vtk(mesh::Mesh,phi::Array{Float64},vtuFile::ASCIIString)
  results_vtu        = open(vtuFile, "w")
  numVertexesPerCell = 4
  VTKCellCode        = 9
  # Write headers
  write(results_vtu, "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"   > \n")
  write(results_vtu, "<UnstructuredGrid> \n");
  write(results_vtu, "<Piece  NumberOfPoints=\"$(mesh.nodeCount)\" NumberOfCells=\"$(mesh.elemCount)\"> \n")

  
  # Write point data
  write(results_vtu, "<Points> \n");
  write(results_vtu, "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\" format=\"ascii\" > \n");
  
  for i=1:mesh.nodeCount
    write(results_vtu, "$(mesh.nodes[1,i])  $(mesh.nodes[2,i]) 0.0\n");
  end
  
  write(results_vtu, "</DataArray> \n");
  write(results_vtu, "</Points> \n");
  
  # Print cells
  write(results_vtu, "<Cells> \n");
  write(results_vtu, "<DataArray  type=\"Int32\"  Name=\"connectivity\"  format=\"ascii\"> \n");
  
  for i=1:mesh.elemCount
    writedlm(results_vtu, transpose(mesh.elements[:,i])-1 );
  end
  
  write(results_vtu, "</DataArray> \n");
  write(results_vtu, "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\"> \n");
  
  offset = 0;
  for i=1:mesh.elemCount
      offset = offset + numVertexesPerCell;
      write(results_vtu, "$offset \n");
  end
  
  write(results_vtu, "</DataArray> \n");
  write(results_vtu, "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\"> \n");
  
  for i=1:mesh.elemCount
     write(results_vtu, "$VTKCellCode \n");
  end
  
  write(results_vtu, "</DataArray> \n");
  write(results_vtu, "</Cells> \n");
  
  write(results_vtu, "<PointData  Scalars=\"U\"> \n");
  write(results_vtu, "<DataArray  type=\"Float64\"  Name=\"U\" NumberOfComponents=\"1\"
          format=\"ascii\"> \n");
  for i=1:length(phi)
    write(results_vtu, "$(phi[i]) \n");
  end
  write(results_vtu, "</DataArray> \n");
  write(results_vtu, "</PointData> \n");
  
  
  write(results_vtu, "</Piece> \n");
  write(results_vtu, "</UnstructuredGrid> \n");
  write(results_vtu, "</VTKFile> \n");
  #
  
  close(results_vtu);
end

function update( mesh::Mesh, materialPoints::Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic})
  elem2MPTemp = Dict{Int64,Array{Int64,1}}()
  for iIndex_MP in 1:length(materialPoints)
    thisMP      = materialPoints[iIndex_MP]
    x           = thisMP.v2Centroid[1]
    y           = thisMP.v2Centroid[2]
    el          = floor(x/mesh.deltaX) + 1 + mesh.elemCountX*floor(y/mesh.deltaY)
    if haskey(elem2MPTemp,el) == false
      elem2MPTemp[el] = [iIndex_MP]
    else
      push!(elem2MPTemp[el],iIndex_MP)
    end
  end
  mesh.elem2MP = elem2MPTemp
end

function findAdjacentPoints(xp::Array{Float64}, mesh::Mesh )
    iBottomLeft_i	= Int64( (floor(xp[1] / mesh.deltaX) + 1.0) )
    iBottomLeft_j	= Int64( (floor(xp[2] / mesh.deltaY) + 1.0) )

    corner1         = (mesh.elemCountX+1)*(iBottomLeft_j-1) + iBottomLeft_i
    corner4         = (mesh.elemCountX+1)*(iBottomLeft_j  ) + iBottomLeft_i


    thisAdjacentGridPoints = [corner1, corner1+1, corner4, corner4+1]

	return((thisAdjacentGridPoints))
end

function plot_mesh(mesh::Mesh)
  ord = [1,2,3,4,1]
  xpt = zeros(5)
  ypt = zeros(5)
  zpt = zeros(5)
  for e=1:mesh.elemCount
    for n=1:5
      xpt[n] = mesh.nodes[1,mesh.elements[ord[n],e]]
      ypt[n] = mesh.nodes[2,mesh.elements[ord[n],e]]
    end
    plot(xpt,ypt, color="blue", linewidth=2.0, linestyle="-")
    hold(true)
  end
end

end
