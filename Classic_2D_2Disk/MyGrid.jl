module moduleGrid
   import moduleMaterialPoint #sina, do not use include here, since you have already included the module in Main.jl

	fTime = 0.0

	function index2DTo1D(i::Int64, j::Int64, nColumns::Int64, nRows::Int64)
		index = nColumns*(j-1) + i

		if(index > nRows*nColumns || index < 1)
			@printf("Index out of bounds: i, j: %d, %d \n", i, j)
		end

		return(Int64(index))
	end

   type mpmGridPoint
	  v2Fixed         :: Array{Bool}
      fMass           :: Float64
      v2Position      :: Array{Float64}
      v2Momentum      :: Array{Float64}
	  v2Force         :: Array{Float64}

      function mpmGridPoint()
         new([false; false],
				0.0,
				zeros(2),
				zeros(2),
				zeros(2))
      end
   end

   #grid container
   type mpmGrid   
      v2Length_Grid        :: Array{Float64}       # length in x and y dirs
      v2Nodes              :: Array{Int64}  
      iNodes               :: Int64                # number of grid nodes
      v2Length_Cell        :: Array{Float64}       # size of each cell/element
      v2Length_CellI       :: Array{Float64}       # inverse of size of each cell/element
      GridPoints           :: Array{mpmGridPoint}  # array of all grid points

      # constructor, GL_x is length of the grid in x dir
      # iN_x: number of nodes in x dir 
      function mpmGrid(fGL_x, fGL_y, iN_x, iN_y)

			v2CL   = zeros(2)
			v2CLI  = zeros(2)
         v2CL[1]   = fGL_x / Float64(iN_x - 1.0)
         v2CL[2]   = fGL_y / Float64(iN_y - 1.0)
         v2CLI[1]  = 1.0 / v2CL[1]
         v2CLI[2]  = 1.0 / v2CL[2]

         thisGridPoint = Array{mpmGridPoint}(iN_x * iN_y)
         for j=1:iN_y
            for i=1:iN_x
               x = (i-1) * v2CL[1]
               y = (j-1) * v2CL[2]
			   index = index2DTo1D(i, j, iN_x, iN_y)

				thisGridPoint[index] = mpmGridPoint() # starts with every member equal to zero
				thisGridPoint[index].v2Fixed = [false; false]
				thisGridPoint[index].v2Position = [x; y]
			end
         end

         new([fGL_x; fGL_y], [iN_x; iN_y], iN_x*iN_y, v2CL, v2CLI, thisGridPoint)
      end
   end

	function getAdjacentGridPoints(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_2D_Classic, thisGrid::mpmGrid)
#		thisAdjacentGridPoints = Array{Int64}(0)

		v2Coordinate   = thisMaterialPoint.v2Centroid
		fLength_Cell_x = thisGrid.v2Length_CellI[1]
		fLength_Cell_y = thisGrid.v2Length_CellI[2]

		iBottomLeft_i  = (floor(v2Coordinate[1] * fLength_Cell_x) + 1.)
		iBottomLeft_j  = (floor(v2Coordinate[2] * fLength_Cell_y) + 1.)

		if(iBottomLeft_j < 1 || iBottomLeft_j > thisGrid.v2Nodes[2])
			@printf("Index out of bounds: j: %d \n", iBottomLeft_j)
			@printf("v2Coordinate[2]: %e \n", v2Coordinate[2])
		end
#
#		for i = iBottomLeft_i:1:iBottomLeft_i+1
#			for j = iBottomLeft_j:1:iBottomLeft_j+1
#				iIndex = index2DTo1D(Int64(i), Int64(j), thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])
#
#				push!(thisAdjacentGridPoints, iIndex)
#			end
#		end
				
        iIndex = index2DTo1D(Int64(iBottomLeft_i), Int64(iBottomLeft_j),   thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])
        jIndex = index2DTo1D(Int64(iBottomLeft_i), Int64(iBottomLeft_j)+1, thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])
        thisAdjacentGridPoints = [iIndex, jIndex, iIndex+1, jIndex+1 ]
		return((thisAdjacentGridPoints))
	end
end
