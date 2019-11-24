module moduleGrid
    import moduleMaterialPoint #sina, do not use include here, since you have already included the module in Main.jl

    export point2ElemIndexIJ, getAdjacentGridPoints

	function index2DTo1D(i::Int64, j::Int64, nColumns::Int64, nRows::Int64)
		index = nColumns*(j-1) + i

		if(index > nRows*nColumns || index < 1)
			@printf("Index out of bounds: i, j: %d, %d \n", i, j)
		end

		return(Int64(index))
	end

   
   ###################################################
   #  mpmGridPoint
   ###################################################

   type mpmGridPoint
	 v2Fixed     :: Array{Bool}               # fix in two directions
     fMass       :: Float64                   # mass of a g rid point
     v2Position  :: Array{Float64}            # coordinates
     v2Velocity  :: Array{Float64}            # velocities
	 v2Momentum  :: Array{Float64}            # momentum
	 v2Force     :: Array{Float64}            # forces

     # constructor to initialise the above data of a grid point
     function mpmGridPoint()
          new([false; false],
	 			0.0,
	 			zeros(2),
	 			zeros(2),
	 			zeros(2),
	 			zeros(2))
      end
   end

   ###################################################
   # mpmGrid, contains an array of mpmGridPoints
   ###################################################

   type mpmGrid   
      v2Length_Grid     :: Array{Float64}    # dimensions of the grid in 2 directions
      v2Nodes           :: Array{Int64}      # number of grid nodes along x and y direction
      iNodes            :: Int64             # total number of grid nodes
      v2Length_Cell     :: Array{Float64}    # Deltax and Deltay, lengths of one cell/element
      GridPoints        :: Array{mpmGridPoint}
  
      elem2MP           :: Dict{Int64,Array{Int64,1}} # cell => MP, for nonlocal damage models

      function mpmGrid(fGL_x, fGL_y, iN_x, iN_y)

			v2CL = zeros(2)
         v2CL[1] = fGL_x / Float64(iN_x - 1.0)
         v2CL[2] = fGL_y / Float64(iN_y - 1.0)

         thisGridPoint = Array{mpmGridPoint}(iN_x * iN_y)
         for j=1:1:iN_y
            for i=1:1:iN_x
              x = (i-1) * v2CL[1]
              y = (j-1) * v2CL[2]
			  index = index2DTo1D(i, j, iN_x, iN_y)

			  bFixed_x = false
			  bFixed_y = false
			  thisGridPoint[index] = mpmGridPoint() # starts with every member equal to zero
			  thisGridPoint[index].v2Fixed = [bFixed_x; bFixed_y]
			  thisGridPoint[index].v2Position = [x; y]
			end
         end
         new([fGL_x; fGL_y], [iN_x; iN_y], iN_x*iN_y, v2CL, thisGridPoint)
      end
   end

	function getAdjacentGridPoints(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_2D_Classic, thisGrid::mpmGrid)
		thisAdjacentGridPoints = Array{Int64}(0)

		v2Coordinate = thisMaterialPoint.v2Centroid

		fLength_Cell_x = thisGrid.v2Length_Cell[1]
		fLength_Cell_y = thisGrid.v2Length_Cell[2]

		iBottomLeft_i	= (floor(v2Coordinate[1] / fLength_Cell_x) + 1.0)
		iBottomLeft_j	= (floor(v2Coordinate[2] / fLength_Cell_y) + 1.0)

		if(iBottomLeft_j < 1 || iBottomLeft_j > thisGrid.v2Nodes[2])
			@printf("Index out of bounds: j: %d \n", iBottomLeft_j)
			@printf("v2Coordinate[2]: %e \n", v2Coordinate[2])
		end

		for i = iBottomLeft_i:1:iBottomLeft_i+1
			for j = iBottomLeft_j:1:iBottomLeft_j+1
				iIndex = index2DTo1D(Int64(i), Int64(j), thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])

				push!(thisAdjacentGridPoints, iIndex)
			end
		end

		return((thisAdjacentGridPoints))
	end
    
    function point2ElemIndexIJ(x::Array{Float64},grid::mpmGrid)
      deltax = grid.v2Length_Cell[1]
      deltay = grid.v2Length_Cell[2]

      indexI = Int64((floor(x[1] / deltax) + 1.0))
      indexJ = Int64((floor(x[2] / deltay) + 1.0))
      
      if indexI <= 0
        indexI = 1
      end
      
      if indexI > grid.v2Nodes[1]
        indexI = grid.v2Nodes[1]
      end
      
      if indexJ <= 0
        indexJ = 1
      end
      
      if indexJ > grid.v2Nodes[2]
        indexJ = grid.v2Nodes[2]
      end

      return (indexI, indexJ)

	end
	
    # find the grid points of elements that contain particles 
    # return Array{Int64,1}
    function getGridPointsForParticles(particles::Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}, thisGrid::mpmGrid)
      set = IntSet() 
      for p=1:length(particles)
        thisParticle = particles[p]
        nearPoints   = getAdjacentGridPoints(thisParticle, thisGrid)
        for np=1:length(nearPoints)
          push!(set, nearPoints[np])
        end
      end

      return collect(set)
    end

    # update the material points inside every elements/cells
    function update( mesh::mpmGrid, materialPoints::Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic})
      elem2MPTemp = Dict{Int64,Array{Int64,1}}()
      for iIndex_MP in 1:length(materialPoints)
        thisMP      = materialPoints[iIndex_MP]
        x           = thisMP.v2Centroid[1]
        y           = thisMP.v2Centroid[2]
        el          = floor(x/mesh.v2Length_Cell[1]) + 1 + (mesh.v2Nodes[1]-1)*floor(y/mesh.v2Length_Cell[2])
        if haskey(elem2MPTemp,el) == false
          elem2MPTemp[el] = [iIndex_MP]
        else
          push!(elem2MPTemp[el],iIndex_MP)
        end
      end
      mesh.elem2MP = elem2MPTemp
    end

    # find the ids of all particles fall within a circle around a given particle
    # @input:
    # xp             = coordinates of the particle of interest
    # radius         = nonlocal radius
    # grid           = the background grid
    # materialPOints = material points 
    # @output: an array of integers which are the indices (in materialPoints array) of sought for particles 
    function getNeibourParticleIDs(xp::Array{Float64}, radius::Float64, grid::mpmGrid, 
                                  materialPoints::Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic})
      xmin = xp[1] - radius
      xmax = xp[1] + radius
      ymin = xp[2] - radius
      ymax = xp[2] + radius

      minI, minJ = point2ElemIndexIJ([xmin;ymin], grid)
      maxI, maxJ = point2ElemIndexIJ([xmax;ymax], grid)

      result     = Array{Int64}(0)
      noElemX    =  grid.v2Nodes[1] - 1
      
      for i = minI:maxI
        for j = minJ:maxJ
          iE               = i + noElemX * (j-1)
          if haskey(grid.elem2MP,iE) == false continue end
          materialPointIds = grid.elem2MP[iE]

          for iP=1:length(materialPointIds)
            mpID   = materialPointIds[iP]
            thisMP = materialPoints[mpID]
            xip    = thisMP.v2Centroid
            if ( (xp[1]-xip[1])^2 + (xp[2]-xip[2])^2 < radius^2)
              push!(result,mpID)
            end
          end
        end
      end

      return result

    end # end of function

end
