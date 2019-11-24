module moduleGrid
   import moduleMaterialPoint #sina, do not use include here, since you have already included the module in Main.jl

	fTime = 0.0

	function index3DTo1D(i::Int, j::Int, k::Int, nColumns::Int, nRows::Int, nLayers::Int)
		index = nColumns*nRows*(k-1) + nColumns*(j-1) + i

		if(index > nRows*nColumns*nLayers)
			@printf("Index out of bounds")
		end

		return(Int(index))
	end

   type mpmGridPoint
		v3Fixed::Array{Bool}
      fMass::Real
      v3Position::Array{Float64}
      v3Momentum::Array{Float64}
		v3Force::Array{Float64}

      function mpmGridPoint()
         new([false; false; false],
				0.0,
				zeros(3),
				zeros(3),
				zeros(3))
      end
   end

   type mpmGrid   #grid container
      v3Length_Grid::Array{Float64}
      v3Nodes::Array{Int}
      iNodes::Int
      v3Length_Cell::Array{Float64}

      GridPoints::Array{mpmGridPoint}

      function mpmGrid(fGL_x, fGL_y, fGL_z, iN_x, iN_y, iN_z)

			v3CL = zeros(3)
         v3CL[1] = fGL_x / (iN_x - 1)
         v3CL[2] = fGL_y / (iN_y - 1)
			v3CL[3] = fGL_z / (iN_z - 1)

         thisGridPoint = Array{mpmGridPoint}(iN_x * iN_y * iN_z)
			for k=1:1:iN_z
	         for j=1:1:iN_y
	            for i=1:1:iN_x
	               x = (i-1) * v3CL[1]
	               y = (j-1) * v3CL[2]
						z = (k-1) * v3CL[3]
						index = index3DTo1D(i, j, k, iN_x, iN_y, iN_z)

						bFixed_x = false
						bFixed_y = false
						bFixed_z = false
						if(x < 0.51)
							bFixed_x = true
							bFixed_y = true
							bFixed_z = true
						end

						thisGridPoint[index] = mpmGridPoint() # starts with every member equal to zero

						thisGridPoint[index].v3Fixed = [bFixed_x; bFixed_y; bFixed_z]

						thisGridPoint[index].v3Position = [x; y; z]
					end
            end
         end

         new([fGL_x; fGL_y; fGL_z], [iN_x; iN_y; iN_z], iN_x*iN_y*iN_z, v3CL, thisGridPoint)
      end
   end

   function isAdjacentGridPoint(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_Tet4, thisGridPoint::mpmGridPoint, thisGrid::mpmGrid)
      fdx = abs(thisMaterialPoint.v3Position[1,1] - thisGridPoint.v3Position[1,1])
      fdy = abs(thisMaterialPoint.v3Position[1,2] - thisGridPoint.v3Position[1,2])
		fdz = abs(thisMaterialPoint.v3Position[1,3] - thisGridPoint.v3Position[1,3])

      if(fdx < thisGrid.v3Length_Cell[1] && fdy < thisGrid.v3Length_Cell[2] && fdz < thisGrid.v3Length_Cell[3])
         return(true)
      else
         return(false)
      end
   end

	function getAdjacentGridPoints(v3Coordinate::Array{Float64}, thisGrid::mpmGrid)
		thisAdjacentGridPoints = Array{Int}(0)

		v3Length_Cell = thisGrid.v3Length_Cell

		iBottomLeft_i	= Int(floor(v3Coordinate[1] / v3Length_Cell[1]) + 1.0)
		iBottomLeft_j	= Int(floor(v3Coordinate[2] / v3Length_Cell[2]) + 1.0)
		iBottomLeft_k	= Int(floor(v3Coordinate[3] / v3Length_Cell[3]) + 1.0)

		for i in iBottomLeft_i:1:iBottomLeft_i+1
			for j in iBottomLeft_j:1:iBottomLeft_j+1
				for k = iBottomLeft_k:1:iBottomLeft_k+1
					iIndex = index3DTo1D(Int(i), Int(j), Int(k), thisGrid.v3Nodes[1], thisGrid.v3Nodes[2], thisGrid.v3Nodes[3])

					push!(thisAdjacentGridPoints, iIndex)
				end
			end
		end

		return(unique(thisAdjacentGridPoints))
	end

	function getAdjacentGridPoints_CPDI(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_Tet4, thisGrid::mpmGrid)
		#sina, be careful, this does not work if the material point is on the grid edge
		#sina, this is not optimum for the CPDI method
		thisAdjacentGridPoints = Array{Int}(0)

		for indexCorner = 1:1:size(thisMaterialPoint.v3Corner,2)
			thisIndex = getAdjacentGridPoints(thisMaterialPoint.v3Corner[:,indexCorner], thisGrid)
			append!(thisAdjacentGridPoints, thisIndex)# sina, push is for single elements, append for arrays
		end

		return(unique(thisAdjacentGridPoints))
	end
end
