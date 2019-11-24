module moduleGrid
   import moduleMath #sina, do not use include here, since you have already included the module in Main.jl
   import moduleMaterialPoint

	fTime = 0.0

	function indexDoubleToSingle(i::Int, j::Int, nRows::Int, nColumns::Int)
		index = nColumns*(i-1) + j

		if(index > nRows*nColumns)
			@printf("Index out of bounds")
		end

		return(Int(index))
	end

   type mpmGridPoint
		bFixed_x::Bool
		bFixed_y::Bool
      fMass::Real
      v2Position::moduleMath.Vector2D
      v2Momentum::moduleMath.Vector2D
		v2Force::moduleMath.Vector2D

      function mpmGridPoint()
         new(false, false, 0.0, moduleMath.Vector2D(0.0, 0.0), moduleMath.Vector2D(0.0, 0.0), moduleMath.Vector2D(0.0, 0.0))
      end
      function mpmGridPoint(bF_x::Bool, bF_y::Bool, fM::Real, v2P::moduleMath.Vector2D, v2M::moduleMath.Vector2D, v2F::moduleMath.Vector2D)
         new(bF_x, bF_y, fM, moduleMath.Vector2D(v2P), moduleMath.Vector2D(v2M), moduleMath.Vector2D(v2F))
      end
   end

   type mpmGrid   #grid container
      fLength_Grid_x::Real
      fLength_Grid_y::Real
      iNodes_x::Int
      iNodes_y::Int
      iNodes::Int
      fLength_Cell_x::Real
      fLength_Cell_y::Real

      GridPoints::Array{mpmGridPoint}

      function mpmGrid(fGL_x, fGL_y, iN_x, iN_y)
         fCL_x = fGL_x / (iN_x - 1)
         fCL_y = fGL_y / (iN_y - 1)

         N = Array{mpmGridPoint}(iN_x * iN_y)
         for i in 1:1:iN_y   #creates row major grid nodes
            for j in 1:1:iN_x
               x = (j-1) * fCL_x
               y = (i-1) * fCL_y
					index = indexDoubleToSingle(i, j, iN_y, iN_x)
               # index = iN_x*(i-1) + j

					bFixed_x = false
					bFixed_y = false
					if(i == 1 || i == iN_y)
						bFixed_y = true
					end
					if(j == 1 || j == iN_x)
						bFixed_x = true
					end
					if(i == iN_y-1)
						bFixed_y = true
					end
					fMass = 0.0
               v2Position = moduleMath.Vector2D(x, y)
               v2Momentum = moduleMath.Vector2D(0.0, 0.0)
					v2Force = moduleMath.Vector2D(0.0, 0.0)
					N[index] = mpmGridPoint(bFixed_x, bFixed_y, fMass, v2Position, v2Momentum, v2Force)
            end
         end

         new(fGL_x, fGL_y, iN_x, iN_y, iN_x*iN_y, fCL_x, fCL_y, N)
      end
   end

   function isAdjacentGridPoint(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::mpmGridPoint, thisGrid::mpmGrid)
      fdx = abs(thisMaterialPoint.v2Position.fx - thisGridPoint.v2Position.fx)
      fdy = abs(thisMaterialPoint.v2Position.fy - thisGridPoint.v2Position.fy)

      if(fdx < thisGrid.fLength_Cell_x && fdy < thisGrid.fLength_Cell_y)
         return(true)
      else
         return(false)
      end
   end

	# function getAdjacentGridPoints(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGrid::mpmGrid)
	# 	#sina, be careful, this does not work if the material point is on the grid edge
	# 	#sina, this is not optimum for the CPDI method
	# 	thisAdjacentGridPoints = Array{Int}(0)
	#
	# 	fParticle_Length_x = thisMaterialPoint.v2Length.fx * thisMaterialPoint.mDeformationGradient[1,1]
	# 	fParticle_Length_y = thisMaterialPoint.v2Length.fy * thisMaterialPoint.mDeformationGradient[2,2]
	#
	# 	iBottomLeft_j	= Int(trunc(floor((thisMaterialPoint.v2Position.fx - 0.5*fParticle_Length_x) / thisGrid.fLength_Cell_x) + 1.0))
	# 	iBottomLeft_i	= Int(trunc(floor((thisMaterialPoint.v2Position.fy - 0.5*fParticle_Length_y) / thisGrid.fLength_Cell_y) + 1.0))
	#
	# 	iTopRight_j		= Int(trunc(ceil((thisMaterialPoint.v2Position.fx + 0.5*fParticle_Length_x) / thisGrid.fLength_Cell_x) + 1.0))
	# 	iTopRight_i		= Int(trunc(ceil((thisMaterialPoint.v2Position.fy + 0.5*fParticle_Length_y) / thisGrid.fLength_Cell_y) + 1.0))
	#
	# 	for i in iBottomLeft_i:1:iTopRight_i
	# 		for j in iBottomLeft_j:1:iTopRight_j
	# 			iIndex = indexDoubleToSingle(Int(i), Int(j), thisGrid.iNodes_y, thisGrid.iNodes_x)
	#
	# 			push!(thisAdjacentGridPoints, iIndex)
	# 		end
	# 	end
	#
	# 	return(thisAdjacentGridPoints)
	# end

	function getAdjacentGridPoints(v2Coordinate::moduleMath.Vector2D, thisGrid::mpmGrid, fTime::Real)
		thisAdjacentGridPoints = Array{Int}(0)

		fLength_Cell_x = thisGrid.fLength_Cell_x
		fLength_Cell_y = thisGrid.fLength_Cell_y

		iBottomLeft_j	= Int(floor(v2Coordinate.fx / fLength_Cell_x) + 1.0)
		iBottomLeft_i	= Int(floor(v2Coordinate.fy / fLength_Cell_y) + 1.0)
		for i in iBottomLeft_i:1:iBottomLeft_i+1
			for j in iBottomLeft_j:1:iBottomLeft_j+1
				iIndex = indexDoubleToSingle(Int(i), Int(j), thisGrid.iNodes_y, thisGrid.iNodes_x)

				push!(thisAdjacentGridPoints, iIndex)
			end
		end

		return(unique(thisAdjacentGridPoints))
	end

	function getAdjacentGridPoints_GIMP(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGrid::mpmGrid, fTime::Real)
		#sina, be careful, this does not work if the material point is on the grid edge
		#sina, this is not optimum for the CPDI method
		thisAdjacentGridPoints = Array{Int}(0)

		fP_x = thisMaterialPoint.v2Position.fx
		fP_y = thisMaterialPoint.v2Position.fy

		mR1 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1
		mR2 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2

		fLength_Cell_x = thisGrid.fLength_Cell_x
		fLength_Cell_y = thisGrid.fLength_Cell_y

		iBottomLeft_j	= Int(floor((fP_x-mR1[1]-mR2[1]) / fLength_Cell_x) + 1.0)
		iBottomLeft_i	= Int(floor((fP_y-mR1[2]-mR2[2]) / fLength_Cell_y) + 1.0)
		iTopRight_j		= Int(ceil((fP_x+mR1[1]+mR2[1]) / fLength_Cell_x) + 1.0)
		iTopRight_i		= Int(ceil((fP_y+mR1[2]+mR2[2]) / fLength_Cell_y) + 1.0)

		for i in iBottomLeft_i:1:iTopRight_i
			for j in iBottomLeft_j:1:iTopRight_j
				iIndex = indexDoubleToSingle(Int(i), Int(j), thisGrid.iNodes_y, thisGrid.iNodes_x)

				push!(thisAdjacentGridPoints, iIndex)
			end
		end

		println("here1")
		println(" mR1: ", mR1, " mR2: ", mR2)

		iBottomLeft_j	= Int(floor((fP_x-mR1[1]-mR2[1]) / fLength_Cell_x) + 1.0)
		iBottomLeft_i	= Int(floor((fP_y-mR1[2]-mR2[2]) / fLength_Cell_y) + 1.0)

		for i in iBottomLeft_i:1:iBottomLeft_i+1
			for j in iBottomLeft_j:1:iBottomLeft_j+1
				iIndex = indexDoubleToSingle(Int(i), Int(j), thisGrid.iNodes_y, thisGrid.iNodes_x)

				push!(thisAdjacentGridPoints, iIndex)
			end
		end

		return(unique(thisAdjacentGridPoints))
	end

	function getAdjacentGridPoints_CPDI(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGrid::mpmGrid, fTime::Real)
		#sina, be careful, this does not work if the material point is on the grid edge
		#sina, this is not optimum for the CPDI method
		thisAdjacentGridPoints = Array{Int}(0)

		fP_x = thisMaterialPoint.v2Position.fx
		fP_y = thisMaterialPoint.v2Position.fy

		mR1 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1
		mR2 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2

		fLength_Cell_x = thisGrid.fLength_Cell_x
		fLength_Cell_y = thisGrid.fLength_Cell_y

		v2Coordinate = moduleMath.Vector2D(0.0, 0.0)

		v2Coordinate.fx = fP_x-mR1[1]-mR2[1]
		v2Coordinate.fy = fP_y-mR1[2]-mR2[2]
		theseIndex = getAdjacentGridPoints(v2Coordinate, thisGrid, fTime)
		append!(thisAdjacentGridPoints, theseIndex)# sina, push is for single elements, append for arrays

		v2Coordinate.fx = fP_x+mR1[1]+mR2[1]
		v2Coordinate.fy = fP_y-mR1[2]-mR2[2]
		theseIndex = getAdjacentGridPoints(v2Coordinate, thisGrid, fTime)
		append!(thisAdjacentGridPoints, theseIndex)# sina, push is for single elements, append for arrays

		v2Coordinate.fx = fP_x+mR1[1]+mR2[1]
		v2Coordinate.fy = fP_y+mR1[2]+mR2[2]
		theseIndex = getAdjacentGridPoints(v2Coordinate, thisGrid, fTime)
		append!(thisAdjacentGridPoints, theseIndex)# sina, push is for single elements, append for arrays

		v2Coordinate.fx = fP_x-mR1[1]-mR2[1]
		v2Coordinate.fy = fP_y+mR1[2]+mR2[2]
		theseIndex = getAdjacentGridPoints(v2Coordinate, thisGrid, fTime)
		append!(thisAdjacentGridPoints, theseIndex)# sina, push is for single elements, append for arrays

		return(unique(thisAdjacentGridPoints))
	end
end
