module moduleMaterialPoint
	type mpmMaterialPoint_Tet4   #material point container
      fMass::Float64
		fVolumeInitial::Float64
		fVolume::Float64

		fElasticModulus::Float64
		fPoissonRatio::Float64

		v3Centroid::Array{Float64}
		v3CentroidIncrement::Array{Float64}
		v3Velocity::Array{Float64}
      v3Momentum::Array{Float64}
		v3ExternalForce::Array{Float64}
		v3Restraint::Array{Float64}	# 0.0=no restraint, 1.0=fully restrained

		v3Corner::Array{Float64} # corner position
		v3Corner_Increment::Array{Float64}

		m33DeformationGradient::Array{Float64}
		m33DeformationGradientIncrement::Array{Float64}

		v6Strain::Array{Float64} # xx, yy, zz, xy, yz, zx
		v6StrainIncrement::Array{Float64}

		v6Stress::Array{Float64}
		v6StressIncrement::Array{Float64}

      function mpmMaterialPoint_Tet4()
         new(1.0, 1.0, 1.0, # Mass, initial volume, volume
				1.0, 0.3, # elastic modulus, poisson ratio
				zeros(3), # centroid position
				zeros(3), # centroid position increment
				zeros(3), # velocity
				zeros(3), # momentum
				zeros(3), # external force
				zeros(3), # restraint
				zeros(3,4), # array of corner positions, 3d coordinates
				zeros(3,4), # array of corner position incretement
				eye(3,3), # deformation gradient
				eye(3,3), # deformation gradient increment
				zeros(6), # strain
				zeros(6), # strain increment
				zeros(6), # stress
				zeros(6) # stress increment
				)
      end
   end

	function getVolume(thisMaterialPoint::mpmMaterialPoint_Tet4)
		x = vec(thisMaterialPoint.v3Corner[1,:])
		y = vec(thisMaterialPoint.v3Corner[2,:])
		z = vec(thisMaterialPoint.v3Corner[3,:])

		x21 = (x[2] - x[1]); x23 = (x[2] - x[3]); x24 = (x[2] - x[4])
		x12 = (x[1] - x[2]); x13 = (x[1] - x[3]); x14 = (x[1] - x[4])
		x31 = (x[3] - x[1]); x32 = (x[3] - x[2]); x34 = (x[3] - x[4])
		x41 = (x[4] - x[1]); x42 = (x[4] - x[2]); x43 = (x[4] - x[3])

		y12 = (y[1] - y[2]); y13 = (y[1] - y[3]); y14 = (y[1] - y[4])
		y21 = (y[2] - y[1]); y23 = (y[2] - y[3]); y24 = (y[2] - y[4])
		y31 = (y[3] - y[1]); y32 = (y[3] - y[2]); y34 = (y[3] - y[4])
		y41 = (y[4] - y[1]); y42 = (y[4] - y[2]); y43 = (y[4] - y[3])

		z12 = (z[1] - z[2]); z13 = (z[1] - z[3]); z14 = (z[1] - z[4])
		z21 = (z[2] - z[1]); z23 = (z[2] - z[3]); z24 = (z[2] - z[4])
		z31 = (z[3] - z[1]); z32 = (z[3] - z[2]); z34 = (z[3] - z[4])
		z41 = (z[4] - z[1]); z42 = (z[4] - z[2]); z43 = (z[4] - z[3])

		fV = 0.0
		fV += x21*(y23*z34 - y34*z23)
		fV += x32*(y34*z12 - y12*z34)
		fV += x43*(y12*z23 - y23*z12)
		fV /= 6.0;

		return(fV)
	end

	function loadMSH(sFile::String)
		@printf("\nLoading .msh file: %s \n", sFile)
		hFile = open(sFile)
		arrayLine = readlines(hFile)

		# read the total number of nodes and elements from file
		nNodes::Int = 0
		nElements::Int = 0
		for indexLine = 1:1:length(arrayLine)
			if(contains(arrayLine[indexLine], "\$Nodes"))
				nNodes = parse(Int, arrayLine[indexLine+1])
			end
			if(contains(arrayLine[indexLine], "\$Elements"))
				nElements = parse(Int, arrayLine[indexLine+1])
			end
		end

		# create arrays to store node data
		arrayNode_ID = Array{Int}(nNodes)
		arrayNode_Coordinate = Array{Float64}(3,nNodes)
		for indexLine = 1:1:length(arrayLine)
			if(contains(arrayLine[indexLine], "\$Nodes"))
				@printf("	Reading %d nodes\n", nNodes)
				for indexNode = 1:1:nNodes
					arrayTemp = Array{Float64}(4)
					arrayTemp = readdlm(IOBuffer(arrayLine[indexLine + 1 + indexNode])) # convert multi number string into array of numbers

					arrayNode_ID[indexNode] = arrayTemp[1]
					arrayNode_Coordinate[:,indexNode] = [arrayTemp[2]; arrayTemp[3]; arrayTemp[4]]

					# @printf("	Node_ID %d: Coordinates (%f,%f,%f)\n", arrayNode_ID[indexNode], arrayNode_Coordinate[indexNode,1], arrayNode_Coordinate[indexNode,2], arrayNode_Coordinate[indexNode,3])
				end
			end
		end

		# create arrays to store element data
		arrayElement_ID = Array{Int}(0)
		arrayElement_Corner = Array{Int}(4,0) # for 4-node tetrahedron
		for indexLine = 1:1:length(arrayLine)
			if(contains(arrayLine[indexLine], "\$Elements"))
				@printf("	Reading %d elements\n", nElements)
				for indexElement = 1:1:nElements
					# from gmsh document -> elm-number elm-type number-of-tags < tag > ... node-number-list
					arrayTemp = readdlm(IOBuffer(arrayLine[indexLine + 1 + indexElement])) # convert multi number string into array of numbers

					if(arrayTemp[2] != 4) # element type = 4-node tetrahedron
						nElements -= 1
					end
					if(arrayTemp[2] == 4) # element type = 4-node tetrahedron
						# arrayElement_ID[indexElement] = arrayTemp[1,1]
						arrayElement_ID = vcat(arrayElement_ID, arrayTemp[1])

						nTags = arrayTemp[3]
						newCorners = [arrayTemp[3+nTags+1]; arrayTemp[3+nTags+2]; arrayTemp[3+nTags+3]; arrayTemp[3+nTags+4]]

						arrayElement_Corner = hcat(arrayElement_Corner, newCorners)
					end
				end
			end
		end
		# for indexElement = 1:1:nElements
		# 	@printf("	Element_ID %d: Corners (%d,%d,%d,%d)\n", arrayElement_ID[indexElement], arrayElement_Corner[indexElement,1], arrayElement_Corner[indexElement,2], arrayElement_Corner[indexElement,3], arrayElement_Corner[indexElement,4])
		# end

		@printf("	nNodes: %d, nElements: %d\n", nNodes, nElements)

		# create material points from node and element data
		thisMaterialDomain = Array{mpmMaterialPoint_Tet4}(nElements)
		@printf("	Creating %d material points...", nElements)

		for indexElement = 1:1:nElements
			thisMaterialPoint = mpmMaterialPoint_Tet4()

			# println(arrayElement_Corner[indexElement,1])
			thisMaterialPoint.v3Corner[:,1] = arrayNode_Coordinate[:,arrayElement_Corner[1,indexElement]]
			thisMaterialPoint.v3Corner[:,2] = arrayNode_Coordinate[:,arrayElement_Corner[2,indexElement]]
			thisMaterialPoint.v3Corner[:,3] = arrayNode_Coordinate[:,arrayElement_Corner[3,indexElement]]
			thisMaterialPoint.v3Corner[:,4] = arrayNode_Coordinate[:,arrayElement_Corner[4,indexElement]]
			# println("createdMaterialPointCorner: ", arrayNode_Coordinate[arrayElement_Corner[indexElement,1]])

			v3Centroid = [0.0; 0.0; 0.0]
			v3Centroid += thisMaterialPoint.v3Corner[:,1]
			v3Centroid += thisMaterialPoint.v3Corner[:,2]
			v3Centroid += thisMaterialPoint.v3Corner[:,3]
			v3Centroid += thisMaterialPoint.v3Corner[:,4]
			v3Centroid *= 0.25

			thisMaterialPoint.v3Centroid = v3Centroid

			# println("createdMaterialPointCorner: ", thisMaterialPoint.v3Corner)

			thisMaterialDomain[indexElement] = thisMaterialPoint
		end
		@printf("done\n")
		# for indexElement = 1:1:1
		# 	println("createdMaterialPointCorner: ", thisMaterialDomain[indexElement].v3Corner)
		# end

		return(thisMaterialDomain)
	end

	# function createMaterialDomain_Circle(fCenter_x::Real, fCenter_y::Real, fRadius::Real, fOffset::Real)
	# 	thisMaterialDomain = Array{mpmMaterialPoint}(0)
	#
	# 	fRadius = floor(fRadius/fOffset) * fOffset	#just in case radius is not a multiple of offset
	#
	# 	for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
	# 		for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
	# 			if(fx^2 + fy^2 < fRadius^2)
	# 				thisMaterialPoint = mpmMaterialPoint()
	# 				thisMaterialPoint.v2Position.fx = fCenter_x + fx
	# 				thisMaterialPoint.v2Position.fy = fCenter_y + fy
	# 				push!(thisMaterialDomain, thisMaterialPoint)
	# 			end
	# 		end
	# 	end
	#
	# 	return(thisMaterialDomain)
	# end
end
