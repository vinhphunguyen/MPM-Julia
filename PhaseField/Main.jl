# MPM implementation for the phase field brittle fracture model.
# Two grids are used: one for the mechanical fields similar to standard MPM for solid mechanics.
# Another grid (called mesh) is used to solve for the phase field. They can be different.
# For the moment, the phase field grid is simply a Cartesian grid.
# 2D plain strain with Amor tension/compression split only. 
#
# Written by:
# Vinh Phu Nguyen, Monash University, Ausgust 2016
# Based on the MPM implementation done by Sina Sinaie, Monash Uni.

# cd("E:\\MyPublications\\MPM_Julia\\Codes\\Classic_2D_SteelAluminum")
import PyPlot

push!(LOAD_PATH, pwd())

#pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time", figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

include("MyMaterialPoint.jl")
include("MyGrid.jl")
include("mesh.jl")
include("FiniteElements.jl")
include("MyBasis.jl")

function mpmMain()
    const nsd      = 2  
	const fGravity = 0.0
	const density  = 1.0
	const young    = 1.0
	const poisson  = 0.3
	const k        = 1.0e-10
    const shear    = young/2./(1.+poisson)
          lambda   = young*poisson/(1.+poisson)/(1.-2.*poisson)
    const bulk     = lambda + 2.*shear/nsd
          P        = [0.5 0.5 0.;-0.5 0.5 0.;0. 0. 1.0]

    const fTimeEnd = 1.0 
          fTime    = 0.

    ###############################################################
    # grid creation, this is for MPM (displacement field)
    ###############################################################
    lx = 1.0                  # length in X direction
    ly = 1.0                  # length in Y direction
    ex = 50                   # number of elements in X direction
    ey = 50                   # number of elements in Y direction
    thisGrid = moduleGrid.mpmGrid(lx, ly, ex+1, ey+1)

    ###############################################################
    # material points, generation
    ###############################################################
	# array holding all material points (these are references to MaterialDomain_01 & 02)
	allMaterialPoint = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)
    fOffset = lx/(ex+1.)/2
	thisMaterialDomain_02 = moduleMaterialPoint.createMaterialDomain_Rectangle([0.5; 0.5],
                                                                            lx, ly, fOffset)
	for iIndex_MP = 1:1:length(thisMaterialDomain_02)
		fVolume   = fOffset*fOffset#3.14159*0.2*0.2/length(thisMaterialDomain_02)
		fMass     = density*fVolume
		thisMaterialDomain_02[iIndex_MP].fMass           = fMass
		thisMaterialDomain_02[iIndex_MP].fVolumeInitial  = fVolume
		thisMaterialDomain_02[iIndex_MP].fVolume         = fVolume
		thisMaterialDomain_02[iIndex_MP].fElasticModulus = young
		thisMaterialDomain_02[iIndex_MP].fPoissonRatio   = poisson
        thisMaterialDomain_02[iIndex_MP].fShear          = shear
        thisMaterialDomain_02[iIndex_MP].fBulk           = bulk
		thisMaterialDomain_02[iIndex_MP].v2Velocity = [0.0; 0.0]
		thisMaterialDomain_02[iIndex_MP].v2Momentum = fMass*thisMaterialDomain_02[iIndex_MP].v2Velocity
		thisMaterialDomain_02[iIndex_MP].v2ExternalForce = [0.0; -fGravity*fMass]

		thisMaterialDomain_02[iIndex_MP].m22DeformationGradient = eye(2,2)
		thisMaterialDomain_02[iIndex_MP].m22DeformationGradientIncrement = eye(2,2)
		
        push!(allMaterialPoint, thisMaterialDomain_02[iIndex_MP])
	end
	#for iIndex_MP = 1:1:length(thisMaterialDomain_02)
	#	push!(allMaterialPoint, thisMaterialDomain_02[iIndex_MP])
	#end

    ###############################################################
    # creating the FE mesh for the phase field equation
    ###############################################################
    exP = 50                   # number of elements in X direction
    eyP = 50                   # number of elements in Y direction
    feMesh      = FeMesh.Mesh(2, lx, ly, exP, eyP) 
    phase_field = zeros(feMesh.nodeCount)
    FeMesh.update(feMesh,allMaterialPoint)     # find which material points locate in which elements
    #solve_fe(phi,feMesh,allMaterialPoint)
    #FeMesh.vtk(feMesh, phi, "phi.vtu" )
    

    #print(feMesh.elem2MP[1])

	# ---------------------------------------------------------------------------
	# information about the created domain
	fMass = 0.0
	for iIndex_MP in 1:1:length(allMaterialPoint)
		fMass += allMaterialPoint[iIndex_MP].fMass
	end
	@printf("Initial configuration: \n")
	@printf("	Single particle Mass: %+.6e \n", allMaterialPoint[1].fMass)
	@printf("	Total mass: %+.6e \n", fMass)
	@printf("	Total number of material points: %d \n", length(allMaterialPoint))
    @printf("	Total number of finite elements: %d \n", feMesh.elemCount)

    #=
	iMaterialPoints = length(allMaterialPoint)
	array_x         = [allMaterialPoint[i].v2Centroid[1] for i in 1:iMaterialPoints]
    array_y         = [allMaterialPoint[i].v2Centroid[2] for i in 1:iMaterialPoints]
    array_color     = Array{Real}(iMaterialPoints, 3)
	array_size      = Array{Real}(iMaterialPoints, 1)
    for iIndex in 1:1:iMaterialPoints
      array_color[iIndex, :] = [1.0, 0.0, 0.0]#[thisGrid.GridPoints[iIndex].fMass/iMaterialPoints, 0.0, 0.0]
	  array_size[iIndex, :] = [5.0]
    end

	pyPlot01 = PyPlot.gca()
	# pyPlot01 = PyPlot.subplot2grid((1,1), (0,0), colspan=1, rowspan=1, aspect="equal")
	PyPlot.scatter(array_x, array_y, c=array_color, lw=0, s=array_size)
    =#
    
    c               = sqrt(young/density)
    criticalTimeInc = thisGrid.v2Length_Cell[1] / c
    fTimeIncrement  = 0.5 * criticalTimeInc
	
    # main analysis loop
	while fTime < fTimeEnd
        #######################################
        #  solve for the phase field
        #######################################
        solve_fe(phase_field, feMesh, allMaterialPoint)

        #######################################
        #  solve for the mechanical fields 
        #######################################

		#reset grid---------------------------------------------------------------------
		for iIndex in 1:1:thisGrid.iNodes
			thisGrid.GridPoints[iIndex].fMass      = 0.0
			thisGrid.GridPoints[iIndex].v2Momentum = [0.0; 0.0]
			thisGrid.GridPoints[iIndex].v2Force    = [0.0; 0.0]
            # Note that Julia stores arrays column wise, so should use column vectors
		end
		# material to grid ------------------------------------------------------------
		for iIndex_MP in 1:length(allMaterialPoint)
			thisMaterialPoint      = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			for iIndex in 1:length(thisAdjacentGridPoints)
				# sina, be careful here, this might not be by reference and might not be good for assignment
                # ???
				thisGridPoint    = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]
				fShapeValue, v2ShapeGradient  = moduleBasis.getShapeAndGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				# mass, momentum, internal force and external force
				thisGridPoint.fMass      += fShapeValue * thisMaterialPoint.fMass
				thisGridPoint.v2Momentum += fShapeValue * thisMaterialPoint.fMass * thisMaterialPoint.v2Velocity
				fVolume                   = thisMaterialPoint.fVolume
				thisGridPoint.v2Force[1] -= fVolume * (v2ShapeGradient[1]*thisMaterialPoint.v3Stress[1] + 
                                                       v2ShapeGradient[2]*thisMaterialPoint.v3Stress[3])
				thisGridPoint.v2Force[2] -= fVolume * (v2ShapeGradient[2]*thisMaterialPoint.v3Stress[2] +
                                                       v2ShapeGradient[1]*thisMaterialPoint.v3Stress[3])
                # unlike Matlab you do not need ... to break codes into lines
				# external forces
				thisGridPoint.v2Force    += fShapeValue*thisMaterialPoint.v2ExternalForce
		   end
		end
		# update grid momentum and apply boundary conditions ---------------------
		for iIndex_GP in 1:1:thisGrid.iNodes
			thisGridPoint = thisGrid.GridPoints[iIndex_GP]
			thisGridPoint.v2Momentum += thisGridPoint.v2Force * fTimeIncrement

			if(thisGridPoint.v2Fixed[1] == true)
				thisGridPoint.v2Momentum[1] = 0.0
				thisGridPoint.v2Force[1]    = 0.0
			end
			if(thisGridPoint.v2Fixed[2] == true)
				thisGridPoint.v2Momentum[2] = 0.0
				thisGridPoint.v2Force[2]    = 0.0
			end
		end

		# ------------------------------------------------------------------------
		# grid to material -------------------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint      = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			v2CentroidIncrement = zeros(2)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				thisGridPoint   = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]
				fShapeValue, v2ShapeGradient  = moduleBasis.getShapeAndGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				v2GridPointVelocity = zeros(2)
				if(thisGridPoint.fMass > 1.0e-8)
				  v2GridPointVelocity = thisGridPoint.v2Momentum / thisGridPoint.fMass
				  thisMaterialPoint.v2Velocity += (fShapeValue * thisGridPoint.v2Force / thisGridPoint.fMass) * fTimeIncrement
				  v2CentroidIncrement += (fShapeValue * thisGridPoint.v2Momentum / thisGridPoint.fMass) * fTimeIncrement
				end

				# from (2011) A convected particle domain interpolation technique to extend ...
                # Note that m22DeformationGradientIncrement was initialised to identity matrix
				thisMaterialPoint.m22DeformationGradientIncrement += v2GridPointVelocity*transpose(v2ShapeGradient)*fTimeIncrement;
		   end
			thisMaterialPoint.v2Centroid += v2CentroidIncrement
			thisMaterialPoint.m22DeformationGradient = thisMaterialPoint.m22DeformationGradientIncrement * thisMaterialPoint.m22DeformationGradient
			v3StrainIncrement = zeros(3)
			v3StrainIncrement[1] = thisMaterialPoint.m22DeformationGradientIncrement[1,1] - 1.0
			v3StrainIncrement[2] = thisMaterialPoint.m22DeformationGradientIncrement[2,2] - 1.0
			v3StrainIncrement[3] = thisMaterialPoint.m22DeformationGradientIncrement[1,2] + 
                                   thisMaterialPoint.m22DeformationGradientIncrement[2,1]
			thisMaterialPoint.m22DeformationGradientIncrement = eye(2,2)

			thisMaterialPoint.v3Strain[1] += v3StrainIncrement[1]
			thisMaterialPoint.v3Strain[2] += v3StrainIncrement[2]
			thisMaterialPoint.v3Strain[3] += v3StrainIncrement[3]

            # get phase  field from mesh to material points
            neighborMeshPnts = FeMesh.findAdjacentPoints(thisMaterialPoint.v2Centroid, feMesh )
            phi = 0.0
            for i=1:length(neighborMeshPnts)
              ni       = neighborMeshPnts[i]
              xi       = feMesh.nodes[:,ni]
              Ni       = getN(thisMaterialPoint.v2Centroid, xi, feMesh.deltaX, feMesh.deltaY)
              phi     += Ni*phase_field[ni]
            end

            stress = zeros(3)
            stress = moduleMaterialPoint.getStress(bulk, shear, phi, k, thisMaterialPoint.v3Strain)
			thisMaterialPoint.v3Stress     = stress

            # update history H
            traceEps   = thisMaterialPoint.v3Strain[1] + thisMaterialPoint.v3Strain[2]
            traceEpsP  = 0.5(traceEps+abs(traceEps))
            newHistory = 0.5 * bulk * traceEpsP^2
            + shear*(thisMaterialPoint.v3Strain'*P*thisMaterialPoint.v3Strain)[1]
            # without ()[1] it is an array of 1 element!!! not scalar


			thisMaterialPoint.fVolume      = det(thisMaterialPoint.m22DeformationGradient) * thisMaterialPoint.fVolumeInitial
			thisMaterialPoint.v2Momentum   = thisMaterialPoint.v2Velocity * thisMaterialPoint.fMass
			thisMaterialPoint.fHistory     = newHistory > thisMaterialPoint.fHistory ? newHistory : thisMaterialPoint.fHistory

            # which material points in which finite elements
            FeMesh.update(feMesh,allMaterialPoint)     

            fTime += fTimeIncrement
		end
  end           # end of 
end             # end of mpmMain()

mpmMain()
