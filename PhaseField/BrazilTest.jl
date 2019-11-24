# MPM implementation for the phase field brittle fracture model.
# Two grids are used: one for the mechanical fields similar to standard MPM for solid mechanics.
# Another grid (called mesh) is used to solve for the phase field. They can be different.
# For the moment, the phase field grid is simply a Cartesian grid.
# 2D plain strain with Amor tension/compression split only.
#
# Brazilliant test.
#
# Written by:
# Vinh Phu Nguyen, Monash University, Ausgust 2016
# Based on the MPM implementation done by Sina Sinaie, Monash Uni.

# cd("E:\\MyPublications\\MPM_Julia\\Codes\\Classic_2D_SteelAluminum")
# In Julia terminal, press ; to enter SHELL commands, very convenient.


import PyPlot

push!(LOAD_PATH, pwd())

#pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time", figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

include("MyMaterialPoint.jl")
include("MyGrid.jl")
include("mesh.jl")
include("FiniteElements.jl")
include("MyBasis.jl")
include("ParticleGeneration.jl")

using moduleGrid
using moduleParticleGen

function info(allDeformMP, allRigidMP, allMP, feMesh)
	# ---------------------------------------------------------------------------
	# information about the created domain
	fMass = 0.0
	for iIndex_MP in 1:1:length(allDeformMP)
		fMass += allDeformMP[iIndex_MP].fMass
	end
	@printf("Initial configuration: \n")
	@printf("	Single particle Mass                : %+.6e \n", allDeformMP[1].fMass)
	@printf("	Total mass                          : %+.6e \n", fMass)
	@printf("	Total number of deformable particles: %d \n", length(allDeformMP))
	@printf("	Total number of rigid particles     : %d \n", length(allRigidMP))
	@printf("	Total number of particles           : %d \n", length(allMP))
    @printf("	Total number of finite elements     : %d \n", feMesh.elemCount)
    @printf("	Dimension of finite elements        : %f %f \n", feMesh.deltaX, feMesh.deltaY)
end


# compute shear and bulk modulus given Young and Poisson
function getShearBulkMod(young, poisson, nsd)
 shear    = young/2./(1.+poisson)
 lambda   = young*poisson/(1.+poisson)/(1.-2.*poisson)
 bulk     = lambda + 2.*shear/nsd

 return shear, bulk
end

# Main function, where the problem is defined and solved
function mpmMain()
    const nsd      = 2
	const fGravity = 0.0
    const density  = 2700.0e-12 # density of concrete [kg/m^3]
	const young    = 50.0e3     # MPa
	const poisson  = 0.25
	const k        = 1.0e-18
    const ft       = 10.         # tensile strength in MPa
    const Gc       = 50.0e-3     # fracture energy [N/mm]
    # length scale in phase field model [mm]
    # from 1D exact solution 
    const l0       = (9/16.)^2*young*Gc/ft^2    # length scale in phase field model [mm]
          P        = [0.5 0.5 0.;-0.5 0.5 0.;0. 0. 1.0]

    shear, bulk    = getShearBulkMod(young , poisson, nsd)

    const smallMass= 1.e-12

    const v0       = 0.1e3     # velocity of the impactor [mm/s]

    const fTimeEnd = 0.02
          fTime    = 0.

          ppc      = [3,3]

          phiFile  = "Brazil"
    const interval = 100       # interval for output

    ###############################################################
    # grid creation, this is for MPM (displacement field)
    ###############################################################
    rad      = 50.                  # radius of the disk [mm]
    thick    = 10                   # thickness of the loading platen
    len      = 20                   # length of the loading platen
    Delta    = 3*2                  # three elements added to the left and right
    lx       = 2rad                 # length in X direction [m]
    ly       = 2rad + 2thick        # length in Y direction
    lxn      = lx + Delta           # length in X direction, with extra layers of elements
    lyn      = ly                   # length in Y direction
    ex       = 50+3*2               # number of elements in X direction, with extra layers left/right
    ey       = 60                   # number of elements in Y direction
    thisGrid = moduleGrid.mpmGrid(lxn, lyn, ex+1, ey+1)

    ###############################################################
    # material points, generation
    ###############################################################
    # array holding all material points (these are references to rollerLeft etc.)
    # and array holding all rigid particles
	allDeformMP        = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)
	allRigidMP         = Array{Any,1}(0)
	allMaterialPoints  = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)

    disk    = moduleParticleGen.createMaterialDomain_Circle([rad+0.5Delta; rad+thick], rad, thisGrid, ppc)
    for iIndex_MP = 1:1:length(disk) 
        disk[iIndex_MP].fMass           = disk[iIndex_MP].fVolume*density 
        disk[iIndex_MP].fElasticModulus = young 
        disk[iIndex_MP].fPoissonRatio   = poisson
        disk[iIndex_MP].fShear          = shear
        disk[iIndex_MP].fBulk           = bulk
		disk[iIndex_MP].v2Velocity      = [0.0; 0.0]
        disk[iIndex_MP].v2Momentum      = disk[iIndex_MP].fMass*disk[iIndex_MP].v2Velocity
		disk[iIndex_MP].fColor          = 1.0
        
        push!(allDeformMP, disk[iIndex_MP])
        push!(allMaterialPoints, disk[iIndex_MP])
	end

    
    # IN Julia you can simply write a=2x, which means a=2*x

    platenB = moduleParticleGen.createMaterialDomain_Rectangle([0.5(lxn-len) 0.0;0.5(lxn+len) thick], thisGrid, ppc)
	for iIndex_MP = 1:length(platenB)
        platenB[iIndex_MP].fMass           = platenB[iIndex_MP].fVolume*density
		platenB[iIndex_MP].fElasticModulus = young
		platenB[iIndex_MP].fPoissonRatio   = poisson
        platenB[iIndex_MP].fShear          = shear
        platenB[iIndex_MP].fBulk           = bulk
		platenB[iIndex_MP].v2Velocity      = [0.0;v0]
        platenB[iIndex_MP].v2Momentum      = platenB[iIndex_MP].fMass*platenB[iIndex_MP].v2Velocity
		platenB[iIndex_MP].fColor          = 2.0

        push!(allMaterialPoints, platenB[iIndex_MP])
	end
    
    platenT = moduleParticleGen.createMaterialDomain_Rectangle([0.5(lxn-len) lyn-thick; 0.5(lxn+len) lyn], thisGrid, ppc)
	for iIndex_MP = 1:length(platenT)
        platenT[iIndex_MP].fMass           = platenT[iIndex_MP].fVolume*density
		platenT[iIndex_MP].fElasticModulus = young
		platenT[iIndex_MP].fPoissonRatio   = poisson
        platenT[iIndex_MP].fShear          = shear
        platenT[iIndex_MP].fBulk           = bulk
		platenT[iIndex_MP].v2Velocity      = [0.0;-v0]
        platenT[iIndex_MP].v2Momentum      = platenT[iIndex_MP].fMass*platenT[iIndex_MP].v2Velocity
		platenT[iIndex_MP].fColor          = 2.0

        push!(allMaterialPoints, platenT[iIndex_MP])
	end
    
    push!(allRigidMP, platenT)
    push!(allRigidMP, platenB)
    
    #=
    @printf("# of partiles of left roller : %d\n", length(rollerLeft))
    @printf("# of partiles of right roller: %d\n", length(rollerRight))
    @printf("# of partiles of mid roller  : %d\n", length(rollerMid))
    @printf("# of partiles of beam        : %d\n", length(beam))
    =#
   

    ###############################################################
    # creating the FE mesh for the phase field equation
    ###############################################################
    exP = 50+3*2              # number of elements in X direction
    eyP = 60                   # number of elements in Y direction
    feMesh      = FeMesh.Mesh(2, lxn, lyn, exP, eyP)
    phase_field = zeros(feMesh.nodeCount)
    FeMesh.update(feMesh,allDeformMP)     # find which material points locate in which elements

    #print(feMesh.elem2MP)

    # summary of problem data
    info(allDeformMP, allRigidMP, allMaterialPoints, feMesh)
    @printf("	length scale                : %6f3 \n", l0)


    
    # PyPlot to plot the particles and the grid
    #=
	iMaterialPoints = length(allDeformMP)
	array_x         = [allDeformMP[i].v2Centroid[1] for i in 1:iMaterialPoints]
    array_y         = [allDeformMP[i].v2Centroid[2] for i in 1:iMaterialPoints]
    array_color     = Array{Real}(iMaterialPoints, 3)
	array_size      = Array{Real}(iMaterialPoints, 1)
    for iIndex in 1:1:iMaterialPoints
      array_color[iIndex, :] = [1.0, 0.0, 0.0]#[thisGrid.GridPoints[iIndex].fMass/iMaterialPoints, 0.0, 0.0]
	  array_size[iIndex, :] = [5.0]
    end

    pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time", figsize=(6., 6.), edgecolor="white", facecolor="white")
	pyPlot01 = PyPlot.gca()
	# pyPlot01 = PyPlot.subplot2grid((1,1), (0,0), colspan=1, rowspan=1, aspect="equal")
	PyPlot.scatter(array_x, array_y, c=array_color, lw=0, s=array_size)

    FeMesh.plot_mesh(feMesh)
    moduleMaterialPoint.VTKParticles(allDeformMP,"ThreePointBending0.vtp")
    =#
    

    c               = sqrt(young/density)
    criticalTimeInc = thisGrid.v2Length_Cell[1] / c
    fTimeIncrement  = 0.2 * criticalTimeInc
    @printf("	Initial time step                : %6f3 \n", fTimeIncrement)
    @printf("	Total time                       : %6f3 \n", fTimeEnd      )
    @printf("	Total number of time steps       : %d \n",   fTimeEnd /fTimeIncrement     )

    iStep           = 0

    # main analysis loop
	while fTime < fTimeEnd
	    @printf("Time step %d                     : %+.6e \n", iStep, fTime)
        #######################################
        #  solve for the phase field
        #######################################

        solve_fe(phase_field, feMesh, allDeformMP, l0, Gc)

        #######################################
        #  solve for the mechanical fields
        #######################################

		#reset grid---------------------------------------------------------------------
		for iIndex in 1:thisGrid.iNodes
			thisGrid.GridPoints[iIndex].fMass      = 0.0
			thisGrid.GridPoints[iIndex].v2Velocity = [0.0; 0.0]  # needed in double mapping technique
			thisGrid.GridPoints[iIndex].v2Momentum = [0.0; 0.0]
			thisGrid.GridPoints[iIndex].v2Force    = [0.0; 0.0]
            # Note that Julia stores arrays column wise, so should use column vectors
		end
		# material to grid ------------------------------------------------------------
		for iIndex_MP in 1:length(allDeformMP)
			thisMaterialPoint      = allDeformMP[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			for iIndex in 1:length(thisAdjacentGridPoints)
				# sina, be careful here, this might not be by reference and might not be good for assignment
                # ???
				thisGridPoint    = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]
				fShapeValue, v2ShapeGradient  = moduleBasis.getShapeAndGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				# mass, momentum, internal force and external force
				thisGridPoint.fMass      += fShapeValue * thisMaterialPoint.fMass
				thisGridPoint.v2Momentum += fShapeValue * thisMaterialPoint.v2Momentum
				fVolume                   = thisMaterialPoint.fVolume
				thisGridPoint.v2Force[1] -= fVolume * (v2ShapeGradient[1]*thisMaterialPoint.v3Stress[1] +
                                                       v2ShapeGradient[2]*thisMaterialPoint.v3Stress[3])
				thisGridPoint.v2Force[2] -= fVolume * (v2ShapeGradient[2]*thisMaterialPoint.v3Stress[2] +
                                                       v2ShapeGradient[1]*thisMaterialPoint.v3Stress[3])
                # unlike Matlab you do not need ... to break codes into lines
				# external forces
				#thisGridPoint.v2Force    += fShapeValue*thisMaterialPoint.v2ExternalForce
		   end
		end
		# update grid momentum and apply boundary conditions ---------------------
		for iIndex_GP in 1:thisGrid.iNodes
			thisGridPoint = thisGrid.GridPoints[iIndex_GP]
			thisGridPoint.v2Momentum += thisGridPoint.v2Force * fTimeIncrement

            #=
			if(thisGridPoint.v2Fixed[1] == true)
				thisGridPoint.v2Momentum[1] = 0.0
				thisGridPoint.v2Force[1]    = 0.0
			end
			if(thisGridPoint.v2Fixed[2] == true)
				thisGridPoint.v2Momentum[2] = 0.0
				thisGridPoint.v2Force[2]    = 0.0
			end
            =#
		end
    
		# update grid momentum and apply boundary conditions ---------------------
        # for rigid material points
        for rb=1:length(allRigidMP)
          rigidParticles = allRigidMP[rb]
          moveNodes = moduleGrid.getGridPointsForParticles(rigidParticles, thisGrid)
          for i=1:length(moveNodes)
            thisGridPoint = thisGrid.GridPoints[moveNodes[i]]
            thisGridPoint.v2Momentum = thisGridPoint.fMass * rigidParticles[1].v2Velocity
		    thisGridPoint.v2Force[2]    = 0.0
          end
        end
		
        #= double mapping technique of Sulsky
        # ------------------------------------------------------------------------
		# grid to material pass 1--update particle velocity only,
		for iIndex_MP in 1:length(allDeformMP)
			thisMaterialPoint = allDeformMP[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			for iIndex in 1:length(thisAdjacentGridPoints)
				thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]
				fShapeValue = moduleBasis.getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
				thisMaterialPoint.v2Velocity += (fShapeValue * thisGridPoint.v2Force / thisGridPoint.fMass) * fTimeIncrement
		   end
		end
		# map particle momenta back to grid------------------------------
		# mass in NOT mapped here
		for iIndex_MP in 1:length(allDeformMP)
			thisMaterialPoint = allDeformMP[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			for iIndex in 1:length(thisAdjacentGridPoints)
				thisGridPoint             = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]
				fShapeValue               = moduleBasis.getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
				thisGridPoint.v2Velocity += fShapeValue * thisMaterialPoint.fMass * thisMaterialPoint.v2Velocity / thisGridPoint.fMass
			end
		end
		#apply boundary conditions velocity------------------------------------
		for iIndex_GP in 1:thisGrid.iNodes
			thisGridPoint = thisGrid.GridPoints[iIndex_GP]
			if(thisGridPoint.v2Fixed[1] == true)
				thisGridPoint.v2Velocity[1] = 0.0
				thisGridPoint.v2Momentum[1] = 0.0
				thisGridPoint.v2Force[1] = 0.0
			end
			if(thisGridPoint.v2Fixed[2] == true)
				thisGridPoint.v2Velocity[2] = 0.0
				thisGridPoint.v2Momentum[2] = 0.0
				thisGridPoint.v2Force[2] = 0.0
			end
			if(thisGridPoint.v2Fixed[2] == true)
				thisGridPoint.v2Velocity[2] = 0.0
				thisGridPoint.v2Momentum[2] = 0.0
				thisGridPoint.v2Force[2] = 0.0
			end
		end
        # apply boundary condiftions velocity for rigid particles-----------------
        for i=1:length(moveNodes)
          thisGridPoint = thisGrid.GridPoints[i]
          thisGridPoint.v2Velocity[1] = 0.
          thisGridPoint.v2Velocity[2] = allRigidMP[1].v2Velocity[2]
        end
        =#

		# ------------------------------------------------------------------------
        # grid to material pass2 (update remaining stuff)-------------------------
		for iIndex_MP in 1:length(allDeformMP)
			thisMaterialPoint      = allDeformMP[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			v2CentroidIncrement = zeros(2)
			for iIndex in 1:length(thisAdjacentGridPoints)
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
				#thisMaterialPoint.m22DeformationGradientIncrement += thisGridPoint.v2Velocity*transpose(v2ShapeGradient)*fTimeIncrement;
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
              if ni > size(feMesh.nodes, 2) 
                @printf("%f %f\n",thisMaterialPoint.v2Centroid[1],thisMaterialPoint.v2Centroid[2]) 
              end
              xi       = feMesh.nodes[:,ni]
              Ni       = getN(thisMaterialPoint.v2Centroid, xi, feMesh.deltaX, feMesh.deltaY)
              phi     += Ni*phase_field[ni]
            end
            # phi = 0.0
			thisMaterialPoint.fPhase  = phi

            stress = zeros(3)
            stress = moduleMaterialPoint.getStress(thisMaterialPoint.fBulk,
                                                   thisMaterialPoint.fShear, phi, k, thisMaterialPoint.v3Strain)
			thisMaterialPoint.v3Stress     = stress
            #=
            v3StressIncrement = moduleMaterialPoint.getStressIncrement_Elastic(thisMaterialPoint.fElasticModulus,
                                       thisMaterialPoint.fPoissonRatio, v3StrainIncrement)
			thisMaterialPoint.v3Stress     += v3StressIncrement
    =#

            # update history H
            traceEps   = thisMaterialPoint.v3Strain[1] + thisMaterialPoint.v3Strain[2]
            traceEpsP  = 0.5(traceEps+abs(traceEps))
            newHistory = 0.5 * thisMaterialPoint.fBulk * traceEpsP^2
            + thisMaterialPoint.fShear*(thisMaterialPoint.v3Strain'*P*thisMaterialPoint.v3Strain)[1]
            # without ()[1] it is an array of 1 element!!! not scalar

			thisMaterialPoint.fVolume      = det(thisMaterialPoint.m22DeformationGradient) * thisMaterialPoint.fVolumeInitial
			thisMaterialPoint.v2Momentum   = thisMaterialPoint.v2Velocity * thisMaterialPoint.fMass
			thisMaterialPoint.fHistory     = newHistory > thisMaterialPoint.fHistory ? newHistory : thisMaterialPoint.fHistory
		end
        
        # update position of rigid particles
        for rb=1:length(allRigidMP)
          rigidParticles = allRigidMP[rb]
		  for iIndex_MP in 1:length(rigidParticles)
	        thisMaterialPoint             = rigidParticles[iIndex_MP]
		    thisMaterialPoint.v2Centroid += fTimeIncrement*thisMaterialPoint.v2Velocity
		  end
        end


		# which material points in which finite elements
		FeMesh.update(feMesh,allDeformMP)

		#@printf("	Initial time step   : %6f3 \n", fTime)

		if ( iStep % interval == 0 )
			FeMesh.vtk(feMesh, phase_field, string(phiFile,"$iStep.vtu" ))
			moduleMaterialPoint.VTKParticles(allMaterialPoints,"Brazil$(iStep).vtp")
		end
		fTime += fTimeIncrement
		iStep += 1
  end           # end of time loop
end             # end of mpmMain()

mpmMain()
