# MPM implementation for the non local Mazars damage model.
#
# Three point bending test with roller supports. And rigid particles are used to model
# imposed velocity.
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
include("MyBasis.jl")
include("ParticleGeneration.jl")

using moduleGrid
using moduleParticleGen

function info(allDeformMP, allRigidMP, allMP)
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
    const density  = 2400.0    # density of concrete [kg/m^3]
	const young    = 70.0e9    # Pa
	const poisson  = 0.3
	const young1   = 70.0e18   # Pa
    const l0       = 0.0005    # length scale in phase field model [m]
    const kappa0   = 0.001
    const kappa1   = 0.01


    shear, bulk    = getShearBulkMod(young , poisson, nsd)
    shear1,bulk1   = getShearBulkMod(young1, poisson, nsd)

    const smallMass= 1.e-12

    const v0       = 1.0     # velocity of the impactor [m/s]

    const fTimeEnd = 0.02
          fTime    = 0.

          ppc      = [3,3]

          phiFile  = "ThreePointNL"
    const interval = 100       # interval for output

    ###############################################################
    # grid creation, this is for MPM (displacement field)
    ###############################################################
    delta    = 0.0/1000             # gap between the impacter and the beam
    Delta    = 3*0.00025            # three elements added to the left and right
    lx       = 0.04                 # length in X direction [m]
    lxn      = 0.04 + 2*Delta       # length in X direction, with extra layers of elements [m]
    ly       = 0.005 + delta        # length in Y direction
    lyn      = ly + Delta           # length in Y direction
    ddd      = 0.0005               # from the left extremity of beam to the left roller
    ey       = 20                   # number of elements in Y direction
    ex       = 160+3*2              # number of elements in X direction
    thisGrid = moduleGrid.mpmGrid(lxn, lyn, ex+1, ey+1)

    fixGridNodes = collect(1:thisGrid.v2Nodes[1])
    
    const nonlocalRad = 2thisGrid.v2Length_Cell[1]
    const c           = 3./(pi*nonlocalRad^2)

    # fix boundary conditions on the grid		
    for iIndex in 1:thisGrid.iNodes
        if ( length(find(fixGridNodes .== iIndex)) != 0 )
          thisGrid.GridPoints[iIndex].v2Fixed = [true;true]
        end
	end

    ###############################################################
    # material points, generation
    ###############################################################
    # array holding all material points (these are references to rollerLeft etc.)
    # and array holding all rigid particles
	allDeformMP        = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)
	allRigidMP         = Array{Any,1}(0)
	allMaterialPoints  = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)

    radius  = 0.001/2.     # radius of rollers

    rollerLeft    = moduleParticleGen.createMaterialDomain_Circle([ddd+Delta; radius], radius, thisGrid, ppc)
    for iIndex_MP = 1:1:length(rollerLeft) 
        rollerLeft[iIndex_MP].fMass           = rollerLeft[iIndex_MP].fVolume*density 
        rollerLeft[iIndex_MP].fElasticModulus = young1 
        rollerLeft[iIndex_MP].fPoissonRatio   = poisson
        rollerLeft[iIndex_MP].fShear          = shear1
        rollerLeft[iIndex_MP].fBulk           = bulk1
		rollerLeft[iIndex_MP].v2Velocity = [0.0; 0.0]
        rollerLeft[iIndex_MP].v2Momentum = rollerLeft[iIndex_MP].fMass*rollerLeft[iIndex_MP].v2Velocity
		rollerLeft[iIndex_MP].fColor          = 1.0
        
        #push!(allDeformMP, rollerLeft[iIndex_MP])
        push!(allMaterialPoints, rollerLeft[iIndex_MP])
	end

    rollerRight = moduleParticleGen.createMaterialDomain_Circle([lx-ddd+Delta; radius], radius, thisGrid, ppc)
	for iIndex_MP = 1:1:length(rollerRight)
        rollerRight[iIndex_MP].fMass           = rollerRight[iIndex_MP].fVolume*density 
		rollerRight[iIndex_MP].fElasticModulus = young1
		rollerRight[iIndex_MP].fPoissonRatio   = poisson
        rollerRight[iIndex_MP].fShear          = shear1
        rollerRight[iIndex_MP].fBulk           = bulk1
		rollerRight[iIndex_MP].v2Velocity = [0.0; 0.0]
        rollerRight[iIndex_MP].v2Momentum = rollerRight[iIndex_MP].fMass*rollerRight[iIndex_MP].v2Velocity
        rollerRight[iIndex_MP].v2ExternalForce = [0.0; 0.0]
		rollerRight[iIndex_MP].fColor          = 1.0

        #push!(allDeformMP, rollerRight[iIndex_MP])
        push!(allMaterialPoints, rollerRight[iIndex_MP])
	end

    
    # IN Julia you can simply write a=2x, which means a=2*x
    # beam with a notch of 2 cells x 3 cells in the middle
    nc   = 1 # => 2 cells with
    #nc   = 0.5 # => 1 cells with

    beam = moduleParticleGen.createMaterialDomain_RectangleWithNotch ([0.+Delta 2radius;lx+Delta 2radius+0.003], 
                                                               [0.5lxn-nc*thisGrid.v2Length_Cell[1] 2radius;
                                                                0.5lxn+nc*thisGrid.v2Length_Cell[1] 2radius+5*thisGrid.v2Length_Cell[2]], thisGrid, ppc)
	for iIndex_MP = 1:length(beam)
        beam[iIndex_MP].fMass           = beam[iIndex_MP].fVolume*density
		beam[iIndex_MP].fElasticModulus = young
		beam[iIndex_MP].fPoissonRatio   = poisson
        beam[iIndex_MP].fShear          = shear
        beam[iIndex_MP].fBulk           = bulk
		beam[iIndex_MP].v2Velocity = [0.0; 0.0]
        beam[iIndex_MP].v2Momentum = beam[iIndex_MP].fMass*beam[iIndex_MP].v2Velocity
		beam[iIndex_MP].fColor         = 3.0

        push!(allDeformMP, beam[iIndex_MP])
        push!(allMaterialPoints, beam[iIndex_MP])
	end
    
    rollerMid = moduleParticleGen.createMaterialDomain_Circle([0.5lx+Delta; 0.003+3*radius+delta], radius, thisGrid, ppc)
	for iIndex_MP = 1:length(rollerMid)
        rollerMid[iIndex_MP].fMass           = rollerMid[iIndex_MP].fVolume*density
		rollerMid[iIndex_MP].fElasticModulus = young1
		rollerMid[iIndex_MP].fPoissonRatio   = poisson
        rollerMid[iIndex_MP].fShear          = shear1
        rollerMid[iIndex_MP].fBulk           = bulk1
		rollerMid[iIndex_MP].v2Velocity = [0.0; -v0]
        rollerMid[iIndex_MP].v2Momentum = rollerMid[iIndex_MP].fMass*rollerMid[iIndex_MP].v2Velocity
		rollerMid[iIndex_MP].fColor         = 2.0

        #push!(allRigidMP, rollerMid[iIndex_MP])
        push!(allMaterialPoints, rollerMid[iIndex_MP])
	end
   
    push!(allRigidMP, rollerLeft )
    push!(allRigidMP, rollerRight)
    push!(allRigidMP, rollerMid  )

    
    #=
    @printf("# of partiles of left roller : %d\n", length(rollerLeft))
    @printf("# of partiles of right roller: %d\n", length(rollerRight))
    @printf("# of partiles of mid roller  : %d\n", length(rollerMid))
    @printf("# of partiles of beam        : %d\n", length(beam))
    =#
   
    # summary of problem data
    info(allDeformMP, allRigidMP, allMaterialPoints)

    moduleGrid.update( thisGrid, allDeformMP )
    
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

			if(thisGridPoint.v2Fixed[1] == true)
				thisGridPoint.v2Momentum[1] = 0.0
				thisGridPoint.v2Force[1]    = 0.0
			end
			if(thisGridPoint.v2Fixed[2] == true)
				thisGridPoint.v2Momentum[2] = 0.0
				thisGridPoint.v2Force[2]    = 0.0
			end
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

           # total stress, (1-D) not taken into account 
           v3StressIncrement = moduleMaterialPoint.getStressIncrement_Elastic(thisMaterialPoint.fElasticModulus,
                                       thisMaterialPoint.fPoissonRatio, v3StrainIncrement)
		   thisMaterialPoint.v3Stress     += v3StressIncrement
		  
           eps = moduleMaterialPoint.getPrincipalStrains(thisMaterialPoint.v3Strain)
           thisMaterialPoint.fEqvStrain   = moduleMaterialPoint.getMazarsStrain(eps)
		   thisMaterialPoint.fVolume      = det(thisMaterialPoint.m22DeformationGradient) * thisMaterialPoint.fVolumeInitial
		   thisMaterialPoint.v2Momentum   = thisMaterialPoint.v2Velocity * thisMaterialPoint.fMass
		   thisMaterialPoint.fVolume      = det(thisMaterialPoint.m22DeformationGradient) * thisMaterialPoint.fVolumeInitial
		end

		# ------------------------------------------------------------------------
        # grid to material pass2 (update stresses)
		for iIndex_MP in 1:length(allDeformMP)
		  thisMaterialPoint = allDeformMP[iIndex_MP]
          xp                = thisMaterialPoint.v2Centroid
          neighbors         = moduleGrid.getNeibourParticleIDs(xp, nonlocalRad, thisGrid, allDeformMP)

          # determine nonlocal equivalent strain
          neqvStrain        = 0.
          weight            = 0.

          for ip=1:length(neighbors)
            id  = neighbors[ip]
            mp  = allDeformMP[id]
            xi  = mp.v2Centroid
            eqv = mp.fEqvStrain
            ri  = sqrt( (xp[1] - xi[1])^2 + (xp[2]-xi[2])^2  )
            wi  = c*(1.-(ri/nonlocalRad)^2)^2 

            neqvStrain += wi * eqv 
            weight     += wi;
          end

          # final nonlocal equivalent strain
          neqvStrain /= weight

          D = (neqvStrain > kappa0) ?  kappa0/neqvStrain*exp((kappa0-neqvStrain)/(kappa1-kappa0)):0.
          D = ( thisMaterialPoint.fPhase < D ) ? D : thisMaterialPoint.fPhase 
          thisMaterialPoint.fPhase      = D
          thisMaterialPoint.v3Stress   *= (1.-D)
		end
        
        # update position of rigid particles
        for rb=1:length(allRigidMP)
          rigidParticles = allRigidMP[rb]
		  for iIndex_MP in 1:length(rigidParticles)
	        thisMaterialPoint             = rigidParticles[iIndex_MP]
		    thisMaterialPoint.v2Centroid += fTimeIncrement*thisMaterialPoint.v2Velocity
		  end
        end

		#@printf("	Initial time step   : %6f3 \n", fTime)

		if ( iStep % interval == 0 )
			moduleMaterialPoint.VTKParticles(allMaterialPoints,"ThreePointNLD$(iStep).vtp")
		end
		fTime += fTimeIncrement
		iStep += 1
  end           # end of time loop
end             # end of mpmMain()

mpmMain()
