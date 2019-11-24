# MPM implementation for the 2disk impact problem
# Sina Sinaie, Monash University
# July, 2016

# cd("E:\\MyPublications\\MPM_Julia\\Codes\\Classic_2D_2Disk")
# import Gadfly
import PyPlot

pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time", figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

include("./MyMaterialPoint.jl")
include("./MyGrid.jl")
include("./MyBasis.jl")

function mpmMain()
    # problem parameters 
	const fGravity      = 0.0
	const density       = 1000.0
	const youngModulus  = 1000.0
	const poissonRatio  = 0.3

	# grid creation, Lx, Ly, # nodes x, # nodes y
	thisGrid = moduleGrid.mpmGrid(1.0, 1.0, 21, 21)

	# array holding all material points (these are references to MaterialDomain_01 & 02)
	allMaterialPoint = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)

    fOffset = 0.2/8.0 # there are 8 material points over the radius (16 MPs)
    # how to calculate fOffset:
    # length of 1 cell = 1/20=0.05 => there are 0.2/0.05=4 cells over the half of the
    # circle. If you want to MPs per cell (in x dir.), then there are 8 MPs there.
	thisMaterialDomain_01 = moduleMaterialPoint.createMaterialDomain_Circle([0.2; 0.2], 0.2, fOffset)
    # for iIndex_MP=1:domain1PCount, where domain1PCount=length(thisMaterialDomain_01)
    # does not significantly perform better!!!, so the following is fine
	for iIndex_MP = 1:length(thisMaterialDomain_01)
		fVolume = fOffset*fOffset#3.14159*0.2*0.2/length(thisMaterialDomain_01)
		fMass   = density*fVolume
		thisMaterialDomain_01[iIndex_MP].fMass           = fMass
		thisMaterialDomain_01[iIndex_MP].fVolumeInitial  = fVolume
		thisMaterialDomain_01[iIndex_MP].fVolume         = fVolume
		thisMaterialDomain_01[iIndex_MP].fElasticModulus = youngModulus
		thisMaterialDomain_01[iIndex_MP].fPoissonRatio   = poissonRatio
		thisMaterialDomain_01[iIndex_MP].v2Velocity      = [0.1; 0.1]
		thisMaterialDomain_01[iIndex_MP].v2Momentum      = fMass*thisMaterialDomain_01[iIndex_MP].v2Velocity
		thisMaterialDomain_01[iIndex_MP].v2ExternalForce = [0.0; -fGravity*fMass]

		thisMaterialDomain_01[iIndex_MP].m22DeformationGradient = eye(2,2)
		thisMaterialDomain_01[iIndex_MP].m22DeformationGradientIncrement = eye(2,2)

		push!(allMaterialPoint, thisMaterialDomain_01[iIndex_MP])
	end

	#for iIndex_MP = 1:1:length(thisMaterialDomain_01)
	#	push!(allMaterialPoint, thisMaterialDomain_01[iIndex_MP])
	#end

	thisMaterialDomain_02 = moduleMaterialPoint.createMaterialDomain_Circle([0.8; 0.8], 0.2, fOffset)
	for iIndex_MP = 1:1:length(thisMaterialDomain_02)
		fVolume = fOffset*fOffset#3.14159*0.2*0.2/length(thisMaterialDomain_02)
		fMass   = density*fVolume
		thisMaterialDomain_02[iIndex_MP].fMass           = fMass
		thisMaterialDomain_02[iIndex_MP].fVolumeInitial  = fVolume
		thisMaterialDomain_02[iIndex_MP].fVolume         = fVolume
		thisMaterialDomain_02[iIndex_MP].fElasticModulus = youngModulus
		thisMaterialDomain_02[iIndex_MP].fPoissonRatio   = poissonRatio
		thisMaterialDomain_02[iIndex_MP].v2Velocity = [-0.1; -0.1]
		thisMaterialDomain_02[iIndex_MP].v2Momentum = fMass*thisMaterialDomain_02[iIndex_MP].v2Velocity
		thisMaterialDomain_02[iIndex_MP].v2ExternalForce = [0.0; -fGravity*fMass]

		thisMaterialDomain_02[iIndex_MP].m22DeformationGradient = eye(2,2)
		thisMaterialDomain_02[iIndex_MP].m22DeformationGradientIncrement = eye(2,2)
		
        push!(allMaterialPoint, thisMaterialDomain_02[iIndex_MP])
	end
	#for iIndex_MP = 1:1:length(thisMaterialDomain_02)
	#	push!(allMaterialPoint, thisMaterialDomain_02[iIndex_MP])
	#end

	# ---------------------------------------------------------------------------
	# information about the created domain
	fMass = 0.0
	for iIndex_MP in 1:1:length(allMaterialPoint)
		fMass += allMaterialPoint[iIndex_MP].fMass
	end
	@printf("Initial configuration: \n")
	@printf("	Mass: %+.6e \n", fMass)
	@printf("Total number of material points: %d \n", length(allMaterialPoint))

	# ---------------------------------------------------------------------------
	# timers
	# ---------------------------------------------------------------------------
	# analysis timer
	fTimeIncrement = 1.0e-3
	fTimeEnd       = 3.5e-0

	# realtime graphics timer
	fPlotTimeInterval = 10.5*fTimeEnd#1000.0*fTimeIncrement
	fPlotTime = 0

	# final results plot timer
	fResultTimeInterval = fTimeIncrement#*fTimeEnd
	fResultTime = 0

	# console output timer
	fConsolTimeInterval = 0.1*fTimeEnd#1000.0*fTimeIncrement
	fConsolTime = fConsolTimeInterval

	# profiler timers
	fProfiler_MainLoop      = 0.0
	fProfiler_Particle2Grid = 0.0
	fProfiler_Grid2Particle = 0.0
	# ---------------------------------------------------------------------------
	# plot arrays
	# ---------------------------------------------------------------------------
	# final results plot holder arrays
	fMarkedParicle_y   = thisMaterialDomain_01[1].v2Centroid[2]
	plot_Time          = Array{Real}(0) #??? Float64??
	plot_Displacement  = Array{Real}(0)
	plot_KineticEnergy = Array{Real}(0)
	plot_StrainEnergy  = Array{Real}(0)

	# main analysis loop
	for fTime in 0.0:fTimeIncrement:fTimeEnd
		tic();
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

		fProfiler_Particle2Grid += toq()

		tic()
		# ------------------------------------------------------------------------
		# grid to material -------------------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint      = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			v2CentroidIncrement = zeros(2)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				thisGridPoint   = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]
				#fShapeValue     = moduleBasis.getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
				#v2ShapeGradient = moduleBasis.getShapeGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
				fShapeValue, v2ShapeGradient  = moduleBasis.getShapeAndGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				v2GridPointVelocity = zeros(2)
				if(thisGridPoint.fMass > 1.0e-8)
					v2GridPointVelocity = thisGridPoint.v2Momentum / thisGridPoint.fMass
					thisMaterialPoint.v2Velocity += (fShapeValue * thisGridPoint.v2Force / thisGridPoint.fMass) * fTimeIncrement
					# thisMaterialPoint.v2Centroid += (fShapeValue * thisGridPoint.v2Momentum / thisGridPoint.fMass) * fTimeIncrement
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

			fE        = thisMaterialPoint.fElasticModulus;
			fNu       = thisMaterialPoint.fPoissonRatio
			fConstant = fE/(1.0 + fNu)/(1.0 - 2.0*fNu)

			thisMaterialPoint.v3Strain[1] += v3StrainIncrement[1]
			thisMaterialPoint.v3Strain[2] += v3StrainIncrement[2]
			thisMaterialPoint.v3Strain[3] += v3StrainIncrement[3]


			thisMaterialPoint.v3Stress[1] += fConstant * ((1.0-fNu)*v3StrainIncrement[1] + fNu*v3StrainIncrement[2])
			thisMaterialPoint.v3Stress[2] += fConstant * ((1.0-fNu)*v3StrainIncrement[2] + fNu*v3StrainIncrement[1])
			thisMaterialPoint.v3Stress[3] += fConstant * ((0.5-fNu)*v3StrainIncrement[3])

			thisMaterialPoint.fVolume      = det(thisMaterialPoint.m22DeformationGradient) * thisMaterialPoint.fVolumeInitial
			thisMaterialPoint.v2Momentum   = thisMaterialPoint.v2Velocity * thisMaterialPoint.fMass
		end

		fProfiler_Grid2Particle += toq()

		# ------------------------------------------------------------------------
		# calculating strain and kinetic energy for final results plot
		# ------------------------------------------------------------------------
		fResultTime += fTimeIncrement
		if(fResultTime > fResultTimeInterval)
			fResultTime = 0.0
			fStrainEnergy = 0.0
			fKineticEnergy = 0.0
			for iIndex_MP in 1:1:length(allMaterialPoint)
				thisMaterialPoint = allMaterialPoint[iIndex_MP]
				fStrainEnergy += 0.5*thisMaterialPoint.v3Strain[1] * thisMaterialPoint.v3Stress[1] * thisMaterialPoint.fVolume
				fStrainEnergy += 0.5*thisMaterialPoint.v3Strain[2] * thisMaterialPoint.v3Stress[2] * thisMaterialPoint.fVolume
				fStrainEnergy += 0.5*thisMaterialPoint.v3Strain[3] * thisMaterialPoint.v3Stress[3] * thisMaterialPoint.fVolume

				fVelocity = thisMaterialPoint.v2Velocity[1]^2 + thisMaterialPoint.v2Velocity[2]^2
				fKineticEnergy += 0.5*fVelocity * thisMaterialPoint.fMass
			end

			#save to plot arrays
			push!(plot_Time, fTime)
			push!(plot_KineticEnergy, fKineticEnergy)
			push!(plot_StrainEnergy, fStrainEnergy)
		end
		# ------------------------------------------------------------------------
		# consol output
		# ------------------------------------------------------------------------
		fConsolTime += fTimeIncrement
		if(fConsolTime > fConsolTimeInterval)
			fConsolTime = 0.0

			fMass = 0.0
			fMomentum_x = 0.0
			for iIndex_MP in 1:1:length(thisMaterialDomain_01)
				fMass += thisMaterialDomain_01[iIndex_MP].fMass

				fMomentum_x += thisMaterialDomain_01[iIndex_MP].v2Momentum[1]
			end
			fProfiler_Total = fProfiler_Particle2Grid + fProfiler_Grid2Particle
			@printf("fTime: %+.3e |", fTime)
			@printf("M_x: %+.3e |", fMomentum_x)
			@printf("(Profiler) Total: %+.3e ", fProfiler_Total)
			@printf("P2G: %+.3e (%+.2f) ", fProfiler_Particle2Grid, fProfiler_Particle2Grid/fProfiler_Total)
			@printf("G2P: %+.3e (%+.2f) \n", fProfiler_Grid2Particle, fProfiler_Grid2Particle/fProfiler_Total)
		end
		# ------------------------------------------------------------------------
		# realtime graphical plotting routines
		# @printf("Plotting...")
		# ------------------------------------------------------------------------
		fPlotTime += fTimeIncrement
		if(fPlotTime > fPlotTimeInterval)
			fPlotTime       = 0.0
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
			pyPlot01[:spines]["top"][:set_color]("gray")
			pyPlot01[:spines]["right"][:set_color]("gray")
			pyPlot01[:spines]["bottom"][:set_color]("gray")
			pyPlot01[:spines]["left"][:set_color]("gray")
			# pyPlot01[:axhline](linewidth=4, color="g")
			# pyPlot01[:axvline](linewidth=4, color="g")
			pyPlot01[:set_xlim](0.0, 1.0)
			pyPlot01[:set_ylim](0.0, 1.0)
			# pyPlot01[:set_xlabel]("")
			# pyPlot01[:set_ylabel]("")
			pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
			pyPlot01[:set_axisbelow](true)
			pyPlot01[:set_xticks]([])# empty to have no major ticks and grids
			pyPlot01[:set_xticks](collect(0.0:0.05:1.0),minor=true)
			pyPlot01[:set_yticks]([])# empty to have no major ticks and grids
			pyPlot01[:set_yticks](collect(0.0:0.05:1.0),minor=true)

			# PyPlot.show()
			# PyPlot.hold(true)

			strFileName = ".\\Figs\\2Disk_$(Int(fTime*1000)).png"
			PyPlot.savefig(strFileName, bbox_inches="tight")
		   # PyPlot.hold(false)
		end
	end

	# final plots
	pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
	PyPlot.clf()
	pyPlot01 = PyPlot.gca()
	PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	pyPlot01[:set_xlim](0.0, 4.0)
	pyPlot01[:set_ylim](0.0, 3.0)
	pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
	pyPlot01[:set_xticks](collect(0.0:1.0:4.0))
	pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	pyPlot01[:set_yticks](collect(0.0:1.0:3.0))
	PyPlot.plot(plot_Time, c="blue", plot_KineticEnergy, "-", label="\$ K \$", linewidth=1.0)
	PyPlot.hold(true)
	PyPlot.plot(plot_Time, c="red", plot_StrainEnergy, "-", label="\$ U \$", linewidth=1.0)
	PyPlot.plot(plot_Time, c="green", plot_KineticEnergy + plot_StrainEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	PyPlot.savefig("..\\..\\Figs\\plot_2Disk_Julia.pdf")
end

mpmMain()
