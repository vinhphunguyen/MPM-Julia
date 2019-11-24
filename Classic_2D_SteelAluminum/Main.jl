# MPM implementation for the steel disk impact problem
# Sina Sinaie, Monash University
# July, 2016

# cd("E:\\MyPublications\\MPM_Julia\\Codes\\Classic_2D_SteelAluminum")
import PyPlot

pyFig_RealTime = PyPlot.figure("MPM Disk impact", figsize=(16/2.54, 16/2.54), edgecolor="white", facecolor="white")

include("./MyMaterialPoint.jl")
include("./MyGrid.jl")
include("./MyBasis.jl")

function mpmMain()
	const fGravity = 0.0
	# grid creation
    # nodes where fixation boundary conditions are hard coded in the following code!!!
	thisGrid = moduleGrid.mpmGrid(60.0, 60.0, 51, 51)

	# array holding all material points (these are references to MaterialDomain_01 & 02)
	allMaterialPoint = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)

	fOffset = 60.0/50/2.0
	thisMaterialDomain_01 = moduleMaterialPoint.createMaterialDomain_Circle([30.0; 50.0], 9.6/2.0, fOffset)
	for iIndex_MP = 1:1:length(thisMaterialDomain_01)
		fVolume = fOffset*fOffset#3.14159*0.2*0.2/length(thisMaterialDomain_01)
		fMass = 7850e-12*fVolume
		thisMaterialDomain_01[iIndex_MP].fMass = fMass
		thisMaterialDomain_01[iIndex_MP].fVolumeInitial = fVolume
		thisMaterialDomain_01[iIndex_MP].fVolume = fVolume
		thisMaterialDomain_01[iIndex_MP].fElasticModulus = 200.0e3
		thisMaterialDomain_01[iIndex_MP].fPoissonRatio = 0.3
		thisMaterialDomain_01[iIndex_MP].fYieldStress = 1.0e24
		thisMaterialDomain_01[iIndex_MP].v2Velocity = [0.0; -1160.0e3]
		thisMaterialDomain_01[iIndex_MP].v2Momentum = fMass*thisMaterialDomain_01[iIndex_MP].v2Velocity
		thisMaterialDomain_01[iIndex_MP].v2ExternalForce = [0.0; -fGravity*fMass]

		thisMaterialDomain_01[iIndex_MP].m22DeformationGradient = eye(2,2)
		thisMaterialDomain_01[iIndex_MP].m22DeformationGradientIncrement = eye(2,2)
	end
	for iIndex_MP = 1:1:length(thisMaterialDomain_01)
		push!(allMaterialPoint, thisMaterialDomain_01[iIndex_MP])
	end

	thisMaterialDomain_02 = moduleMaterialPoint.createMaterialDomain_Rectangle([30.0; 20.0], 60.0, 40.6, fOffset)
	for iIndex_MP = 1:1:length(thisMaterialDomain_02)
		fVolume = fOffset*fOffset#3.14159*0.2*0.2/length(thisMaterialDomain_02)
		fMass = 2700e-12*fVolume
		thisMaterialDomain_02[iIndex_MP].fMass = fMass
		thisMaterialDomain_02[iIndex_MP].fVolumeInitial = fVolume
		thisMaterialDomain_02[iIndex_MP].fVolume = fVolume
		thisMaterialDomain_02[iIndex_MP].fElasticModulus = 78.2e3
		thisMaterialDomain_02[iIndex_MP].fPoissonRatio = 0.3
		thisMaterialDomain_02[iIndex_MP].fYieldStress = 300.0
		thisMaterialDomain_02[iIndex_MP].v2Velocity = [0.0; 0.0]
		thisMaterialDomain_02[iIndex_MP].v2Momentum = fMass*thisMaterialDomain_02[iIndex_MP].v2Velocity
		thisMaterialDomain_02[iIndex_MP].v2ExternalForce = [0.0; -fGravity*fMass]

		thisMaterialDomain_02[iIndex_MP].m22DeformationGradient = eye(2,2)
		thisMaterialDomain_02[iIndex_MP].m22DeformationGradientIncrement = eye(2,2)
	end
	for iIndex_MP = 1:1:length(thisMaterialDomain_02)
		push!(allMaterialPoint, thisMaterialDomain_02[iIndex_MP])
	end

	# ---------------------------------------------------------------------------
	# information about the created domain
	fMass = 0.0
	for iIndex_MP in 1:1:length(allMaterialPoint)
		fMass += allMaterialPoint[iIndex_MP].fMass
	end
	@printf("Initial configuration: \n")
	@printf("	Single particle Mass: %+.6e \n", allMaterialPoint[1].fMass)
	@printf("	Total mass: %+.6e \n", fMass)

	@printf("	Disk, number of material points: %d \n", length(thisMaterialDomain_01))
	@printf("	Target, number of material points: %d \n", length(thisMaterialDomain_02))
	@printf("	Total number of material points: %d \n", length(allMaterialPoint))

	# ---------------------------------------------------------------------------
	# timers
	# ---------------------------------------------------------------------------
	# analysis timer
	fTimeIncrement = 1.0e-8
	fTimeEnd = 1.0e-4
	iTimeCycle = 0

	# realtime graphics timer
	fPlotTimeInterval = 100.0*fTimeIncrement
	fPlotTime = fPlotTimeInterval

	# final results plot timer
	fResultTimeInterval = fTimeIncrement#*fTimeEnd
	fResultTime = 0

	# console output timer
	fConsolTimeInterval = 0.1*fTimeEnd#1000.0*fTimeIncrement
	fConsolTime = fConsolTimeInterval

	# profiler timers
	fProfiler_MainLoop = 0.0
	fProfiler_Particle2Grid = 0.0
	fProfiler_Grid2Particle = 0.0
	# ---------------------------------------------------------------------------
	# plot arrays
	# ---------------------------------------------------------------------------
	# final results plot holder arrays
	fMarkedParicle_y = thisMaterialDomain_01[1].v2Centroid[2]
	plot_Time = Array{Real}(0)
	plot_Displacement = Array{Real}(0)
	plot_KineticEnergy = Array{Real}(0)
	plot_StrainEnergy = Array{Real}(0)

	# main analysis loop
	for fTime in 0.0:fTimeIncrement:fTimeEnd
		iTimeCycle += 1
		# ------------------------------------------------------------------------
		# realtime graphical plotting routines
		# @printf("Plotting...")
		# ------------------------------------------------------------------------
		fPlotTime += fTimeIncrement
		if(fPlotTime > fPlotTimeInterval)
			fPlotTime = 0.0

			iMaterialPoints = length(allMaterialPoint)
			array_x = [allMaterialPoint[i].v2Centroid[1] for i in 1:iMaterialPoints]
		   array_y = [allMaterialPoint[i].v2Centroid[2] for i in 1:iMaterialPoints]
		   array_color = Array{Real}(iMaterialPoints, 3)
			array_size = Array{Real}(iMaterialPoints, 1)
		   for iIndex in 1:1:iMaterialPoints
				thisColor = allMaterialPoint[iIndex].fAlpha
				thisColor /= (allMaterialPoint[iIndex].fYieldStress/allMaterialPoint[iIndex].fElasticModulus)*500

				if(thisColor > 1.0)
					# @printf("thiscolor %f\n", thisColor)
					thisColor = 1.0
				end

				if(allMaterialPoint[iIndex].fYieldStress > 1.0e3)
					array_color[iIndex, :] = [1.0, 0.0, 0.0]
				else
		      	array_color[iIndex, :] = [thisColor, 0.5*thisColor, 1.0-thisColor]#[thisGrid.GridPoints[iIndex].fMass/iMaterialPoints, 0.0, 0.0]
				end
				array_size[iIndex, :] = [4.0]
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
			pyPlot01[:set_xlim](0.0, 60.0)
			pyPlot01[:set_ylim](0.0, 60.0)
			# pyPlot01[:set_xlabel]("")
			# pyPlot01[:set_ylabel]("")
			pyPlot01[:grid](b=true, which="both", color="white", linestyle="-", linewidth=0.2)
			pyPlot01[:set_axisbelow](true)
			pyPlot01[:set_xticks]([])# empty to have no major ticks and grids
			pyPlot01[:set_xticks](collect(0.0:1.2:60.0),minor=true)
			pyPlot01[:set_yticks]([])# empty to have no major ticks and grids
			pyPlot01[:set_yticks](collect(0.0:1.2:60.0),minor=true)

			# PyPlot.show()
			# PyPlot.hold(true)

			strFileName = ".\\Figs\\DiskImpact_$(iTimeCycle).png"
			PyPlot.savefig(strFileName, bbox_inches="tight")
		   PyPlot.hold(false)
		end
		tic();
		#reset grid------------------------------------
		for iIndex in 1:1:thisGrid.iNodes
			thisGrid.GridPoints[iIndex].fMass = 0.0
			thisGrid.GridPoints[iIndex].v2Velocity = [0.0; 0.0]
			thisGrid.GridPoints[iIndex].v2Momentum = [0.0; 0.0]
			thisGrid.GridPoints[iIndex].v2Force = [0.0; 0.0]
		end
		# material to grid -------------------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(allMaterialPoint[iIndex_MP], thisGrid)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				# sina, be careful here, this might not be by reference and might not be good for assignment
				thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]

				fShapeValue = moduleBasis.getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				v2ShapeGradient = moduleBasis.getShapeGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				# mass
				thisGridPoint.fMass += fShapeValue * thisMaterialPoint.fMass
				# momentum
				thisGridPoint.v2Momentum += fShapeValue * thisMaterialPoint.fMass * thisMaterialPoint.v2Velocity
				# internal forces
				fVolume = thisMaterialPoint.fVolume
				thisGridPoint.v2Force[1] += -fVolume * (v2ShapeGradient[1]*thisMaterialPoint.v3Stress[1] + v2ShapeGradient[2]*thisMaterialPoint.v3Stress[3])
				thisGridPoint.v2Force[2] += -fVolume * (v2ShapeGradient[2]*thisMaterialPoint.v3Stress[2] + v2ShapeGradient[1]*thisMaterialPoint.v3Stress[3])
				# external forces
				thisGridPoint.v2Force += fShapeValue*thisMaterialPoint.v2ExternalForce
		   end
		end
		# update grid momentum and apply boundary conditions ---------------------
		for iIndex_GP in 1:1:thisGrid.iNodes
			thisGridPoint = thisGrid.GridPoints[iIndex_GP]

			thisGridPoint.v2Momentum += thisGridPoint.v2Force * fTimeIncrement

			if(thisGridPoint.v2Fixed[1] == true)
				thisGridPoint.v2Momentum[1] = 0.0
				thisGridPoint.v2Force[1] = 0.0
			end
			if(thisGridPoint.v2Fixed[2] == true)
				thisGridPoint.v2Momentum[2] = 0.0
				thisGridPoint.v2Force[2] = 0.0
			end
		end

		fProfiler_Particle2Grid += toq()

		tic()
		# ------------------------------------------------------------------------
		# grid to material pass 1-------------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			v2CentroidIncrement = zeros(2)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]

				fShapeValue = moduleBasis.getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
				v2ShapeGradient = moduleBasis.getShapeGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				thisMaterialPoint.v2Velocity += (fShapeValue * thisGridPoint.v2Force / thisGridPoint.fMass) * fTimeIncrement
		   end
		end
		# reset grid momenta -----------------------------------------------------
		# for iIndex_GP in 1:1:thisGrid.iNodes
		# 	thisGrid.GridPoints[iIndex_GP].v2Momentum = [0.0; 0.0]
		# end
		# map particle momenta back to grid------------------------------
		# mass in NOT mapped here
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]

				fShapeValue = moduleBasis.getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				# thisGridPoint.v2Momentum += fShapeValue * thisMaterialPoint.fMass * thisMaterialPoint.v2Velocity
				thisGridPoint.v2Velocity += fShapeValue * thisMaterialPoint.fMass * thisMaterialPoint.v2Velocity / thisGridPoint.fMass
			end
		end
		#apply boundary conditions velocity------------------------------------
		for iIndex_GP in 1:1:thisGrid.iNodes
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
		end
		# ------------------------------------------------------------------------
		# grid to material pass 2 ------------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint, thisGrid)
			v2CentroidIncrement = zeros(2)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]

				fShapeValue = moduleBasis.getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
				v2ShapeGradient = moduleBasis.getShapeGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

				v2GridPointVelocity = thisGridPoint.v2Velocity#v2Momentum / thisGridPoint.fMass

				v2CentroidIncrement += (fShapeValue * thisGridPoint.v2Momentum / thisGridPoint.fMass) * fTimeIncrement

				# from (2011) A convected particle domain interpolation technique to extend ...
				thisMaterialPoint.m22DeformationGradientIncrement += v2GridPointVelocity*transpose(v2ShapeGradient)*fTimeIncrement;
		   end
			thisMaterialPoint.v2Centroid += v2CentroidIncrement
			thisMaterialPoint.m22DeformationGradient = thisMaterialPoint.m22DeformationGradientIncrement * thisMaterialPoint.m22DeformationGradient
			v3StrainIncrement = zeros(3)
			v3StrainIncrement[1] = thisMaterialPoint.m22DeformationGradientIncrement[1,1] - 1.0
			v3StrainIncrement[2] = thisMaterialPoint.m22DeformationGradientIncrement[2,2] - 1.0
			v3StrainIncrement[3] = thisMaterialPoint.m22DeformationGradientIncrement[1,2] + thisMaterialPoint.m22DeformationGradientIncrement[2,1]
			thisMaterialPoint.m22DeformationGradientIncrement = eye(2,2)
			thisMaterialPoint.v3Strain[1] += v3StrainIncrement[1]
			thisMaterialPoint.v3Strain[2] += v3StrainIncrement[2]
			thisMaterialPoint.v3Strain[3] += v3StrainIncrement[3]

			fE = thisMaterialPoint.fElasticModulus;
			fNu = thisMaterialPoint.fPoissonRatio
			fYield = thisMaterialPoint.fYieldStress

			v3StrainCurrent = thisMaterialPoint.v3Strain
			v3StressCurrent = thisMaterialPoint.v3Stress
			v3PlasticStrainCurrent = thisMaterialPoint.v3PlasticStrain
			fAlphaCurrent = thisMaterialPoint.fAlpha

			v32Result = zeros(3,2)
			v32Result = moduleMaterialPoint.getIncrement_Plastic(fE, fNu, fYield, fAlphaCurrent, v3StressCurrent, v3StrainCurrent, v3PlasticStrainCurrent, v3StrainIncrement)

			v3StressIncrement = v32Result[:, 1]
			v3PlasticStrainIncrement = v32Result[:, 2]
			fAlphaIncrement = v32Result[1, 3]

			thisMaterialPoint.v3Stress += v3StressIncrement
			thisMaterialPoint.v3PlasticStrain += v3PlasticStrainIncrement
			thisMaterialPoint.fAlpha += fAlphaIncrement

			thisMaterialPoint.fVolume = det(thisMaterialPoint.m22DeformationGradient) * thisMaterialPoint.fVolumeInitial

			thisMaterialPoint.v2Momentum = thisMaterialPoint.v2Velocity * thisMaterialPoint.fMass
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
	end

	# final plots
	# pyFig_RealTime = PyPlot.figure("2Disk Impact FinalPlot", figsize=(8/2.54, 4/2.54))
	# PyPlot.clf()
	# pyPlot01 = PyPlot.gca()
	# PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	# pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	# pyPlot01[:set_axisbelow](true)
	# pyPlot01[:set_xlim](0.0, 4.0)
	# pyPlot01[:set_ylim](0.0, 3.0)
	# pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	# pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
	# pyPlot01[:set_xticks](collect(0.0:1.0:4.0))
	# pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	# pyPlot01[:set_yticks](collect(0.0:1.0:3.0))
	# PyPlot.plot(plot_Time, c="blue", plot_KineticEnergy, "-", label="\$ K \$", linewidth=1.0)
	# PyPlot.hold(true)
	# PyPlot.plot(plot_Time, c="red", plot_StrainEnergy, "-", label="\$ U \$", linewidth=1.0)
	# PyPlot.plot(plot_Time, c="green", plot_KineticEnergy + plot_StrainEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	# PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	# PyPlot.savefig("..\\..\\Figs\\plot_2Disk_Julia.pdf")
end

mpmMain()
