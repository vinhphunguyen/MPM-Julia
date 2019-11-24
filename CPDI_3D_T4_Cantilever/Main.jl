# MPM implementation for the cantilever problem
# Sina Sinaie, Monash University
# July, 2016

# cd("E:\\MyPublications\\MPM_Julia\\Codes\\CPDI_3D_T4")
# import Gadfly
import PyPlot

pyFig_RealTime = PyPlot.figure("MPM 3D", figsize=(18/2.54, 18/2.54), edgecolor="white", facecolor="white")

include("./MyMaterialPoint.jl")
include("./MyGrid.jl")
include("./MyBasis.jl")

function mpmMain()
	const fGravity = 10.0
	# grid creation
	thisGrid = moduleGrid.mpmGrid(8.0, 2.0, 8.0, 17, 3, 17)

	# array holding all material points (these are references to MaterialDomain_01 & 02)
	allMaterialPoint = Array{moduleMaterialPoint.mpmMaterialPoint_Tet4}(0)

	thisMaterialDomain = moduleMaterialPoint.loadMSH("beam.msh")
	for iIndex_MP = 1:1:length(thisMaterialDomain)
		fVolume = moduleMaterialPoint.getVolume(thisMaterialDomain[iIndex_MP])
		fMass = 1050.0*fVolume
		# fMass = 1050.0*1.0*1.0*4.0/length(thisMaterialDomain)
		# fVolume = 1.0*1.0*4.0/length(thisMaterialDomain)
		thisMaterialDomain[iIndex_MP].fMass = fMass
		thisMaterialDomain[iIndex_MP].fVolumeInitial = fVolume
		thisMaterialDomain[iIndex_MP].fVolume = fVolume
		thisMaterialDomain[iIndex_MP].fElasticModulus = 1.0e6
		thisMaterialDomain[iIndex_MP].fPoissonRatio = 0.3
		thisMaterialDomain[iIndex_MP].v3Centroid += [0.5; 0.5; 3.5]
		for index_Corner = 1:1:4
			thisMaterialDomain[iIndex_MP].v3Corner[:,index_Corner] += [0.5; 0.5; 3.5]
		end
		thisMaterialDomain[iIndex_MP].v3Velocity = [0.0; 0.0; 0.0]
		thisMaterialDomain[iIndex_MP].v3Momentum = fMass*thisMaterialDomain[iIndex_MP].v3Velocity
		thisMaterialDomain[iIndex_MP].v3ExternalForce = [0.0; 0.0; -fGravity*fMass]

		thisMaterialDomain[iIndex_MP].m33DeformationGradient = eye(3,3)
		thisMaterialDomain[iIndex_MP].m33DeformationGradientIncrement = eye(3,3)
	end
	for iIndex_MP = 1:1:length(thisMaterialDomain)
		# if(thisMaterialDomain[iIndex_MP].v3Centroid[3] > 4.0 && thisMaterialDomain[iIndex_MP].v3Centroid[2] > 1.0)
			push!(allMaterialPoint, thisMaterialDomain[iIndex_MP])
		# end
	end
	for iIndex_MP = 1:1:length(allMaterialPoint)
		for iIndex_Corner = 1:1:4
			v3Coordinate = allMaterialPoint[iIndex_MP].v3Corner[:, iIndex_Corner]
			if(abs(v3Coordinate[1]-4.5) < 1e-10 && abs(v3Coordinate[2]-0.5) < 1e-10 && abs(v3Coordinate[3]-3.5) < 1e-10)
				@printf("Corner Mark: %d, %d\n", iIndex_MP, iIndex_Corner)
			end
		end
	end
	# ---------------------------------------------------------------------------
	# information about the created domain
	fMass = 0.0
	for iIndex_MP in 1:1:length(allMaterialPoint)
		fMass += allMaterialPoint[iIndex_MP].fMass
	end
	@printf("Initial configuration: \n")
	@printf("	Element count: %d \n", length(allMaterialPoint))
	@printf("	Mass: %+.6e \n", fMass)

	# @printf("Total number of material points: %d \n", length(allMaterialPoint))

	# ---------------------------------------------------------------------------
	# timers
	# ---------------------------------------------------------------------------
	# analysis timer
	fTimeIncrement = 1.0e-3
	fTimeEnd = 3.0e-0

	# realtime graphics timer
	fPlotTimeInterval = 10.6#*fTimeEnd#1000.0*fTimeIncrement
	fPlotTime = fPlotTimeInterval

	# final results plot timer
	fResultTimeInterval = fTimeIncrement
	fResultTime = 0

	# console output timer
	fConsolTimeInterval = 0.1*fTimeEnd
	fConsolTime = fConsolTimeInterval

	# profiler timers
	fProfiler_MainLoop = 0.0
	fProfiler_Particle2Grid = 0.0
	fProfiler_Grid2Particle = 0.0
	fProfiler_Corners = 0.0
	# ---------------------------------------------------------------------------
	# plot arrays
	# ---------------------------------------------------------------------------
	# final results plot holder arrays
	fMarkedParicle_y = thisMaterialDomain[39].v3Centroid[3]
	plot_Time = Array{Real}(0)
	plot_Displacement = Array{Real}(0)
	plot_KineticEnergy = Array{Real}(0)
	plot_StrainEnergy = Array{Real}(0)

	# main analysis loop
	for fTime in 0.0:fTimeIncrement:fTimeEnd
		#save to plot arrays
		push!(plot_Time, fTime)
		push!(plot_Displacement, allMaterialPoint[918].v3Corner[3,4] - 3.5)#2.4433-1.5)#sina, remember 39

		tic();
		#reset grid------------------------------------
		for iIndex in 1:1:thisGrid.iNodes
			thisGrid.GridPoints[iIndex].fMass = 0.0
			thisGrid.GridPoints[iIndex].v3Momentum = [0.0; 0.0; 0.0]
			thisGrid.GridPoints[iIndex].v3Force = [0.0; 0.0; 0.0]
		end
		# material to grid -------------------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints_CPDI(allMaterialPoint[iIndex_MP], thisGrid)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				# sina, be careful here, this might not be by reference and might not be good for assignment
				thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]

				fShapeValue = moduleBasis.getShapeValue_CPDI_Tet4(thisMaterialPoint, thisGridPoint, thisGrid)
				v3ShapeGradient = moduleBasis.getShapeGradient_CPDI_Tet4(thisMaterialPoint, thisGridPoint, thisGrid)

				# mass
				thisGridPoint.fMass += fShapeValue * thisMaterialPoint.fMass
				# momentum
				thisGridPoint.v3Momentum += fShapeValue * thisMaterialPoint.v3Momentum
				# internal forces
				fVolume = thisMaterialPoint.fVolume
				thisGridPoint.v3Force[1] += -fVolume * (v3ShapeGradient[1]*thisMaterialPoint.v6Stress[1] + v3ShapeGradient[2]*thisMaterialPoint.v6Stress[4] + v3ShapeGradient[3]*thisMaterialPoint.v6Stress[6])
				thisGridPoint.v3Force[2] += -fVolume * (v3ShapeGradient[2]*thisMaterialPoint.v6Stress[2] + v3ShapeGradient[1]*thisMaterialPoint.v6Stress[4] + v3ShapeGradient[3]*thisMaterialPoint.v6Stress[5])
				thisGridPoint.v3Force[3] += -fVolume * (v3ShapeGradient[3]*thisMaterialPoint.v6Stress[3] + v3ShapeGradient[1]*thisMaterialPoint.v6Stress[6] + v3ShapeGradient[2]*thisMaterialPoint.v6Stress[5])
				# # external forces
				thisGridPoint.v3Force += fShapeValue*thisMaterialPoint.v3ExternalForce
		   end
		end
		# update grid momentum and apply boundary conditions ---------------------
		for iIndex_GP in 1:1:thisGrid.iNodes
			thisGridPoint = thisGrid.GridPoints[iIndex_GP]

			thisGridPoint.v3Momentum += thisGridPoint.v3Force * fTimeIncrement

			if(thisGridPoint.v3Fixed[1] == true)
				thisGridPoint.v3Momentum[1] = 0.0
				thisGridPoint.v3Force[1] = 0.0
			end
			if(thisGridPoint.v3Fixed[2] == true)
				thisGridPoint.v3Momentum[2] = 0.0
				thisGridPoint.v3Force[2] = 0.0
			end
			if(thisGridPoint.v3Fixed[3] == true)
				thisGridPoint.v3Momentum[3] = 0.0
				thisGridPoint.v3Force[3] = 0.0
			end
		end

		fProfiler_Particle2Grid += toq()

		tic()
		# ------------------------------------------------------------------------
		# grid to material -------------------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint = allMaterialPoint[iIndex_MP]
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints_CPDI(thisMaterialPoint, thisGrid)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]

				fShapeValue = moduleBasis.getShapeValue_CPDI_Tet4(thisMaterialPoint, thisGridPoint, thisGrid)
				v3ShapeGradient = moduleBasis.getShapeGradient_CPDI_Tet4(thisMaterialPoint, thisGridPoint, thisGrid)

				v3GridPointVelocity = zeros(3)
				if(thisGridPoint.fMass > 1.0e-16)
					v3GridPointVelocity = thisGridPoint.v3Momentum / thisGridPoint.fMass

					thisMaterialPoint.v3Velocity += (fShapeValue * thisGridPoint.v3Force / thisGridPoint.fMass) * fTimeIncrement
				end

				thisMaterialPoint.m33DeformationGradientIncrement += v3GridPointVelocity*transpose(v3ShapeGradient)*fTimeIncrement;
		   end

			thisMaterialPoint.m33DeformationGradient = thisMaterialPoint.m33DeformationGradientIncrement * thisMaterialPoint.m33DeformationGradient
			thisMaterialPoint.m33DeformationGradientIncrement = eye(3,3)

			fE = thisMaterialPoint.fElasticModulus;
			fNu = thisMaterialPoint.fPoissonRatio

			# for a Neo_Hookean model, (2011) A convected particle domain interpolation technique to extend
			fLame1 = 0.5*fE / (1.0 + fNu)
			fLame2 = fNu*fE / ((1+fNu)*(1.0-2.0*fNu))

			mF = thisMaterialPoint.m33DeformationGradient
			fJ = det(mF)

			mP = fLame2*log(fJ)/fJ * eye(3) + fLame1/fJ*(mF*mF'-eye(3))
			thisMaterialPoint.v6Stress[1] = mP[1,1]
			thisMaterialPoint.v6Stress[2] = mP[2,2]
			thisMaterialPoint.v6Stress[3] = mP[3,3]
			thisMaterialPoint.v6Stress[4] = mP[1,2]
			thisMaterialPoint.v6Stress[5] = mP[2,3]
			thisMaterialPoint.v6Stress[6] = mP[3,1]

			thisMaterialPoint.v3Momentum = thisMaterialPoint.v3Velocity * thisMaterialPoint.fMass
		end
		fProfiler_Grid2Particle += toq()

		tic()
		# calculate corner increments---------------------------------------------
		for iIndex_MP in 1:1:length(allMaterialPoint)
			thisMaterialPoint = allMaterialPoint[iIndex_MP]
			thisMaterialPoint.v3Centroid = zeros(3) # reset, to be calculated shortly
			for indexCorner = 1:1:size(thisMaterialPoint.v3Corner, 2)
				# v3CornerCoordinate = thisMaterialPoint.v3Corner[:, indexCorner]
				thisCorner = thisMaterialPoint.v3Corner[:, indexCorner]
				v3CornerIncrement = zeros(3)

				thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisCorner, thisGrid)
				for iIndex in 1:1:length(thisAdjacentGridPoints)
					thisGridPoint = thisGrid.GridPoints[thisAdjacentGridPoints[iIndex]]

					v3Distance = thisCorner - thisGridPoint.v3Position
					v3CellLength = thisGrid.v3Length_Cell

					fShapeValue = moduleBasis.getShapeValue_Classic(v3Distance, v3CellLength)

					v3GridPointVelocity = zeros(3)
					if(thisGridPoint.fMass > 1.0e-16)
						v3GridPointVelocity = thisGridPoint.v3Momentum / thisGridPoint.fMass
					end

					v3CornerIncrement += fShapeValue * v3GridPointVelocity * fTimeIncrement
					# thisMaterialPoint.v3Corner[:, indexCorner] += fShapeValue * v3GridPointVelocity * fTimeIncrement
				end
				thisMaterialPoint.v3Corner[:, indexCorner] += v3CornerIncrement

				thisMaterialPoint.v3Centroid += 1.0/4.0 * thisMaterialPoint.v3Corner[:, indexCorner] # 4.0 for 4 corners, the sum will the average of all corners
			end

			thisMaterialPoint.fVolume = det(thisMaterialPoint.m33DeformationGradient) * thisMaterialPoint.fVolumeInitial
		end

		fProfiler_Corners += toq()

		# ------------------------------------------------------------------------
		# consol output
		# ------------------------------------------------------------------------
		fConsolTime += fTimeIncrement
		if(fConsolTime > fConsolTimeInterval)
			fConsolTime = 0.0

			fProfiler_Total = fProfiler_Particle2Grid + fProfiler_Grid2Particle + fProfiler_Corners
			@printf("fTime: %+.3e |", fTime)
			@printf("Displacement: %+.3e ", thisMaterialDomain[39].v3Centroid[3] - fMarkedParicle_y)
			@printf("(Profiler) Total: %+.3e ", fProfiler_Total)
			# @printf("P2G: %+.3e (%+.2f) ", fProfiler_Particle2Grid, fProfiler_Particle2Grid/fProfiler_Total)
			# @printf("G2P: %+.3e (%+.2f) \n", fProfiler_Grid2Particle, fProfiler_Grid2Particle/fProfiler_Total)
			# @printf("Cor: %+.3e (%+.2f) \n", fProfiler_Corners, fProfiler_Corners/fProfiler_Total)
			@printf(" \n")
		end
		# ------------------------------------------------------------------------
		# realtime graphical plotting routines
		# ------------------------------------------------------------------------
		fPlotTime += fTimeIncrement
		if(fPlotTime > fPlotTimeInterval)
			fPlotTime = 0.0

			# array_x = [thisGrid.GridPoints[i].v3Position[1] for i in 1:thisGrid.iNodes]
		   # array_y = [thisGrid.GridPoints[i].v3Position[2] for i in 1:thisGrid.iNodes]
			# array_z = [thisGrid.GridPoints[i].v3Position[3] for i in 1:thisGrid.iNodes]
			array_x = Array{Float64}(0)
		   array_y = Array{Float64}(0)
			array_z = Array{Float64}(0)
		   array_color = Array{Float64}(0,3)
			array_size = Array{Float64}(0)
		   for iIndex in 1:1:thisGrid.iNodes
		      # array_color[iIndex, :] = [0.0, 0.5, 0.0]#[thisGrid.GridPoints[iIndex].fMass/iMaterialPoints, 0.0, 0.0]
				# if(iIndex == 89 || iIndex == 90 || iIndex == 110 || iIndex == 111)
				# 	array_color[iIndex, :] = [0.0, 1.0, 0.0]
				# end
				if(thisGrid.GridPoints[iIndex].v3Fixed[1] == true || thisGrid.GridPoints[iIndex].v3Fixed[2] == true || thisGrid.GridPoints[iIndex].v3Fixed[3] == true)
					array_color = vcat(array_color, [0.0 0.0 0.0])
					push!(array_size, 20.0)
					push!(array_x, thisGrid.GridPoints[iIndex].v3Position[1])
					push!(array_y, thisGrid.GridPoints[iIndex].v3Position[2])
					push!(array_z, thisGrid.GridPoints[iIndex].v3Position[3])
				end
				# println(array_color)
		   end

			pyPlot01 = PyPlot.gca(projection="3d")
			# pyPlot01 = PyPlot.subplot2grid((5,5), (0,0), colspan=4, rowspan=5, projection="3d")
			scat = PyPlot.scatter3D(array_x, array_y, array_z, lw=0, c=array_color, s=array_size)
			PyPlot.hold(true)
			# pyPlot01[:set_dist](1000.0)
			pyPlot01[:w_xaxis][:line][:set_color]((1.0, 1.0, 1.0, 0.0))
			pyPlot01[:w_yaxis][:line][:set_color]((1.0, 1.0, 1.0, 0.0))
			pyPlot01[:w_zaxis][:line][:set_color]((1.0, 1.0, 1.0, 0.0))
			pyPlot01[:w_xaxis][:gridlines][:set_color]([1.0, 0.0, 0.0, 1.0])
			pyPlot01[:w_yaxis][:gridlines][:set_color]((1.0, 0.0, 0.0, 1.0))
			pyPlot01[:w_zaxis][:gridlines][:set_color]((1.0, 0.0, 0.0, 1.0))
			# pyPlot01[:w_xaxis][:gridlines][:set_lw](10.0)
			pyPlot01[:set_xlim](0.0, 5.0)
			pyPlot01[:set_ylim](-3.0, 5.0)
			pyPlot01[:set_zlim](1.5, 6.5)
			# pyPlot01[:set_xlabel]("X")
			# pyPlot01[:set_ylabel]("Y")
			# pyPlot01[:set_zlabel]("Z")
			pyPlot01[:grid](b=true, which="major", color=[1.0; 0.0; 0.0], linestyle="--", linewidth=10.5)
			# pyPlot01[:set_axis_off]()
			# pyPlot01[:set_xticks](collect(0.0:thisGrid.v3Length_Cell[1]:thisGrid.v3Length_Grid[1]))# empty to have no major ticks and grids
			# pyPlot01[:set_yticks](collect(0.0:thisGrid.v3Length_Cell[2]:thisGrid.v3Length_Grid[2]))# empty to have no major ticks and grids
			# pyPlot01[:set_zticks](collect(0.0:thisGrid.v3Length_Cell[3]:thisGrid.v3Length_Grid[3]))# empty to have no major ticks and grids
			pyPlot01[:set_xticks]([])# empty to have no major ticks and grids
			pyPlot01[:set_yticks]([])# empty to have no major ticks and grids
			pyPlot01[:set_zticks]([])# empty to have no major ticks and grids
			pyPlot01[:set_xticklabels]([])# empty to have no major ticks and grids
			pyPlot01[:set_yticklabels]([])# empty to have no major ticks and grids
			pyPlot01[:set_zticklabels]([])# empty to have no major ticks and grids
			# pyPlot01[:tick_params](axis='x',length=0)
			# pyPlot01[:set_view](distance=1000.0)
			# pyPlot01[:view_init](elev=30.0, azim=-60.0)
			pyPlot01[:view_init](30,-60.0)
			# pyPlot01[:view_init](10.0, -70.0)

			# plot all material points (3d tet4)
			# @printf("Plotting...")
			for iIndex_MP in 1:1:size(allMaterialPoint, 1)
				thisColor = [1.0, 0.0, 0.0]
				thisSize = 1.0

				# if(iIndex_MP == 2)
				# 	thisColor = [0.0, 0.0, 0.0]
				# 	thisSize = 50.0
				# end

				nCorners = size(allMaterialPoint[iIndex_MP].v3Corner, 2)
				# println("Corner: ", thisMaterialDomain[iIndex_MP].v3Corner[1,:])
				# thisMarker = Array{Real}(nCorners+4,3) # the marker is draw over 7 lines
				thisMarker_x = zeros(nCorners+4)#Array{Real}(nCorners+4,3) # the marker is draw over 7 lines
				thisMarker_y = zeros(nCorners+4)
				thisMarker_z = zeros(nCorners+4)

				thisMarker_x[1] = allMaterialPoint[iIndex_MP].v3Corner[1,1]; thisMarker_y[1] = allMaterialPoint[iIndex_MP].v3Corner[2,1]; thisMarker_z[1] = allMaterialPoint[iIndex_MP].v3Corner[3,1]
				thisMarker_x[2] = allMaterialPoint[iIndex_MP].v3Corner[1,2]; thisMarker_y[2] = allMaterialPoint[iIndex_MP].v3Corner[2,2]; thisMarker_z[2] = allMaterialPoint[iIndex_MP].v3Corner[3,2]
				thisMarker_x[3] = allMaterialPoint[iIndex_MP].v3Corner[1,3]; thisMarker_y[3] = allMaterialPoint[iIndex_MP].v3Corner[2,3]; thisMarker_z[3] = allMaterialPoint[iIndex_MP].v3Corner[3,3]
				thisMarker_x[4] = allMaterialPoint[iIndex_MP].v3Corner[1,1]; thisMarker_y[4] = allMaterialPoint[iIndex_MP].v3Corner[2,1]; thisMarker_z[4] = allMaterialPoint[iIndex_MP].v3Corner[3,1]
				thisMarker_x[5] = allMaterialPoint[iIndex_MP].v3Corner[1,4]; thisMarker_y[5] = allMaterialPoint[iIndex_MP].v3Corner[2,4]; thisMarker_z[5] = allMaterialPoint[iIndex_MP].v3Corner[3,4]
				thisMarker_x[6] = allMaterialPoint[iIndex_MP].v3Corner[1,3]; thisMarker_y[6] = allMaterialPoint[iIndex_MP].v3Corner[2,3]; thisMarker_z[6] = allMaterialPoint[iIndex_MP].v3Corner[3,3]
				thisMarker_x[7] = allMaterialPoint[iIndex_MP].v3Corner[1,2]; thisMarker_y[7] = allMaterialPoint[iIndex_MP].v3Corner[2,2]; thisMarker_z[7] = allMaterialPoint[iIndex_MP].v3Corner[3,2]
				thisMarker_x[8] = allMaterialPoint[iIndex_MP].v3Corner[1,4]; thisMarker_y[8] = allMaterialPoint[iIndex_MP].v3Corner[2,4]; thisMarker_z[8] = allMaterialPoint[iIndex_MP].v3Corner[3,4]

				PyPlot.plot(thisMarker_x, thisMarker_y, thisMarker_z, c=thisColor, linewidth=0.1)
				# PyPlot.scatter3D([allMaterialPoint[iIndex_MP].v3Centroid[1]], [allMaterialPoint[iIndex_MP].v3Centroid[2]], [allMaterialPoint[iIndex_MP].v3Centroid[3]], c=thisColor, lw = 0, s=thisSize, marker="o")
			end
			# marked corner
			# PyPlot.scatter3D([allMaterialPoint[918].v3Corner[1,4]], [allMaterialPoint[918].v3Corner[2,4]], [allMaterialPoint[918].v3Corner[3,4]], c=[0.0, 1.0, 0.0], lw = 0, s=40, marker="o")
			PyPlot.hold(false)

			strFileName = ".\\Figs\\Cantelever_$(Int(fTime*1000)).png"
			# PyPlot.tight_layout(pad=0)
			PyPlot.savefig(strFileName, bbox_inches="tight")

			# @printf("done\n")

			# pyPlot02 = PyPlot.subplot2grid((5,5), (0,4), colspan=1, rowspan=2)
			# PyPlot.plot(plot_Time, plot_Displacement, "-")
			# PyPlot.grid("on")
			#
		   # PyPlot.pause(0.01)
		   # PyPlot.draw()
		   # PyPlot.hold(false)
		end
	end

	# final plots
	# pyFig_FinalResults = PyPlot.figure("Time-Displacement", figsize=(8.0/2.54, 6/2.54))
	# PyPlot.clf()
	# pyPlot01 = PyPlot.gca()
	# PyPlot.subplots_adjust(left=0.2, bottom=0.15, right=0.9)
	# pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	# pyPlot01[:set_axisbelow](true)
	# pyPlot01[:set_xlim](0.0, 3.0)
	# pyPlot01[:set_ylim](-3.5, 0.0)
	# pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	# pyPlot01[:set_ylabel]("displacement (m)", fontsize=8)
	# pyPlot01[:set_xticks](collect(0.0:0.5:3.0))
	# pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	# pyPlot01[:set_yticks](collect(0.0:-0.5:-3.5))
	#
	# PyPlot.hold(true)
	# PyPlot.plot(plot_Time, plot_Displacement, "-", linewidth=1.0)
	#
	# PyPlot.savefig("..\\..\\Figs\\plot_Cantilever_Julia.pdf")
end

mpmMain()
