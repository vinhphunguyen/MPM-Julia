# cd("Julia/MPM")
# cd("E:\\MyCodes\\Julia\\CPDI_T3")
# import Gadfly
import PyPlot

pyFig_RealTime = PyPlot.figure("MPM Two-Disk Impact", figsize=(8, 8))

include("./MyMath.jl")
include("./MyMaterialPoint.jl")
include("./MyGrid.jl")
include("./MyBasis.jl")

function mpmMain()
	const fGravity = 1000.0
	# grid creation
	thisGrid = moduleGrid.mpmGrid(4.0, 10.0, 9, 21)

	# array holding all material points (these are references to MaterialDomain_01 & 02)
	thisMaterialPoint = Array{moduleMaterialPoint.mpmMaterialPoint}(0)

	fOffset = 1.0/6.0
	thisMaterialDomain_03 = moduleMaterialPoint.createMaterialDomain_Rectangle("triangle", 2.0, 9.0, 1.0, 1.0, fOffset)
	for iIndex_MP in 1:1:length(thisMaterialDomain_03)
		fMass = 1050.0*1.0*1.0/length(thisMaterialDomain_03)
		fVolume = 1.0*1.0/length(thisMaterialDomain_03)
		thisMaterialDomain_03[iIndex_MP].fMass = fMass
		thisMaterialDomain_03[iIndex_MP].fVolumeInitial = fVolume
		thisMaterialDomain_03[iIndex_MP].fVolume = fVolume
		thisMaterialDomain_03[iIndex_MP].v2Length.fx = fOffset
		thisMaterialDomain_03[iIndex_MP].v2Length.fy = fOffset
		thisMaterialDomain_03[iIndex_MP].fElasticModulus = 1.0e6
		thisMaterialDomain_03[iIndex_MP].fPoissonRatio = 0.3
		thisMaterialDomain_03[iIndex_MP].v2Velocity.fx = 0.0
		thisMaterialDomain_03[iIndex_MP].v2Velocity.fy = 0.0
		thisMaterialDomain_03[iIndex_MP].v2Momentum.fx = thisMaterialDomain_03[iIndex_MP].v2Velocity.fx*fMass
		thisMaterialDomain_03[iIndex_MP].v2Momentum.fy = thisMaterialDomain_03[iIndex_MP].v2Velocity.fy*fMass
		thisMaterialDomain_03[iIndex_MP].v2ExternalForce.fx = 0.0
		thisMaterialDomain_03[iIndex_MP].v2ExternalForce.fy = -fGravity * fMass

		thisMaterialDomain_03[iIndex_MP].mRadial1 = [0.5*fOffset, 0.0]
		thisMaterialDomain_03[iIndex_MP].mRadial2 = [0.0, 0.5*fOffset]
		thisMaterialDomain_03[iIndex_MP].mDeformationGradient = eye(2,2)
		thisMaterialDomain_03[iIndex_MP].mDeformationGradientIncrement = eye(2,2)

		push!(thisMaterialPoint, thisMaterialDomain_03[iIndex_MP])
	end

	fMass = 0.0
	for iIndex_MP in 1:1:length(thisMaterialPoint)
		fMass += thisMaterialPoint[iIndex_MP].fMass
	end
	@printf("(Initial configuration) Mass: %+.6e \t", fMass)

	iMaterialPoints = length(thisMaterialPoint)
	println("Total number of material points: ", iMaterialPoints)

	# analysis timer
	fTimeIncrement = 6.0e-5
	fTimeEnd = 25.0e-2

	# realtime graphics timer
	fPlotTimeInterval = 100.0*fTimeIncrement
	fPlotTime = fPlotTimeInterval

	# final results plot timer
	fResultTimeInterval = 0.05
	fResultTime = fResultTimeInterval

	# console output timer
	fConsolTimeInterval = 100.0*fTimeIncrement
	fConsolTime = fConsolTimeInterval

	# final results plot holder arrays
	plot_Time = Array{Real}(0)
	plot_Displacement = Array{Real}(0)
	plot_KineticEnergy = Array{Real}(0)
	plot_StrainEnergy = Array{Real}(0)

	# main analysis loop
	for fTime in 0.0:fTimeIncrement:fTimeEnd
		#reset grid------------------------------------
		for iIndex in 1:1:thisGrid.iNodes
			thisGrid.GridPoints[iIndex].fMass = 0.0
			thisGrid.GridPoints[iIndex].v2Momentum.fx = 0.0
			thisGrid.GridPoints[iIndex].v2Momentum.fy = 0.0
			thisGrid.GridPoints[iIndex].v2Force.fx = 0.0
			thisGrid.GridPoints[iIndex].v2Force.fy = 0.0
		end
		# material to grid (mass, momentum, force)------------------------------
		# println("material to grid")
		for iIndex_MP in 1:1:iMaterialPoints
			# thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint[iIndex_MP], thisGrid)
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints_CPDI(thisMaterialPoint[iIndex_MP], thisGrid, fTime)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				iIndex_GP = thisAdjacentGridPoints[iIndex]
				thisGridPoint = thisGrid.GridPoints[iIndex_GP]

				fShapeValue = moduleBasis.getShapeValue(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)

				fShapeGradient_x = moduleBasis.getShapeGradient_x(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)
				fShapeGradient_y = moduleBasis.getShapeGradient_y(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)

				thisGrid.GridPoints[iIndex_GP].fMass += fShapeValue * thisMaterialPoint[iIndex_MP].fMass
				thisGrid.GridPoints[iIndex_GP].v2Momentum.fx += fShapeValue * thisMaterialPoint[iIndex_MP].v2Momentum.fx
	         thisGrid.GridPoints[iIndex_GP].v2Momentum.fy += fShapeValue * thisMaterialPoint[iIndex_MP].v2Momentum.fy
				# internal forces
				thisGrid.GridPoints[iIndex_GP].v2Force.fx += -thisMaterialPoint[iIndex_MP].fVolume * (fShapeGradient_x*thisMaterialPoint[iIndex_MP].v3Stress.f1 + fShapeGradient_y*thisMaterialPoint[iIndex_MP].v3Stress.f3)
				thisGrid.GridPoints[iIndex_GP].v2Force.fy += -thisMaterialPoint[iIndex_MP].fVolume * (fShapeGradient_y*thisMaterialPoint[iIndex_MP].v3Stress.f2 + fShapeGradient_x*thisMaterialPoint[iIndex_MP].v3Stress.f3)
				# external forces
				thisGrid.GridPoints[iIndex_GP].v2Force.fx += fShapeValue*thisMaterialPoint[iIndex_MP].v2ExternalForce.fx
				thisGrid.GridPoints[iIndex_GP].v2Force.fy += fShapeValue*thisMaterialPoint[iIndex_MP].v2ExternalForce.fy
		   end
		end
		# apply boundary conditions to grid points
		# println("boundary conditions")
		for iIndex_GP in 1:1:thisGrid.iNodes
			thisGridPoint = thisGrid.GridPoints[iIndex_GP]
			if(thisGridPoint.bFixed_x == true)
				thisGridPoint.v2Momentum.fx = 0.0
				thisGridPoint.v2Force.fx = 0.0
			end
			if(thisGridPoint.bFixed_y == true)
				thisGridPoint.v2Momentum.fy = 0.0
				thisGridPoint.v2Force.fy = 0.0
			end
		end
	 	# update grid momentum------------------------------------
		# println("update grid momenta")
		for iIndex in 1:1:thisGrid.iNodes
	      thisGridPoint = thisGrid.GridPoints[iIndex]

			thisGrid.GridPoints[iIndex].v2Momentum.fx += thisGrid.GridPoints[iIndex].v2Force.fx * fTimeIncrement
			thisGrid.GridPoints[iIndex].v2Momentum.fy += thisGrid.GridPoints[iIndex].v2Force.fy * fTimeIncrement
		end
	 	# grid to material (position increments and velocity updates)------------------------------
		# material point positions are not updated here, since they are still to be used
		# position increments are calculated to eventually be added to the positions at the final step of the loop
		# println("grid to material")
		for iIndex_MP in 1:1:iMaterialPoints
			# thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint[iIndex_MP], thisGrid)
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints_CPDI(thisMaterialPoint[iIndex_MP], thisGrid, fTime)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				iIndex_GP = thisAdjacentGridPoints[iIndex]
				thisGridPoint = thisGrid.GridPoints[iIndex_GP]

				fShapeValue = moduleBasis.getShapeValue(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)

				thisGridPointVelocity_x = 0.0
				thisGridPointVelocity_y = 0.0
				if(thisGrid.GridPoints[iIndex_GP].fMass > 1.0e-32)
					thisGridPointVelocity_x = thisGrid.GridPoints[iIndex_GP].v2Momentum.fx / thisGrid.GridPoints[iIndex_GP].fMass
					thisGridPointVelocity_y = thisGrid.GridPoints[iIndex_GP].v2Momentum.fy / thisGrid.GridPoints[iIndex_GP].fMass

					thisMaterialPoint[iIndex_MP].v2Velocity.fx += (fShapeValue * thisGrid.GridPoints[iIndex_GP].v2Force.fx / thisGrid.GridPoints[iIndex_GP].fMass) * fTimeIncrement
					thisMaterialPoint[iIndex_MP].v2Velocity.fy += (fShapeValue * thisGrid.GridPoints[iIndex_GP].v2Force.fy / thisGrid.GridPoints[iIndex_GP].fMass) * fTimeIncrement
				end

				thisMaterialPoint[iIndex_MP].v2PositionIncrement.fx += fShapeValue*thisGridPointVelocity_x * fTimeIncrement
				thisMaterialPoint[iIndex_MP].v2PositionIncrement.fy += fShapeValue*thisGridPointVelocity_y * fTimeIncrement
		   end
		end
		#reset grid momenta, to be used for re-mapping the material point to grid momenta (velocities)------------------------------------
		# println("reset grid momenta")
		for iIndex in 1:1:thisGrid.iNodes
			thisGrid.GridPoints[iIndex].v2Momentum.fx = 0.0
			thisGrid.GridPoints[iIndex].v2Momentum.fy = 0.0
		end
		# map particle momenta back to grid------------------------------
		# mass in NOT mapped here
		# println("remapping to grid")
		for iIndex_MP in 1:1:iMaterialPoints
			# thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint[iIndex_MP], thisGrid)
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints_CPDI(thisMaterialPoint[iIndex_MP], thisGrid, fTime)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				iIndex_GP = thisAdjacentGridPoints[iIndex]
				thisGridPoint = thisGrid.GridPoints[iIndex_GP]

				fShapeValue = moduleBasis.getShapeValue(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)

				thisGrid.GridPoints[iIndex_GP].v2Momentum.fx += fShapeValue * thisMaterialPoint[iIndex_MP].fMass * thisMaterialPoint[iIndex_MP].v2Velocity.fx
				thisGrid.GridPoints[iIndex_GP].v2Momentum.fy += fShapeValue * thisMaterialPoint[iIndex_MP].fMass * thisMaterialPoint[iIndex_MP].v2Velocity.fy
		   end
		end
		# apply boundary conditions to grid points
		# println("boundary conditions2")
		for iIndex_GP in 1:1:thisGrid.iNodes
			thisGridPoint = thisGrid.GridPoints[iIndex_GP]
			if(thisGridPoint.bFixed_x == true)
				# thisGridPoint.fMass = 0.0
				thisGridPoint.v2Momentum.fx = 0.0
				thisGridPoint.v2Force.fx = 0.0
			end
			if(thisGridPoint.bFixed_y == true)
				# thisGridPoint.fMass = 0.0
				thisGridPoint.v2Momentum.fy = 0.0
				thisGridPoint.v2Force.fy = 0.0
			end
		end
		#strain calculations------------------------------
		# println("strain calculations")
		for iIndex_MP in 1:1:iMaterialPoints
			# thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisMaterialPoint[iIndex_MP], thisGrid)
			thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints_CPDI(thisMaterialPoint[iIndex_MP], thisGrid, fTime)
			for iIndex in 1:1:length(thisAdjacentGridPoints)
				iIndex_GP = thisAdjacentGridPoints[iIndex]
				thisGridPoint = thisGrid.GridPoints[iIndex_GP]

				fShapeValue = moduleBasis.getShapeValue(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)

				fShapeGradient_x = moduleBasis.getShapeGradient_x(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)
				fShapeGradient_y = moduleBasis.getShapeGradient_y(thisMaterialPoint[iIndex_MP], thisGridPoint, thisGrid)

				thisGridPointVelocity_x = 0.0
				thisGridPointVelocity_y = 0.0
				if(thisGrid.GridPoints[iIndex_GP].fMass > 1.0e-32)
					thisGridPointVelocity_x = thisGrid.GridPoints[iIndex_GP].v2Momentum.fx / thisGrid.GridPoints[iIndex_GP].fMass
					thisGridPointVelocity_y = thisGrid.GridPoints[iIndex_GP].v2Momentum.fy / thisGrid.GridPoints[iIndex_GP].fMass
				end

				thisMaterialPoint[iIndex_MP].mDeformationGradientIncrement[1,1] += fShapeGradient_x * thisGridPointVelocity_x * fTimeIncrement
				thisMaterialPoint[iIndex_MP].mDeformationGradientIncrement[2,2] += fShapeGradient_y * thisGridPointVelocity_y * fTimeIncrement
				thisMaterialPoint[iIndex_MP].mDeformationGradientIncrement[1,2] += fShapeGradient_y * thisGridPointVelocity_x * fTimeIncrement
				thisMaterialPoint[iIndex_MP].mDeformationGradientIncrement[2,1] += fShapeGradient_x * thisGridPointVelocity_y * fTimeIncrement
		   end
		end
		# calculate corner increments-------------------------------------------------
		# println("corner increments")
		for iIndex_MP in 1:1:iMaterialPoints
			for iIndex_Corner = 1:1:size(thisMaterialPoint[iIndex_MP].mCorner, 1)
				thisCorner = moduleMath.Vector2D(0.0, 0.0)
				thisCorner.fx = thisMaterialPoint[iIndex_MP].mCorner[iIndex_Corner, 1]
				thisCorner.fy = thisMaterialPoint[iIndex_MP].mCorner[iIndex_Corner, 2]
				thisAdjacentGridPoints = moduleGrid.getAdjacentGridPoints(thisCorner, thisGrid, fTime)
				for iIndex in 1:1:length(thisAdjacentGridPoints)
					iIndex_GP = thisAdjacentGridPoints[iIndex]
					thisGridPoint = thisGrid.GridPoints[iIndex_GP]

					fDistance_x = thisCorner.fx - thisGridPoint.v2Position.fx
					fDistance_y = thisCorner.fy - thisGridPoint.v2Position.fy
					fShapeValue = moduleBasis.getShapeValue_Classic(fDistance_x, fDistance_y, thisGrid.fLength_Cell_x, thisGrid.fLength_Cell_y)
					thisGridPointVelocity_x = 0.0
					thisGridPointVelocity_y = 0.0

					if(thisGrid.GridPoints[iIndex_GP].fMass > 1.0e-32)
						thisGridPointVelocity_x = thisGridPoint.v2Momentum.fx / thisGridPoint.fMass
						thisGridPointVelocity_y = thisGridPoint.v2Momentum.fy / thisGridPoint.fMass
					end

					thisMaterialPoint[iIndex_MP].mCorner_Increment[iIndex_Corner, 1] += fShapeValue * thisGridPointVelocity_x * fTimeIncrement
					thisMaterialPoint[iIndex_MP].mCorner_Increment[iIndex_Corner, 2] += fShapeValue * thisGridPointVelocity_y * fTimeIncrement
		# println("mCorner_Increment: ", thisMaterialPoint[iIndex_MP].mCorner_Increment[iIndex_Corner, :])
				end
			end
		end
		# update particle corners, strains, deformation gradient, volume-------------------
		# println("update particles")
		for iIndex_MP in 1:1:iMaterialPoints
			thisMaterialPoint[iIndex_MP].v3Strain.f1 += thisMaterialPoint[iIndex_MP].v3StrainIncrement.f1
			thisMaterialPoint[iIndex_MP].v3Strain.f2 += thisMaterialPoint[iIndex_MP].v3StrainIncrement.f2
			thisMaterialPoint[iIndex_MP].v3Strain.f3 += thisMaterialPoint[iIndex_MP].v3StrainIncrement.f3

			thisMaterialPoint[iIndex_MP].mDeformationGradient = thisMaterialPoint[iIndex_MP].mDeformationGradientIncrement * thisMaterialPoint[iIndex_MP].mDeformationGradient
			thisMaterialPoint[iIndex_MP].mDeformationGradientIncrement = eye(2,2)

			for iIndex_Corner = 1:1:size(thisMaterialPoint[iIndex_MP].mCorner, 1)
				thisMaterialPoint[iIndex_MP].mCorner[iIndex_Corner, :] += thisMaterialPoint[iIndex_MP].mCorner_Increment[iIndex_Corner, :]
				thisMaterialPoint[iIndex_MP].mCorner_Increment[iIndex_Corner, :] = [0.0, 0.0]
			end

			thisMaterialPoint[iIndex_MP].v3StrainIncrement.f1 = 0.0
			thisMaterialPoint[iIndex_MP].v3StrainIncrement.f2 = 0.0
			thisMaterialPoint[iIndex_MP].v3StrainIncrement.f3 = 0.0

			thisMaterialPoint[iIndex_MP].fVolume = det(thisMaterialPoint[iIndex_MP].mDeformationGradient) * thisMaterialPoint[iIndex_MP].fVolumeInitial
		end
		#stress calculations------------------------------------------------------
		# println("stress calculations")
		for iIndex_MP in 1:1:iMaterialPoints
			fE = thisMaterialPoint[iIndex_MP].fElasticModulus;
			fNu = thisMaterialPoint[iIndex_MP].fPoissonRatio

			# for a Neo_Hookean model, (2011) A convected particle domain interpolation technique to extend
			fLame1 = 0.5*fE / (1.0 + fNu)
			fLame2 = fNu*fE / ((1+fNu)*(1.0-2.0*fNu))

			mF = thisMaterialPoint[iIndex_MP].mDeformationGradient
			fJ = det(mF)

			mP = fLame2*log(fJ)/fJ * eye(2) + fLame1/fJ*(mF*mF'-eye(2))
			thisMaterialPoint[iIndex_MP].v3Stress.f1 = mP[1,1]
			thisMaterialPoint[iIndex_MP].v3Stress.f2 = mP[2,2]
			thisMaterialPoint[iIndex_MP].v3Stress.f3 = mP[1,2]
		end
		# update particle positions, momenta----------------
		# println("update particles2")
		fStrainEnergy = 0.0
		fKineticEnergy = 0.0
		for iIndex_MP in 1:1:iMaterialPoints
			thisMaterialPoint[iIndex_MP].v2Position.fx += thisMaterialPoint[iIndex_MP].v2PositionIncrement.fx
			thisMaterialPoint[iIndex_MP].v2Position.fy += thisMaterialPoint[iIndex_MP].v2PositionIncrement.fy

			thisMaterialPoint[iIndex_MP].v2PositionIncrement.fx = 0.0
			thisMaterialPoint[iIndex_MP].v2PositionIncrement.fy = 0.0

			thisMaterialPoint[iIndex_MP].v2Momentum.fx = thisMaterialPoint[iIndex_MP].v2Velocity.fx * thisMaterialPoint[iIndex_MP].fMass
			thisMaterialPoint[iIndex_MP].v2Momentum.fy = thisMaterialPoint[iIndex_MP].v2Velocity.fy * thisMaterialPoint[iIndex_MP].fMass
		end

		# calculating strain and kinetic energy for final results plot
		fStrainEnergy = 0.0
		fKineticEnergy = 0.0
		for iIndex_MP in 1:1:iMaterialPoints
			fStrainEnergy += 0.5*thisMaterialPoint[iIndex_MP].v3Strain.f1 * thisMaterialPoint[iIndex_MP].v3Stress.f1 * thisMaterialPoint[iIndex_MP].fVolume
			fStrainEnergy += 0.5*thisMaterialPoint[iIndex_MP].v3Strain.f2 * thisMaterialPoint[iIndex_MP].v3Stress.f2 * thisMaterialPoint[iIndex_MP].fVolume
			fStrainEnergy += 0.25*thisMaterialPoint[iIndex_MP].v3Strain.f3 * thisMaterialPoint[iIndex_MP].v3Stress.f3 * thisMaterialPoint[iIndex_MP].fVolume

			fVelocity2 = thisMaterialPoint[iIndex_MP].v2Velocity.fx^2 + thisMaterialPoint[iIndex_MP].v2Velocity.fy^2
			fKineticEnergy += 0.5*fVelocity2 * thisMaterialPoint[iIndex_MP].fMass
		end

		#save to plot arrays
		push!(plot_Time, fTime)
		push!(plot_Displacement, thisMaterialPoint[3].v2Position.fy - 8.58333333)
		push!(plot_KineticEnergy, fKineticEnergy)
		push!(plot_StrainEnergy, fStrainEnergy)

		# consol output
		fConsolTime += fTimeIncrement
		if(fConsolTime > fConsolTimeInterval)
			fConsolTime = 0.0

			fMass = 0.0
			fMomentum_x = 0.0
			fKineticEnergy = 0.0
			for iIndex_MP in 1:1:length(thisMaterialPoint)
				fMass += thisMaterialPoint[iIndex_MP].fMass

				fMomentum_x += thisMaterialPoint[iIndex_MP].v2Momentum.fx

				fVelocity2 = thisMaterialPoint[iIndex_MP].v2Velocity.fx^2 + thisMaterialPoint[iIndex_MP].v2Velocity.fy^2
				fKineticEnergy += 0.5*fVelocity2 * thisMaterialPoint[iIndex_MP].fMass
			end
			@printf("fTime: %+.6e \t fMomentum: %+.6e \t fKineticEnergy: : %+.6e \n", fTime, fMomentum_x, fKineticEnergy)
		end

		# realtime graphical plotting routines
		fPlotTime += fTimeIncrement
		if(fPlotTime > fPlotTimeInterval)
			fPlotTime = 0.0

			array_x = [thisGrid.GridPoints[i].v2Position.fx for i in 1:thisGrid.iNodes]
		   array_y = [thisGrid.GridPoints[i].v2Position.fy for i in 1:thisGrid.iNodes]
		   array_color = Array{Real}(thisGrid.iNodes, 3)
			array_size = Array{Real}(thisGrid.iNodes, 1)
		   for iIndex in 1:1:thisGrid.iNodes
		      array_color[iIndex, :] = [0.0, 0.5, 0.0]#[thisGrid.GridPoints[iIndex].fMass/iMaterialPoints, 0.0, 0.0]
				# if(iIndex == 89 || iIndex == 90 || iIndex == 110 || iIndex == 111)
				# 	array_color[iIndex, :] = [0.0, 1.0, 0.0]
				# end
				if(thisGrid.GridPoints[iIndex].bFixed_x == true || thisGrid.GridPoints[iIndex].bFixed_y == true)
					array_color[iIndex, :] = [0.0, 0.0, 0.0]
					array_size[iIndex, :] = [10]
				else
					array_size[iIndex, :] = [1+thisGrid.GridPoints[iIndex].fMass*0.1]
				end
		   end

			pyPlot01 = PyPlot.subplot2grid((5,5), (0,0), colspan=2, rowspan=5)
			PyPlot.scatter(array_x, array_y, c=array_color, lw=0, s=array_size)
			PyPlot.hold(true)

			fStrainEnergy = 0.0
			for iIndex_MP in 1:1:iMaterialPoints
				fStrainEnergy = 0.0
				fStrainEnergy += 0.5*thisMaterialPoint[iIndex_MP].v3Strain.f1 * thisMaterialPoint[iIndex_MP].v3Stress.f1 * thisMaterialPoint[iIndex_MP].fVolume
				fStrainEnergy += 0.5*thisMaterialPoint[iIndex_MP].v3Strain.f2 * thisMaterialPoint[iIndex_MP].v3Stress.f2 * thisMaterialPoint[iIndex_MP].fVolume
				fStrainEnergy += 0.25*thisMaterialPoint[iIndex_MP].v3Strain.f3 * thisMaterialPoint[iIndex_MP].v3Stress.f3 * thisMaterialPoint[iIndex_MP].fVolume
				fStrainEnergy *= 10.0
				if(fStrainEnergy > 1.0)
					fStrainEnergy = 1.0
				end

				thisColor = [fStrainEnergy, 0.0, 0.0]
				if(iIndex_MP == 3)
					thisColor = [0.0, 1.0, 0.0]
				end

				thisSize = 5 + fStrainEnergy * 10.0

				nCorners = size(thisMaterialPoint[iIndex_MP].mCorner, 1)
				thisMarker = Array{Real}(nCorners+1,2)

				for iIndex_Corner = 1:1:nCorners
					thisMarker[iIndex_Corner,:] = thisMaterialPoint[iIndex_MP].mCorner[iIndex_Corner,:]
				end
				thisMarker[nCorners+1,:] = thisMaterialPoint[iIndex_MP].mCorner[1,:]
				# thisMarker[1,:] = thisMaterialPoint[iIndex_MP].mCorner[1,:]
				# thisRectangle[2,:] = thisMaterialPoint[iIndex_MP].mCorner[2,:]
				# thisRectangle[3,:] = thisMaterialPoint[iIndex_MP].mCorner[3,:]
				# thisRectangle[4,:] = thisMaterialPoint[iIndex_MP].mCorner[4,:]
				# thisRectangle[5,:] = thisMaterialPoint[iIndex_MP].mCorner[1,:]

				PyPlot.plot(thisMarker[:,1], thisMarker[:,2], c=[1.0, 0.0, 0.0])
				PyPlot.scatter([thisMaterialPoint[iIndex_MP].v2Position.fx], [thisMaterialPoint[iIndex_MP].v2Position.fy], c=thisColor, lw = 0, s=thisSize, marker="o")
			end
			PyPlot.hold(false)

			pyPlot02 = PyPlot.subplot2grid((5,5), (0,3), colspan=2, rowspan=2)
			PyPlot.plot(plot_Time, plot_Displacement, "-")
			# PyPlot.plot(plot_Time, plot_KineticEnergy, "-")
			# PyPlot.hold(true)
			# PyPlot.plot(plot_Time, plot_StrainEnergy, "-")
			# PyPlot.plot(plot_Time, plot_KineticEnergy + plot_StrainEnergy, "-")
			PyPlot.grid("on")

		   PyPlot.pause(0.01)
		   PyPlot.draw()
		   PyPlot.hold(false)
		end
	end

	# final plots
	pyFig_FinalResults = PyPlot.figure("Time-Displacement", figsize=(6, 6))
	PyPlot.plot(plot_Time, plot_Displacement, "-")

	# PyPlot.plot(plot_Time, plot_KineticEnergy, "-")
	# PyPlot.hold(true)
	# PyPlot.plot(plot_Time, plot_StrainEnergy, "-")
	# PyPlot.plot(plot_Time, plot_KineticEnergy + plot_StrainEnergy, "-")
end

mpmMain()
