module moduleMaterialPoint
   import moduleMath #sina, do not use include here, since you have already included the module in Main.jl

   type mpmMaterialPoint   #material point container
      fMass::Real
		fVolumeInitial::Real
		fVolume::Real
		v2Length::moduleMath.Vector2D

		fElasticModulus::Real
		fPoissonRatio::Real

		v2Position::moduleMath.Vector2D
		v2PositionIncrement::moduleMath.Vector2D
		v2Velocity::moduleMath.Vector2D
      v2Momentum::moduleMath.Vector2D
		v2ExternalForce::moduleMath.Vector2D
		v2Restraint::moduleMath.Vector2D	# 0.0=no restraint, 1.0=fully restrained

		mCorner::Array{Real}
		mCorner_Increment::Array{Real}

		mRadial1::Array{Real}
		mRadial2::Array{Real}

		mDeformationGradient::Array{Real}
		mDeformationGradientIncrement::Array{Real}

		v3Strain::moduleMath.Vector3D
		v3Stress::moduleMath.Vector3D
		v3StrainIncrement::moduleMath.Vector3D
		v3StressIncrement::moduleMath.Vector3D

      function mpmMaterialPoint()
         new(1.0, 1.0, 1.0,
				moduleMath.Vector2D(0.0, 0.0),
				1.0, 0.3,
				moduleMath.Vector2D(0.0, 0.0),
				moduleMath.Vector2D(0.0, 0.0),
				moduleMath.Vector2D(0.0, 0.0),
				moduleMath.Vector2D(0.0, 0.0),
				moduleMath.Vector2D(0.0, 0.0),
				moduleMath.Vector2D(0.0, 0.0),
				zeros(0,2),
				zeros(0,2),
				ones(2),
				ones(2),
				eye(2,2),
				eye(2,2),
				moduleMath.Vector3D(0.0, 0.0, 0.0),
				moduleMath.Vector3D(0.0, 0.0, 0.0),
				moduleMath.Vector3D(0.0, 0.0, 0.0),
				moduleMath.Vector3D(0.0, 0.0, 0.0))
      end
      function mpmMaterialPoint(fM::Real, fV::Real, fEM::Real, fPR::Real, v2P::moduleMath.Vector2D, v2V::moduleMath.Vector2D, v2M::moduleMath.Vector2D, v2ExternalForce::moduleMath.Vector2D, v3Strain::moduleMath.Vector3D, v3Stress::moduleMath.Vector3D)
         new(fM, fV, fv,
				moduleMath.Vector2D(0.0, 0.0),
				fEM, fPR,
				moduleMath.Vector2D(v2P),
				moduleMath.Vector2D(0.0, 0.0),
				moduleMath.Vector2D(v2V),
				moduleMath.Vector2D(v2M),
				moduleMath.Vector2D(v2ExternalForce),
				moduleMath.Vector2D(0.0, 0.0),
				zeros(0,2),
				zeros(0,2),
				ones(2),
				ones(2),
				eye(2,2),
				eye(2,2),
				moduleMath.Vector3D(0.0, 0.0, 0.0),
				moduleMath.Vector3D(0.0, 0.0, 0.0),
				moduleMath.Vector3D(0.0, 0.0, 0.0),
				moduleMath.Vector3D(0.0, 0.0, 0.0))
      end
   end

	function createMaterialDomain_Rectangle(sParticleShape::String, fCenter_x::Real, fCenter_y::Real, fWidth::Real, fHeight::Real, fOffset::Real)
		thisMaterialDomain = Array{mpmMaterialPoint}(0)

		if(sParticleShape == "triangle")
			# left triangles
			for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
				for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
					v2Corner = Array{Real}(3,2)
					v2Corner[1,1] = fBaseCorner_x + 0.0
					v2Corner[1,2] = fBaseCorner_y + 0.0

					v2Corner[2,1] = fBaseCorner_x + 0.5*fOffset
					v2Corner[2,2] = fBaseCorner_y + 0.5*fOffset

					v2Corner[3,1] = fBaseCorner_x + 0.0
					v2Corner[3,2] = fBaseCorner_y + 1.0*fOffset

					v2Centroid = [0.0 0.0]

					fWeight = 1.0 / size(v2Corner,1)
					for iIndex_Corner = 1:1:size(v2Corner,1)
						v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
						v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
					end

					if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
						if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
							thisMaterialPoint = mpmMaterialPoint()
							thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
							thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

							for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
								newCorner = [0.0 0.0]
								newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
								newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
								thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

								thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
							end

							push!(thisMaterialDomain, thisMaterialPoint)
						end
					end
				end
			end
			# bottom triangles
			for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
				for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
					v2Corner = Array{Real}(3,2)
					v2Corner[1,1] = fBaseCorner_x + 0.0
					v2Corner[1,2] = fBaseCorner_y + 0.0

					v2Corner[2,1] = fBaseCorner_x + 1.0*fOffset
					v2Corner[2,2] = fBaseCorner_y + 0.0

					v2Corner[3,1] = fBaseCorner_x + 0.5*fOffset
					v2Corner[3,2] = fBaseCorner_y + 0.5*fOffset

					v2Centroid = [0.0 0.0]

					fWeight = 1.0 / size(v2Corner,1)
					for iIndex_Corner = 1:1:size(v2Corner,1)
						v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
						v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
					end

					if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
						if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
							thisMaterialPoint = mpmMaterialPoint()
							thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
							thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

							for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
								newCorner = [0.0 0.0]
								newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
								newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
								thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

								thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
							end

							push!(thisMaterialDomain, thisMaterialPoint)
						end
					end
				end
			end
			# right triangles
			for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
				for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
					v2Corner = Array{Real}(3,2)
					v2Corner[1,1] = fBaseCorner_x + 1.0*fOffset
					v2Corner[1,2] = fBaseCorner_y + 0.0

					v2Corner[2,1] = fBaseCorner_x + 1.0*fOffset
					v2Corner[2,2] = fBaseCorner_y + 1.0*fOffset

					v2Corner[3,1] = fBaseCorner_x + 0.5*fOffset
					v2Corner[3,2] = fBaseCorner_y + 0.5*fOffset

					v2Centroid = [0.0 0.0]

					fWeight = 1.0 / size(v2Corner,1)
					for iIndex_Corner = 1:1:size(v2Corner,1)
						v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
						v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
					end

					if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
						if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
							thisMaterialPoint = mpmMaterialPoint()
							thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
							thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

							for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
								newCorner = [0.0 0.0]
								newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
								newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
								thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

								thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
							end

							push!(thisMaterialDomain, thisMaterialPoint)
						end
					end
				end
			end
			# top triangles
			for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
				for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
					v2Corner = Array{Real}(3,2)
					v2Corner[1,1] = fBaseCorner_x + 1.0*fOffset
					v2Corner[1,2] = fBaseCorner_y + 1.0*fOffset

					v2Corner[2,1] = fBaseCorner_x + 0.0
					v2Corner[2,2] = fBaseCorner_y + 1.0*fOffset

					v2Corner[3,1] = fBaseCorner_x + 0.5*fOffset
					v2Corner[3,2] = fBaseCorner_y + 0.5*fOffset

					v2Centroid = [0.0 0.0]

					fWeight = 1.0 / size(v2Corner,1)
					for iIndex_Corner = 1:1:size(v2Corner,1)
						v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
						v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
					end

					if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
						if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
							thisMaterialPoint = mpmMaterialPoint()
							thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
							thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

							for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
								newCorner = [0.0 0.0]
								newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
								newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
								thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

								thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
							end

							push!(thisMaterialDomain, thisMaterialPoint)
						end
					end
				end
			end
		end

		if(sParticleShape == "rectangle")
			for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
				for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
					v2Corner = Array{Real}(4,2)
					v2Corner[1,1] = fBaseCorner_x + 0.0
					v2Corner[1,2] = fBaseCorner_y + 0.0

					v2Corner[2,1] = fBaseCorner_x + fOffset
					v2Corner[2,2] = fBaseCorner_y + 0.0

					v2Corner[3,1] = fBaseCorner_x + fOffset
					v2Corner[3,2] = fBaseCorner_y + fOffset

					v2Corner[4,1] = fBaseCorner_x + 0.0
					v2Corner[4,2] = fBaseCorner_y + fOffset

					v2Centroid = [0.0 0.0]

					for iIndex_Corner = 1:1:size(v2Corner,1)
						v2Centroid[1,1] += 0.25 * v2Corner[iIndex_Corner,1]
						v2Centroid[1,2] += 0.25 * v2Corner[iIndex_Corner,2]
					end

					if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
						if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
							thisMaterialPoint = mpmMaterialPoint()
							thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
							thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

							for iIndex_Corner = 1:1:size(v2Corner,1)
								newCorner = [0.0 0.0]
								newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
								newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
								thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

								thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
							end

							push!(thisMaterialDomain, thisMaterialPoint)
						end
					end
				end
			end
		end

		return(thisMaterialDomain)
	end

	function createMaterialDomain_Rectangle(fCenter_x::Real, fCenter_y::Real, fWidth::Real, fHeight::Real, fOffset::Real)
		thisMaterialDomain = Array{mpmMaterialPoint}(0)

		fWidth	= floor(fWidth/fOffset) * fOffset	#just in case width is not a multiple of offset
		fHeight	= floor(fHeight/fOffset) * fOffset	#just in case height is not a multiple of offset

		for fy in -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
			for fx in -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
				thisMaterialPoint = mpmMaterialPoint()
				thisMaterialPoint.v2Position.fx = fCenter_x + fx
				thisMaterialPoint.v2Position.fy = fCenter_y + fy
				push!(thisMaterialDomain, thisMaterialPoint)
			end
		end

		return(thisMaterialDomain)
	end

	function createMaterialDomain_Circle(fCenter_x::Real, fCenter_y::Real, fRadius::Real, fOffset::Real)
		thisMaterialDomain = Array{mpmMaterialPoint}(0)

		fRadius = floor(fRadius/fOffset) * fOffset	#just in case radius is not a multiple of offset

		for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
			for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
				if(fx^2 + fy^2 < fRadius^2)
					thisMaterialPoint = mpmMaterialPoint()
					thisMaterialPoint.v2Position.fx = fCenter_x + fx
					thisMaterialPoint.v2Position.fy = fCenter_y + fy
					push!(thisMaterialDomain, thisMaterialPoint)
				end
			end
		end

		return(thisMaterialDomain)
	end
end
