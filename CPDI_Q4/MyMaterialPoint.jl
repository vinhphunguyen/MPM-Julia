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
				zeros(4,2),
				zeros(4,2),
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
				zeros(4,2),
				zeros(4,2),
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
