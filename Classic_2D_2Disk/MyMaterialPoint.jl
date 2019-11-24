module moduleMaterialPoint
    # material point container
	type mpmMaterialPoint_2D_Classic   
        fMass            :: Float64
		fVolumeInitial   :: Float64
		fVolume          :: Float64

		fElasticModulus  :: Float64
		fPoissonRatio    :: Float64

		v2Centroid       :: Array{Float64}  # position
		v2Velocity       :: Array{Float64}
        v2Momentum       :: Array{Float64}
		v2ExternalForce  :: Array{Float64}
		v2Restraint      :: Array{Float64}	# 0.0=no restraint, 1.0=fully restrained

		v2Corner         :: Array{Float64} # corner position

		m22DeformationGradient          :: Array{Float64}
		m22DeformationGradientIncrement :: Array{Float64}

		v3Strain         :: Array{Float64} # xx, yy, zz, xy, yz, zx
		v3Stress         :: Array{Float64}

      function mpmMaterialPoint_2D_Classic()
         new(1.0, 1.0, 1.0, # Mass, initial volume, volume
				1.0, 0.3, # elastic modulus, poisson ratio
				zeros(2), # centroid position
				zeros(2), # velocity
				zeros(2), # momentum
				zeros(2), # external force
				zeros(2), # restraint
				zeros(2,4), # array of corner positions, 3d coordinates
				eye(2,2), # deformation gradient
				eye(2,2), # deformation gradient increment
				zeros(3), # strain
				zeros(3), # stress
				)
      end
   end

	function createMaterialDomain_Circle(fCenter::Array{Float64}, fRadius::Float64, fOffset::Float64)
		thisMaterialDomain = Array{mpmMaterialPoint_2D_Classic}(0)

		fRadius = floor(fRadius/fOffset) * fOffset	#just in case radius is not a multiple of offset

		for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
			for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
				if(fx^2 + fy^2 < fRadius^2)
					thisMaterialPoint = mpmMaterialPoint_2D_Classic()
					thisMaterialPoint.v2Centroid = [fCenter[1] + fx; fCenter[2] + fy]
					push!(thisMaterialDomain, thisMaterialPoint)
				end
			end
		end

		return(thisMaterialDomain)
	end
end
