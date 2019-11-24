module moduleMaterialPoint
	type mpmMaterialPoint_2D_Classic   #material point container
      fMass::Float64
		fVolumeInitial::Float64
		fVolume::Float64

		fElasticModulus::Float64
		fPoissonRatio::Float64
		fYieldStress::Float64

		v2Centroid::Array{Float64}
		v2Velocity::Array{Float64}
      v2Momentum::Array{Float64}
		v2ExternalForce::Array{Float64}
		v2Restraint::Array{Float64}	# 0.0=no restraint, 1.0=fully restrained

		v2Corner::Array{Float64} # corner position

		m22DeformationGradient::Array{Float64}
		m22DeformationGradientIncrement::Array{Float64}

		v3Strain::Array{Float64} # xx, yy, zz, xy, yz, zx
		v3PlasticStrain::Array{Float64} # xx, yy, zz, xy, yz, zx
		fAlpha::Float64 # equivalent plastic strain

		v3Stress::Array{Float64}

      function mpmMaterialPoint_2D_Classic()
         new(1.0, 1.0, 1.0, # Mass, initial volume, volume
				1.0, 0.3, 1.0e24, # elastic modulus, poisson ratio
				zeros(2), # centroid position
				zeros(2), # velocity
				zeros(2), # momentum
				zeros(2), # external force
				zeros(2), # restraint
				zeros(2,4), # array of corner positions, 3d coordinates
				eye(2,2), # deformation gradient
				eye(2,2), # deformation gradient increment
				zeros(3), # strain
				zeros(3), # plastic strain
				0.0, # equivalent plastic strain
				zeros(3), # stress
				)
      end
   end

	function createMaterialDomain_Circle(fCenter::Array{Float64}, fRadius::Float64, fOffset::Float64)
		thisMaterialDomain = Array{mpmMaterialPoint_2D_Classic}(0)

		fRadius = floor(fRadius/fOffset) * fOffset	#just in case radius is not a multiple of offset

		for fy = -fRadius-0.5*fOffset:fOffset:+fRadius+0.5*fOffset
			for fx = -fRadius-0.5*fOffset:fOffset:+fRadius+0.5*fOffset
				if(fx^2 + fy^2 <= fRadius^2)
					thisMaterialPoint = mpmMaterialPoint_2D_Classic()
					thisMaterialPoint.v2Centroid = [fCenter[1] + fx; fCenter[2] + fy]
					push!(thisMaterialDomain, thisMaterialPoint)
				end
			end
		end

		return(thisMaterialDomain)
	end

	function createMaterialDomain_Rectangle(v2Center::Array{Float64}, fWidth::Real, fHeight::Real, fOffset::Real)
		thisMaterialDomain = Array{mpmMaterialPoint_2D_Classic}(0)

		fWidth	= floor(fWidth/fOffset) * fOffset	#just in case width is not a multiple of offset
		fHeight	= floor(fHeight/fOffset) * fOffset	#just in case height is not a multiple of offset

		for fy = -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
			for fx = -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
				thisMaterialPoint = mpmMaterialPoint_2D_Classic()
				thisMaterialPoint.v2Centroid = v2Center + [fx; fy]
				push!(thisMaterialDomain, thisMaterialPoint)
			end
		end

		return(thisMaterialDomain)
	end

	function getStressIncrement_Elastic(fE::Float64, fNu::Float64, v3StrainIncrement::Array{Float64})
		v3Result = zeros(3)

		fConstant = fE/(1.0 + fNu)/(1.0 - 2.0*fNu)

		v3Result[1] = fConstant * ((1.0-fNu)*v3StrainIncrement[1] + fNu*v3StrainIncrement[2])
		v3Result[2] = fConstant * ((1.0-fNu)*v3StrainIncrement[2] + fNu*v3StrainIncrement[1])
		v3Result[3] = fConstant * ((0.5-fNu)*v3StrainIncrement[3])

		return(v3Result)
	end
	function getIncrement_Plastic(fE::Float64, fNu::Float64, fK0::Float64, fAlphaCurrent::Float64, v3StressCurrent::Array{Float64}, v3StrainCurrent::Array{Float64}, v3PlasticStrainCurrent::Array{Float64}, v3StrainIncrement::Array{Float64})
		v32Result = zeros(3,2)

		v3StressIncrement = zeros(3)
		v3PlasticStrainIncrement = zeros(3)
		fAlphaIncrement = 0.0

		mu      = fE/2.0/(1.0+fNu);    # shear modulus
		lambda  = fE*fNu/((1.0+fNu)*(1.0-2.0*fNu));
		kappa   = lambda + mu;

		eye2   = [1.0; 1.0; 0.0];
		eye2x2 = [1.0 1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 0.0];
		I_dev  = eye(3) - 0.5*eye2x2;
		I      = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.5];#0.5 to make engineering strain to physical one
		Iinv   = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0];

		# Compute trial stress
		epsilon_dev  = I_dev * v3StrainCurrent;
		s_trial      = 2.0 * mu * I * (epsilon_dev - v3PlasticStrainCurrent);
		norm_s_trial = sqrt(s_trial[1]^2 + s_trial[2]^2 + 2*s_trial[3]^2);
		sigma_trial  = kappa*sum(v3StrainCurrent[1] + v3StrainCurrent[2])*eye2 + s_trial;

		# Check yield condition
		# alpha0 = 0.0 # sina, perfect plasticity for now
		k1 = 0.0 #sina, perfect plasticity for now
		f_trial = norm_s_trial - (k1*fAlphaCurrent + fK0);

		if f_trial <= 0.0 # elastic update
		   fAlphaIncrement = 0.0
			v3PlasticStrainIncrement = zeros(3)
			v3StressIncrement = sigma_trial - v3StressCurrent
		else # plastic step
		    normal = s_trial/norm_s_trial;
		    lambda = (norm_s_trial - k1*fAlphaCurrent - fK0)/(2.0*mu + k1);
			 fAlphaIncrement = lambda
		    # Update plastic strain and stress
			 v3PlasticStrainIncrement = lambda*Iinv*normal;
			 v3StressIncrement = kappa*sum(v3StrainCurrent[1] + v3StrainCurrent[2])*eye2 + s_trial - 2.0*mu*lambda*normal;
			 v3StressIncrement -= v3StressCurrent
		end

		v32Result = hcat(v3StressIncrement, v3PlasticStrainIncrement)
		v32Result = hcat(v32Result, [fAlphaIncrement; fAlphaIncrement; fAlphaIncrement])
		return(v32Result)
	end
end
