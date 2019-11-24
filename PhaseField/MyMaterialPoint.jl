module moduleMaterialPoint

    export mpmMaterialPoint_2D_Classic
    
	type mpmMaterialPoint_2D_Classic   #material point container
      fMass           :: Float64
	  fVolumeInitial  :: Float64
	  fVolume         :: Float64

	  fElasticModulus :: Float64         # Young modulus    
	  fPoissonRatio   :: Float64         # Poisson ratio
      fShear          :: Float64         # shear modulus
      fBulk           :: Float64         # bulk modulus

	  v2Centroid      :: Array{Float64}
	  v2Velocity      :: Array{Float64}
      v2Momentum      :: Array{Float64}
	  v2ExternalForce :: Array{Float64}
	  v2Restraint     :: Array{Float64}	# 0.0=no restraint, 1.0=fully restrained

	  m22DeformationGradient::Array{Float64}
	  m22DeformationGradientIncrement::Array{Float64}

	  v3Strain        :: Array{Float64} # xx, yy, zz, xy, yz, zx
	  v3Stress        :: Array{Float64}

      fHistory        :: Float64   # local history field H of Miehe
      fColor          :: Float64   # color for visualisation
      fPhase          :: Float64   # phase field
      fEqvStrain      :: Float64   # local equivalent strain

      function mpmMaterialPoint_2D_Classic()
         new(1.0, 1.0, 1.0, # Mass, initial volume, volume
			 1.0, 0.3,0.,0.,  # elastic modulus, poisson ratio, bulk, shear
			 zeros(2), # centroid position
			 zeros(2), # velocity
			 zeros(2), # momentum
			 zeros(2), # external force
			 zeros(2), # restraint
			 eye(2,2), # deformation gradient
			 eye(2,2), # deformation gradient increment
			 zeros(3), # strain
			 zeros(3), # stress
             10.0, 0.0, 0., # local history, color, phase field
             0.0
			)
      end
   end

    # compute stress using Amor's tension/compression split  
    # phi: phase field
    function getStress(bulk, shear, phi, k, strain::Array{Float64})
		v3Result     = zeros(3)

        traceEps     = strain[1] + strain[2]
        traceEpsP    = 0.5 * (traceEps + abs(traceEps))
        traceEpsM    = 0.5 * (traceEps - abs(traceEps))

        degrad       = (1.0-phi)^2 + k   

        epsXX        = strain[1]
        epsYY        = strain[2]
        epsXY        = strain[3]

        v3Result[1] = degrad * ( bulk*traceEpsP + 2. * shear * (epsXX-epsYY) ) + bulk * traceEpsM 
        v3Result[2] = degrad * ( bulk*traceEpsP - 2. * shear * (epsXX-epsYY) ) + bulk * traceEpsM 
        v3Result[3] = degrad * (                  2. * shear * epsXY ) 

		return(v3Result)
	end
	
    function getStressIncrement_Elastic(fE::Float64, fNu::Float64, v3StrainIncrement::Array{Float64})
		v3Result = zeros(3)

		fConstant = fE/(1.0 + fNu)/(1.0 - 2.0*fNu)

		v3Result[1] = fConstant * ((1.0-fNu)*v3StrainIncrement[1] + fNu*v3StrainIncrement[2])
		v3Result[2] = fConstant * ((1.0-fNu)*v3StrainIncrement[2] + fNu*v3StrainIncrement[1])
		v3Result[3] = fConstant * ((0.5-fNu)*v3StrainIncrement[3])

		return(v3Result)
	end

    # compute plane strain principal strains
    function getPrincipalStrains(v3Strain::Array{Float64})

      exx =       v3Strain[1]
      eyy =       v3Strain[2]
      exy = 0.5 * v3Strain[3]

      # principal strains are root of a quadratic equation with det = d

      d   = ( exx - eyy ) * ( exx - eyy ) + 4.0 * exy * exy
      d   = sqrt( d )
		
      v3Result    = zeros(3)
      v3Result[1] = 0.5 * ( exx + eyy + d )
      v3Result[2] = 0.5 * ( exx + eyy - d )

      return v3Result
    end

    # compute Mazars equivalent strain (scalar) from the principal strains
    function getMazarsStrain(v3PStrain::Array{Float64})
      result = 0.

      for i=1:3
        xi      = v3PStrain[i]
        result += (xi > 0.) ? xi*xi : 0
      end

      return sqrt(result)
    end

    # write particle data to vtp files for postprocessing in Paraview or VisIt
    function VTKParticles(materialPoints,vtuFile)
      results_vtu = open(vtuFile, "w");
      
      # Write headers
      write(results_vtu, "<VTKFile type=\"PolyData\"  version=\"0.1\"   > \n");
      write(results_vtu, "<PolyData> \n");
      write(results_vtu, "<Piece  NumberOfPoints=\"$(length(materialPoints))\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\"> \n") 
      
      # Write point coordinates
      write(results_vtu, "<Points> \n");
      write(results_vtu, "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\" format=\"ascii\" > \n");
      
      for i=1:length(materialPoints)
          mp = materialPoints[i]
          write(results_vtu, "$(mp.v2Centroid[1]) $(mp.v2Centroid[2]) 0.0 \n")
      end
      
      write(results_vtu, "</DataArray> \n");
      write(results_vtu, "</Points> \n");
      
      # write point data
      write(results_vtu, "<PointData  Scalars=\"vonMises\" Vectors=\"sigma\"> \n");
      write(results_vtu, "<DataArray  type=\"Float64\"  Name=\"color\"  format=\"ascii\"> \n");
      for i=1:length(materialPoints)
         mp = materialPoints[i]
         write(results_vtu, "$(mp.fColor) \n");
      end
      write(results_vtu, "</DataArray> \n");
      
      write(results_vtu, "<DataArray  type=\"Float64\"  Name=\"phase field\"  format=\"ascii\"> \n");
      for i=1:length(materialPoints)
         mp = materialPoints[i]
         write(results_vtu, "$(mp.fPhase) \n");
      end
      write(results_vtu, "</DataArray> \n");
      
      write(results_vtu, "<DataArray  type=\"Float64\"  Name=\"history\"  format=\"ascii\"> \n");
      for i=1:length(materialPoints)
         mp = materialPoints[i]
         write(results_vtu, "$(mp.fHistory) \n");
      end
      write(results_vtu, "</DataArray> \n");
      
      write(results_vtu, "<DataArray  type=\"Float64\"  Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\"> \n");
      for i=1:length(materialPoints)
         mp = materialPoints[i]
         write(results_vtu, "$(mp.v2Velocity[1]) $(mp.v2Velocity[2]) 0.0\n");
      end
      write(results_vtu, "</DataArray> \n");
      
      write(results_vtu, "<DataArray  type=\"Float64\"  Name=\"sigma-xx\"  format=\"ascii\"> \n");
      for i=1:length(materialPoints)
         mp = materialPoints[i]
         write(results_vtu, "$(mp.v3Stress[1]) \n");
      end
      write(results_vtu, "</DataArray> \n");

      write(results_vtu, "</PointData> \n");
      
      # end of VTK file
      write(results_vtu, "</Piece> \n");
      write(results_vtu, "</PolyData> \n");
      write(results_vtu, "</VTKFile> \n");
      
      # close file
      close(results_vtu);
    end

end
