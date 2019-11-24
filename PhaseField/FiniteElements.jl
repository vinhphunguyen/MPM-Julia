
import FeMesh 
import moduleMaterialPoint

function getN(x::Array{Float64}, xI::Array{Float64}, dx::Float64, dy::Float64)
		fShapeValue   = 0.0
		v2Distance    = x - xI
		v2ShapeValue1 = 1.0 - abs(v2Distance[1]) / dx
		v2ShapeValue2 = 1.0 - abs(v2Distance[2]) / dy

		if(v2ShapeValue1 < 0.0)
			# @printf("Negative shape value!!! %e \n", v2Distance[1])
			v2ShapeValue1 = 0.0
		end
		if(v2ShapeValue2 < 0.0)
			# @printf("Negative shape value!!! \n")
			v2ShapeValue2 = 0.0
		end

		fShapeValue = v2ShapeValue1 * v2ShapeValue2
		#
		return(fShapeValue)
end

function getdNdx(x::Array{Float64}, xI::Array{Float64}, dx::Float64, dy::Float64)
		v2Result = zeros(2)

		v2Distance   = x - xI
		v2ShapeValue1 = 1.0 - abs(v2Distance[1]) / dx
		v2ShapeValue2 = 1.0 - abs(v2Distance[2]) / dy

		if(v2ShapeValue1 < 0.0)
			v2ShapeValue1 = 0.0
		end
		if(v2ShapeValue2 < 0.0)
			v2ShapeValue2 = 0.0
		end

		v2Result[1] = -v2ShapeValue2*sign(v2Distance[1]) / dx
		v2Result[2] = -v2ShapeValue1*sign(v2Distance[2]) / dy

		return(v2Result)
end

function getNdNdx(x::Array{Float64}, xI::Array{Float64}, dx::Float64, dy::Float64)
		v2Result = zeros(2)

		v2Distance   = x - xI
		v2ShapeValue1 = 1.0 - abs(v2Distance[1]) / dx
		v2ShapeValue2 = 1.0 - abs(v2Distance[2]) / dy

		if(v2ShapeValue1 < 0.0)
			v2ShapeValue1 = 0.0
		end
		if(v2ShapeValue2 < 0.0)
			v2ShapeValue2 = 0.0
		end

		v2Result[1] = -v2ShapeValue2*sign(v2Distance[1]) / dx
		v2Result[2] = -v2ShapeValue1*sign(v2Distance[2]) / dy

		return (v2ShapeValue1*v2ShapeValue2, v2Result)
end

function fe_matrices(mesh::FeMesh.Mesh,
                     materialPoints::Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic},
                    F, l0, gc)
  iu  = zeros(Int,20*mesh.nodeCount)
  ju  = zeros(Int,20*mesh.nodeCount)
  vu  = zeros(    20*mesh.nodeCount)

  ke  = zeros(4,4)
  fe  = zeros(4)
  sparseCount = 0

  real_dofs = IntSet()

  for iE=1:mesh.elemCount
    if haskey(mesh.elem2MP,iE) == false continue end  # skip elements with no material points
    materialPointIds = mesh.elem2MP[iE]
    fill!(ke,0.0)
    fill!(fe,0.0)
    #print(materialPointIds)
    #@printf("\n")
    for iP=1:length(materialPointIds)
      thisMP = materialPoints[materialPointIds[iP]]
      xp     = thisMP.v2Centroid
      H      = thisMP.fHistory
      Vp     = thisMP.fVolume
      for i=1:4
          ni       = mesh.elements[i,iE]
          xi       = mesh.nodes[:,ni]
          Ni,dNi   = getNdNdx(xp, xi, mesh.deltaX, mesh.deltaY)
          fe[i]   += 2. * Ni * H * Vp
          push!(real_dofs, ni)
        for j=1:4
          nj       = mesh.elements[j,iE]
          xj       = mesh.nodes[:,nj]
          Nj,dNj   = getNdNdx(xp, xj,mesh.deltaX,mesh.deltaY)
          ke[i,j] += Vp * ( gc*l0*dot(dNi,dNj) + gc/l0*Ni*Nj )
        end
      end
    end

    for i=1:4
      ni       = mesh.elements[i,iE]
      F[ni]   += fe[i]
      for j=1:4
          nj   = mesh.elements[j,iE]
          sparseCount += 1
          try
            iu[sparseCount] = ni;
            ju[sparseCount] = nj;
            vu[sparseCount] = ke[i,j];
          catch BoundsError()
            iu = vcat(iu,ni)
            ju = vcat(ju,nj)
            vu = vcat(vu,ke[i,j])
          end
       end
    end
  end

  iu = iu[1:sparseCount]
  ju = ju[1:sparseCount]
  vu = vu[1:sparseCount]

  # convert from integer set to array to be used in sparse matrix later
  dofs = collect(real_dofs)

  return iu,ju,vu,dofs
end



function solve_fe(phi::Array{Float64}, mesh::FeMesh.Mesh,
                  materialPoints::Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic},
                 l0, Gc)
 
  f                  = zeros(mesh.nodeCount) 
  iu,ju,vu,d_free    = fe_matrices(mesh,materialPoints,f, l0, Gc)
  K                  = sparse(iu,ju,vu,mesh.nodeCount,mesh.nodeCount)

  #@printf("sum of force: %.6e \n",sum(f))
  #= applying boundary conditions phi=1
  udofs    = collect( Int64((mesh.elemCountX+1)*(mesh.elemCountY/2)+1):
                      Int64((mesh.elemCountX+1)*Int64((mesh.elemCountY/2))+1+mesh.elemCountX/2) )
  uFix     = ones(length(udofs))

  bcwt     = mean(diag(K))

  f=f-K[:,udofs]*uFix
  K[udofs,:]=0.
  K[:,udofs]=0.
  K[udofs,udofs]=bcwt*speye(length(udofs))
  f[udofs]=bcwt*speye(length(udofs))*uFix
  
  U     = K\f;
  =#

  u_free        = K[d_free, d_free] \ f[d_free]
  phi[d_free]   = u_free 

  #print(phi)
end


