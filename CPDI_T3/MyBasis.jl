module moduleBasis
   import moduleMath #sina, do not use include here, since you have already included the module in Main.jl
	import moduleGrid
	import moduleMaterialPoint

	# -------------------------------------------------------------
	# Template functions-------------------------------------------
	function getShapeValue(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		# fShapeValue = getShapeValue_cpGIMP(thisMaterialPoint, thisGridPoint, thisGrid)
		# fShapeValue = getShapeValue_CPDI(thisMaterialPoint, thisGridPoint, thisGrid)
		# fShapeValue = getShapeValue_CPDI_Q4(thisMaterialPoint, thisGridPoint, thisGrid)
		fShapeValue = getShapeValue_CPDI_T3(thisMaterialPoint, thisGridPoint, thisGrid)
		return(fShapeValue)
	end
	function getShapeGradient_x(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		# fShapeGradient_x = getShapeGradient_cpGIMP_x(thisMaterialPoint, thisGridPoint, thisGrid)
		# fShapeGradient_x = getShapeGradient_CPDI_x(thisMaterialPoint, thisGridPoint, thisGrid)
		# fShapeGradient_x = getShapeGradient_CPDI_Q4_x(thisMaterialPoint, thisGridPoint, thisGrid)
		fShapeGradient_x = getShapeGradient_CPDI_T3_x(thisMaterialPoint, thisGridPoint, thisGrid)
		return(fShapeGradient_x)
	end
	function getShapeGradient_y(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		# fShapeGradient_y = getShapeGradient_cpGIMP_y(thisMaterialPoint, thisGridPoint, thisGrid)
		# fShapeGradient_y = getShapeGradient_CPDI_y(thisMaterialPoint, thisGridPoint, thisGrid)
		# fShapeGradient_y = getShapeGradient_CPDI_Q4_y(thisMaterialPoint, thisGridPoint, thisGrid)
		fShapeGradient_y = getShapeGradient_CPDI_T3_y(thisMaterialPoint, thisGridPoint, thisGrid)
		return(fShapeGradient_y)
	end
	# -------------------------------------------------------------
	# Classic functions--------------------------------------------
	function getShapeValue_Classic(fDistance_x::Real, fDistance_y::Real, fCell_Length_x::Real, fCell_Length_y::Real)
		fShapeValue_x = 1.0 - abs(fDistance_x) / fCell_Length_x
		if(fShapeValue_x < 0.0)
			fShapeValue_x = 0.0
		end

		fShapeValue_y = 1.0 - abs(fDistance_y) / fCell_Length_y
		if(fShapeValue_y < 0.0)
			fShapeValue_y = 0.0
		end

		fShapeValue = fShapeValue_x * fShapeValue_y

		return(fShapeValue)
	end

	function getShapeGradient_Classic_x(fDistance_x::Real, fDistance_y::Real, fCell_Length_x::Real, fCell_Length_y::Real)
		fShapeValue_y = 1.0 - abs(fDistance_y) / fCell_Length_y
		if(fShapeValue_y < 0.0)
			fShapeValue_y = 0.0
		end

		fShapeGradient_x = -fShapeValue_y*sign(fDistance_x) / fCell_Length_x

		return(fShapeGradient_x)
	end

	function getShapeGradient_Classic_y(fDistance_x::Real, fDistance_y::Real, fCell_Length_x::Real, fCell_Length_y::Real)
		fShapeValue_x = 1.0 - abs(fDistance_x) / fCell_Length_x
		if(fShapeValue_x < 0.0)
			fShapeValue_x = 0.0
		end

		fShapeGradient_y = -fShapeValue_x*sign(fDistance_y) / fCell_Length_y

		return(fShapeGradient_y)
	end
	# cpGimp functions-----------------------------------------------
	# note: input  distance values should not be absolute value
	function getShapeValue_cpGIMP_x(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fDistance_x = thisMaterialPoint.v2Position.fx - thisGridPoint.v2Position.fx
		fDistance_y = thisMaterialPoint.v2Position.fy - thisGridPoint.v2Position.fy

		fParticle_Length_x = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1)[1]
		fParticle_Length_y = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2)[2]

		fCell_Length_x = thisGrid.fLength_Cell_x
		fCell_Length_y = thisGrid.fLength_Cell_y

		if(abs(fDistance_x) < 0.5*fParticle_Length_x)
			fResult = 1.0 - (4.0*fDistance_x^2 + fParticle_Length_x^2) / (4.0*fCell_Length_x*fParticle_Length_x)
		elseif(0.5*fParticle_Length_x <= abs(fDistance_x) <= fCell_Length_x - 0.5*fParticle_Length_x)
			fResult = 1.0 - abs(fDistance_x)/fCell_Length_x
		elseif(fCell_Length_x - 0.5*fParticle_Length_x <= abs(fDistance_x) < fCell_Length_x + 0.5*fParticle_Length_x)
			fResult = ((fCell_Length_x + 0.5*fParticle_Length_x - abs(fDistance_x))^2) / (2.0*fParticle_Length_x*fCell_Length_x)
		else
			fResult = 0.0
		end

		return(fResult)
	end

	function getShapeValue_cpGIMP_y(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fDistance_x = thisMaterialPoint.v2Position.fx - thisGridPoint.v2Position.fx
		fDistance_y = thisMaterialPoint.v2Position.fy - thisGridPoint.v2Position.fy

		fParticle_Length_x = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1)[1]
		fParticle_Length_y = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2)[2]

		fCell_Length_x = thisGrid.fLength_Cell_x
		fCell_Length_y = thisGrid.fLength_Cell_y

		if(abs(fDistance_y) < 0.5*fParticle_Length_y)
			fResult = 1.0 - (4.0*fDistance_y^2 + fParticle_Length_y^2) / (4.0*fCell_Length_y*fParticle_Length_y)
		elseif(0.5*fParticle_Length_y <= abs(fDistance_y) <= fCell_Length_y - 0.5*fParticle_Length_y)
			fResult = 1.0 - abs(fDistance_y)/fCell_Length_y
		elseif(fCell_Length_y - 0.5*fParticle_Length_y <= abs(fDistance_y) < fCell_Length_y + 0.5*fParticle_Length_y)
			fResult = ((fCell_Length_y + 0.5*fParticle_Length_y - abs(fDistance_y))^2) / (2.0*fParticle_Length_y*fCell_Length_y)
		else
			fResult = 0.0
		end

		return(fResult)
	end

	function getShapeValue_cpGIMP(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fShapeValue_x = getShapeValue_cpGIMP_x(thisMaterialPoint, thisGridPoint, thisGrid)
		fShapeValue_y = getShapeValue_cpGIMP_y(thisMaterialPoint, thisGridPoint, thisGrid)

		return(fShapeValue_x*fShapeValue_y)
	end

	function getShapeGradient_cpGIMP_x(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fDistance_x = thisMaterialPoint.v2Position.fx - thisGridPoint.v2Position.fx
		fDistance_y = thisMaterialPoint.v2Position.fy - thisGridPoint.v2Position.fy

		fParticle_Length_x = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1)[1]
		fParticle_Length_y = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2)[2]

		fCell_Length_x = thisGrid.fLength_Cell_x
		fCell_Length_y = thisGrid.fLength_Cell_y

		fShapeValue_y = getShapeValue_cpGIMP_y(thisMaterialPoint, thisGridPoint, thisGrid)

		if(abs(fDistance_x) < 0.5*fParticle_Length_x)
			fResult = (-8.0*fDistance_x)/(4.0*fCell_Length_x*fParticle_Length_x)
		elseif(0.5*fParticle_Length_x <= abs(fDistance_x) <= fCell_Length_x - 0.5*fParticle_Length_x)
			fResult = (-1.0/fCell_Length_x)*sign(fDistance_x)
		elseif(fCell_Length_x - 0.5*fParticle_Length_x <= abs(fDistance_x) < fCell_Length_x + 0.5*fParticle_Length_x)
			fResult = -sign(fDistance_x)*(fCell_Length_x+0.5*fParticle_Length_x-abs(fDistance_x))/(fCell_Length_x*fParticle_Length_x)
		else
			fResult = 0.0
		end

		fResult *= fShapeValue_y

		return(fResult)
	end

	function getShapeGradient_cpGIMP_y(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fDistance_x = thisMaterialPoint.v2Position.fx - thisGridPoint.v2Position.fx
		fDistance_y = thisMaterialPoint.v2Position.fy - thisGridPoint.v2Position.fy

		fParticle_Length_x = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1)[1]
		fParticle_Length_y = 2.0*(thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2)[2]

		fCell_Length_x = thisGrid.fLength_Cell_x
		fCell_Length_y = thisGrid.fLength_Cell_y

		fShapeValue_x = getShapeValue_cpGIMP_x(thisMaterialPoint, thisGridPoint, thisGrid)

		if(abs(fDistance_y) < 0.5*fParticle_Length_y)
			fResult = (-8.0*fDistance_y)/(4.0*fCell_Length_y*fParticle_Length_y)
		elseif(0.5*fParticle_Length_y <= abs(fDistance_y) < fCell_Length_y - 0.5*fParticle_Length_y)
			fResult = (-1.0/fCell_Length_y)*sign(fDistance_y)
		elseif(fCell_Length_y - 0.5*fParticle_Length_y <= abs(fDistance_y) < fCell_Length_y + 0.5*fParticle_Length_y)
			fResult = -sign(fDistance_y)*(fCell_Length_y+0.5*fParticle_Length_y-abs(fDistance_y))/(fCell_Length_y*fParticle_Length_y)
		else
			fResult = 0.0
		end

		fResult *= fShapeValue_x

		return(fResult)
	end
	# -------------------------------------------------------------
	# CPDI_R4 functions-----------------------------------------------
	function getShapeValue_CPDI(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fP_x = thisMaterialPoint.v2Position.fx
		fP_y = thisMaterialPoint.v2Position.fy

		fN_x = thisGridPoint.v2Position.fx
		fN_y = thisGridPoint.v2Position.fy

		mR1 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1
		mR2 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2

		fCell_Length_x = thisGrid.fLength_Cell_x
		fCell_Length_y = thisGrid.fLength_Cell_y

		# lower left corner
		fS_1 = getShapeValue_Classic((fP_x-mR1[1]-mR2[1])-fN_x, (fP_y-mR1[2]-mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# lower right corner
		fS_2 = getShapeValue_Classic((fP_x+mR1[1]-mR2[1])-fN_x, (fP_y+mR1[2]-mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# upper right corner
		fS_3 = getShapeValue_Classic((fP_x+mR1[1]+mR2[1])-fN_x, (fP_y+mR1[2]+mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# upper left corner
		fS_4 = getShapeValue_Classic((fP_x-mR1[1]+mR2[1])-fN_x, (fP_y-mR1[2]+mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)

		fResult = 0.25 * (fS_1 + fS_2 + fS_3 + fS_4)

		return(fResult)
	end

	function getShapeGradient_CPDI_x(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fShapeGradient_x = 0.0

		fP_x = thisMaterialPoint.v2Position.fx
		fP_y = thisMaterialPoint.v2Position.fy

		fN_x = thisGridPoint.v2Position.fx
		fN_y = thisGridPoint.v2Position.fy

		mR1 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1
		mR2 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2

		fCell_Length_x = thisGrid.fLength_Cell_x
		fCell_Length_y = thisGrid.fLength_Cell_y

		# lower left corner
		fS_1 = getShapeValue_Classic((fP_x-mR1[1]-mR2[1])-fN_x, (fP_y-mR1[2]-mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# lower right corner
		fS_2 = getShapeValue_Classic((fP_x+mR1[1]-mR2[1])-fN_x, (fP_y+mR1[2]-mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# upper right corner
		fS_3 = getShapeValue_Classic((fP_x+mR1[1]+mR2[1])-fN_x, (fP_y+mR1[2]+mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# upper left corner
		fS_4 = getShapeValue_Classic((fP_x-mR1[1]+mR2[1])-fN_x, (fP_y-mR1[2]+mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)

		fShapeGradient_x = ((fS_1 - fS_3)*(mR1[2] - mR2[2]) + (fS_2 - fS_4)*(mR1[2] + mR2[2])) / (1.0*thisMaterialPoint.fVolume)

		return(fShapeGradient_x)
	end

	function getShapeGradient_CPDI_y(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fShapeGradient_y = 0.0

		fP_x = thisMaterialPoint.v2Position.fx
		fP_y = thisMaterialPoint.v2Position.fy

		fN_x = thisGridPoint.v2Position.fx
		fN_y = thisGridPoint.v2Position.fy

		mR1 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial1
		mR2 = thisMaterialPoint.mDeformationGradient * thisMaterialPoint.mRadial2

		fCell_Length_x = thisGrid.fLength_Cell_x
		fCell_Length_y = thisGrid.fLength_Cell_y

		# lower left corner
		fS_1 = getShapeValue_Classic((fP_x-mR1[1]-mR2[1])-fN_x, (fP_y-mR1[2]-mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# lower right corner
		fS_2 = getShapeValue_Classic((fP_x+mR1[1]-mR2[1])-fN_x, (fP_y+mR1[2]-mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# upper right corner
		fS_3 = getShapeValue_Classic((fP_x+mR1[1]+mR2[1])-fN_x, (fP_y+mR1[2]+mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)
		# upper left corner
		fS_4 = getShapeValue_Classic((fP_x-mR1[1]+mR2[1])-fN_x, (fP_y-mR1[2]+mR2[2])-fN_y, fCell_Length_x, fCell_Length_y)

		fShapeGradient_y = ((fS_1 - fS_3)*(mR2[1] - mR1[1]) + (fS_2 - fS_4)*(-mR1[1] - mR2[1])) / (1.0*thisMaterialPoint.fVolume)

		return(fShapeGradient_y)
	end
	# -------------------------------------------------------------
	# CPDI_Q4 functions-----------------------------------------------
	function getShapeValue_CPDI_Q4(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fa = (thisMaterialPoint.mCorner[4,1] - thisMaterialPoint.mCorner[1,1])*(thisMaterialPoint.mCorner[2,2] - thisMaterialPoint.mCorner[3,2]) - (thisMaterialPoint.mCorner[2,1] - thisMaterialPoint.mCorner[3,1])*(thisMaterialPoint.mCorner[4,2] - thisMaterialPoint.mCorner[1,2])
		fb = (thisMaterialPoint.mCorner[3,1] - thisMaterialPoint.mCorner[4,1])*(thisMaterialPoint.mCorner[1,2] - thisMaterialPoint.mCorner[2,2]) - (thisMaterialPoint.mCorner[1,1] - thisMaterialPoint.mCorner[2,1])*(thisMaterialPoint.mCorner[3,2] - thisMaterialPoint.mCorner[4,2])
		fV = 0.0
		fV += (thisMaterialPoint.mCorner[1,1]*thisMaterialPoint.mCorner[2,2] - thisMaterialPoint.mCorner[2,1]*thisMaterialPoint.mCorner[1,2])
		fV += (thisMaterialPoint.mCorner[2,1]*thisMaterialPoint.mCorner[3,2] - thisMaterialPoint.mCorner[3,1]*thisMaterialPoint.mCorner[2,2])
		fV += (thisMaterialPoint.mCorner[3,1]*thisMaterialPoint.mCorner[4,2] - thisMaterialPoint.mCorner[4,1]*thisMaterialPoint.mCorner[3,2])
		fV += (thisMaterialPoint.mCorner[4,1]*thisMaterialPoint.mCorner[1,2] - thisMaterialPoint.mCorner[1,1]*thisMaterialPoint.mCorner[4,2])
		fV *= 0.5

		fWeight = Array{Real}(size(thisMaterialPoint.mCorner, 1))
		fWeight[1] = 6.0*fV - fa - fb
		fWeight[2] = 6.0*fV - fa + fb
		fWeight[3] = 6.0*fV + fa + fb
		fWeight[4] = 6.0*fV + fa - fb

		fShapeValue = Array{Real}(size(thisMaterialPoint.mCorner, 1))

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			thisCorner = moduleMath.Vector2D(0.0, 0.0)
			thisCorner.fx = thisMaterialPoint.mCorner[iIndex_Corner, 1]
			thisCorner.fy = thisMaterialPoint.mCorner[iIndex_Corner, 2]

			fDistance_x = thisCorner.fx - thisGridPoint.v2Position.fx
			fDistance_y = thisCorner.fy - thisGridPoint.v2Position.fy

			fShapeValue[iIndex_Corner] = moduleBasis.getShapeValue_Classic(fDistance_x, fDistance_y, thisGrid.fLength_Cell_x, thisGrid.fLength_Cell_y)
		end

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			fResult += fWeight[iIndex_Corner]*fShapeValue[iIndex_Corner]
		end

		fResult /= (24.0*fV)

		return(fResult)
	end

	function getShapeGradient_CPDI_Q4_x(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fV = 0.0
		fV += (thisMaterialPoint.mCorner[1,1]*thisMaterialPoint.mCorner[2,2] - thisMaterialPoint.mCorner[2,1]*thisMaterialPoint.mCorner[1,2])
		fV += (thisMaterialPoint.mCorner[2,1]*thisMaterialPoint.mCorner[3,2] - thisMaterialPoint.mCorner[3,1]*thisMaterialPoint.mCorner[2,2])
		fV += (thisMaterialPoint.mCorner[3,1]*thisMaterialPoint.mCorner[4,2] - thisMaterialPoint.mCorner[4,1]*thisMaterialPoint.mCorner[3,2])
		fV += (thisMaterialPoint.mCorner[4,1]*thisMaterialPoint.mCorner[1,2] - thisMaterialPoint.mCorner[1,1]*thisMaterialPoint.mCorner[4,2])
		fV *= 0.5

		fWeight = Array{Real}(size(thisMaterialPoint.mCorner, 1))
		fWeight[1] = thisMaterialPoint.mCorner[2,2] - thisMaterialPoint.mCorner[4,2]
		fWeight[2] = thisMaterialPoint.mCorner[3,2] - thisMaterialPoint.mCorner[1,2]
		fWeight[3] = thisMaterialPoint.mCorner[4,2] - thisMaterialPoint.mCorner[2,2]
		fWeight[4] = thisMaterialPoint.mCorner[1,2] - thisMaterialPoint.mCorner[3,2]

		fShapeValue = Array{Real}(size(thisMaterialPoint.mCorner, 1))

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			thisCorner = moduleMath.Vector2D(0.0, 0.0)
			thisCorner.fx = thisMaterialPoint.mCorner[iIndex_Corner, 1]
			thisCorner.fy = thisMaterialPoint.mCorner[iIndex_Corner, 2]

			fDistance_x = thisCorner.fx - thisGridPoint.v2Position.fx
			fDistance_y = thisCorner.fy - thisGridPoint.v2Position.fy

			fShapeValue[iIndex_Corner] = moduleBasis.getShapeValue_Classic(fDistance_x, fDistance_y, thisGrid.fLength_Cell_x, thisGrid.fLength_Cell_y)
		end

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			fResult += fWeight[iIndex_Corner]*fShapeValue[iIndex_Corner]
		end

		fResult /= (2.0*fV)

		return(fResult)
	end

	function getShapeGradient_CPDI_Q4_y(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fV = 0.0
		fV += (thisMaterialPoint.mCorner[1,1]*thisMaterialPoint.mCorner[2,2] - thisMaterialPoint.mCorner[2,1]*thisMaterialPoint.mCorner[1,2])
		fV += (thisMaterialPoint.mCorner[2,1]*thisMaterialPoint.mCorner[3,2] - thisMaterialPoint.mCorner[3,1]*thisMaterialPoint.mCorner[2,2])
		fV += (thisMaterialPoint.mCorner[3,1]*thisMaterialPoint.mCorner[4,2] - thisMaterialPoint.mCorner[4,1]*thisMaterialPoint.mCorner[3,2])
		fV += (thisMaterialPoint.mCorner[4,1]*thisMaterialPoint.mCorner[1,2] - thisMaterialPoint.mCorner[1,1]*thisMaterialPoint.mCorner[4,2])
		fV *= 0.5

		fWeight = Array{Real}(size(thisMaterialPoint.mCorner, 1))
		fWeight[1] = thisMaterialPoint.mCorner[4,1] - thisMaterialPoint.mCorner[2,1]
		fWeight[2] = thisMaterialPoint.mCorner[1,1] - thisMaterialPoint.mCorner[3,1]
		fWeight[3] = thisMaterialPoint.mCorner[2,1] - thisMaterialPoint.mCorner[4,1]
		fWeight[4] = thisMaterialPoint.mCorner[3,1] - thisMaterialPoint.mCorner[1,1]

		fShapeValue = Array{Real}(size(thisMaterialPoint.mCorner, 1))

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			thisCorner = moduleMath.Vector2D(0.0, 0.0)
			thisCorner.fx = thisMaterialPoint.mCorner[iIndex_Corner, 1]
			thisCorner.fy = thisMaterialPoint.mCorner[iIndex_Corner, 2]

			fDistance_x = thisCorner.fx - thisGridPoint.v2Position.fx
			fDistance_y = thisCorner.fy - thisGridPoint.v2Position.fy

			fShapeValue[iIndex_Corner] = moduleBasis.getShapeValue_Classic(fDistance_x, fDistance_y, thisGrid.fLength_Cell_x, thisGrid.fLength_Cell_y)
		end

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			fResult += fWeight[iIndex_Corner]*fShapeValue[iIndex_Corner]
		end

		fResult /= (2.0*fV)

		return(fResult)
	end
	# -------------------------------------------------------------
	# CPDI_T3 functions-----------------------------------------------
	function getShapeValue_CPDI_T3(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fWeight = Array{Real}(size(thisMaterialPoint.mCorner, 1))
		fWeight[1] = 1.0/3.0
		fWeight[2] = 1.0/3.0
		fWeight[3] = 1.0/3.0

		fShapeValue = Array{Real}(size(thisMaterialPoint.mCorner, 1))

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			thisCorner = moduleMath.Vector2D(0.0, 0.0)
			thisCorner.fx = thisMaterialPoint.mCorner[iIndex_Corner, 1]
			thisCorner.fy = thisMaterialPoint.mCorner[iIndex_Corner, 2]

			fDistance_x = thisCorner.fx - thisGridPoint.v2Position.fx
			fDistance_y = thisCorner.fy - thisGridPoint.v2Position.fy

			fShapeValue[iIndex_Corner] = moduleBasis.getShapeValue_Classic(fDistance_x, fDistance_y, thisGrid.fLength_Cell_x, thisGrid.fLength_Cell_y)
		end

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			fResult += fWeight[iIndex_Corner]*fShapeValue[iIndex_Corner]
		end

		return(fResult)
	end

	function getShapeGradient_CPDI_T3_x(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fV = 0.0
		fV += (thisMaterialPoint.mCorner[1,1]-thisMaterialPoint.mCorner[3,1]) * (thisMaterialPoint.mCorner[2,2]-thisMaterialPoint.mCorner[1,2])
		fV -= (thisMaterialPoint.mCorner[1,1]-thisMaterialPoint.mCorner[2,1]) * (thisMaterialPoint.mCorner[3,2]-thisMaterialPoint.mCorner[1,2])
		fV *= 0.5

		fWeight = Array{Real}(size(thisMaterialPoint.mCorner, 1))
		fWeight[1] = thisMaterialPoint.mCorner[2,2] - thisMaterialPoint.mCorner[3,2]
		fWeight[2] = thisMaterialPoint.mCorner[3,2] - thisMaterialPoint.mCorner[1,2]
		fWeight[3] = thisMaterialPoint.mCorner[1,2] - thisMaterialPoint.mCorner[2,2]

		fShapeValue = Array{Real}(size(thisMaterialPoint.mCorner, 1))

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			thisCorner = moduleMath.Vector2D(0.0, 0.0)
			thisCorner.fx = thisMaterialPoint.mCorner[iIndex_Corner, 1]
			thisCorner.fy = thisMaterialPoint.mCorner[iIndex_Corner, 2]

			fDistance_x = thisCorner.fx - thisGridPoint.v2Position.fx
			fDistance_y = thisCorner.fy - thisGridPoint.v2Position.fy

			fShapeValue[iIndex_Corner] = moduleBasis.getShapeValue_Classic(fDistance_x, fDistance_y, thisGrid.fLength_Cell_x, thisGrid.fLength_Cell_y)
		end

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			fResult += fWeight[iIndex_Corner]*fShapeValue[iIndex_Corner]
		end

		fResult /= (2.0*fV)

		return(fResult)
	end

	function getShapeGradient_CPDI_T3_y(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fV = 0.0
		fV += (thisMaterialPoint.mCorner[1,1]-thisMaterialPoint.mCorner[3,1]) * (thisMaterialPoint.mCorner[2,2]-thisMaterialPoint.mCorner[1,2])
		fV -= (thisMaterialPoint.mCorner[1,1]-thisMaterialPoint.mCorner[2,1]) * (thisMaterialPoint.mCorner[3,2]-thisMaterialPoint.mCorner[1,2])
		fV *= 0.5

		fWeight = Array{Real}(size(thisMaterialPoint.mCorner, 1))
		fWeight[1] = thisMaterialPoint.mCorner[3,1] - thisMaterialPoint.mCorner[2,1]
		fWeight[2] = thisMaterialPoint.mCorner[1,1] - thisMaterialPoint.mCorner[3,1]
		fWeight[3] = thisMaterialPoint.mCorner[2,1] - thisMaterialPoint.mCorner[1,1]

		fShapeValue = Array{Real}(size(thisMaterialPoint.mCorner, 1))

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			thisCorner = moduleMath.Vector2D(0.0, 0.0)
			thisCorner.fx = thisMaterialPoint.mCorner[iIndex_Corner, 1]
			thisCorner.fy = thisMaterialPoint.mCorner[iIndex_Corner, 2]

			fDistance_x = thisCorner.fx - thisGridPoint.v2Position.fx
			fDistance_y = thisCorner.fy - thisGridPoint.v2Position.fy

			fShapeValue[iIndex_Corner] = moduleBasis.getShapeValue_Classic(fDistance_x, fDistance_y, thisGrid.fLength_Cell_x, thisGrid.fLength_Cell_y)
		end

		for iIndex_Corner = 1:1:size(thisMaterialPoint.mCorner, 1)
			fResult += fWeight[iIndex_Corner]*fShapeValue[iIndex_Corner]
		end

		fResult /= (2.0*fV)

		return(fResult)
	end
end
