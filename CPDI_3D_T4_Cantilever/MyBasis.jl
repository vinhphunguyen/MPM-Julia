module moduleBasis
	import moduleGrid #sina, do not use include here, since you have already included the module in Main.jl
	import moduleMaterialPoint

	# -------------------------------------------------------------
	# Classic functions--------------------------------------------
	function getShapeValue_Classic(v3Distance::Array{Float64,1}, v3CellLength::Array{Float64,1})
		fShapeValue = 0.0

		# v3ShapeValue = ones(3) - abs(v3Distance) ./ v3CellLength

		v3ShapeValue = zeros(3)
		v3ShapeValue[1] = 1.0 - abs(v3Distance[1]) / v3CellLength[1]
		v3ShapeValue[2] = 1.0 - abs(v3Distance[2]) / v3CellLength[2]
		v3ShapeValue[3] = 1.0 - abs(v3Distance[3]) / v3CellLength[3]

		if(v3ShapeValue[1] < 0.0)
			v3ShapeValue[1] = 0.0
		end
		if(v3ShapeValue[2] < 0.0)
			v3ShapeValue[2] = 0.0
		end
		if(v3ShapeValue[3] < 0.0)
			v3ShapeValue[3] = 0.0
		end

		fShapeValue = v3ShapeValue[1] * v3ShapeValue[2] * v3ShapeValue[3]
		#
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
	# -------------------------------------------------------------
	# CPDI_Tet4 functions-----------------------------------------------
	function getShapeValue_CPDI_Tet4(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_Tet4, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = 0.0

		fWeight = [0.25 0.25 0.25 0.25]

		fShapeValue = [0.0 0.0 0.0 0.0]

		for indexCorner = 1:1:4
			v3CornerCoordinate = thisMaterialPoint.v3Corner[:, indexCorner]
			v3Distance = v3CornerCoordinate - thisGridPoint.v3Position
			v3CellLength = thisGrid.v3Length_Cell

			fShapeValue[indexCorner] = moduleBasis.getShapeValue_Classic(v3Distance, v3CellLength)
		end

		for indexCorner = 1:1:4
			fResult += fWeight[indexCorner]*fShapeValue[indexCorner]
		end

		return(fResult)
	end
	function getShapeGradient_CPDI_Tet4(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_Tet4, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fResult = zeros(3)

		x = vec(thisMaterialPoint.v3Corner[1,:])
		y = vec(thisMaterialPoint.v3Corner[2,:])
		z = vec(thisMaterialPoint.v3Corner[3,:])

		x21 = (x[2] - x[1]); x23 = (x[2] - x[3]); x24 = (x[2] - x[4])
		x12 = (x[1] - x[2]); x13 = (x[1] - x[3]); x14 = (x[1] - x[4])
		x31 = (x[3] - x[1]); x32 = (x[3] - x[2]); x34 = (x[3] - x[4])
		x41 = (x[4] - x[1]); x42 = (x[4] - x[2]); x43 = (x[4] - x[3])

		y12 = (y[1] - y[2]); y13 = (y[1] - y[3]); y14 = (y[1] - y[4])
		y21 = (y[2] - y[1]); y23 = (y[2] - y[3]); y24 = (y[2] - y[4])
		y31 = (y[3] - y[1]); y32 = (y[3] - y[2]); y34 = (y[3] - y[4])
		y41 = (y[4] - y[1]); y42 = (y[4] - y[2]); y43 = (y[4] - y[3])

		z12 = (z[1] - z[2]); z13 = (z[1] - z[3]); z14 = (z[1] - z[4])
		z21 = (z[2] - z[1]); z23 = (z[2] - z[3]); z24 = (z[2] - z[4])
		z31 = (z[3] - z[1]); z32 = (z[3] - z[2]); z34 = (z[3] - z[4])
		z41 = (z[4] - z[1]); z42 = (z[4] - z[2]); z43 = (z[4] - z[3])

		fV = 0.0
		fV += x21*(y23*z34 - y34*z23)
		fV += x32*(y34*z12 - y12*z34)
		fV += x43*(y12*z23 - y23*z12)
		fV *= 1.0;

		fa = zeros(4)
		fa[1] = y42*z32 - y32*z42
		fa[2] = y31*z43 - y34*z13
		fa[3] = y24*z14 - y14*z24
		fa[4] = y13*z21 - y12*z31

		fb = zeros(4)
		fb[1] = x32*z42 - x42*z32
		fb[2] = x43*z31 - x13*z34
		fb[3] = x14*z24 - x24*z14
		fb[4] = x21*z13 - x31*z12

		fc = zeros(4)
		fc[1] = x42*y32 - x32*y42
		fc[2] = x31*y43 - x34*y13
		fc[3] = x24*y14 - x14*y24
		fc[4] = x13*y21 - x12*y31

		v3Weight = zeros(3,4)
		v3Weight[1,:] = fa
		v3Weight[2,:] = fb
		v3Weight[3,:] = fc

		fShapeValue = zeros(4)

		for indexCorner = 1:1:4
			v3CornerCoordinate = thisMaterialPoint.v3Corner[:, indexCorner]
			v3Distance = v3CornerCoordinate - thisGridPoint.v3Position
			v3CellLength = thisGrid.v3Length_Cell

			fShapeValue[indexCorner] = moduleBasis.getShapeValue_Classic(v3Distance, v3CellLength)
		end

		for indexCorner = 1:1:4
			fResult += fShapeValue[indexCorner] * v3Weight[:,indexCorner]
			# fResult += fWeight[iIndex_Corner]*fShapeValue[iIndex_Corner]
		end

		fResult /= (1.0*fV)

		return(fResult)
	end
end
