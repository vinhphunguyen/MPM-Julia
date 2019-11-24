module moduleBasis
	import moduleGrid #sina, do not use include here, since you have already included the module in Main.jl
	import moduleMaterialPoint

	# -------------------------------------------------------------
	# Classic functions--------------------------------------------
	function getShapeValue_Classic(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_2D_Classic, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		fShapeValue = 0.0

		v2Distance = thisMaterialPoint.v2Centroid - thisGridPoint.v2Position
		v2CellLength = thisGrid.v2Length_Cell;

		v2ShapeValue = zeros(2)
		v2ShapeValue[1] = 1.0 - abs(v2Distance[1]) / v2CellLength[1]
		v2ShapeValue[2] = 1.0 - abs(v2Distance[2]) / v2CellLength[2]

		if(v2ShapeValue[1] < 0.0)
			# @printf("Negative shape value!!! %e \n", v2Distance[1])
			v2ShapeValue[1] = 0.0
		end
		if(v2ShapeValue[2] < 0.0)
			# @printf("Negative shape value!!! \n")
			v2ShapeValue[2] = 0.0
		end

		fShapeValue = v2ShapeValue[1] * v2ShapeValue[2]
		#
		return(fShapeValue)
	end
	
    function getShapeGradient_Classic(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_2D_Classic, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		v2Result = zeros(2)

		v2Distance = thisMaterialPoint.v2Centroid - thisGridPoint.v2Position
		v2CellLength = thisGrid.v2Length_Cell;

		v2ShapeValue = zeros(2)
		v2ShapeValue[1] = 1.0 - abs(v2Distance[1]) / v2CellLength[1]
		v2ShapeValue[2] = 1.0 - abs(v2Distance[2]) / v2CellLength[2]

		if(v2ShapeValue[1] < 0.0)
			v2ShapeValue[1] = 0.0
		end
		if(v2ShapeValue[2] < 0.0)
			v2ShapeValue[2] = 0.0
		end

		v2Result[1] = -v2ShapeValue[2]*sign(v2Distance[1]) / v2CellLength[1]
		v2Result[2] = -v2ShapeValue[1]*sign(v2Distance[2]) / v2CellLength[2]

		return(v2Result)
	end

	function getShapeAndGradient_Classic(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_2D_Classic, thisGridPoint::moduleGrid.mpmGridPoint, thisGrid::moduleGrid.mpmGrid)
		v2Result = zeros(2)

		v2Distance = thisMaterialPoint.v2Centroid - thisGridPoint.v2Position
		v2CellLength = thisGrid.v2Length_Cell;

		v2ShapeValue = zeros(2)
		v2ShapeValue[1] = 1.0 - abs(v2Distance[1]) / v2CellLength[1]
		v2ShapeValue[2] = 1.0 - abs(v2Distance[2]) / v2CellLength[2]

		if(v2ShapeValue[1] < 0.0)
			v2ShapeValue[1] = 0.0
		end
		if(v2ShapeValue[2] < 0.0)
			v2ShapeValue[2] = 0.0
		end

		v2Result[1] = -v2ShapeValue[2]*sign(v2Distance[1]) / v2CellLength[1]
		v2Result[2] = -v2ShapeValue[1]*sign(v2Distance[2]) / v2CellLength[2]

        return(v2ShapeValue[1]*v2ShapeValue[2],v2Result)
	end
end
