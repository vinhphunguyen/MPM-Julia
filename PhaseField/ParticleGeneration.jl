module moduleParticleGen

   using  moduleGrid
   using  moduleMaterialPoint
    
   #export createMaterialDomain_Rectangle, createMaterialDomain_Circle


   function createMaterialDomain_Circle(fCenter::Array{Float64}, 
                                        fRadius::Float64, 
                                        grid,
                                        ppc::Array{Int64})

     xc   = fCenter[1]
     yc   = fCenter[2]

     xmin = xc - fRadius
     xmax = xc + fRadius
     ymin = yc - fRadius
     ymax = yc + fRadius

     minI, minJ = moduleGrid.point2ElemIndexIJ([xmin;ymin], grid)
     maxI, maxJ = moduleGrid.point2ElemIndexIJ([xmax;ymax], grid)

     deltaX     = grid.v2Length_Cell[1]
     deltaY     = grid.v2Length_Cell[2]

     dx         = deltaX/ppc[1]
     dy         = deltaY/ppc[2]

     thisMaterialDomain = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)
     
     for i = minI:maxI
       for j = minJ:maxJ
         x1  = (i-1)*deltaX   # first node of cell (i,j), x-coord
         y1  = (j-1)*deltaY   # first node of cell (i,j), y-coord
         for ip=1:ppc[1]
           for jp=1:ppc[2]
              x = x1 + 0.5dx + (jp-1)*dx
              y = y1 + 0.5dy + (ip-1)*dy
              if ( (x-xc)^2 + (y-yc)^2 < fRadius^2  )
                thisMaterialPoint = mpmMaterialPoint_2D_Classic()
                thisMaterialPoint.v2Centroid     = [x; y]
                thisMaterialPoint.fVolumeInitial = dx*dy
                thisMaterialPoint.fVolume        = dx*dy
                push!(thisMaterialDomain, thisMaterialPoint)
             end
           end
        end
       end
     end
     
     return(thisMaterialDomain)
    end
   
    function createMaterialDomain_Rectangle(corners::Array{Float64}, 
                                           grid,
                                           ppc::Array{Int64})

     xmin = corners[1,1]
     ymin = corners[1,2]
     xmax = corners[2,1]
     ymax = corners[2,2]

     minI, minJ = point2ElemIndexIJ([xmin;ymin], grid)
     maxI, maxJ = point2ElemIndexIJ([xmax;ymax], grid)

     deltaX     = grid.v2Length_Cell[1]
     deltaY     = grid.v2Length_Cell[2]

     dx         = deltaX/ppc[1]
     dy         = deltaY/ppc[2]

     thisMaterialDomain = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)
     
     for i = minI:maxI
       for j = minJ:maxJ
         x1  = (i-1)*deltaX
         y1  = (j-1)*deltaY
         for ip=1:ppc[1]
           for jp=1:ppc[2]
              x = x1 + 0.5dx + (jp-1)*dx
              y = y1 + 0.5dy + (ip-1)*dy
              if ( x > xmin && x < xmax && y > ymin && y < ymax  )
                thisMaterialPoint = mpmMaterialPoint_2D_Classic()
                thisMaterialPoint.v2Centroid = [x; y]
                thisMaterialPoint.fVolumeInitial = dx*dy
                thisMaterialPoint.fVolume        = dx*dy
                push!(thisMaterialDomain, thisMaterialPoint)
             end
           end
        end
       end
     end
     
     return(thisMaterialDomain)
    end

   function createMaterialDomain_RectangleWithNotch(corners::Array{Float64}, 
                                                    notch::Array{Float64},
                                                    grid,
                                                    ppc::Array{Int64})

     xmin = corners[1,1]
     ymin = corners[1,2]
     xmax = corners[2,1]
     ymax = corners[2,2]
     
     xnmin = notch[1,1]
     ynmin = notch[1,2]
     xnmax = notch[2,1]
     ynmax = notch[2,2]

     minI, minJ = point2ElemIndexIJ([xmin;ymin], grid)
     maxI, maxJ = point2ElemIndexIJ([xmax;ymax], grid)

     deltaX     = grid.v2Length_Cell[1]
     deltaY     = grid.v2Length_Cell[2]

     dx         = deltaX/ppc[1]
     dy         = deltaY/ppc[2]

     thisMaterialDomain = Array{moduleMaterialPoint.mpmMaterialPoint_2D_Classic}(0)
     
     for i = minI:maxI
       for j = minJ:maxJ
         x1  = (i-1)*deltaX
         y1  = (j-1)*deltaY
         for ip=1:ppc[1]
           for jp=1:ppc[2]
              x = x1 + 0.5dx + (jp-1)*dx
              y = y1 + 0.5dy + (ip-1)*dy

              if ( x > xnmin && x < xnmax && y > ynmin && y < ynmax  )
                continue
              end

              if ( x > xmin && x < xmax && y > ymin && y < ymax  )
                thisMaterialPoint = mpmMaterialPoint_2D_Classic()
                thisMaterialPoint.v2Centroid = [x; y]
                thisMaterialPoint.fVolumeInitial = dx*dy
                thisMaterialPoint.fVolume        = dx*dy
                push!(thisMaterialDomain, thisMaterialPoint)
             end
           end
        end
       end
     end
     
     return(thisMaterialDomain)
    end

end
