# cd("Julia/MPM")
# cd("E:\\MyPublications\\MPM_Julia\\Codes\\SpeedTest")

# simple implementation to test the processing speed of iterative for loops and vectors

import PyPlot

function myMain()
	iIterations = 1000
	# Data Structure-------------------------------------------------------------
	# ---------------------------------------------------------------------------
	@printf("Creating arrays...")
	tic()
	thisArray01 = Array{Float64}(10,1000000)
	thisArray02 = Array{Float64}(10,1000000)
	thisArray03 = Array{Float64}(10,1000000)
	@printf("...done ")
	toc()
	@printf("\n")

	@printf("Assigning arrays...")
	tic()
	for i = 1:1:size(thisArray01,1)
		for j = 1:1:size(thisArray01,2)
			thisArray01[i,j] = 1.0*i
			thisArray02[i,j] = -1.0*i
			thisArray03[i,j] = 0
		end
	end
	@printf("...done ")
	toc()
	@printf("\n")

	@printf("Manipulating arrays through loop ...")
	tic()
	ni = size(thisArray01,1)
	nj = size(thisArray01,2)
	for count in 1:1:iIterations
		for i = 1:1:ni
		for j = 1:1:nj
				thisArray03[i,j] = thisArray01[i,j] + thisArray02[i,j]
			end
		end
	end
	@printf("...done ")
	toc()
	@printf("\n")
	# println("Array03: ", thisArray03)

	@printf("Manipulating arrays through vectorization ...")
	tic()
	for count in 1:1:iIterations
		thisArray03 = thisArray01 + thisArray02
	end
	@printf("...done ")
	toc()
	@printf("\n")

	# println("Array03: ", thisArray03)
end

myMain()
