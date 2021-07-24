#2.25D First Scalar Operator Bound

#Import python packages
import numpy as np
import matplotlib.pyplot as plt
# Imports the package
import bootstrap
# The conformal blocks needed for a given run are calculated as a sum over poles.
# Demand that poles with small residues are approximated by poles with large ones.
bootstrap.cutoff = 1e-10

def cprint(message):
    print("\033[94m" + message + "\033[0m")

# Spatial dimension.
dim = 2.25
# Dictates the number of poles to keep and therefore the accuracy of a conformal block.
k_max = 20
# Says that conformal blocks for spin-0 up to and including spin-14 should be computed.
l_max = 14
# Conformal blocks are functions of (a, b) and many derivatives of each should be kept for strong bounds.
# This says to keep derivatives up to fourth order in b.
n_max = 4
# For a given n, this states how many a derivatives should be included beyond 2 * (n - n_max).
m_max = 2


dim_phi= [0.13, 0.15, 0.16, 0.17, 0.18, 0.19, 0.193, 0.196, 0.198, 0.2, 0.202, 0.204, 0.206, 0.208, 0.30, 0.302, 0.305, 0.31, 0.32, 0.33, 0.34, 0.35, 0.4]
#dim_phi= np.linspace(0, 0.4, 41)
length= len(dim_phi)

# Generates the table.
table1 = bootstrap.ConformalBlockTable(dim, k_max, l_max, m_max, n_max)
# Computes the convolution.
table2 = bootstrap.ConvolvedBlockTable(table1)

# We think it is perfectly fine for all internal scalars coupling to our external one to have dimension above 1.7.
lower = 0.25
# Conversely, we think it is a problem for crossing symmetry if they all have dimension above 5.7
upper = 4
# The boundary between these regions will be found within an error of 0.01.
tol = 0.01
# The 0.7 and 1.7 are our guesses for scalars, not some other type of operator.
channel = 0


#store result as pairs
result = [0] * length

#saving to a file
f = open("JNrychkov2p25D.dat", "wb")

for i in range(0, length):
	# Sets up a semidefinite program that we can use to study this.
	sdp = bootstrap.SDP(dim_phi[i], table2)	
	# Calls SDPB to compute the bound.
    	result[i] = sdp.bisect(lower, upper, tol, channel)
	f.write("%.4f \t %.6f \n" % (dim_phi[i], result[i]))
	#update lower and upper
f.close()

for i in range(0, length):
	print dim_phi[i], "\t", result[i]	


#plotting the figure
plt.plot(dim_phi, result, label = "D=2.25")
# naming the x axis
plt.xlabel('d')
# naming the y axis
plt.ylabel('f(d)')
# giving a title to my graph
plt.title('Lowest scalar')
# show a legend on the plot
plt.legend()  
# function to show the plot
plt.show()

