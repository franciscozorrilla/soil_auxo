import glob
import cobra
import reframed
from reframed import load_cbmodel
from reframed.cobra import auxotrophy
from reframed import Environment

for file in glob.iglob(r'data/models/*.xml'):

	# Use try to handle exceptions with reading SMBL model
	try:

		#load model
		model = load_cbmodel(file)

		#load env media
		env = Environment.from_compounds(["fru","h2o","met__L","cys__L","arg__L","ala__L","asn__L","asp__L","glu__L","gln__L","gly","his__L","ile__L","leu__L","lys__L","phe__L","pro__L","ser__L","thr__L","trp__L","tyr__L","val__L","ribflv","cbl1","thm","pydxn","btn","ade","ura","csn","thym","k","pi","h2","na1","cl","nh4","mg2","so4","ca2","fe","fe2","fe3","cobalt2","cu2","mn2","mobd","ni2","zn2","o2"])

		# Open text file with .aux extension
		text_file = open(file+".aux", "w")

		# Use model.medium COBRA function to extract active fluxes from model
		n = text_file.write(str(auxotrophy.auxotrophies(model,min_rel_growth=0.975,min_abs_growth=0.1)))

		# Close file
		text_file.close()

	# Catch cobra.io.sbml.CobraSBMLError and continue loop    
	except cobra.io.sbml.CobraSBMLError:

		# Print message with model ID for log
		print("The model",file,"raised an SBML error and could not be read in by COBRA")

		#Nothing to see here fellas
		pass
