##This script takes a path with PDB files, an output path and a csv file as inputs, and creates a high quality static image
#of each uniprot in the csv file in the right orientation, coloured by secondary structure
#The CSV file contains a column with uniprot ids, a column with the pdb file that contains the best chain for that uniprot, and the chain ID of that chain within the file
#(this csv file was generated with the script "generate_uniprot_beststructure_list.py")

##Environment variables for the three paths (reference file path, pdb files path and output path) must be set before running this script.

##MODULES IMPORT
import csv
import os

##FUNCTION DECLARATION
def create_picture_chain(uniprot_id,filepath,chain,out_path, count):
###This function takes a PDB filename and a chain identifier and generates previews of that chain in the best possible orientation
###The previews are stored in the out_path provided as input
	
	cmd.load(filepath)
	cmd.hide("everything")
	cmd.show("cartoon", "chain " + chain)
	cmd.orient("chain " + chain)
	cmd.color("red", "ss h")
	cmd.color("yellow", "ss s")
	cmd.color("green", "ss l")
	cmd.ray(1920,1080)
	cmd.png(out_path + uniprot_id, dpi="300")
	count = count + 1
	cmd.delete("all")
	count = count + 1
	return count
####End of function


##SCRIPT SETTINGS
reference_path = os.environ['PERSONAL_DATA'] + "/uniprot2bestfile.csv"
files_path = os.environ['DATA'] + "/pdb_interactome3d_mod/"
out_path = "/data/protein_previews/"
count = 0

##GENERAL PYMOL SETTINGS
cmd.bg_color("white")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_fancy_sheets", 1)
cmd.set("cartoon_highlight_color", "lightorange")
cmd.set("cartoon_discrete_colors", 1)
cmd.set("ray_trace_mode", 0)


#PARSING OF THE CSV INPUT FILE
ref_file = open(reference_path, "r")
reader = csv.reader(ref_file)

for row in reader:
	uniprot_id = row[0]
	filepath = files_path + row[1]
	chain = row[2]
	count = create_picture_chain(uniprot_id,filepath,chain,out_path,count)

