# Script which creates tables and then populates all

## DATABASE SETTINGS: ##

import MySQLdb
import os
import re
import getpass

# Database connection object
# Will prompt for 
db = MySQLdb.connect(host=raw_input('MySQL Host: '),
					user=raw_input('MySQL Username: '),
					passwd=getpass.getpass('MySQL Password: '),
					db=raw_input('MySQL DB: '))



##METHOD DEFINITION:

def create_tables():
	print("Creating tables")

	cur = db.cursor()
	# Get path for table containing file
	sql_path = "database_creation.sql"
	# Open file and execute sql
	with open(sql_path) as f:
		table_creation = f.read()
		cur.execute(table_creation)

	print("Tables created")
	

def uniprot_import():
	print("Populating uniprot tables")
	cur = db.cursor()
	#Script to populate uniprot table
	
	fasta_path = os.environ['SNV_DATA'] + "/FASTA/Fasta"

	filenames = os.listdir(fasta_path)

	# For files in directory
	# Find ones ending in '.txt' i.e. FASTA sequences
	# Extract uniprot id and sequence
	# Insert into uniprot table
	for filename in filenames:
		if filename[-4:] == ".txt":
			uniprot = filename[:-4]
			fasta_file = open(fasta_path+"/"+filename)
			seq = ""
			for line in fasta_file:
				if line[0] != ">":
					seq += line.rstrip()
			try:
				cur.execute('INSERT INTO uniprot VALUES (%s,%s)',(uniprot,seq))
				db.commit()
			except:
				db.rollback()
	cur.close()

	print("Populated uniprot table")

	

def snv_type_import():
	
	print("Populating snv_type table")
	#This script populates the table "snv_type" using data from the "humsavar.txt" file.
	#Environment variables with the path to this file must be changed if this script is run on another machine.

	cur = db.cursor()

	in_file_path = os.environ['SNV_DATA'] + "/humsavar.txt"
	in_file = open(in_file_path, "r")

	types_list = []

	for line in in_file:
		if not line.startswith('#'):
			splitline = line.split(None,6)

			type = splitline[4]

			if type not in types_list:
				types_list.append(type)


	for element in types_list:
		try:
			cur.execute('INSERT INTO snv_type VALUES (%s)', element)
			db.commit()

		except:
			db.rollback()

	cur.close()
	print("Populated snv_type table")


def amino_acid_import():
	# Import amino_acid data into db (from scratch, no other files used for this)

	print("Populating amino_acid table")

	cur = db.cursor()

	data = [
		["A","Ala","Alanine"],
		["I","Ile","Isoleucine"],
		["L","Leu","Leucine"],
		["V","Val","Valine"],
		["F","Phe","Phenylalanine"],
		["W","Trp","Tryptophan"],
		["Y","Tyr","Tyrosine"],
		["N","Asn","Asparagine"],
		["C","Cys","Cysteine"],
		["Q","Gln","Glutamine"],
		["M","Met","Methionine"],
		["S","Ser","Serine"],
		["T","Thr","Threonine"],
		["D","Asp","Aspartic Acid"],
		["E","Glu","Glutamic Acid"],
		["H","His","Histidine"],
		["K","Lys","Lysine"],
		["G","Gly","Glycine"],
		["P","Pro","Proline"],
		["R","Arg","Arginine"],
		["B","Asx","Aspartic Acid or Asparagine"],
		["X","Unk","Any amino acid"],
		["Z","Glx","Glutamic Acid or Glutamine"],
		["U","Sec","Selenocysteine"]
		]

	for element in data:
		try:
			cur.execute('INSERT INTO amino_acid VALUES (%s,%s,%s)',(element[0],element[1],element[2]))
			db.commit()
		except:
			db.rollback()

	cur.close()
	
	print("Populated amino_acid table")

def uniprot_residue_import():
	
	#This populates the "uniprot_residue" table using the data in the "uniprot" table and the "amino_acid" table.
	#No interaction with the data files is required if these tables are already populated.

	print("Populating uniprot_residue table")

	cur = db.cursor()
	
	cur.execute('SELECT acc_number, sequence FROM  uniprot')

	uniprots = cur.fetchall()


	for uniprot in uniprots:
			# Accession number of uniprot
			acc_number = uniprot[0] 
			
			# Sequence of uniprot
			seq = uniprot[1]

			for i in range(len(seq)):
				amino_acid = seq[i] #get amino_acid
				uniprot_position = i + 1 #get uniprot_position

				try:
					#Add the complete row into the table
					cur.execute('INSERT INTO uniprot_residue (uniprot_acc_number,uniprot_position,amino_acid) VALUES (%s,%s,%s)', (acc_number, uniprot_position, amino_acid)) 

					db.commit()

				except:
					db.rollback()

	cur.close()

	print("Populated uniprot_residue table")

		
def disease_import():
	
	##This script populates the "disease" table using the information in the file "humsavar.txt".

	print("Populating disease table")

	cur = db.cursor()

	
	in_data_path = os.environ['SHARED'] + "/snv/data/humsavar.txt" #This passes the information to the script from the corrected copy of humsavar in the shared folder
	in_file = open(in_data_path, "r")

	mim_regex = re.compile('^.*\[MIM:(\d*)')
	##REGEX for the OMIM code: Start_line + any characters + "[MIM:" + (any digits)

	entries = []

	for line in in_file:
		if not line.startswith("#"):
			splitline = line.split(None,6)
			splitline[5] = splitline[5].rstrip()
			type = splitline[4]

			if type == "Disease":

				data = splitline[6]
				
				if  data != "-":
							
					name = data.split('[',1)[0].rstrip() ## split once by "[", take the first element, delete the whitespaces and newlines in the end

					try:
						mim_match = re.match(mim_regex, data) ## apply REGEX to the disease column
						mim = mim_match.group(1)
					except:
						mim = None
					
					entry = [mim,name]

					if entry not in entries:
						entries.append(entry) ## Find out the unique set of diseases/OMIM codes

	for element in entries:
		try:
			cur.execute('INSERT INTO disease VALUES (%s,%s)', element)
			db.commit()

		except:
			db.rollback()



	cur.close()
	in_file.close()

	print("Populated disease table")

		
				
def snv_import():
	
	#This script populates the "snv" table using the information in the "humsavar.txt" file and the information in the "amino_acid" table.
	# Also populates the snv_uniprot_residue table for snvs which map to uniprots we have represented

	#The environment variables to open the file must be changed accordingly if this script is being run on another machine.
	#The tables "snv_type" and "amino_acid" must be populated before this script is run.

	print("Populating snv table")

	cur = db.cursor()
		
	in_file_path = os.environ['SHARED'] + "/snv/data/humsavar.txt"

	in_file = open(in_file_path, 'r')


	wtallele_regex = re.compile('^p\.(\D*)\d*')
	##REGEX for the ancestral allele: Start_line + "p." + (any non-digits) + any digits
	mallele_regex = re.compile('^p\.\D*\d*(\D*)$')
	##REGEX for the mutated allele: Start_line + "p." + any non-digits + any digits + (any non-digits) + end_line

	position_regex = re.compile('^p\.\D*(\d*)')
	##REGEX for the mutated position (Uniprot): Start_line + "p." + any non-digits + (any digits)

	for line in in_file:
		if not line.startswith('#'):
			splitline = line.split(None,6)

			#Extraction of the FT_id:
			ftid = splitline[2]

			#Extraction of the type:
			type = splitline[4]

			#Extraction of the wt_allele:
			wtallele_match = re.match(wtallele_regex, splitline[3])
			wt_allele = wtallele_match.group(1)

			#Extraction of the mutant_allele:
			mallele_match = re.match(mallele_regex, splitline[3])
			m_allele = mallele_match.group(1)

			#Extraction of the uniprot_acc_number:
			uniprot_acc_number = splitline[1]

			#Extraction of the uniprot_position:
			position_match = re.match(position_regex, splitline[3])
			uniprot_position = position_match.group(1)

			#Extraction of the gene_code:
			gene_code_input = splitline[0]

			if gene_code_input == "-":
				gene_code = None
			else:
				gene_code = gene_code_input

			#Extraction of the SNV_id in dbSNP:
			db_SNP_data = splitline[5]

			if db_SNP_data == "-":
				db_SNP = None
			else:
				db_SNP = db_SNP_data.rstrip()

			# Get one letter code from three letter code in file
			# wt_allele
			cur.execute('SELECT one_letter_code FROM amino_acid	WHERE three_letter_code=%s', (wt_allele))
			converted_aa = cur.fetchone()
			wt_allele = converted_aa[0]
			# m_allele
			cur.execute('SELECT one_letter_code FROM amino_acid	WHERE three_letter_code=%s', (m_allele))
			converted_aa = cur.fetchone()
			m_allele = converted_aa[0]

			# Print warning if amino acids not one letter
			if len(wt_allele)>1:
				print "excessive length of wt_aa. Variant code: ", ftid
			elif len(m_allele)>1:
				print "excessive length of m_aa. Variant code: ", ftid
			# Print warning if over length db_snp
			if db_SNP != None:
				if len(db_SNP)>15:
					print "excessive length of db_SNP. Variant code: ", ftid, len(db_SNP)

			#Final compilation as a list
			data_list = [ftid, type, wt_allele, m_allele, uniprot_acc_number, uniprot_position, gene_code, db_SNP]

			try:
				cur.execute('INSERT INTO snv (ft_id,type,wt_aa,mutant_aa,uniprot_acc_number,uniprot_position,gene_code,db_snp) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)', data_list)
				db.commit()

			except MySQLdb.IntegrityError:
				db.rollback()
			# Get uniprot_residue_id
			cur.execute('SELECT id FROM uniprot_residue WHERE uniprot_acc_number=%s AND uniprot_position=%s',(uniprot_acc_number,uniprot_position))
			result = cur.fetchone()
			# If uniprot residue exists insert value into mapping table
			if result != None:
				uniprot_res_id = result[0]
				try:
					cur.execute('INSERT INTO snv_uniprot_residue (ft_id,uniprot_residue_id) VALUES (%s,%s)',(ftid,uniprot_res_id))
					db.commit()
				except MySQLdb.IntegrityError:
					db.rollback()
	cur.close()

	print("Populated snv table")


	

def snv_disease_import():
	
	## This script populates the "snv_disease" table using the information in the "humsavar.txt" file.
	
	print("Populating snv_disease table")

	cur = db.cursor()
	
	in_file_path = os.environ['SHARED'] + "/snv/data/humsavar.txt"
	in_file = open(in_file_path,"r");

	mim_regex = re.compile('^.*\[MIM:(\d*)')
	##REGEX for the OMIM code: Start_line + any word characters + "[MIM:" + (any digits)

	for line in in_file:
		if not line.startswith("#"):
			splitline = line.split(None,6)

			ftid = splitline[2]

			# Check if snv is in database


			mim_data = splitline[6]
			
			if not mim_data == "-":
				mim_match = re.match(mim_regex, mim_data)

				if mim_match != None:
					mim = mim_match.group(1)

					try:
						cur.execute('INSERT INTO snv_disease (ft_id,mim) VALUES(%s,%s)', [ftid,mim])
						db.commit()

					except:
						db.rollback()
							
	cur.close()
	in_file.close()

	print("Populated snv_disease table")
	


def chain_interaction_interaction_type_import():
	# Parse interactome file v3 

	# Version 3 - Not using peewee


	# TODO
	# For each row in file pull both chains and save
	# Know the pk id's of those two 
	# Check they are unique in db

	# File Format
	# 21 or 22 white-space delimited columns
	# DOMAIN2 sometimes blank and therefore 21 rather than 22
	# PROT1	0
	#PROT2	1
	#RANK_MAJOR	2
	#RANK_MINOR	3
	#TYPE	4
	#PDB_ID	5
	#BIO_UNIT	6
	#CHAIN1	7
	#MODEL1	8
	#SEQ_IDENT1	9
	#COVERAGE1	10
	#SEQ_BEGIN1	11
	#SEQ_END1	12
	#DOMAIN1	13
	#CHAIN2	14
	#MODEL2	15
	#SEQ_IDENT2	16
	#COVERAGE2	17
	#SEQ_BEGIN2	18
	#SEQ_END2	19
	#DOMAIN2	20
	#FILENAME	21

	print("Populating chain, interaction and interaction_type tables")

	cur = db.cursor()

	snv_path = os.environ['SNV_DATA']

	flatfile = open(snv_path+'/NON_CANCELLARE/interactome3d_filtered_V1_NUOVO1.txt')

	# Setup interaction_types list
	cur.execute('SELECT type FROM interaction_type')
	inter_types = list(cur.fetchall())

	for line in flatfile:
		# If line not header
		if line[:5] != "PROT1":
			# Split line on whitespace
			columns = line.split()
			# Setup pdb
			pdb_id = columns[5]
			if len(pdb_id) != 4:
				print(pdb_id)
			# Create chain tuples
			# New Order: pdb, pdb_chain, pdb_model, seq_identity, coverage, seq_start, seq_end, uniprot_acc
			chain_1 = (pdb_id,columns[7],columns[8],columns[9],columns[10],columns[11],columns[12],columns[0])
			chain_2 = (pdb_id,columns[14],columns[15],columns[16],columns[17],columns[18],columns[19],columns[1])
			# Insert chains then add id
			pks = []
			try:
				# Chain 1
				try:
					cur.execute('INSERT INTO chain (pdb_id,pdb_chain,pdb_model,seq_identity,coverage,seq_start,seq_end,uniprot_acc_number) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)', chain_1)
					db.commit()
					# Append last modified row id i.e. chain 1
					pks.append(cur.lastrowid)
				except MySQLdb.IntegrityError:
					cur.execute('SELECT id FROM chain WHERE pdb_id=%s AND pdb_chain=%s AND pdb_model=%s AND uniprot_acc_number=%s',(pdb_id,columns[7],columns[8],columns[0]))
					pks.append(cur.fetchone()[0])
				# Chain 2
				try:
					cur.execute('INSERT INTO chain (pdb_id,pdb_chain,pdb_model,seq_identity,coverage,seq_start,seq_end,uniprot_acc_number) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)', chain_2)
					db.commit()
					# Append last modified row id i.e. chain 2
					pks.append(cur.lastrowid)
				except MySQLdb.IntegrityError:
					cur.execute('SELECT id FROM chain WHERE pdb_id=%s AND pdb_chain=%s AND pdb_model=%s  AND uniprot_acc_number=%s',(pdb_id,columns[14],columns[15],columns[1]))
					pks.append(cur.fetchone()[0])
				# Interactions
				inter_type = columns[4]
				if (inter_type,) not in inter_types:
					cur.execute('INSERT INTO interaction_type VALUES (%s)',(inter_type))
					db.commit()
					inter_types.append((inter_type,))
				# New Order: chain_1_id, chain_2_id, type, filename
				interaction = (pks[0],pks[1],inter_type,columns[-1])
				cur.execute('INSERT INTO interaction (chain_1_id,chain_2_id,type,filename) VALUES (%s,%s,%s,%s)', interaction)
				db.commit()
			except:
				db.rollback()
	cur.close()
	print("Populated chain, interactionn and interaction_type tables")
			

def chain_residue_position_mapping_import():
	#Script to populate chain_residue and position_mapping tables
	#(formerly "import_mapping")

	## Parse Mapping Files

	## TODO
	# Populate chain residues using mapping file 
	# Each interaction unique to filename
	# Then use A or B to define whether it is chain 1 or chain 2
	# Then write residue at a time the chain residue and mapping table

	print("Populating chain_residue and position_mapping tables")

	cur = db.cursor()

	snv_path = os.environ['SNV_DATA']

	mapping_path = snv_path+"/alignFASTA/"

	files = os.listdir(mapping_path)

	for filename in files:
		if filename[-4:] == ".out":
			# Partner_id is either A or B depending on if it is chain 1 or chain 2
			partner_id = filename[13]
			# Name of pdb file extracted from clustal files
			pdb_filename = filename[16:-12]
			# Get chain_id of the chain being mapped
			try:
				if partner_id == "A":
					cur.execute('SELECT chain_1_id FROM interaction WHERE filename=%s',(pdb_filename))
					chain_id = cur.fetchone()[0]
				elif partner_id == "B":
					cur.execute('SELECT chain_2_id FROM interaction WHERE filename=%s',(pdb_filename))
					chain_id = cur.fetchone()[0]
			except TypeError:
				print(pdb_filename)
				raise
			# Check uniprot being mapped matches 
			file_uniprot = filename[:6]
			cur.execute('SELECT uniprot_acc_number FROM chain WHERE id=%s',(chain_id))
			chain_uniprot = cur.fetchone()[0]
			if file_uniprot != chain_uniprot:
				print("NON MATCHING:")
				print(file_uniprot)
				print(chain_uniprot)
				print(filename)
			# Open mapping file and iterate over each line
			mapping_file = open(mapping_path+filename)
			for line in mapping_file:
				columns = line.split()
				# If line is data line not header
				if len(columns) == 5 or len(columns) == 6:
					pdb_residue = columns[1]
					uni_residue = columns[2]
					pdb_res_position = columns[-2]
					uni_res_position = columns[-1]

					if pdb_residue != "-":
						try:
							cur.execute('INSERT INTO chain_residue (chain_id, chain_position, amino_acid) VALUES (%s,%s,%s)', (chain_id, pdb_res_position, pdb_residue))
							db.commit()
							# Get id of inserted chain residue
							chain_res_id = cur.lastrowid
							if uni_residue != "-":
								# Get id of uniprot residue to be mapped
								cur.execute('SELECT id FROM uniprot_residue WHERE uniprot_acc_number=%s AND uniprot_position=%s',(chain_uniprot,uni_res_position))
								uniprot_res_id = cur.fetchone()[0]
								# Insert position mapping of residue ids
								cur.execute('INSERT INTO position_mapping (uniprot_residue_id,chain_residue_id) VALUES (%s,%s)',(uniprot_res_id,chain_res_id))
								db.commit()
						except MySQLdb.IntegrityError:
							db.rollback()
			mapping_file.close()

	cur.close()
	print("Populated chain_residue and position_mapping tables")
	

def interface_residue_import():
	#Script to populate interface_residue table

	##This script populates the interface_residue table, using the information in the pdbinter unique files
	print("Populating interface_residue table")

	cur = db.cursor()

	pdbinter_path = os.environ['SNV_DATA'] + "/pdbinter_results/"

	files = os.listdir(pdbinter_path)

	for filename in files:

		if filename[-4:] == ".uni":
		#Check for file extension

			pdb_filename = filename[:-4]

			cur.execute('SELECT id,chain_1_id,chain_2_id FROM interaction WHERE filename=%s', (pdb_filename))
			interaction = cur.fetchone()
			#Get information from the DB: chain IDs and interaction ID


			interaction_id = interaction[0]
			# Get interaction id
			
			in_file = open(pdbinter_path + filename,"r")

			for line in in_file:
				splitline = line.split()

				chain_position = splitline[0]
				partner_id = splitline[2]
				#Parse lines to get the "partner ID" (the A/B parameter in the pdbinter file, from which we will get the chain id by querying the DB) and the chain_position of each residue

				#This is the if construction to get the chain_ID from the "parther_ID" that we parsed from the file
				if partner_id == "A":
					chain_id = interaction[1]
				
				elif partner_id == "B":
					chain_id = interaction[2]
				
				# Get chain_residue_id for this position

				cur.execute('SELECT id FROM chain_residue WHERE chain_id=%s AND chain_position=%s',(chain_id,chain_position))
				chain_residue = cur.fetchone()
				# Some chains have no residues due to malformed pdbs preventing alignments
				# Discard interfaces as no way of mapping
				if chain_residue != None:
					chain_residue_id = chain_residue[0]
					try:
						cur.execute('INSERT INTO interface_residue (chain_residue_id, interaction_id) VALUES(%s,%s)', (chain_residue_id,interaction_id))
						db.commit()
					except MySQLdb.IntegrityError:
						db.rollback()

			in_file.close()
	cur.close()

	print("Populated interface_residue table")

def accessibility_import():
	#Script to populate accesibility table
	#(formerly "import accessibility")

	## Accessibility Parser

	## TODO ##
	# Use part of filename to identify the interaction
	# Iterate over each row of file
	# Extract bound, unbound and disulphide bridge from row

	## FILE FORMAT ##
	# Split on whitespace
	## Name - List index ##
	# Data - 0
	# Chain - 1
	# Residue no - 2
	# Amino acid - 3
	# Area changed between bound and unbound - 4[0]
	# In disulphide bridge - 4[1]
	# Disulphide bridge no - 5
	# Tot area res (Unbound) - 6
	# Tot area main chain (Unbound) - 7
	# Tot area side chain (Unbound) - 8
	# Tot RELA area residue (Unbound) - 9
	# Tot RELA area main chain (Unbound) - 10
	# Tot RELA area side chain (Unbound) - 11
	# Tot area res (Bound) - 12
	# Tot area main chain (Bound) - 13
	# Tot area side chain (Bound) - 14
	# Tot RELA area residue (Bound) - 15
	# Tot RELA area main chain (Bound) - 16
	# Tot RELA area side chain (Bound) - 17

	# Columns to extract: 9,15,5

	## FILE FORMAT ##
	# Use splices 
	# DATA [0:4]
	# Chain [5]
	# Residue No [6:10]
	# Amino Acid [11:15]
	# Disulphide Bridge [19]
	# RELA Area Residue (Unbound) [33:37]
	# RELA Area Residue (Bound) [57:61]

	def is_numeric(x):
		try:
			float(x)
			return True
		except:
			return False

	print("Populating accessibility table")

	cur = db.cursor()

	snv_path = os.environ['SNV_DATA']

	access_path = snv_path+"/rarfiles/"

	files = os.listdir(access_path)

	# Build Amino Acid Dictionary
	cur.execute('SELECT three_letter_code,one_letter_code FROM amino_acid')
	amino_acids = cur.fetchall()
	three2one = {}
	for aa in amino_acids:
		three2one[aa[0].upper()] = aa[1]

	for filename in files:
		if filename[-4:] == ".rar":
			# Partner_id is either A or B depending on if it is chain 1 or chain 2
			partner_id = filename[13]
			# Name of pdb file extracted from rar filenames
			pdb_filename = filename[16:-4]
			# Get interaction id
			cur.execute('SELECT id,chain_1_id,chain_2_id FROM interaction WHERE filename=%s', (pdb_filename))
			interaction = cur.fetchone()
			interaction_id = interaction[0]
			# Set chain 
			if partner_id == "A":
				chain_id = interaction[1]
			elif partner_id == "B":
				chain_id = interaction[2]
			# Uniprot id from filename
			file_uniprot = filename[:6]
			cur.execute('SELECT uniprot_acc_number FROM chain WHERE id=%s',(chain_id))
			chain_uniprot = cur.fetchone()[0]
			if file_uniprot != chain_uniprot:
				print("NON MATCHING:")
				print(file_uniprot)
				print(chain_uniprot)
				print(filename)
			# Open Accessibility file
			access_file = open(access_path+filename)
			for line in access_file:
				# Check it is a DATA line
				if line[:4] == "DATA":
					# Extract variables from line
					amino_acid = three2one[line[11:15].strip().upper()]
					pdb_position = line[6:11].strip()
					apo_rel_area = line[33:37].strip()
					complex_rel_area = line[57:61].strip()
					disulphide_bridge = line[19]
					# Get chain_residue_id
					cur.execute('SELECT id FROM chain_residue WHERE chain_id=%s AND chain_position=%s',(chain_id,pdb_position))
					chain_residue = cur.fetchone()
					# Some chains have no residues due to malformed pdbs preventing alignments
					# Discard accessibilities if this is the case
					if chain_residue != None:
						chain_residue_id = chain_residue[0]
						# Attempt to insert or rollback
						try:
							cur.execute('''INSERT INTO accessibility 
								(interaction_id, 
								chain_residue_id,
								bound_acc,
								unbound_acc,
								disulphide_bridge_no)
								VALUES (%s,%s,%s,%s,%s)''',
								(interaction_id,chain_residue_id,apo_rel_area,complex_rel_area,disulphide_bridge))
							db.commit()
						except MySQLdb.IntegrityError:
							db.rollback()
	cur.close()

	print("Populated accessibility table")
	

##METHOD CALLS:

create_tables()

uniprot_import()

snv_type_import()

amino_acid_import()

uniprot_residue_import()

disease_import()

snv_import()

snv_disease_import()

chain_interaction_interaction_type_import()

chain_residue_position_mapping_import()

interface_residue_import()

accessibility_import()


db.close()
