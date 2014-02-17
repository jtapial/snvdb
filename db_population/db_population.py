# Script which creates tables and then populates all

## DATABASE SETTINGS: ##

import MySQLdb
import os
import re
import getpass
from decimal import *

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
	

def uniprot_import_from_fasta():
	#This populates the uniprot table in two stages from the data in the provided FASTA files
		
	print("Populating uniprot table from FASTA files")
	cur = db.cursor()
		
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
				print "Error adding uniprot from FASTA file. Accession number: ", uniprot
				db.rollback()
				continue

	cur.close()
	print("Added Uniprot entries from FASTA files")

	

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


def uniprot_gene_code_check(tuples):
	#This method takes a list of tuples (uniprot ID, gene code), then queries Uniprot with the batch of IDs and corrects the gene codes.
	#The output is a list of corrected gene codes
	
	print "Initiating retrieval from Uniprot"
	code_list = [x[1] for x in tuples]
	uniprot_ids = [x[0] for x in tuples]
	
	uniprot_string = ' '.join(uniprot_ids)
	
	#Generate url POST request for Uniprot batch retrieval

	url = "http://www.uniprot.org/batch/"

	params = {
		'format': 'txt',
		'query': uniprot_string
	}

	data = urllib.urlencode(params)
	request = urllib2.Request(url, data)
	contact = "javier.tapial-rodriguez13@imperial.ac.uk"
	request.add_header('User-Agent', 'Python %s' % contact)

	#Open and read the request (text flatfile)
	
	response = urllib2.urlopen(request)
	
	retrieval_file = response.read()

	#Split the file by protein ("//"" Characters at the end of a line)

	splitpage = retrieval_file.split('//\n')

	#Regex compilation
	uniprot_name_regex = re.compile('RecName:.*?Full=(.*?);')
	gene_code_regex = re.compile('Name=(.*?);')

	#Initialise variables for the dictionary creation (will use this for the final error report):
	no_uniprot_found = []

	#Creation of a dictionary (the keys are the uniprot accession numbers on the entries (maybe many per entry), and the values are the text entries for each protein
	
	uniprot_dict = {}
	
	for protein in splitpage:
		extracted_uniprot_candidates = []
	
		for line in protein.split('\n'):
			if line.startswith("AC"):
				for item in line.split()[1:]:
					extracted_uniprot_candidates.append(item.rstrip(";"))
		
		for item in extracted_uniprot_candidates:
			uniprot_dict[item] = splitpage.index(protein)
	
	#Parsing of each uniprot entry

	#Initialise general variables for the parsing (for the final error report)
		
		code_mismatch_list = []
		code_mismatch_candidates_list = []
		no_uniprot_name = []

		#Initialise output list
		correct_codes = []

	for element in tuples:
		u_id = element[0]
		code = element[1]

		#Initialise variables for each parsing
		
		extracted_code_candidates = []
		extracted_uniprot_name = None
		extracted_code = None

		#Look for the dictionary value for each uniprot
		try:
			protein = splitpage[uniprot_dict[u_id]]
		
		#If no uniprot found in the dictionary, add this case to the final error report
		except:
			no_uniprot_found.append(u_id)
			continue
						
		for line in protein.split('\n'):

			#Look for gene codes and append them to the candidate list
			if line.startswith('GN') and gene_code_regex.search(line):
				for item in gene_code_regex.findall(line):
					extracted_code_candidates.append(item) 
			
			#Look for uniprot full names and take the first one (we do not have any previous information about this for the edge cases)
			elif line.startswith("DE") and uniprot_name_regex.search(line):
				extracted_uniprot_name = uniprot_name_regex.search(line).group(1)
			

			#Try to fit the previous gene code to any of the candidates. If not possible, then pick the first one add this case to the final error report
			for element in extracted_code_candidates:
				if element.startswith(code):
					extracted_code = element
					break
			
		if extracted_code == None:
			extracted_code = extracted_code_candidates[0]
			code_mismatch_list.append(code)
			code_mismatch_candidates_list.append(extracted_code_candidates)
			
		#Append the extracted code to the output
		correct_codes.append(extracted_code)


		#If no uniprot full name found, then add the uniprot accession number to the final error report
		if extracted_uniprot_name == None:
			no_uniprot_name.append(u_id)
			
	#Final message
	print "Retrieval from uniprot complete. Fixed codes for: ", len(correct_codes), "entries"

	##Final report of failures:

	print "No uniprot found for: ", len(no_uniprot_found), "entries"
	print no_uniprot_found
	
	print "Found: ", len(code_mismatch_list), "gene code mismatches"
	for element1, element2 in zip(code_mismatch_list, code_mismatch_candidates_list):
		print element1
		print element2

	print "No uniprot name found for: ", len(no_uniprot_name), "uniprots"
	for element in no_uniprot_name:
		print element
	
	#Return of results
	return correct_codes
	

def idretrieval(codes):
#Method to retrieve the ids using a list of gene codes as input
	print "Retrieving GenBank ID from the local download for :", len(codes), "gene codes"

	data_file = open('/data/geneDB/ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession', 'r')
	id_dict = {}
		
	dataline = data_file.readline()
	while dataline:
		splitdataline = [x.rstrip() for x in dataline.split('\t')]
		if splitdataline[0] == "9606":
			for code in codes:
				if splitdataline[-1] == code:
					id_dict[code] = splitdataline[1]
					codes.remove(code)
					break
					
		
		dataline = data_file.readline()

	data_file.close()
	
	print "GenBank ID retrieval complete. Found ID for: ", len(id_dict.items()), "entries"
	return id_dict.items()
	


def gene_import():
	#Method to populate the gene table
	#This method takes information from the humsavar file only. No interaction with the database.
	#Code extraction from the humsavar file

	print "Populating gene table"
	
	in_file_path = os.environ['SHARED_DATA'] + '/humsavar.txt'
	in_file = open(in_file_path, 'r')
	in_file_lines = in_file.readlines()

	short_codes_raw = []
	long_codes_tuples = []

	for line in in_file_lines:
		if not line.startswith('#'):
			splitline = line.split(None,6)
			code = splitline[0]
			if len(code) < 9 :
				short_codes_raw.append(code)
			else:
				uniprot_id = splitline[1]
				long_codes_tuples.append((uniprot_id,code))
			

	short_codes = list(set(short_codes_raw))

	long_codes_unique_tuples = list(set(long_codes_tuples))

	print "Gene code extraction complete."
	print "Total number of codes: ", len(short_codes + long_codes_unique_tuples)
	print len(short_codes), "codes extracted normally"
	print len(long_codes_unique_tuples), "codes are likely to be truncated and will be checked online with Uniprot"

	#Method calling:
	correct_long_codes = uniprot_gene_code_check(long_codes_unique_tuples)

	final_codes = short_codes + correct_long_codes

	output = idretrieval(final_codes)


	#Write results into the database:


	for entry in output:
		try:
			cur.execute('INSERT INTO gene VALUES (%s,%s)', entry)
			db.commit()

		except:
			print "Database error for entry: ", entry
			db.rollback()
			continue

	print "Populated gene table"
	in_file.close()


def snv_import():

#This script populates the "snv" table using the information in the "humsavar.txt" file.
#It also calls the methods:
# - uniprot_import_from_uniprot
# - check_uniprot_online
# - add_this_snv

#The environment variables to open the file must be changed accordingly if this script is being run on another machine.
#The tables "snv_type" and "amino_acid" must be populated before this script is run.
#The uniprot table must have been populated with the import_uniprot_from_fasta method


	print("Populating snv table. Uniprots not present in the FASTA files will also be added")

	cur = db.cursor()
		
	in_file_path = os.environ['SHARED'] + "/snv/data/humsavar.txt"

	in_file = open(in_file_path, 'r')


	wtallele_regex = re.compile('^p\.(\D*)\d*')
	##REGEX for the wt allele: Start_line + "p." + (any non-digits) + any digits
	mallele_regex = re.compile('^p\.\D*\d*(\D*)$')
	##REGEX for the mutant allele: Start_line + "p." + any non-digits + any digits + (any non-digits) + end_line

	position_regex = re.compile('^p\.\D*(\d*)')
	##REGEX for the mutant position (Uniprot): Start_line + "p." + any non-digits + (any digits)

	for line in in_file:
		if not line.startswith('#'):
			splitline = line.split(None,6)

			#Extraction of the FT_id:
			ftid = splitline[2]

			#Extraction of the uniprot_acc_number:
			uniprot_acc_number = splitline[1]

			#Extraction of the uniprot_position:
			position_match = re.match(position_regex, splitline[3])
			uniprot_position = position_match.group(1)

			#Extraction of the type:
			snv_type = splitline[4]

			#Extraction of the wt_allele:
			wtallele_match = re.match(wtallele_regex, splitline[3])
			wt_allele = wtallele_match.group(1)

			#Extraction of the mutant_allele:
			mallele_match = re.match(mallele_regex, splitline[3])
			m_allele = mallele_match.group(1)

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
			data_list = [ftid, snv_type, wt_allele, m_allele, uniprot_acc_number, uniprot_position, gene_code, db_SNP]
			

			#Check if the uniprot exists in the database
			cur.execute('SELECT EXISTS(SELECT 1 FROM uniprot WHERE acc_number=%s)', (uniprot_acc_number))
			check_database = cur.fetchone()[0]

			#If it exists, add the snv normally calling the method "add_this_snv"
			if check_database == 1:
				add_this_snv(data_list)

			#If it does not exists, check the snv with uniprot online:
			elif check_database == 0:

				#If the uniprot checker detects a sequence on Uniprot online and adds it to the local database, extract more data for the snv and add the snv
				if check_uniprot_online(uniprot_acc_number) == True:

					#Addition of the snv (call method)
					add_this_snv(data_list)
				
				#If the uniprot checker does not detect a valid sequence, print a message and do not add the snv (nor the uniprot)
				else:
					print "Skipped SNV, FTID: ", ftid, ", Uniprot accession number: ", uniprot_acc_number, " Could not detect a valid uniprot entry"
					continue
	
	cur.close()
	db.close()
	in_file.close()

def check_uniprot_online(uniprot_acc_number):

	#This method checks online if there is a valid sequence for a given uniprot accession number.
	#If there is, it adds the Uniprot entry to the local database and returns true
	#If there is not, it returns false
	
	cur = db.cursor()

	webpage = urllib.urlopen("http://www.uniprot.org/uniprot/" + uniprot_acc_number + ".fasta")
	lines = webpage.readlines()

	for i in range(len(lines)):
		lines[i] = lines[i].rstrip()
	
	seq_lines = lines[1:]
	seq_lines = ''.join(seq_lines)

	#Check if the sequence is empty.
	if seq_lines != "":
		try:
			cur.execute('INSERT INTO uniprot VALUES (%s, %s)', (uniprot_acc_number,seq_lines))
			return True
		except:
			print seq_lines
			db.rollback()
			print "ERROR adding Uniprot from the online database. Uniprot accession number: ", uniprot_acc_number
			return False
	else:
		return False


def add_this_snv(data_list):
	#This method adds a snv to the local database. It requires a uniprot entry with the adequate accession number alredy in the database
	# It also requires the amino_acid table to have been populated
	# The information from the snvs comes from the humsavar.txt file opened by the snv_import method


	#Check that the line on the humsavar.txt file is still open and stored in splitline

	cur = db.cursor()

	#Addition to the database
	try:
		cur.execute('INSERT INTO snv (ft_id,type,wt_aa,mutant_aa,uniprot_acc_number,uniprot_position,gene_code,db_snp) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)', data_list)
		db.commit()

	except MySQLdb.IntegrityError:
		#Check if the SNV is already in the database (THE FILE HAS DUPLICATE ENTRIES FOR A SNV IF HAS MORE THAN ONE DISEASE)
		cur.execute('SELECT EXISTS(SELECT 1 FROM snv WHERE ft_id=%s)', data_list[0])
		check_already_added = cur.fetchone()[0]
		if check_already_added == 1:
			db.rollback()
			pass
		
		else:
			db.rollback()
			print "Integrity error adding SNV FTID: ", data_list
			pass


def uniprot_residue_import():
	
	#This populates the "uniprot_residue" table using the data in the "uniprot" table and the "amino_acid" table.
	#No interaction with the data files is required if these tables are already populated.

	print("Populating uniprot_residue table")

	cur = db.cursor()
	
	cur.execute('SELECT acc_number, sequence FROM uniprot')

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
					print "Error importing uniprot residue. Uniprot accession number: ", acc_number, " Position: ", uniprot_position
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

def pfam_import():

	print("Populating pfam_hmm, uniprot_pfam_mapping and active_site_residue tables")

	cur = db.cursor()

	pfam_file = open(os.environ['SHARED']+'/snv/data/Pfam_results.txt')

	for line in pfam_file:
		# Split on whitespace
		values = line.split()
		# Check not header line
		if values[0][0] != '#':
			# hmm values - to be inserted into pfam_hmm
			hmm_acc = values[5]
			hmm_name = values[6]
			hmm_type = values[7]
			hmm_length = int(values[10])
			hmm_clan = values[12]
			# Insert into pfam_hmm
			try:
				cur.execute('INSERT INTO pfam_hmm (hmm_acc,name,type,length,clan) VALUES (%s,%s,%s,%s,%s)',(hmm_acc,hmm_name,hmm_type,hmm_length,hmm_clan))
				db.commit()
			except MySQLdb.IntegrityError:
				pass
			# mapping values - insert into uniprot_pfam_mapping
			uniprot_acc_number = values[0][3:9]
			ur_ids = []
			# get all uniprot_residue ids for alignment and envelope starts and ends
			for i in range(1,5):
				cur.execute('SELECT id FROM uniprot_residue WHERE uniprot_acc_number=%s AND uniprot_position=%s',(uniprot_acc_number,values[i]))
				ur_ids.append(cur.fetchone())
			alignment_start_ur_id = ur_ids[0][0]
			alignment_end_ur_id = ur_ids[1][0]
			envelope_start_ur_id = ur_ids[2][0]
			envelope_end_ur_id = ur_ids[3][0]
			hmm_start = int(values[8])
			hmm_end = int(values[9])
			bit_score = Decimal(values[11])
			e_value = values[12]
			significance = int(values[13])
			mapping_insert = (uniprot_acc_number,hmm_acc,alignment_start_ur_id,alignment_end_ur_id,envelope_start_ur_id,envelope_end_ur_id,hmm_start,hmm_end,bit_score,e_value,significance)
			mapping_id = ""
			try:
				cur.execute('''INSERT INTO uniprot_pfam_mapping 
					(uniprot_acc_number,
					hmm_acc,
					alignment_start_ur_id,
					alignment_end_ur_id,
					envelope_start_ur_id,
					envelope_end_ur_id,
					hmm_start,
					hmm_end,
					bit_score,
					e_value,
					significance)
					VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''',mapping_insert)
				mapping_id = int(cur.lastrowid)
			except MySQLdb.IntegrityError:
				cur.execute('SELECT id FROM uniprot_pfam_mapping WHERE uniprot_acc_number=%s AND hmm_acc=%s AND alignment_start_ur_id=%s AND alignment_end_ur_id=%s',(
					uniprot_acc_number,
					hmm_acc,
					alignment_start_ur_id,
					alignment_end_ur_id
					))
				mapping_id = cur.fetchone()[0]
			
			# active site insert
			if len(values) == 16:
				active_site = values[15]
				match = re.search('\[(.*?)\]',active_site)
				if match:
					active_site_residues = match.group(1)
				for position in set(active_site_residues.split(',')):
					try:
						cur.execute('SELECT id FROM uniprot_residue WHERE uniprot_acc_number=%s AND uniprot_position=%s',(uniprot_acc_number,position))
						active_site_ur_id = cur.fetchone()[0]
						cur.execute('INSERT INTO active_site_residue (up_mapping_id,active_site_ur_id) VALUES (%s,%s)',(mapping_id,active_site_ur_id))
						db.commit()
					except MySQLdb.IntegrityError:
						pass
	print("Populated pfam_hmm, uniprot_pfam_mapping and active_site_residue tables")
	cur.close()
	
##METHOD CALLS:

create_tables()

uniprot_import_from_fasta()

snv_type_import()

amino_acid_import()

gene_import()

snv_import()

uniprot_residue_import()

disease_import()

snv_disease_import()

chain_interaction_interaction_type_import()

chain_residue_position_mapping_import()

interface_residue_import()

accessibility_import()

pfam_import()

db.close()
