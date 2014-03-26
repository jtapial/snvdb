# Script which creates tables and then populates all

## DATABASE SETTINGS: ##

import MySQLdb
import os
import re
import getpass
from decimal import *
import urllib
import urllib2
from Bio import Entrez
import fnmatch
from itertools import izip
# Database connection object
# Will prompt for 
db = MySQLdb.connect(host=raw_input('MySQL Host: '),
user=raw_input('MySQL Username: '),
passwd=getpass.getpass('MySQL Password: '),
db=raw_input('MySQL DB: '))

## HELPER METHODS ##

def is_numeric(x):
		try:
			float(x)
			return True
		except:
			return False

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
	

def uniprot_extract_from_fasta():
	#This performs the first stage of the uniprot table population (from the data in the provided fasta files). The second stage is done along with the snv table population
		
	print("Extracting uniprots from FASTA files")
	cur = db.cursor()
		
	fasta_path = os.environ['SNV_DATA'] + "/FASTA/Fasta"

	filenames = os.listdir(fasta_path)
	seq_dict = {}
	
	# For files in directory
	# Find ones ending in '.txt' i.e. FASTA sequences
	# Extract uniprot id and sequence and make a dictionary with them
	
	for filename in filenames:
		if filename[-4:] == ".txt":
			uniprot = filename[:-4]
			fasta_file = open(fasta_path+"/"+filename)
			seq = ""
			for line in fasta_file:
				if line[0] != ">":
					seq += line.rstrip()
			
			seq_dict[uniprot] = seq
			
					
	uniprots_list = seq_dict.keys()

	print "Extracted:", len(uniprots_list), "uniprots from FASTA files"
	print "Extracted:", len(seq_dict), "uniprot sequences from FASTA files"
	
	return(uniprots_list, seq_dict)
	

def uniprot_extract_from_humsavar():
#This method takes information from the humsavar file only. No interaction with the database.
#Code extraction from the humsavar file
	print "Extracting uniprots from humsavar.txt"
	in_file_path = os.environ['SHARED'] + '/snv/data/humsavar.txt'
	in_file = open(in_file_path, 'r')
	in_file_lines = in_file.readlines()

	initial_codes_dict = {}
	extracted_uniprots_list = []
	no_code_list = []

	for line in in_file_lines:
		if not line.startswith('#'):
			splitline = line.split(None,6)
			code = splitline[0]
			uniprot_id = splitline[1]
			if code == "-":
				code = None
				no_code_list.append(code)

			extracted_uniprots_list.append(uniprot_id)
			if uniprot_id not in initial_codes_dict:
				initial_codes_dict[uniprot_id] = code
			

	print "Extraction from humsavar.txt complete."
	print "Total number of uniprots from humsavar: ", len(extracted_uniprots_list)
	print "Extracted: ", len(initial_codes_dict), "gene codes, of which", len(no_code_list), "were found to be NULL"

	in_file.close()

	return (extracted_uniprots_list, initial_codes_dict)
			
def uniprot_check (uniprots_list, seq_dict, initial_codes_dict):
	#Batch request to uniprot. The dataset is divided in chunks of 200 to avoid server overload:
	chunks = [uniprots_list[i:i+200] for i in range(0,len(uniprots_list),200)]
		
	uniprot_output_list = [None]*len(chunks)
	count = 0
	for element in chunks:
		uniprot_string = ' '.join(element)
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
		uniprot_output_list[count] = response.read()
		count += 1
	#Join the results in a single string
	uniprot_output = ''.join(uniprot_output_list)	
	
	#Split the file by protein ("//" Characters at the end of a line)
	splitpage = uniprot_output.split('//\n')
	print "Total number of uniprot entries found: ", len(splitpage)

	#Creation of an entries dictionary (the keys are the uniprot accession numbers on the entries (maybe many per entry), and the values are the text entries for each protein)
	entries_dict = {}
	for entry in splitpage:
		extracted_uniprot_candidates = []
	
		for line in entry.split('\n'):
			if line.startswith("AC"):
				for item in line.split()[1:]:
					extracted_uniprot_candidates.append(item.rstrip(";"))
		
		for item in extracted_uniprot_candidates:
			entries_dict[item] = splitpage.index(entry)
	
	print "Length of uniprot entries dict: ", len(entries_dict)
	
	#Check if there are any uniprots missing and generate a list:
	no_uniprot_found = [x for x in uniprots_list if x not in entries_dict]
	print len(no_uniprot_found), "entries were not found in the batch retrieval"
	if no_uniprot_found:
		print no_uniprot_found
	
	#Update uniprots_list removing the non existent uniprots
	uniprots_list = [x for x in uniprots_list if x not in no_uniprot_found]
	
	
	#PARSING OF EACH UNIPROT ENTRY:
	#Initialise general variables for the parsing (for the final error report)
	no_code_list = []
	no_name_list = []
	
	#Initialise output dicts
	codes_dict = {}
	synonyms_dict = {}
	names_dict = {}
	
	#Regex compilation
	gene_code_regex = re.compile('Name=(.*?);')
	gene_synonym_regex = re.compile('Synonyms=(.*?);')
	uniprot_name_regex = re.compile('RecName:.*?Full=(.*?);')

	for uniprot in uniprots_list:
			
		#Initialise variables for each parsing
		
		extracted_code_candidates = []
		extracted_code = None
		extracted_uniprot_name = None

		#Look for the entries dictionary value for each uniprot
		try:
			entry = splitpage[entries_dict[uniprot]]
	
		#If no uniprot found in the dictionary, add this case to the final error report
		except:
			print "ERROR FOR UNIPROT: ", uniprot, ". NO UNIPROT ENTRY FOUND. CHECK MANUALLY"
			continue
						
		entrylines = entry.split('\n')
		
		for line in entrylines:

			#Look for gene codes and append them to the candidate list
			if line.startswith('GN'):
				if gene_code_regex.search(line):
					for item in gene_code_regex.findall(line):
						extracted_code_candidates.append(item)
					
				#Include synonym codes too	
				synonyms_match = gene_synonym_regex.search(line)
				
				if synonyms_match:
					synonyms = synonyms_match.group(1)
					synonyms_list = [x.strip() for x in synonyms.split(',')]
					for item in synonyms_list:
						extracted_code_candidates.append(item)
				
				continue
			
			#Look for uniprot full names and take the first one (we do not have any previous information about this)
			elif line.startswith("DE") and uniprot_name_regex.search(line):
				extracted_uniprot_name = uniprot_name_regex.search(line).group(1)
				continue
			
		#Pull the sequence if no previous sequence had been found (that is, the uniprot comes from the humsavar file and not from FASTA)
		
			elif line.startswith("SQ"):
				index = entrylines.index(line)
				rawseq = entrylines[(index+1):-1]
				stripseq_1 = [x.strip() for x in rawseq]
				seq_draft = ''.join(stripseq_1)
				stripseq_2 = seq_draft.split()
				seq = ''.join(stripseq_2)					
							
					
		#Data analysis of the results from the just-parsed entry:
		
		#Pick a gene code (preferably the most similar to the code provided, if the uniprot comes from the humsavar file).
		#If not possible, add this case to the final error report
		
		extracted_code = None
		
		if extracted_code_candidates:
			if uniprot in initial_codes_dict:
				for code in extracted_code_candidates:
					if code.startswith(initial_codes_dict[uniprot]):
						extracted_code = code
						extracted_code_candidates.remove(code)
						break
				
				if extracted_code == None:
					extracted_code = extracted_code_candidates.pop(0)
			
			else:
				extracted_code = extracted_code_candidates.pop(0)
		
		else:
			extracted_code = None
			no_code_list.append(uniprot)
					
		#Append the extracted code to the output
		codes_dict[uniprot] = extracted_code
				
		#Take the rest of the codes and add them to the synonyms dict:
		synonyms_dict[uniprot] = [x for x in extracted_code_candidates if x != extracted_code]
		
		#Add the name to the names dict:
		names_dict[uniprot] = extracted_uniprot_name
		if not extracted_uniprot_name:
			no_name_list.append(uniprot)
		
		#Add the sequence to the seq dict:
		if uniprot not in seq_dict:
			seq_dict[uniprot] = seq
		
	#Pull a list of uniprots without seq and uniprots without synonyms:
	no_seq_list = [x for x in uniprots_list if x not in seq_dict]
	no_synonyms_list = [x for x in uniprots_list if synonyms_dict[x] == []]

				
	#Final message
	print "Retrieval from uniprot complete."

	##Final reports:
	print "Retrieved: ", len(codes_dict), "gene codes, of which ", len(no_code_list), "were found to be NULL"
	print no_code_list
	
	print "Retrieved synonyms for:", len(synonyms_dict), "uniprots, of which", len(no_synonyms_list), "were found to have no synonyms"

	print "Retrieved: ", len(names_dict), "protein names, of which ", len(no_name_list), "were found to be NULL"
	print no_name_list
	
	print "Retrieved: ", len(seq_dict), "protein sequences. No sequence found for: ", len(no_seq_list), "proteins"	
		
	#Return the information
	return (uniprots_list, names_dict, seq_dict, codes_dict, synonyms_dict)

	

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
		["A","Ala","Alanine","Aliphatic"],
		["I","Ile","Isoleucine","Aliphatic"],
		["L","Leu","Leucine","Aliphatic"],
		["V","Val","Valine","Aliphatic"],
		["F","Phe","Phenylalanine","Aromatic"],
		["W","Trp","Tryptophan","Aromatic"],
		["Y","Tyr","Tyrosine","Aromatic"],
		["N","Asn","Asparagine","Acidic/Amidic"],
		["C","Cys","Cysteine","Sulphur-containing"],
		["Q","Gln","Glutamine","Acidic/Amidic"],
		["M","Met","Methionine","Sulphur-containing"],
		["S","Ser","Serine","Hydroxylic"],
		["T","Thr","Threonine","Hydroxylic"],
		["D","Asp","Aspartic Acid","Acidic/Amidic"],
		["E","Glu","Glutamic Acid","Acidic/Amidic"],
		["H","His","Histidine","Basic"],
		["K","Lys","Lysine","Basic"],
		["G","Gly","Glycine","Aliphatic"],
		["P","Pro","Proline","Aliphatic"],
		["R","Arg","Arginine","Basic"],
		["B","Asx","Aspartic Acid or Asparagine","Acidic/Amidic"],
		["X","Unk","Any amino acid",'None'],
		["Z","Glx","Glutamic Acid or Glutamine","Acidic/Amidic"],
		["U","Sec","Selenocysteine","Selenium-containing"]
		]

	for element in data:
		try:
			cur.execute('INSERT INTO amino_acid_group (name) VALUES (%s)',element[3])
			db.commit()
			group_id = cur.lastrowid
		except MySQLdb.IntegrityError:
			cur.execute('SELECT id FROM amino_acid_group WHERE name=%s',element[3])
			group_id = cur.fetchone()[0]
		try:
			cur.execute('INSERT INTO amino_acid (one_letter_code,three_letter_code,name,amino_acid_group_id) VALUES (%s,%s,%s,%s)',(element[0],element[1],element[2],group_id))
			db.commit()
		except MySQLdb.IntegrityError:
			db.rollback()

	cur.close()
	
	print("Populated amino_acid table")


def id_retrieval(uniprots_list, codes_dict, synonyms_dict):
#Method to retrieve the genbank ids using a list of uniprots, a uniprot:code dict, and an uniprot:[synonyms] dict as input
#Returns an uniprot:GenbankID dict

	print "Retrieving GenBank ID from the local download for :", len(uniprots_list), "uniprots"

	genepath = os.environ['GENE_DATA']

	data_file = open(genepath+'/gene2accession', 'r')
	
	remaining_uniprots = list(uniprots_list)
		
	id_dict = {}
	
	
	#First, try to find ids using the original code. The dictionary will maintain only the keys for which no id was found.
	found_uniprots = []	
	dataline = data_file.readline()
	while dataline:
		splitdataline = [x.rstrip() for x in dataline.split('\t')]
		if splitdataline[0] == "9606":
			for uniprot in remaining_uniprots:
				if splitdataline[-1] == codes_dict[uniprot] or splitdataline[-1] in synonyms_dict[uniprot]:
					id_dict[uniprot] = splitdataline[1]
					found_uniprots.append(uniprot)
					remaining_uniprots.remove(uniprot)
					break
					
		dataline = data_file.readline()
	

	print "First round completed using gene codes and synonyms. Found: ", len(id_dict)
	print "Found in this round: ", len(found_uniprots)
	print "Left: ", len(remaining_uniprots)
	print remaining_uniprots
	
	#Finally, go to the gene_refseq_uniprot_collab file and look up uniprots whose gene_codes nor synonyms were found on the gene2accession file
	gene_uniprot_file = open(genepath+'/gene_refseq_uniprotkb_collab', 'r')
	mapping_dict = {}
	no_mapseq_uniprots = remaining_uniprots
	
	gene_uniprot_line = gene_uniprot_file.readline()
	while gene_uniprot_line:
		split_gene_uniprot_line = [x.rstrip() for x in gene_uniprot_line.split('\t')]
		if split_gene_uniprot_line[-1] in no_mapseq_uniprots:
			found_uniprots.append(split_gene_uniprot_line[-1])
			mapping_dict[split_gene_uniprot_line[-1]] = split_gene_uniprot_line[0]
			no_mapseq_uniprots.remove(split_gene_uniprot_line[-1])
		
		gene_uniprot_line = gene_uniprot_file.readline()
		
	print "Remaining uniprots: ", len(remaining_uniprots)
	print "Found mapseq for:", len(mapping_dict), "uniprots"
	print "Could not find mapseqs for:", len(no_mapseq_uniprots)
			
	found_uniprots = []
	remaining_uniprots_with_mapseq = mapping_dict.keys()
	
	data_file.seek(0)
	dataline = data_file.readline()
		
	while dataline:
		splitdataline = [x.rstrip() for x in dataline.split('\t')]
		if splitdataline[0] == "9606":
			for uniprot in remaining_uniprots_with_mapseq:
				if any(item.startswith(mapping_dict[uniprot]) for item in splitdataline): 
					id_dict[uniprot] = splitdataline[1]
					found_uniprots.append(uniprot)
					remaining_uniprots_with_mapseq.remove(uniprot)
					break
		dataline = data_file.readline()
		
	for item in remaining_uniprots:
		if item in found_uniprots:
			remaining_uniprots.remove(item)
			
		
	print "Second round completed using the uniprot-refseq mapping file. Found: ", len(id_dict)
	print "Left: ", len(remaining_uniprots)
	print remaining_uniprots
		
	
	data_file.close()
	gene_uniprot_file.close()
	
	##Third round: query Gene online with the remaining genes
	Entrez.email = "javier.tapial-rodriguez13@imperial.ac.uk"
	found_uniprots = []
	for uniprot in remaining_uniprots:
		if codes_dict[uniprot] != None:
			handle = Entrez.esearch(db="gene", term="%s AND homo sapiens[Organism]" % codes_dict[uniprot])
			record = Entrez.read(handle)
						
			if record["IdList"]:
				first_id = record["IdList"][0]
				id_dict[uniprot] = first_id
				found_uniprots.append(uniprot)
								
			else:
				print "Nothing found online for the uniprot", uniprot
	
	for item in remaining_uniprots:
		if item in found_uniprots:
			remaining_uniprots.remove(item)

	#This is the only entry for which the algorithm cannot find an ID while it is known that there is one.
	#The result is then hard-coded

	id_dict["Q6GMX4"] = "4505"
	found_uniprots.append("Q6GMX4")
	remaining_uniprots.remove("Q6GMX4")

	#Final report		
	print "Third round completed querying Entrez online. Found: ", len(id_dict)
	print "Left: ", len(remaining_uniprots)
	print remaining_uniprots	
	
	print "GenBank ID retrieval complete. Found ID for: ", len(id_dict), "entries"
	print "No ID found for: ", len(remaining_uniprots), "entries"
	
	if remaining_uniprots:
		print "Report of unfound entries (uniprot/code/synonyms):"
		for item in remaining_uniprots:
			print item, codes_dict[item], synonyms_dict[item]

	print "The reported uniprots without a genbank ID will have it set to NULL"
	for uniprot in remaining_uniprots:
			if uniprot not in id_dict:
				id_dict[uniprot] = None
	
	
	return id_dict

def uniprot_import():
	print "Populating uniprot table"

	#Call the method to extract uniprots and seqs from FASTA files
	fasta_uniprots, fasta_seq_dict = uniprot_extract_from_fasta()

	#Call the method to extract uniprots and code candidates from humsavar.txt
	humsavar_uniprots, humsavar_codes_dict = uniprot_extract_from_humsavar()

	#Merge the list of uniprots
	all_uniprots = list(set(fasta_uniprots + humsavar_uniprots))

	print "Total number of uniprots: ", len(all_uniprots)

	#Pass the information to the uniprot_check method to check that the uniprots exist, and get the rest of the sequences, final codes, synonyms and protein names
	uniprots_list, names_dict, seq_dict, codes_dict, synonyms_dict = uniprot_check(all_uniprots, fasta_seq_dict, humsavar_codes_dict)


	#Pass the information from the uniprot_check method to the id_retrieval method to retrieve GenBank IDs
	id_dict = id_retrieval(uniprots_list, codes_dict, synonyms_dict)

	cur = db.cursor()

	for bit in uniprots_list:
		a = names_dict.get(bit, None)
		b = seq_dict.get(bit, None)
		c = codes_dict.get(bit, None)
		d = id_dict.get(bit, None)
		
		try:
			cur.execute("INSERT INTO uniprot VALUES(%s,%s,%s,%s,%s)", (bit,a,b,c,d))
			db.commit()
		except:
			db.rollback()
			print "DATABASE ERROR FOR UNIPROT: ", bit
			print bit, a, b, c, d
	
	print "population of Uniprot table complete"

def snv_import():

#This script populates the "snv" table using the information in the "humsavar.txt" file.

#The environment variables to open the file must be changed accordingly if this script is being run on another machine.
#The tables "snv_type" and "amino_acid" must be populated before this script is run.
#The uniprot table must have been populated


	print("Populating snv table.")

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

			#Check if the uniprot exists in the database
			cur.execute('SELECT EXISTS(SELECT 1 FROM uniprot WHERE acc_number=%s)', (uniprot_acc_number))
			check_database = cur.fetchone()[0]

			#If it exists, add the snv normally calling the method "add_this_snv"
			if check_database == 0:
				print "Skipped SNV, FTID: ", ftid, ", Uniprot accession number: ", uniprot_acc_number, " Could not detect a valid uniprot entry"
				continue
				

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

		
			# Get uniprot_residue_id
			cur.execute('SELECT id FROM uniprot_residue WHERE uniprot_acc_number=%s AND uniprot_position=%s',(uniprot_acc_number,uniprot_position))
			try:
				uniprot_residue_id = cur.fetchone()[0]
			except TypeError:
				print("This uniprot does not exist")
				print(uniprot_acc_number,uniprot_position)
				uniprot_residue_id=None

			#Final compilation as a list

			data_list = [ftid, snv_type, wt_allele, m_allele, uniprot_acc_number, uniprot_residue_id, db_SNP]

			#If it exists, add the snv normally calling the method "add_this_snv"
			if check_database == 1:
				
				add_this_snv(data_list)
	
	cur.close()
	in_file.close()


def add_this_snv(data_list):
	#This method adds a snv to the local database. It requires a uniprot entry with the adequate accession number alredy in the database
	# It also requires the amino_acid table to have been populated
	# The information from the snvs comes from the humsavar.txt file opened by the snv_import method


	#Check that the line on the humsavar.txt file is still open and stored in splitline

	cur = db.cursor()

	#Addition to the database
	try:
		cur.execute('INSERT INTO snv (ft_id,type,wt_aa,mutant_aa,uniprot_acc_number,uniprot_residue_id,db_snp) VALUES (%s,%s,%s,%s,%s,%s,%s)', data_list)
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

	cur.close()

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

			inter_type = columns[4]
			if (inter_type,) not in inter_types:
				cur.execute('INSERT INTO interaction_type VALUES (%s)',(inter_type))
				db.commit()
				inter_types.append((inter_type,))
			biological_unit =  int(columns[6])
			# Create chain tuples
			# New Order: pdb, pdb_chain, pdb_model, seq_identity, coverage, seq_start, seq_end, uniprot_acc, type, biological_unit
			chain_1 = (pdb_id,columns[7],columns[8],columns[9],columns[10],columns[11],columns[12],columns[0],inter_type,biological_unit)
			chain_2 = (pdb_id,columns[14],columns[15],columns[16],columns[17],columns[18],columns[19],columns[1],inter_type,biological_unit)
			# Insert chains then add id
			pks = []
			# Chain 1
			try:
				cur.execute('INSERT INTO chain (pdb_id,pdb_chain,pdb_model,seq_identity,coverage,seq_start,seq_end,uniprot_acc_number,type,biological_unit) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', chain_1)
				db.commit()
				# Append last modified row id i.e. chain 1
				pks.append(cur.lastrowid)
			except MySQLdb.IntegrityError:
				cur.execute('SELECT id FROM chain WHERE pdb_id=%s AND pdb_chain=%s AND pdb_model=%s AND uniprot_acc_number=%s AND type=%s AND biological_unit=%s',(pdb_id,columns[7],columns[8],columns[0],inter_type,biological_unit))
				pks.append(cur.fetchone()[0])
			# Chain 2
			try:
				cur.execute('INSERT INTO chain (pdb_id,pdb_chain,pdb_model,seq_identity,coverage,seq_start,seq_end,uniprot_acc_number,type,biological_unit) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', chain_2)
				db.commit()
				# Append last modified row id i.e. chain 2
				pks.append(cur.lastrowid)
			except MySQLdb.IntegrityError:
				cur.execute('SELECT id FROM chain WHERE pdb_id=%s AND pdb_chain=%s AND pdb_model=%s  AND uniprot_acc_number=%s AND type=%s AND biological_unit=%s',(pdb_id,columns[14],columns[15],columns[1],inter_type,biological_unit))
				pks.append(cur.fetchone()[0])
			# Interactions
			
			# New Order: chain_1_id, chain_2_id, type, filename
			interaction = (pks[0],pks[1],inter_type,columns[-1])
			cur.execute('INSERT INTO interaction (chain_1_id,chain_2_id,type,filename) VALUES (%s,%s,%s,%s)', interaction)
			db.commit()
	cur.close()
	print("Populated chain, interaction and interaction_type tables")
			

def chain_residue_position_mapping_import():
	#Script to populate chain_residue and position_mapping tables
	#(formerly "import_mapping")

	# Populate chain_residue, position_mapping and position_transform tables
	# Position transform populated with adjustment values for chains with different labeled but identically mapped sequences.
	#

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
					cur.execute('SELECT chain_1_id,id FROM interaction WHERE filename=%s',(pdb_filename))
					results = cur.fetchone()
					chain_id = results[0]
				elif partner_id == "B":
					cur.execute('SELECT chain_2_id,id FROM interaction WHERE filename=%s',(pdb_filename))
					results = cur.fetchone()
					chain_id = results[0]
				interaction_id = results[1]
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
			
			# Check whether there are already chain residues for this chain
			if cur.execute('SELECT id FROM chain_residue WHERE chain_id=%s',(chain_id)) == 0:
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
				cur.execute('INSERT INTO position_transform (interaction_id,chain_id,value) VALUES (%s,%s,0)',(interaction_id,chain_id))
				db.commit()

			else:
				# If there are chain residues
				# Check that sequences match
				# Get sequence from input file
				# Attempt to open input file

				# Get a filename for which the transform value is 0
				# Therefore is the reference numbering system
				cur.execute('SELECT interaction_id FROM position_transform WHERE value=0 AND chain_id=%s',(chain_id))
				ref_transform = cur.fetchone()
				ref_interaction_id = ref_transform[0]
				# Find filename and chain ids for that interaction
				cur.execute('SELECT filename,chain_1_id,chain_2_id FROM interaction WHERE id=%s',(ref_interaction_id))
				ref_interaction = cur.fetchone()
				# Set A or B dependent on part of interaction
				ref_base_filename = ref_interaction[0]
				if ref_interaction[1] == chain_id:
					ref_partner = "A"
				elif ref_interaction[2] == chain_id:
					ref_partner = "B"
				else:
					print("No Partner")
					print(ref_base_filename)
					print(ref_interaction)
				# Get ref full file
				ref_match = fnmatch.filter(files,"*"+ref_partner+"__"+ref_base_filename+"*.out")
				if len(ref_match) == 1:
					ref_file = open(mapping_path+ref_match[0])
				else:
					print('More than one filename match')
					print(ref_match)
					print(ref_base_filename)
					print(ref_partner)
				trans_file = open(mapping_path+filename)

				# Check whether every line in file is identical
				diff_count = 0
				for ref,trans in izip(ref_file,trans_file):
					if ref != trans:
						diff_count += 1
				# If files have no differences insert transform of 0
				if diff_count == 0:
					cur.execute('INSERT INTO position_transform (chain_id,interaction_id,value) VALUES (%s,%s,0)',(chain_id,interaction_id))
					db.commit()
				else:
					# Setup lists to hold sequences and positions
					ref_sequence = []
					ref_positions = []
					trans_sequence = []
					trans_positions = []

					ref_file2 = open(mapping_path+ref_match[0])
					trans_file2 = open(mapping_path+filename)

					# Get all amino acid 1 letter codes
					#cur.execute('SELECT one_letter_code FROM amino_acid')
					#one_letter_codes = cur.fetchall()

					for ref,trans in izip(ref_file2,trans_file2):

						# Get positions
						ref_position = ref[19:23].strip()
						trans_position = trans[19:23].strip()
						# Get amino acids
						ref_aa = ref[7]
						trans_aa = trans[7]							
						if ref_aa != "-" and ref_aa != "l":
							if is_numeric(ref_position):
								ref_positions.append(int(ref_position))
							else:
								ref_positions.append(ref_position)
								print("Non-numeric position")
								print(ref_position)
							ref_sequence.append(ref_aa)
						if trans_aa != '-' and ref_aa != 'l':
							if is_numeric(trans_position):
								trans_positions.append(int(trans_position))
							else:
								trans_positions.append(trans_position)
								print("Non-numeric position")
								print(trans_position)
							trans_sequence.append(trans_aa)
					if ref_sequence == trans_sequence:
						if ref_positions != trans_positions:
							value = ref_positions[0] - trans_positions[0]
							amended_positions = []
							for p in trans_positions:
								amended_positions.append(p+value)
							if amended_positions == ref_positions:
								cur.execute('INSERT INTO position_transform (chain_id,interaction_id,value) VALUES (%s,%s,%s)',(chain_id,interaction_id,value))
								db.commit()
							else:
								print('Positions do not match')
								print(filename)
								print(ref_positions)
								print(trans_positions)
								print(amended_positions)
						else:
							print('Something else is different between files')
							print(ref_filename)
							print(filename)
					else:
						print("Non-matching sequences")
						print(ref_sequence)
						print(trans_sequence)
							
				


				

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
			#Get information from the DB: chain IDs and interaction ID
			cur.execute('SELECT id,chain_1_id,chain_2_id FROM interaction WHERE filename=%s', (pdb_filename))
			interaction = cur.fetchone()
			
			# Get interaction id
			interaction_id = interaction[0]

			# Get transform values and map them to chain ids
			cur.execute('SELECT chain_id,value FROM position_transform WHERE interaction_id=%s',(interaction_id))
			chain2value = {}
			for record in cur.fetchall():
				chain2value[int(record[0])] = int(record[1])
			
			in_file = open(pdbinter_path + filename,"r")

			for line in in_file:
				splitline = line.split()
				#Parse lines to get the "partner ID" (the A/B parameter in the pdbinter file, from which we will get the chain id by querying the DB) and the chain_position of each residue

				chain_position = splitline[0]
				partner_id = splitline[2]



				#This is the if construction to get the chain_ID from the "parther_ID" that we parsed from the file
				if partner_id == "A":
					chain_id = interaction[1]
				
				elif partner_id == "B":
					chain_id = interaction[2]
				

				# Adjust chain_position using transform table
				try:
					if chain2value[int(chain_id)] != 0:
						transformed_position = int(chain_position) + chain2value[int(chain_id)]
					else:
						transformed_position = chain_position
					# Get chain_residue_id for this position

					cur.execute('SELECT id FROM chain_residue WHERE chain_id=%s AND chain_position=%s',(chain_id,transformed_position))
					chain_residue = cur.fetchone()
				# KeyError shows there are no chain residues
				except KeyError:
					chain_residue = None
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


			# Check chain has any residues
			residues = cur.execute('SELECT id FROM chain_residue WHERE chain_id=%s',(chain_id))
			if residues != 0:
				# Get transform value
				cur.execute('SELECT value FROM position_transform WHERE chain_id=%s AND interaction_id=%s',(chain_id,interaction_id))
				transform_value = int(cur.fetchone()[0])

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
						# Transform position value
						if transform_value != 0:
							transformed_position = int(pdb_position) + transform_value
						else:
							transformed_position = pdb_position
						# Get chain_residue_id
						results_length = cur.execute('SELECT id FROM chain_residue WHERE chain_id=%s AND chain_position=%s',(chain_id,transformed_position))
						chain_residue = cur.fetchone()
						# Some chains have no residues due to malformed pdbs preventing alignments
						# Discard accessibilities if this is the case
						if results_length != 0:
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
						else:
							print('Unable to find matching residue')
							print(chain_id)
							print(pdb_position)
							print(transformed_position)
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

#create_tables()

#uniprot_import()

#snv_type_import()

amino_acid_import()

#uniprot_residue_import()

#snv_import()

#disease_import()

#snv_disease_import()

#chain_interaction_interaction_type_import()

#chain_residue_position_mapping_import()

#interface_residue_import()

#accessibility_import()

#pfam_import()

db.close()
