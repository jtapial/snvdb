CREATE TABLE uniprot (
acc_number CHAR(6) NOT NULL,
name VARCHAR(250),
sequence TEXT,
gene_code CHAR(15),
genbank_id CHAR(15),
PRIMARY KEY (acc_number)
);
CREATE TABLE amino_acid_group (
id INT NOT NULL AUTO_INCREMENT,
name VARCHAR(255),
PRIMARY KEY (id),
UNIQUE (name)
);
CREATE TABLE amino_acid (
one_letter_code CHAR(1) NOT NULL,
three_letter_code CHAR(3) NOT NULL,
name VARCHAR(255) NOT NULL,
amino_acid_group_id INT NOT NULL,
PRIMARY KEY (one_letter_code),
FOREIGN KEY(amino_acid_group_id) REFERENCES amino_acid_group (id)
);
CREATE TABLE uniprot_residue (
id INT NOT NULL AUTO_INCREMENT,
uniprot_acc_number CHAR(6) NOT NULL,
uniprot_position INT NOT NULL,
amino_acid CHAR(1) NOT NULL,
PRIMARY KEY (id),
UNIQUE (uniprot_acc_number,uniprot_position),
FOREIGN KEY (uniprot_acc_number) REFERENCES uniprot (acc_number),
FOREIGN KEY (amino_acid) REFERENCES amino_acid (one_letter_code)
);
CREATE TABLE chain (
id INT NOT NULL AUTO_INCREMENT,
pdb_id CHAR(4),
pdb_chain VARCHAR(10) NOT NULL,
pdb_model INT NOT NULL,
seq_identity VARCHAR(10) NOT NULL,
coverage VARCHAR(10) NOT NULL,
seq_start VARCHAR(10) NOT NULL,
seq_end VARCHAR(10) NOT NULL,
type VARCHAR(20) NOT NULL,
biological_unit INT NOT NULL,
uniprot_acc_number CHAR(6) NOT NULL,
PRIMARY KEY (id),
FOREIGN KEY (uniprot_acc_number) REFERENCES uniprot (acc_number),
FOREIGN KEY (type) REFERENCES interaction_type (type),
UNIQUE (pdb_id,pdb_chain,uniprot_acc_number,pdb_model,type,biological_unit)
);
CREATE TABLE interaction_type (
type VARCHAR(20) NOT NULL,
PRIMARY KEY (type)
);
CREATE TABLE interaction (
id INT NOT NULL AUTO_INCREMENT,
chain_1_id INT NOT NULL,
chain_2_id INT NOT NULL,
type VARCHAR(20) NOT NULL,
filename VARCHAR(100) NOT NULL,
crossing_number INT,
PRIMARY KEY (id),
FOREIGN KEY (chain_1_id) REFERENCES chain (id),
FOREIGN KEY (chain_2_id) REFERENCES chain (id),
FOREIGN KEY (type) REFERENCES interaction_type (type),
UNIQUE (chain_1_id,chain_2_id)
);
CREATE TABLE chain_residue (
id INT NOT NULL AUTO_INCREMENT,
chain_id INT NOT NULL,
chain_position VARCHAR(10) NOT NULL,
amino_acid CHAR(1) NOT NULL,
PRIMARY KEY (id),
UNIQUE (chain_id, chain_position),
FOREIGN KEY (amino_acid) REFERENCES amino_acid (one_letter_code)
);
CREATE TABLE interface_residue (
id INT NOT NULL AUTO_INCREMENT,
chain_residue_id INT NOT NULL,
interaction_id INT NOT NULL,
res_order INT,
PRIMARY KEY (id),
UNIQUE (chain_residue_id, interaction_id),
FOREIGN KEY (chain_residue_id) REFERENCES chain_residue (id),
FOREIGN KEY (interaction_id) REFERENCES interaction (id)
);
CREATE TABLE accessibility (
id INT NOT NULL AUTO_INCREMENT,
interaction_id INT NOT NULL,
chain_residue_id INT NOT NULL,
bound_acc VARCHAR(10) NOT NULL,
unbound_acc VARCHAR(10) NOT NULL,
disulphide_bridge_no INT NOT NULL,
PRIMARY KEY (id),
UNIQUE (interaction_id,chain_residue_id),
FOREIGN KEY (interaction_id) REFERENCES interaction (id),
FOREIGN KEY (chain_residue_id) REFERENCES chain_residue (id)
);
CREATE TABLE position_mapping (
id INT NOT NULL AUTO_INCREMENT,
uniprot_residue_id INT NOT NULL,
chain_residue_id INT NOT NULL,
PRIMARY KEY (id),
UNIQUE (uniprot_residue_id,chain_residue_id),
FOREIGN KEY (uniprot_residue_id) REFERENCES uniprot_residue (id),
FOREIGN KEY (chain_residue_id) REFERENCES chain_residue (id)
);
CREATE TABLE snv_type (
type VARCHAR(20) NOT NULL,
PRIMARY KEY (type)
);
CREATE TABLE disease (
mim INT(6) NOT NULL,
name VARCHAR(255),
PRIMARY KEY (mim)
);
CREATE TABLE snv (
ft_id VARCHAR(10) NOT NULL,
type VARCHAR(20) NOT NULL,
wt_aa CHAR(1) NOT NULL,
mutant_aa CHAR(1) NOT NULL,
uniprot_acc_number CHAR(6) NOT NULL,
uniprot_residue_id INT NOT NULL,
db_snp VARCHAR(15),
PRIMARY KEY (ft_id),
FOREIGN KEY (type) REFERENCES snv_type (type),
FOREIGN KEY (wt_aa) REFERENCES amino_acid (one_letter_code),
FOREIGN KEY (mutant_aa) REFERENCES amino_acid (one_letter_code),
FOREIGN KEY (uniprot_acc_number) REFERENCES uniprot (acc_number),
FOREIGN KEY (uniprot_residue_id) REFERENCES uniprot_residue (id)
);
CREATE TABLE snv_disease (
id INT NOT NULL AUTO_INCREMENT,
ft_id VARCHAR(10) NOT NULL,
mim INT(6) NOT NULL,
PRIMARY KEY (id),
UNIQUE (ft_id, mim),
FOREIGN KEY (ft_id) REFERENCES snv (ft_id),
FOREIGN KEY (mim) REFERENCES disease (mim)
);
CREATE TABLE pfam_hmm (
hmm_acc VARCHAR(12) NOT NULL,
name VARCHAR(50) NOT NULL,
type VARCHAR(50) NOT NULL,
length INT NOT NULL,
clan VARCHAR(50) NOT NULL,
PRIMARY KEY (hmm_acc)
);
CREATE TABLE uniprot_pfam_mapping (
id INT AUTO_INCREMENT NOT NULL,
uniprot_acc_number VARCHAR(6) NOT NULL,
hmm_acc VARCHAR(12) NOT NULL,
alignment_start_ur_id INT NOT NULL,
alignment_end_ur_id INT NOT NULL,
envelope_start_ur_id INT NOT NULL,
envelope_end_ur_id INT NOT NULL,
hmm_start INT NOT NULL,
hmm_end INT NOT NULL,
bit_score DECIMAL(5,1) NOT NULL,
e_value VARCHAR(255) NOT NULL,
significance BOOL NOT NULL,
PRIMARY KEY (id),
UNIQUE(uniprot_acc_number,hmm_acc,alignment_start_ur_id,alignment_end_ur_id),
FOREIGN KEY (uniprot_acc_number) REFERENCES uniprot (acc_number),
FOREIGN KEY (hmm_acc) REFERENCES pfam_hmm (hmm_acc),
FOREIGN KEY (alignment_start_ur_id) REFERENCES uniprot_residue (id),
FOREIGN KEY (alignment_end_ur_id) REFERENCES uniprot_residue (id),
FOREIGN KEY (envelope_start_ur_id) REFERENCES uniprot_residue (id),
FOREIGN KEY (envelope_end_ur_id) REFERENCES uniprot_residue (id)
);
CREATE TABLE active_site_residue (
id INT AUTO_INCREMENT NOT NULL,
up_mapping_id INT NOT NULL,
active_site_ur_id INT NOT NULL,
PRIMARY KEY (id),
FOREIGN KEY (up_mapping_id) REFERENCES uniprot_pfam_mapping (id),
FOREIGN KEY (active_site_ur_id) REFERENCES uniprot_residue (id),
UNIQUE(up_mapping_id,active_site_ur_id)
);
CREATE TABLE combined_uniprot_mapping (
id INT AUTO_INCREMENT NOT NULL,
ref_chain_id INT NOT NULL,
target_chain_id INT NOT NULL,
target_chain_letter CHAR(1) NOT NULL,
target_interaction_id INT NOT NULL,
PRIMARY KEY (id),
FOREIGN KEY (ref_chain_id) REFERENCES chain (id),
FOREIGN KEY (target_chain_id) REFERENCES chain (id),
FOREIGN KEY (target_interaction_id) REFERENCES interaction (id),
UNIQUE (ref_chain_id,target_chain_id,target_chain_letter)
);
CREATE TABLE position_transform (
id INT AUTO_INCREMENT NOT NULL,
interaction_id INT NOT NULL,
chain_id INT NOT NULL,
value INT NOT NULL,
PRIMARY KEY (id),
FOREIGN KEY (interaction_id) REFERENCES interaction (id),
FOREIGN KEY (chain_id) REFERENCES chain (id),
UNIQUE (interaction_id,chain_id)
);
CREATE TABLE interface_atom (
id INT AUTO_INCREMENT NOT NULL,
interface_residue_id INT NOT NULL,
label CHAR(4) NOT NULL,
PRIMARY KEY (id),
FOREIGN KEY (interface_residue_id) REFERENCES interface_residue (id),
UNIQUE (interface_residue_id,label)
);
CREATE TABLE interface_atom_interaction (
id INT AUTO_INCREMENT NOT NULL,
interface_atom_id_l INT NOT NULL,
interface_atom_id_r INT NOT NULL,
type CHAR(3) NOT NULL,
length DECIMAL(4,3) NOT NULL,
PRIMARY KEY (id),
FOREIGN KEY (interface_atom_id_l) REFERENCES interface_atom (id),
FOREIGN KEY (interface_atom_id_r) REFERENCES interface_atom (id),
UNIQUE (interface_atom_id_l,interface_atom_id_r)
);
