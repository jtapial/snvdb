# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#     * Rearrange models' order
#     * Make sure each model has one field with primary_key=True
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin.py sqlcustom [appname]'
# into your database.
from __future__ import unicode_literals

from django.db import models

from django.core.urlresolvers import reverse

from Bio import Entrez

##############################
#     GLOBAL METHOD
##############################

#This method gets a list of SNVs and return appropriate color code to display
def get_color(snv_list):
	if len(snv_list)==1:
		if snv_list[0].type.type =='Disease':
			return  'red'
		elif snv_list[0].type.type =='Polymorphism':	
			return 'green'
		elif snv_list[0].type.type =='Unclassified':	
			return 'blue'
	return 'purple'			


##############################
#     CLASSES
############################## 
class Uniprot(models.Model):
	acc_number = models.CharField(max_length=6L, primary_key=True)
	sequence = models.TextField(blank=True)
	class Meta:
		db_table = 'uniprot'

	def get_Snv(self):	#Return a list of all related snvs and their statistics
		snvlist=[]
		for item in self.residues.all():					
			snvlist.extend(item.snvs.all())		
		di = po = un = 0		
		for item in snvlist:
			if item.type.type == 'Disease':
				di+=1
			elif item.type.type == 'Polymorphism':
				po+=1
			else:
				un+=1

		return [snvlist,[len(snvlist),di,po,un]]

	def get_seq(self):
		mod_seq = ''
		for curr_pos in range(len(self.sequence)):
			if((curr_pos)%10 ==0):
				mod_seq = mod_seq + ' '			
			mod_seq = mod_seq+self.sequence[curr_pos]
		return mod_seq

	###############################################
	#This method returns html code for a sequence with mapped snvs 	
	def get_mapping_seq(self):
		seq=self.sequence
		mod_seq=''		
		snvlist=self.get_Snv()
		mod = {}
		for snv in snvlist:
			if snv.uniprot_position in mod:
				mod[snv.uniprot_position].append(snv)
			else:
				mod[snv.uniprot_position]=[snv]
		
		if len(mod)==0:
			return '<p>No SNV Found</p>'		
		sorted_mod = sorted(mod.items(), key=lambda t: t[0]) #list of snv [pos0,pos1] pos0 = residue number, pos1 = snv obj 
		pos = sorted_mod.pop(0)		
		for curr_pos in range(len(seq)):
			if((curr_pos)%10 ==0):
				mod_seq = mod_seq + ' '		
			if(curr_pos<pos[0]-1):
				mod_seq = mod_seq + seq[curr_pos]
			else: 
				pos_mute = ''			
				color = get_color(pos[1])
				if color == 'purple':#multiple mutation case
					html = ''
					for item in pos[1]:					
						pos_mute=pos_mute + item.mutant_aa.one_letter_code +'('+item.type.type+') '	
						html = html+"<p><a href = '"+item.get_absolute_url()+"' target='_blank'>ID: "+item.ft_id+", Mutation: "+item.mutant_aa.one_letter_code +" ("+item.type.type+") </a></p>"
					mod_seq = mod_seq + '<mark style="background-color:'+color+'"><a href="#" style="color:white;" class="open-AddBookDialog" data-id="'+html+'" data-toggle="modal" data-target=".bs-modal-sm"  title="Position '+str(curr_pos+1)+', Mutation: ' +pos_mute+'">'+seq[curr_pos]+'</a></mark>'	
				else:#One mutation case
					url = pos[1][0].get_absolute_url()
					pos_mute=pos_mute + pos[1][0].mutant_aa.one_letter_code +'('+pos[1][0].type.type+') '	
					mod_seq = mod_seq + '<mark style="background-color:'+color+'"><a href="'+url+'" target="_blank"style="color:white;" title="Position '+str(curr_pos+1)+', Mutation: ' +pos_mute+'">'+seq[curr_pos]+'</a></mark>'	

					
				
				if(len(sorted_mod)>0):
					pos = sorted_mod.pop(0)				
				else:				
					while curr_pos <len(seq)-1:
						curr_pos+=1
						if((curr_pos)%10 ==0):
							mod_seq = mod_seq + '  '							
						mod_seq = mod_seq+seq[curr_pos]
						
					break		
		return mod_seq
	###########################################################
	def get_color(self):
		if self.type.type =='Disease':
			return  set('red')
		elif self.type.type =='Polymorphism':	
			return 'green'
		elif self.type.type =='Unclassified':	
			return 'blue'
		return set('purple')	
	


	def interacting_partners(self):
		chains = self.chains.all()
		partners = []
		for chain in chains:
			p = chain.binding_partners()
			for partner in p:
				partners.append(partner.uniprot)
		# Return list of Uniprot objects
		return set(partners)

	def interactions(self):
		chains = self.chains.all()
		output = []
		for chain in chains:
			i1 = chain.interactions_1.all()
			i2 = chain.interactions_2.all()
			for interaction in i1:
				output.append(interaction)
			for interaction in i2:
				output.append(interaction)
		#Return a list of Interaction objects
		return set(output)

	def get_absolute_url(self):
		return reverse('uniprot-view', kwargs={'pk': self.acc_number})

class AminoAcid(models.Model):
    one_letter_code = models.CharField(max_length=1L, primary_key=True)
    three_letter_code = models.CharField(max_length=3L)
    name = models.CharField(max_length=255L)
    class Meta:
        db_table = 'amino_acid'

class UniprotResidue(models.Model):
    id = models.IntegerField(primary_key=True)
    uniprot = models.ForeignKey(Uniprot,db_column='uniprot_acc_number',related_name="residues")
    position = models.IntegerField(db_column='uniprot_position')
    amino_acid = models.ForeignKey(AminoAcid,db_column='amino_acid',related_name="+")
    class Meta:
        db_table = 'uniprot_residue'

class Chain(models.Model):
    id = models.IntegerField(primary_key=True)
    pdb_id = models.CharField(max_length=4L, blank=True)
    pdb_chain = models.CharField(max_length=10L)
    pdb_model = models.IntegerField()
    seq_identity = models.CharField(max_length=10L)
    coverage = models.CharField(max_length=10L)
    seq_start = models.CharField(max_length=10L)
    seq_end = models.CharField(max_length=10L)
    uniprot = models.ForeignKey(Uniprot,db_column='uniprot_acc_number',related_name='chains')
    class Meta:
        db_table = 'chain'

    def binding_partners(self):
        partners = []
        i1 = self.interactions_1.all()
        i2 = self.interactions_2.all()

        for interaction in i1:
            partners.append(interaction.chain_2)
        for interaction in i2:
            partners.append(interaction.chain_1)

        # Returns list of chains
        return(partners)



class InteractionType(models.Model):
    inter_type = models.CharField(db_column="type",max_length=20L, primary_key=True)
    class Meta:
        db_table = 'interaction_type'

class Interaction(models.Model):
    id = models.IntegerField(primary_key=True)
    chain_1 = models.ForeignKey(Chain,db_column='chain_1_id',related_name='interactions_1')
    chain_2 = models.ForeignKey(Chain,db_column='chain_2_id',related_name='interactions_2')
    inter_type = models.ForeignKey(InteractionType,db_column='type',related_name='+')
    filename = models.CharField(max_length=100L)
    class Meta:
        db_table = 'interaction'

    def get_snvs(self):
        cr_withsnv_partner1 = []
        for chain_residue in self.chain_1.residues.all():
            for ur in chain_residue.uniprot_residue.all():
                for snv in ur.snvs.all():
                    cr_withsnv_partner1.append(chain_residue)
        cr_withsnv_partner2 = []
        for chain_residue in self.chain_2.residues.all():
            for ur in chain_residue.uniprot_residue.all():
                for snv in ur.snvs.all():
                    cr_withsnv_partner2.append(chain_residue)

        return [set(cr_withsnv_partner1),set(cr_withsnv_partner2)]

class ChainResidue(models.Model):
    id = models.IntegerField(primary_key=True)
    chain = models.ForeignKey(Chain,db_column='chain_id',related_name='residues')
    position = models.CharField(max_length=10L,db_column='chain_position')
    amino_acid = models.ForeignKey(AminoAcid,related_name="+",db_column='amino_acid')
    uniprot_residue = models.ManyToManyField(UniprotResidue,through='PositionMapping',related_name='chain_residues')
    class Meta:
        db_table = 'chain_residue'

class PositionMapping(models.Model):
    id = models.IntegerField(primary_key=True)
    uniprot_residue = models.ForeignKey(UniprotResidue,db_column='uniprot_residue_id',related_name='+')
    chain_residue = models.ForeignKey(ChainResidue,db_column='chain_residue_id',related_name='+')
    class Meta:
        db_table = 'position_mapping'

class Accessibility(models.Model):
    id = models.IntegerField(primary_key=True)
    interaction = models.ForeignKey(Interaction,db_column='interaction_id',related_name='+')
    chain_residues = models.ForeignKey(ChainResidue,db_column='chain_residue_id',related_name='accessibilities')
    bound_acc = models.CharField(max_length=10L)
    unbound_acc = models.CharField(max_length=10L)
    disulphide_bridge_no = models.IntegerField()
    class Meta:
        db_table = 'accessibility'

class InterfaceResidue(models.Model):
    id = models.IntegerField(primary_key=True)
    chain_residue = models.ForeignKey(ChainResidue,db_column='chain_residue_id',related_name='interface_residues')
    interaction = models.ForeignKey(Interaction,db_column='interaction_id',related_name='interface_residues')
    class Meta:
        db_table = 'interface_residue'

class SnvType(models.Model):
    type = models.CharField(max_length=20L, primary_key=True)
    class Meta:
        db_table = 'snv_type'


class Snv(models.Model):
	ft_id = models.CharField(max_length=10L, primary_key=True)
	type = models.ForeignKey(SnvType,db_column='type',related_name='snvs')
	wt_aa = models.ForeignKey(AminoAcid,related_name="+",db_column="wt_aa")
	mutant_aa = models.ForeignKey(AminoAcid,related_name="+",db_column="mutant_aa")
	uniprot = models.CharField(max_length=6L,db_column="uniprot_acc_number")
	uniprot_position = models.IntegerField()
	gene_code = models.CharField(max_length=10L, blank=True)
	db_snp = models.CharField(max_length=15L, blank=True)
	uniprot_residue = models.ManyToManyField(UniprotResidue,through='SnvUniprotResidue',related_name="snvs") #similar to wt_aa, but it will be available only when their uniprots exist in database 

	class Meta:
		db_table = 'snv'
	
	def get_marked_seq(self):
		uniobj = Uniprot.objects.filter(acc_number=self.uniprot) #get a list of object (it should return a list with one element inside)
		mod_seq = ''		
		if(len(uniobj)!=0): #if the object is available in our database
			pre_seq =  uniobj[0].sequence #get sequence data
			for curr_pos in range(len(pre_seq)): #run through all position
				if curr_pos%10 ==0:
					mod_seq = mod_seq + '  '
				if curr_pos==self.uniprot_position-1:
					color = get_color([self])
					mod_seq =  mod_seq+'<mark style="background-color:'+color+'"><a style="color:white;" title="Position:'+str(self.uniprot_position)+', Original Amino Acid:'+self.wt_aa.one_letter_code+' ">'+self.mutant_aa.one_letter_code+'</a></mark>'				 		
				else:				
					mod_seq = mod_seq + pre_seq[curr_pos]	
		return mod_seq

	def get_similar_snv(self):
		snvlist = list(set(Snv.objects.filter(uniprot=self.uniprot).filter(uniprot_position=self.uniprot_position))-set([self]))
		return snvlist
	
	def get_absolute_url(self):
		return reverse('snv-view', kwargs={'pk': self.ft_id})

	def get_genbank_id(self):
		handle = Entrez.esearch(db="gene", term="%s[Gene Name] AND homo sapiens[Organism]" % self.gene_code)
		record = Entrez.read(handle)
		gene_id = record["IdList"][0]
		return gene_id

class Disease(models.Model):
	mim = models.IntegerField(primary_key=True)
	name = models.CharField(max_length=255L, blank=True)
	snvs = models.ManyToManyField(Snv,through='SnvDisease',related_name="diseases")

	class Meta:
		db_table = 'disease'

	def get_Snv(self):		#Return a list of related snvs and diseases

		snvlist=self.snvs.all()
		relatedsnv = []
		disease_sum = []
		relateduni = []
		for item in snvlist:				
			#find all diseases from all snvs that share the same uniprot_acc_number		
			related_snvs = Snv.objects.filter(uniprot=item.uniprot)
			related_disease = []			
			for	snv_item in related_snvs:
				for disease_item in snv_item.diseases.all():
					related_disease.append(disease_item)
					disease_sum.append(disease_item)
			related_disease_set = list(set(related_disease)-set(item.diseases.all()))
			if [item.uniprot,list(set(related_disease))] not in relateduni:
				relateduni.append([item.uniprot,list(set(related_disease))])
			relatedsnv.append({'obj':item, 'disease_set':related_disease_set})
		final_related_uniprot = relateduni
		final_disease_sum = list(set(disease_sum)-set([self]))
		return [relatedsnv,final_disease_sum,final_related_uniprot]


	def get_absolute_url(self):
		return reverse('disease-view', kwargs={'pk': self.mim})


class SnvDisease(models.Model):
    id = models.IntegerField(primary_key=True)
    snv = models.ForeignKey(Snv,db_column='ft_id')
    disease = models.ForeignKey(Disease,db_column='mim')
    class Meta:
        db_table = 'snv_disease'

class SnvUniprotResidue(models.Model):
    id = models.IntegerField(primary_key=True)
    snv = models.ForeignKey(Snv,db_column='ft_id')
    uniprot_residue = models.ForeignKey(UniprotResidue, db_column='uniprot_residue_id') #db_column='uniprot_residue_id'
    class Meta:
        db_table = 'snv_uniprot_residue'


class PfamHmm(models.Model):
    acc = models.CharField(max_length=12L, primary_key=True,db_column="hmm_acc")
    name = models.CharField(max_length=50L)
    type = models.CharField(max_length=50L)
    length = models.IntegerField()
    clan = models.CharField(max_length=50L)
    uniprots =  models.ManyToManyField(Uniprot,through='UniprotPfamMapping',related_name='pfam_hmms')
    class Meta:
        db_table = 'pfam_hmm'

class UniprotPfamMapping(models.Model):
    id = models.IntegerField(primary_key=True)
    uniprot = models.ForeignKey(Uniprot,db_column='uniprot_acc_number',related_name='pfam_mappings')
    hmm = models.ForeignKey(PfamHmm,db_column='hmm_acc')
    alignment_start_residue = models.ForeignKey(UniprotResidue,db_column='alignment_start_ur_id',related_name='+')
    alignment_end_residue = models.ForeignKey(UniprotResidue,db_column='alignment_end_ur_id',related_name='+')
    envelope_start_residue = models.ForeignKey(UniprotResidue,db_column='envelope_start_ur_id',related_name='+')
    envelope_end_residue = models.ForeignKey(UniprotResidue,db_column='envelope_end_ur_id',related_name='+')
    hmm_start = models.IntegerField()
    hmm_end = models.IntegerField()
    bit_score = models.DecimalField(max_digits=7, decimal_places=1)
    e_value = models.CharField(max_length=255L)
    significance = models.IntegerField()
    class Meta:
        db_table = 'uniprot_pfam_mapping'

class ActiveSiteResidue(models.Model):
    id = models.IntegerField(primary_key=True)
    mapping = models.ForeignKey(UniprotPfamMapping,db_column='up_mapping_id',related_name="active_site_residues")
    residue = models.ForeignKey(UniprotResidue,db_column='active_site_ur_id',related_name='active_site_residues')
    class Meta:
        db_table = 'active_site_residue'

class SuperpositionMapping(models.Model):
    id = models.IntegerField(primary_key=True)
    ref_chain = models.ForeignKey(Chain,db_column='ref_chain_id',related_name='superpositions_where_ref')
    target_chain = models.ForeignKey(Chain,db_column='target_chain_id',related_name='superpositions_where_target')
    target_chain_letter = models.CharField(max_length=1L)
    class Meta:
        db_table = 'combined_uniprot_mapping'


















