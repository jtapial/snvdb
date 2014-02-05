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


class Uniprot(models.Model):
    acc_number = models.CharField(max_length=6L, primary_key=True)
    sequence = models.TextField(blank=True)
    class Meta:
        db_table = 'uniprot'
        
    def get_absolute_url(self):
        return reverse('uniprot-view', kwargs={'pk': self.acc_number})

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
    amino_acid = models.ForeignKey(AminoAcid,related_name="+")
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
    type = models.CharField(max_length=20L, primary_key=True)
    class Meta:
        db_table = 'interaction_type'

class Interaction(models.Model):
    id = models.IntegerField(primary_key=True)
    chain_1 = models.ForeignKey(Chain,db_column='chain_1_id',related_name='interactions_1')
    chain_2 = models.ForeignKey(Chain,db_column='chain_2_id',related_name='interactions_2')
    type = models.ForeignKey(InteractionType,db_column='type',related_name='+')
    filename = models.CharField(max_length=100L)
    class Meta:
        db_table = 'interaction'

class ChainResidue(models.Model):
    id = models.IntegerField(primary_key=True)
    chain = models.ForeignKey(Chain,db_column='chain_id',related_name='residues')
    position = models.CharField(max_length=10L,db_column='chain_position')
    amino_acid = models.ForeignKey(AminoAcid,related_name="+")
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
    uniprot = models.ForeignKey(Uniprot,db_column="uniprot_acc_number",related_name="snvs")
    uniprot_position = models.IntegerField()
    gene_code = models.CharField(max_length=10L, blank=True)
    db_snp = models.CharField(max_length=15L, blank=True)

    uniprot_residue = models.ManyToManyField(UniprotResidue,through='SnvUniprotResidue',related_name="snvs")
    class Meta:
        db_table = 'snv'
    def get_absolute_url(self):
        return reverse('snv-view', kwargs={'pk': self.ft_id})

class Disease(models.Model):
    mim = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=255L, blank=True)
    snvs = models.ManyToManyField(Snv,through='SnvDisease',related_name="diseases")
    class Meta:
        db_table = 'disease'
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
    uniprot_residue = models.ForeignKey(UniprotResidue,db_column='uniprot_residue_id')
    class Meta:
        db_table = 'snv_uniprot_residue'
# AUTO




















