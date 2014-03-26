# NEW APP NEW APP NEW APP NEW APP
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django.core.urlresolvers import reverse
from django.shortcuts import render, get_object_or_404

# Class Import
from snv.models import *

#View Import
#for listing
from django.views.generic import ListView, View, DetailView
#for creating
from django.core.urlresolvers import reverse #for edit & create
#from django.views.generic import CreateView
#for editing
#from django.views.generic import UpdateView
#for deleting
#from django.views.generic import DeleteView




class UniprotList(ListView):
	model = Uniprot
	template_name = 'Uniprot_list.html'

class UniprotView(DetailView):
	model = Uniprot
	template_name = 'Uniprot_view.html'
	
	def get_context_data(self, **kwargs):
		# Call the base implementation first to get a context
		data = super(UniprotView, self).get_context_data(**kwargs)
		data['chains']  =self.object.get_pdb_align()
		interact_snvs = self.object.interactions() #also initialise data for marking
		data['interactions'] = interact_snvs[0]
		data['mapped_seq']   = interact_snvs[1]
 		snvs_data   = self.object.get_Snv()		
		data['snv_list']= snvs_data[0]
		data['snv_statistic']=snvs_data[1]
		
		data['marked_reg'] = interact_snvs[2]
		data['annotation'] = self.object.annotation()
		return data
	
class DiseaseView(DetailView):
	model = Disease
	template_name = 'Disease_view.html'   

	def get_context_data(self, **kwargs):
		# Call the base implementation first to get a context
		snv = super(DiseaseView, self).get_context_data(**kwargs)
		# Add in a QuerySet of all the ralated snv
		obtained_data = self.object.get_Snv()		
		snv['snv_list']= obtained_data[0]
		snv['diseases']= obtained_data[1]
		snv['uniprots']= obtained_data[2]
		return snv


class SnvView(DetailView):
	model = Snv
	template_name = 'Snv_view.html'

	def get_context_data(self, **kwargs):	
		response = super(SnvView, self).get_context_data(**kwargs)
		response['marked_seq']= self.object.get_marked_seq()
		response['sim_snv']=self.object.get_similar_snv()
		return response

class InteractionView(DetailView):
	model = Interaction
	template_name = 'Interaction_view.html'

	def get_context_data(self, **kwargs):

		interaction = super(InteractionView, self).get_context_data(**kwargs)
		# Get snvs
		snvs = self.object.get_snv_chain_residues()
		interaction["snvdict1"] = snvs[2]
		interaction["snvdict2"] = snvs[3]
		interaction["sortedsnvdict1"] = snvs[4]
		interaction["sortedsnvdict2"] = snvs[5]
		# Get converted positions for snvs and add to list
		chain1_snv_positions = []
		for cr in snvs[0]:
			chain1_snv_positions.append(cr.get_transformed_position(self.object))
		chain2_snv_positions = []
		for cr in snvs[1]:
			chain2_snv_positions.append(cr.get_transformed_position(self.object))
		# Pass snv position lists to interaction dictionary
		interaction['chain1_snv_positions'] = chain1_snv_positions
		interaction['chain2_snv_positions'] = chain2_snv_positions


		# Will be 2 dictionaries with this format
		# [{mapping:[start_pdb_position,end_pdb_position]} for chain 1, {mapping:[start_pdb_position,end_pdb_position]} for chain 2]
		pfam_positions = self.object.get_pfam_mapping_positions()
		interaction['chain1_pfam_positions'] = pfam_positions[0]
		interaction['chain2_pfam_positions'] = pfam_positions[1]

		#Get pre-computed contacts and add to list
		contacts = {}
		for contact in StoredContact.objects.filter(interaction=self.object):
			try:
				contacts[contact.bond_type].append(contact)
			except KeyError:
				contacts[contact.bond_type] = [contact]

		interaction['contacts'] = contacts
		#interaction['contacts'] = self.object.get_contacts()


		return interaction


class InterfaceView(View):
	#model = Interaction
	#template_name = 'Interface_view.html'
	
	def get(self,request,pk):
		interaction = Interaction.objects.get(id=pk)

		left_interface_residues = interaction.chain_1.get_interface_residues(interaction)
		right_interface_residues = interaction.chain_2.get_interface_residues(interaction)
		interface_residues = left_interface_residues + right_interface_residues
		
		uniprot_positions = {}
		amino_acids = {}
		group2num = {}
		count = 0
		for ir in interface_residues:
			try:
				ur = ir.chain_residue.uniprot_residue.all()[0]
				uniprot_positions[ir] = ur.position
			except IndexError:
				uniprot_positions[ir] = "None"
			group = ur.amino_acid.group
			try:
				group_no = group2num[group]
			except KeyError:
				group2num[group] = count
				group_no = count
				count += 1
			amino_acids[ir] = [ur.amino_acid.three_letter_code,group_no]
			

		# Get orders
		orders = {}
		if left_interface_residues[0].res_order == None:
			print(orders)
			l_order = 1
			for ir in left_interface_residues:
				orders[ir] = l_order
				l_order += 1
		else:
			for ir in left_interface_residues:
				orders[ir] = ir.res_order
		if right_interface_residues[0].res_order == None:
			r_order = 1
			for ir in right_interface_residues:
				orders[ir] = r_order
				r_order += 1
		else:
			for ir in right_interface_residues:
				orders[ir] = ir.res_order
		# JS Object format
		#{ id : , chain : "A", order : , aa : , uniprot_pos : , chain_pos : , aa_group : }
		nodes = []
		for ir in interface_residues:
			print(ir)
			if ir in left_interface_residues:
				chain_id = "A"
			else:
				chain_id = "B"
			empty_js_object = " id : '{id}' , chain : '{chain}', order : {order}, aa : '{aa}', uniprot_pos : '{uni}', chain_pos : '{chain_pos}', aa_group_no : '{aa_group_no}' "
			js_object = "{" + empty_js_object.format(
										id=ir.id,
										chain=chain_id,
										order=orders[ir],
										aa=amino_acids[ir][0],
										uni=uniprot_positions[ir],
										chain_pos=ir.get_position(),
										aa_group_no=amino_acids[ir][1]
										) + "}"
			nodes.append(js_object)
		# Get edges
		# Left irs mapped to right
		edges = {}
		bond_types = []
		for ir in left_interface_residues:
			inter_dict = ir.get_interacting_residues()
			edges[ir] = inter_dict
			for bt in inter_dict.values():
				bond_types.append(bt)

		bond_types = list(set(bond_types))


		# Return dictionary
		interface = {
			'nodes' : nodes,
			'edges' : edges,
			'bond_types' : bond_types,
			'group2num' : group2num
		}

		return render(request,'Interface_view.html',interface)



class SuperpositionView(View):

	def get(self,request,uniprot_acc):
		uniprot = Uniprot.objects.get(acc_number=uniprot_acc)
		ref_chain = uniprot.get_ref_chain()
		assignments = SuperpositionMapping.objects.filter(ref_chain=ref_chain)
		filename = uniprot.acc_number+"-interactions.pdb"
		return render(request,'Superposition_view.html', {'uniprot':uniprot,'ref_chain':ref_chain,'assignments':assignments,'filename':filename})



###################### BASIC VIEW ##################################

def home(request):
	return render(request, 'home.html')

def info(request):
	return render(request, 'info.html')

def about(request):
	return render(request, 'about.html')

def gettingstarted(request):
	return render(request, 'getting-started.html')

#####################  SEARCH #####################################
from django.template.response import TemplateResponse
from django.shortcuts import render
from django.db.models import Q

def search_form(request):
	return render(request, 'search_form.html')

def search(request):
	if 'q' in request.GET and request.GET['q']:
		q = request.GET['q'].strip(' \t\n\r')
		option = request.GET['search_select']
		if option =='1':
			uniprot = list(set(Uniprot.objects.filter(Q(acc_number__icontains=q)|Q(name__icontains=q)|Q(gene_code__icontains=q)|Q(genbank_id__icontains=q)|Q(chains__pdb_id__icontains=q))))
			if len(uniprot) == 1:
				return HttpResponseRedirect(reverse('uniprot-view', kwargs={'pk': uniprot[0].acc_number})) 
			else:
				return render(request, 'Uniprot_search_results.html',
		                    {'uniprot': uniprot, 'query': q})
		elif option =='2':
			snv = list(set(Snv.objects.filter(Q(ft_id__icontains=q)|Q(db_snp__icontains=q)|Q(diseases__name__icontains=q))))
			if len(snv) == 1:
				return HttpResponseRedirect(reverse('snv-view', kwargs={'pk': snv[0].ft_id})) 
			else:
				return render(request, 'Snv_search_results.html',
		                    {'snves': snv, 'query': q})
		else:
			disease = Disease.objects.filter(Q(name__icontains=q)|Q(mim__icontains=q))
			if disease.count() == 1:
				return HttpResponseRedirect(reverse('disease-view', kwargs={'pk': disease[0].mim})) 
			else:			
				return render(request, 'Disease_search_results.html',
		                    {'diseases': disease, 'query': q})
	else:
		return TemplateResponse(request, 'home.html', {'error': 'Please input your query'})

######################################################################
