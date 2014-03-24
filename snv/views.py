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
		data['markscript'] = interact_snvs[3]
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

		interaction['contacts'] = self.object.get_contacts()	

		return interaction


class InterfaceView(DetailView):
	model = Interaction
	template_name = 'Interface_view.html'

	def get_context_data(self, **kwargs):

		interface = super(InterfaceView, self).get_context_data(**kwargs)

		return interface



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
			uniprot = Uniprot.objects.filter(Q(acc_number__icontains=q)|Q(name__icontains=q)|Q(gene_code__icontains=q)|Q(genbank_id__icontains=q))
			if uniprot.count() == 1:
				return HttpResponseRedirect(reverse('uniprot-view', kwargs={'pk': uniprot[0].acc_number})) 
			else:
				return render(request, 'Uniprot_search_results.html',
		                    {'uniprot': uniprot, 'query': q})
		elif option =='2':
			snv = Snv.objects.filter(Q(ft_id__icontains=q)|Q(db_snp__icontains=q))
			if snv.count() == 1:
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
