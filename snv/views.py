# NEW APP NEW APP NEW APP NEW APP
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django.core.urlresolvers import reverse

# Class Import
from snv.models import Uniprot, Snv, Disease

#View Import
#for listing
from django.views.generic import ListView
#for creating
from django.core.urlresolvers import reverse #for edit & create
#from django.views.generic import CreateView
#for editing
#from django.views.generic import UpdateView
#for deleting
#from django.views.generic import DeleteView
#for detail view
from django.views.generic import DetailView


class UniprotList(ListView):
	model = Uniprot
	template_name = 'Uniprot_list.html'

class UniprotView(DetailView):
	model = Uniprot
	template_name = 'Uniprot_view.html'

class DiseaseView(DetailView):
	model = Disease
	template_name = 'Disease_view.html'

class SnvView(DetailView):
	model = Snv
	template_name = 'Snv_view.html'

class home(ListView):
	#return render(request, 'home.html')
	model = Uniprot
	template_name = 'home.html'

#####################  SEARCH #####################################
from django.template.response import TemplateResponse
from django.shortcuts import render
from django.db.models import Q

def search_form(request):
	return render(request, 'search_form.html')

def search(request):
	if 'q' in request.GET and request.GET['q']:
		q = request.GET['q']
		option = request.GET['search_select']
		if option =='1':
			uniprot = Uniprot.objects.filter(acc_number__icontains=q)
			return render(request, 'Uniprot_search_results.html',
		                    {'uniprot': uniprot, 'query': q})
		elif option =='2':
			snv = Snv.objects.filter(Q(ft_id__icontains=q)|Q(gene_code__icontains=q)|Q(db_snp__icontains=q))
			return render(request, 'Snv_search_results.html',
		                    {'snves': snv, 'query': q})
		else:
			disease = Disease.objects.filter(Q(name__icontains=q)|Q(mim__icontains=q))
			return render(request, 'Disease_search_results.html',
		                    {'diseases': disease, 'query': q})
	else:
		return TemplateResponse(request, 'home.html', {'error': 'Please input your query'})
		#return ValidationError('Invalid value')
		#return HttpResponse('Please submit a search term.')
		#url = reverse('home', kwargs={'error': 'Please input your query'})

		#return HttpResponseRedirect(url)

######################################################################
