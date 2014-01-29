# NEW APP NEW APP NEW APP NEW APP
from django.http import HttpResponse 
from django.shortcuts import render

# Class Import
from snv.models import Uniprot

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

class home(ListView):
	#return render(request, 'home.html')
	model = Uniprot
	template_name = 'home.html'

#####################  SEARCH #####################################
from django.shortcuts import render

def search_form(request):
	return render(request, 'search_form.html')

def search(request):
	if 'q' in request.GET and request.GET['q']:
		q = request.GET['q']
		uniprot = Uniprot.objects.filter(acc_number__icontains=q)
		return render(request, 'Uniprot_search_results.html',
                        {'uniprot': uniprot, 'query': q})
	else:
		return HttpResponse('Please submit a search term.')

######################################################################
