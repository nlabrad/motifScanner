from django.core.files.storage import FileSystemStorage
from django.shortcuts import render
from django.views.generic import FormView
from .utils import motifscanner
from .forms import motifForm


# Create your views here.

def upload(request):
    if request.method == 'POST':
        uploaded_file = request.FILES['fasta']
        fs = FileSystemStorage()
        fs.save(uploaded_file.name, uploaded_file)
    return render(request, 'upload.html')


def find(request):
    motifs = "none"
    if request.method == "POST":
        # Get the posted form
        genbankForm = motifForm(request.POST)
        if genbankForm.is_valid():
            genbankId = genbankForm.cleaned_data['genbank_id']
            genbankSeq = str(motifscanner.fetch(genbankId))
            itsseq = motifscanner.extractITS(genbankSeq)
            motifs = motifscanner.findMotifs(itsseq[1])
            return render(request, 'results.html', {'motifs': motifs})
    else:
        genbankForm = motifForm()
    return render(request, 'find.html')


class GenbankSearchView(FormView):
    template_name = "motifsearch.html"
    form_class = motifForm
    success_url = "/results/"

    def form_valid(self, form):
        form.findMotif()
        return super().form_valid(form)

