from django import forms

class motifForm(forms.Form):
    genbank_id = forms.CharField(label='genbank_id', max_length=800)
    def findMotif(self):
        pass