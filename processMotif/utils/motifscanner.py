import re
from Bio import Entrez, __version__, SeqIO, Seq  # import the required libraries


def fetch(accession):  # fetch sequence and taxonomy with accession#
    Entrez.email = 'nlab@fastmail.com'  # email reported to entrez to associate with the query
    with Entrez.efetch(db='nucleotide', id=accession,
                       rettype='gb') as handle:  # id is what genbank ID to query, type is genbank
        record = SeqIO.read(handle, "gb")  # get sequence
        genbank_organism_taxonomy = record.annotations["taxonomy"]
        genbank_organism_sequence = record.seq
        # print ("Taxonomy: ", genbank_organism_taxonomy)
        # print ("Sequence:", genbank_organism_sequence)
    return genbank_organism_sequence


def get_tax_id(organism):  # pulls organism taxaID
    Entrez.email = 'nlab@fastmail.com'  # email reported to entrez to associate with the query
    organism = organism.replace(" ", "+").strip()  # replace whitespaces with +
    handle = Entrez.esearch(term=organism, db="taxonomy", retmode="xml")  # search for the organism
    record = Entrez.read(handle)
    return record['IdList'][0]  # IdList[0] contains the taxID


def get_taxonomy(tax_id):  # uses TaxaID to get the taxa
    with Entrez.efetch(db='taxonomy', id=tax_id, retmode='xml') as handle:
        record = Entrez.read(handle, validate=False)
        print(record)
        return record


# input: fasta file or search genbank
# output: motifs
# 1 Find the ITS region first. 16S ends with CCTCCTT. Ends in CAGAC: 16S...CCTCCTT-ITS...CAGAC-23S...
# 2 Find motifs in ITS region


# Extraction of ITS region from a 16S-23S region.
# Find CCTCCTT, get rid of what's before that.
def extractITS(geneRegion):
    if (re.search("CCTCCTT", geneRegion) == None):  # If It cannot find the start of the ITS region4
        print("ITS Region not found in this sequence")
        return None
    else:
        its_start_position = re.search("CCTCCTT", geneRegion).start() + 7
        itsseq = geneRegion[its_start_position:]
        return its_start_position, itsseq


def findMotifs(its_region):  # Find the motifs

    motifs = {"ITS": its_region}  # dictionary. motifs[motif-name]=[start-position,end-position]
    d1d1_search = re.search(r"GACC[CTG](.*?)AGGTC",
                           its_region)  # find text between basal clamps, starting with GACCT/C to the first AGGTC (*? is lazy search)
    if d1d1_search == None:  # if the above pattern is not found, try the one below
        d1d1_search = re.search(r"GACCA(.*?)[TA]GGTC",
                               its_region)  # find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if d1d1_search == None:  # If that one is not found either, then d1d1 cannot be found
        motifs["d1d1"] = [None, None]
    else:
        # motifs["leader"] = [its_region[0:d1d1_search.start()], len(its_region[0:d1d1_search.start():])]
        motifs["leader"] = [0, d1d1_search.start()]
        motifs["d1d1"] = [d1d1_search.start(), d1d1_search.end()]  # store start/end positions of the d1d1 region
        index_shift = d1d1_search.end()  # count of processed characters. Used to keep track of absolute position in the index.
        its_region = its_region[d1d1_search.end():]  # trim processed its region to avoid overlapping searches.

    # spacer-D2-spacer D3-spacer region

    # find tRNA-Ile(1)
    t_rna1_search = re.search(r"GGGCTATTA(.*?)GGCCCA", its_region)  # find text between basal clamps
    if t_rna1_search is None:
        t_rna1_search = re.search(r"GACCA(.*?)[TA]GGTC",
                                its_region)  # find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if t_rna1_search is None:
        motifs["tRNA1"] = None
    else:
        motifs["d2d3"] = [d1d1_search.end(), t_rna1_search.start() + index_shift]
        motifs["tRNA1"] = [t_rna1_search.start() + index_shift,
                           t_rna1_search.end() + index_shift]  # inded_shift accounts for the trimmed portion of the string.
        index_shift = index_shift + t_rna1_search.end()  # add to the index_shift what will be trimmed in the next line
        its_region = its_region[t_rna1_search.end():]  # trim processed its region

    # find tRNA-Ala(2)
    t_rna2_search = re.search(r"GGGG(.*?)CCTCCA", its_region)  # find text between basal clamps
    if t_rna2_search is None:
        t_rna2_search = re.search(r"GACCA(.*?)[TA]GGTC",
                                its_region)  # find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if t_rna2_search is None:
        motifs["tRNA2"] = None
    else:
        motifs["tRNA2"] = [t_rna2_search.start() + index_shift, t_rna2_search.end() + index_shift]
        index_shift = index_shift + t_rna2_search.end()
        its_region = its_region[t_rna2_search.end():]  # trim processed its region

    # find BoxB
    box_b_search = re.search(r"AGCA(.*?)TGCT", its_region)  # find text between basal clamps
    if box_b_search is None:
        box_b_search = re.search(r"GACCA(.*?)[TA]GGTC",
                               its_region)  # find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if box_b_search is None:
        motifs["BoxB"] = None
    else:
        motifs["BoxB"] = [box_b_search.start() + index_shift, box_b_search.end() + index_shift]
        index_shift = index_shift + box_b_search.end()
        its_region = its_region[box_b_search.end():]  # trim remaining its region

    # find Box-A
    box_a_search = re.search(r"G[CA]A(.*?)GAAA", its_region)  # find text between basal clamps
    if box_a_search is None:
        box_a_search = re.search(r"GACCA(.*?)[TA]GGTC",
                               its_region)  # find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if box_a_search is None:
        motifs["BoxA"] = None
    else:
        motifs["BoxA"] = [box_a_search.start() + index_shift, box_a_search.end() + index_shift]
        index_shift = index_shift + box_a_search.end()
        its_region = its_region[box_a_search.end():]  # trim remaining its region

    # find D4
    d4_search = re.search(r"ACT(.*?)AAA", its_region)
    if d4_search is None:
        d4_search = re.search(r"GACCA(.*?)[TA]GGTC",
                             its_region)  # find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if d4_search is None:
        motifs["D4"] = None
    else:
        motifs["D4"] = [d4_search.start() + index_shift, d4_search.end() + index_shift]
        index_shift = index_shift + d4_search.end()
        its_region = its_region[d4_search.end():]  # trim remaining its region

    return motifs
    # end of extractMotif


# Put your sequence here:

# accessionNum = input("Enter Acession#: ")
# genbankSeq = str(fetch(accessionNum))
# itsseq = extractITS(genbankSeq)
# extractedMotifs = findMotifs(itsseq[1])

# print("\nITS Region: ", extractedMotifs["ITS"])
# if (extractedMotifs['leader'] is not None):
#     print("\nLeader: ", extractedMotifs['ITS'][extractedMotifs['leader'][0]:extractedMotifs['leader'][1]])
#
# if (extractedMotifs['d1d1'] is not None):
#     print("\nD1D1: ", extractedMotifs['ITS'][extractedMotifs['d1d1'][0]:extractedMotifs['d1d1'][1]])
# else:
#     print("Warning! d1d1 was not found, please verify manually")
#
# if (extractedMotifs['d2d3'] is not None):
#     print("\nSP-D2-D3-SP: ", extractedMotifs['ITS'][extractedMotifs['d2d3'][0]:extractedMotifs['d2d3'][1]])
#
# if (extractedMotifs['tRNA1'] is not None):
#     print("\ntRNA1: ", extractedMotifs['ITS'][extractedMotifs['tRNA1'][0]:extractedMotifs['tRNA1'][1]])
# else:
#     print("Warning! tRNA1 was not found, please verify manually")
#
# if (extractedMotifs['tRNA2'] is not None):
#     print("\ntRNA2: ", extractedMotifs['ITS'][extractedMotifs['tRNA2'][0]:extractedMotifs['tRNA2'][1]])
# else:
#     print("Warning! tRNA2 was not found, please verify manually")
#
# if (extractedMotifs['BoxB'] is not None):
#     print("\nBoxB: ", extractedMotifs['ITS'][extractedMotifs['BoxB'][0]:extractedMotifs['BoxB'][1]])
# else:
#     print("Warning! BoxB was not found, please verify manually")
#
# if (extractedMotifs['BoxA'] is not None):
#     print("\nBoxA: ", extractedMotifs['ITS'][extractedMotifs['BoxA'][0]:extractedMotifs['BoxA'][1]])
# else:
#     print("Warning! BoxA was not found, please verify manually")
#
# if (extractedMotifs['D4'] is not None):
#     print("\nD4: ", extractedMotifs['ITS'][extractedMotifs['D4'][0]:extractedMotifs['D4'][1]])
# else:
#     print("Warning! D4 was not found, please verify manually")

# print(extractedMotifs['d1d1'][1])





