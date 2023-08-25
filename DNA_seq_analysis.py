
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# make a codon dictionary for plants

""" plant_archea_bacteria_codons.txt:
        AGT S
        AGC S
        AGA R
        AGG R
        GGT G
        GGC G
        GGA G
        GGG G
"""
# DNA_codons_plants_archea_bacteria = {}
# with open("plant_archea_bacteria_codons.txt", "r") as file:
#     for line in file:
#         codon, aa = line.strip().split()
#         DNA_codons_plants_archea_bacteria[codon] = aa

DNA_codons_plants_archea_bacteria = {
        # 'L, I, M, V' - START, '*' - STOP
        'TTT': 'F', 'TTC': 'F', 
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 
        'CAT': 'H', 'CAC': 'H', 
        'CAA': 'Q', 'CAG': 'Q', 
        'AAT': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K', 
        'GAT': 'D', 'GAC': 'D', 
        'GAA': 'E', 'GAG': 'E', 
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 
        'TGG': 'W', 
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
        'AGT': 'S', 'AGC': 'S', 
        'AGA': 'R', 'AGG': 'R', 
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

# fetch and parse DNA sequence
Entrez.email = "sviatoslavkharuk@gmail.com" 

with Entrez.efetch(db = "nucleotide", 
                   id = "EU490707", 
                   rettype = "fasta",
                   retmode = "text") as handle:
    seq_record = SeqIO.read(handle, "fasta")
    print(f"Sample ID: {seq_record.id}")
    seq_parse = repr(seq_record.seq)
    print(f"The length of sequence is {(len(seq_record.seq))}")
    print(f"Sequence object: {seq_parse}")

# calculate the GC content
my_seq = Seq(seq_parse)
print(f"The GC content of sequence is {gc_fraction(my_seq)*100:.2f}%")

# function for Translation the sequence 
def translateSeq(seq, init_pos = 0):
     """Translates the DNA to protein's amino acid sequence"""
     return [DNA_codons_plants_archea_bacteria[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]

# function for ORFs generation
def generetesORF(seq):
     """Generates Open Reading Frames (ORF)"""
     frame_list = []
     frame_list.append(translateSeq(seq, 0))
     frame_list.append(translateSeq(seq, 1))
     frame_list.append(translateSeq(seq, 2))
     return frame_list

# function for generating the proteins from ORFs
def proteinsfromORF(aa_seq):
     current_prot = []
     proteins = []
          ## якщо є стоп кодон, тоді нічого не роби
     for aa in aa_seq:
        if aa == "*":
            ## якщо стоп кодон, тоді все що є in current_proteins add to protein 
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
     ## якщо є старт кодон, тоді почни брати і робити білок
        else:     
            if aa == "L" or aa == "I" or aa == "M" or aa == "V":
              current_prot.append("")
            ## it will add an amino acid to the end of the i (our protein string)
            for i in range(len(current_prot)):
              current_prot[i] += aa
     return proteins


# generate ORFs
aa_prot = generetesORF(seq_record.seq)

# unpack the lists of ORF into one list
new_prot_list = [i for sublist in aa_prot for i in sublist]

# get proteins from the ORF
# print(proteinsfromORF(new_prot_list))

# first protein from ORF. Just copy and paste
first_prot = "IFYEPVEIFGYDNKSSLVLVKRLITRMYQQNFLISSVNDSNQKGFWGHKHFFSSHFSSQMVSEGFGVILEIPFSSQLVSSLEEKKIPKYQNLRSIHSIFPFLEDKFLHLNYVSDLLIPHPIHLEILVQILQCRIKDVPSLHLLRLLFHEYHNLNSLITSKKFIYAFSKRKKRFLWLLYNSYVYECEYLFQFLRKQSSYLRSTSSGVFLERTHLYVKIEHLLVVCCNSFQRILCFLKDPFMHYVRYQGKAILASKGTLILMKKWKFHLVNFWQSYFHFWSQPYRIHIKQLSNYSFSFLGYFSSVLENHLVVRNQMLENSFIINLLTKKFDTIAPVISLIGSLSKAQFCTVLGHPISKPIWTDFSDSDILDRFCRICRNLCRYHSGSSKKQVLYRIKYILRLSCARTLARKHKSTVRTFMRRLGSGLLEEFFMEEEFFTNLWKFLVMTINLV"

# BLASTp first protein from ORF
result_BLAST = NCBIWWW.qblast("blastp", "nr", first_prot)
with open("new_blast.xml", "w+") as out_handle:
    out_handle.write(result_BLAST.read())

result_BLAST.close()

# parsing BLAST sequence
result_BLAST = open("new_blast.xml")
blast_records = NCBIXML.read(result_BLAST)
e_value_thresh = 0.04
for alignment in blast_records.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < e_value_thresh:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
        

