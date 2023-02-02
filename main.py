from pyfaidx import Fasta
import csv

# reading of the file with amiloids motifs
waltz_seq = []       # list for short motifs of amiloid
waltz_dict = {}      # dictionry for sequences from WaltzDB
with open('WALTZ_DB_amiloid_seq') as waltz_db:
    for line in waltz_db:
        line = line.strip()
        waltz_seq.append(line)
    for number, amiloid_motif in enumerate(waltz_seq):
        if amiloid_motif.strip().isalpha() == False:    # deleting empty lines
            break
        else:
            waltz_dict[number] = amiloid_motif


proteome_seq = {}  # dict for all sequences which are shorter than 1000 residues
with Fasta('uniprot-proteome.fasta') as genes:
    for number_of_peptide in range((len(genes.keys()))):
        proteome_seq[genes[number_of_peptide][:].name] = genes[number_of_peptide][:]

'''
# reading tsv file, which contains UniProt names for all peptides from DisProt. They will be used in final fasta
tsv_list = {}  # dict for pairs "UniProt_name - Sequence"
with open('DisProt release_2022_06.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    for rows in file:
        while '' in rows:     # deleting empty elements from tsv
            rows.remove('')
        tsv_list[rows[len(rows)-1]] = rows[0]
'''

# finding peptides with amiloid motifs
amiloid_motifs_list_for_header = {}  # remembering amiloid motifs for each peptide, which have amiloid-forming region
filter_seq = {}  # dict for seq that have amiloid-forming regions
for key, seq in proteome_seq.items():
    amiloid_motifs_list = []   # updated list for all amiloid motifs in one peptide
    for val in waltz_dict.values():
        if val in str(seq):  #
            filter_seq[key] = seq
            pos = str(seq).find(val)+1
            val = str(pos) + val
            amiloid_motifs_list.append(val)
        amiloid_motifs_list_for_header[key] = amiloid_motifs_list

discriptor = []
with Fasta('uniprot-proteome.fasta') as genes:
    for record in genes:
        line = record.long_name
        discriptor.append(str(line))

amiloid_seq_no_amil_disc = {}
for key, val in filter_seq.items():
    for disc in discriptor:
        if key in disc:
            amiloid_seq_no_amil_disc[disc] = val

amil_seq = {}
fasta_uniprot = {}
for key, val in amiloid_seq_no_amil_disc.items():
    head_for_fasta = str(key)
    for key_amil, val_amil in amiloid_motifs_list_for_header.items():
        if key_amil in key:
            for i in val_amil:
                if ' amyloid_seq=' not in head_for_fasta:
                    head_for_fasta += ' amyloid_seq=' + i
                else:
                    head_for_fasta += '|' + i
    amil_seq[head_for_fasta] = val

table = []
with open('DisProt release_2022_06.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    count = 0
    for rows in file:
        if count == 0:
            table.append(rows)
            count += 1
        for pep, uni in new_tsv.items():
            if uni in rows and pep in rows:
                table.append(rows)

'''
with open('table_amiloids.tsv', 'w') as tab:
    for elem in table:
        for words in elem:
            tab.write(words + ', ')
        tab.write('\n')
'''

with open('Amyloid+Disprot_full_corr', 'w') as A_D:
    for key, val in amil_seq.items():
        A_D.write(str(key) + '\n')
        A_D.write(str(val) + '\n')


with open('UniProt+DisProt_amyloids.tsv', 'w') as tab:
    for elem in table_short:
        for words in elem:
            tab.write(words + ', ')
        tab.write('\n')