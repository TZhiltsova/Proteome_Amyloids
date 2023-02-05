from pyfaidx import Fasta
import csv

# reading of the file with amyloids motifs
waltz_seq = []       # list for short motifs of amyloid
waltz_dict = {}      # dictionary for sequences from WaltzDB
with open('WALTZ_DB_amyloid_seq') as waltz_db:
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

# finding peptides with amyloid motifs
new_csv = []
amyloid_motifs_list_for_header = {}  # remembering amyloid motifs for each peptide, which have amyloid-forming region
filter_seq = {}  # dict for seq that have amyloid-forming regions
for key, seq in proteome_seq.items():
    amyloid_motifs_list = []   # updated list for all amyloid motifs in one peptide
    for val in waltz_dict.values():
        if val in str(seq):
            filter_seq[key] = seq
            pos = str(seq).find(val)+1
            val = str(pos) + val
            amyloid_motifs_list.append(val)
            name = key.split('|')
            if name[1] not in new_csv:
                new_csv.append(name[1])
        amyloid_motifs_list_for_header[key] = amyloid_motifs_list

discriptor = []
with Fasta('uniprot-proteome.fasta') as genes:
    for record in genes:
        line = record.long_name
        discriptor.append(str(line))

amyloid_seq_no_amil_disc = {}
for key, val in filter_seq.items():
    for disc in discriptor:
        if key in disc:
            amyloid_seq_no_amil_disc[disc] = val

amil_seq = {}
fasta_uniprot = {}
for key, val in amyloid_seq_no_amil_disc.items():
    head_for_fasta = str(key)
    for key_amil, val_amil in amyloid_motifs_list_for_header.items():
        if key_amil in key:
            for i in val_amil:
                if ' amyloid_seq=' not in head_for_fasta:
                    head_for_fasta += ' amyloid_seq=' + i
                else:
                    head_for_fasta += '|' + i
    amil_seq[head_for_fasta] = val

table = []
with open('uniprot-proteome.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    count = 0
    for rows in file:
        if count == 0:
            table.append(rows)
            count += 1
        for code in new_csv:
            if code in rows:
                table.append(rows)

table_short = []
t = 1
for line in table:
    entry = line.index('Entry')
    entry_name = line.index('Entry Name')
    reviewed = line.index('Reviewed')
    protein_names = line.index('Protein names')
    gene_names = line.index('Gene Names')
    organism = line.index('Organism')
    length = line.index('Length')
    break

for elem in table:
    table_short_words = []
    table_short_words.append(elem[entry])
    table_short_words.append(elem[entry_name])
    table_short_words.append(elem[reviewed])
    table_short_words.append(elem[protein_names])
    table_short_words.append(elem[gene_names])
    table_short_words.append(elem[organism])
    table_short_words.append(elem[length])
    table_short.append(table_short_words)

with open('Amyloid_Proteome', 'w') as am_seq:
    for key, val in amil_seq.items():
        am_seq.write(str(key) + '\n')
        am_seq.write(str(val) + '\n')


with open('table_prot_amyl.tsv', 'w') as tab:
    for elem in table_short:
        for words in elem:
            tab.write(words + '\t')
        tab.write('\n')
