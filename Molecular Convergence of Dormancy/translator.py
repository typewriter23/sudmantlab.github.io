def transcribe(sequence):
    rna_seq = sequence.replace('T', 'U')
    return(rna_seq)

def translate_rna(sequence):
    codon2aa = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N", 
                "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T", 
                "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S", 
                "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I", 

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H", 
                "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P", 
                "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R", 
                "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L", 

                "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D", 
                "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A", 
                "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G", 
                "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V", 

                "UAA":"_", "UAC":"Y", "UAG":"_", "UAU":"T", 
                "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S", 
                "UGA":"_", "UGC":"C", "UGG":"W", "UGU":"C", 
                "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}

    protein_seq = ''
    for n in range(0, len(sequence), 3):
        if sequence[n:n+3] in codon2aa:
            protein_seq += codon2aa[sequence[n:n+3]]
    return protein_seq


ids, sequences = [], []
n = -1
with open('dna.fasta') as fh:
    for line in fh:
        line = line.strip()
        if line[0] == '>':
            id = line.split()[0][1:]
            ids.append(id)
            sequences.append('')
            n += 1
        else:
            sequences[n] += line


for id, seq in zip(ids, sequences):
    print('>'+id)
    rna = transcribe(seq)
    protein = translate_rna(rna)
    print(protein)
    
    