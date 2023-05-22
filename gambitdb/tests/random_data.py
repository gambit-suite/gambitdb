import random
num_sequences = 10  # Number of sequences in the test set

with open('test_sequences.fasta', 'w') as file:
    for i in range(num_sequences):
        sequence = 'ATGAC' + ''.join(random.choice('ACGT') for _ in range(6)) + '\n'
        file.write(f'>sequence{i+1}\n{sequence}')
