# Find lines that contain “five_prime_UTR”
with open('/home/qurat-ul-ain/Desktop/Homo_sapiens.GRCh38.109.chromosome.21.gff3', 'r') as f:
 for line in f:
 if 'five_prime_UTR' in line:
 print(line)

# Find RNAs that contain more than one 5’UTR
filename = "/home/qurat-ul-ain/Desktop/Homo_sapiens.GRCh38.109.chromosome.21.gff3"
utr_counts = {}
with open(filename, 'r') as f:
 for line in f:
 if not line.startswith('#'):
 fields = line.split('\t')
 if fields[2] == 'five_prime_UTR':
 attributes = fields[8].split(';')
 for attr in attributes:
 if attr.startswith('Parent='):
 rna_id = attr.split('=')[1]
 if rna_id not in utr_counts:
 utr_counts[rna_id] = 1
 else:
 utr_counts[rna_id] += 1
for rna_id, count in utr_counts.items():
 if count > 1:
 print(f"RNA {rna_id} has {count} 5'UTRs")

# Find 5’ UTR sizes
filename = "/home/qurat-ul-ain/Desktop/Homo_sapiens.GRCh38.109.chromosome.21.gff3"
utr_sizes = {}
with open(filename, 'r') as f:
 for line in f:
 if not line.startswith('#'):
 fields = line.split('\t')
 if fields[2] == 'five_prime_UTR':
 attributes = fields[8].split(';')
 for attr in attributes:
 if attr.startswith('Parent='):
 rna_id = attr.split('=')[1]
 if rna_id not in utr_sizes:
 utr_sizes[rna_id] = []
 utr_sizes[rna_id].append(int(fields[4]) - int(fields[3]) + 1)
for rna_id, sizes in utr_sizes.items():
 if len(sizes) > 1:
 print(f"RNA {rna_id} has {len(sizes)} 5'UTRs of sizes: {', '.join(map(str, sizes))}")

# Find introns and exons sizes
filename = "/home/qurat-ul-ain/Desktop/Homo_sapiens.GRCh38.109.chromosome.21.gff3"
rna_features = {}
with open(filename, 'r') as f:
 for line in f:
 if not line.startswith('#'):
 fields = line.strip().split('\t')
 feature_type = fields[2]
 if feature_type in ['exon', 'five_prime_UTR']:
 attributes = fields[8].split(';')
 rna_id = None
 for attr in attributes:
 if attr.startswith('Parent='):
 rna_id = attr.split('=')[1]
 break
 if rna_id is not None:
 if rna_id not in rna_features:
 rna_features[rna_id] = {'exons': [], 'five_prime_utrs': []}
 if feature_type == 'exon':
 rna_features[rna_id]['exons'].append((int(fields[3]), int(fields[4])))
 elif feature_type == 'five_prime_UTR':
 rna_features[rna_id]['five_prime_utrs'].append((int(fields[3]), int(fields[4])))
for rna_id, features in rna_features.items():
 if len(features['five_prime_utrs']) > 1:
 exons = sorted(features['exons'], key=lambda x: x[0])
 utrs = sorted(features['five_prime_utrs'], key=lambda x: x[0])
 introns = []
 for i in range(1, len(exons)):
 start = exons[i-1][1] + 1
 end = exons[i][0] - 1
 introns.append((start, end))
 exon_sizes = [exons[i][1] - exons[i][0] + 1 for i in range(len(exons))]
 intron_sizes = [introns[i][1] - introns[i][0] + 1 for i in range(len(introns))]
 print(f"RNA {rna_id} has {len(exon_sizes)} exons of sizes: {', '.join(map(str, exon_sizes))}")
 print(f"RNA {rna_id} has {len(intron_sizes)} introns of sizes: {', '.join(map(str, 
intron_sizes))}")
 print("\n \n")
