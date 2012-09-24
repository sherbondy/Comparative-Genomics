# Yeast motif-finding code
# Problem 3

import re
from collections import Counter
from operator import itemgetter

def find_motifs():
  inter_data = ""
  motifs = Counter()
  conserved_motifs = Counter()
  conservation_ratios = {}

  with open("genes/allinter") as f:
    inter_data = f.read()

  motifs.update([inter_data[i:(i+6)] for i in range(len(inter_data)-6)])

  print "Total number of hexamers: {0}".format(len(motifs))
  print "50 most frequent motifs:"
  print motifs.most_common(50)
    
  with open("genes/allintercons") as f:
    cons_data = f.read()
    # use finditer
    for match in re.finditer('\s(\*{6})\s', cons_data):
      motif = inter_data[match.start(1):match.end(1)]
      conserved_motifs.update([motif])

    for motif in conserved_motifs:
      conservation_ratios[motif] = float(conserved_motifs[motif]) / motifs[motif]

    print "Total number of conserved hexamers: {0}".format(
      sum(conserved_motifs.values()))
    print "50 most *frequent* conserved motifs:"
    print conserved_motifs.most_common(50)

    print "Top 50 most *conserved* motifs (conserved / total):"
    print sorted(conservation_ratios.items(), 
                  key=itemgetter(1), reverse=True)[:50]

def main():
  find_motifs()

if __name__ == "__main__":
  main()