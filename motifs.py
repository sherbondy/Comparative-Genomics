# Yeast motif-finding code
# Problem 3

import re
from collections import Counter

def find_motifs():
  inter_data = ""
  motifs = Counter()
  conserved_motifs = Counter()

  with open("genes/allinter") as f:
    inter_data = f.read()

  motifs.update([inter_data[i:(i+6)] for i in range(len(inter_data)-6)])

  print "Total number of motifs: {0}".format(len(motifs))
  print "Most frequent motifs:"
  print motifs.most_common(50)
    
  with open("genes/allintercons") as f:
    cons_data = f.read()
    # use finditer
    for match in re.finditer('\s(\*{6})\s', cons_data):
      motif = inter_data[match.start(1):match.end(1)]
      conserved_motifs.update([motif])

    print "Total number of hexamers: {0}".format(sum(conserved_motifs.values()))
    print "Most frequent conserved motifs:"
    return conserved_motifs.most_common(50)

def main():
  print "testing"
  print find_motifs()

if __name__ == "__main__":
  main()