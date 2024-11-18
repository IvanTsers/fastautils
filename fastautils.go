package fastautils

import (
	"github.com/evolbioinf/fasta"
)

// Function Clean removes non-canonical nucleotides from a Sequence (that is, keeps only ATGC/atgc). The function updates the input sequence in place.
func Clean(s *fasta.Sequence) {
	d := s.Data()
	i := 0
	for _, c := range d {
		if c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
			c == 'a' || c == 'c' || c == 'g' || c == 't' {
			d[i] = c
			i++
		}
	}
	d = d[:i]
	*s = *fasta.NewSequence(s.Header(), d)
}
