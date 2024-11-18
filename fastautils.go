package fastautils

import (
	"fmt"
	"github.com/evolbioinf/fasta"
	"os"
)

// Function Clean removes non-canonical nucleotides from a Sequence (that is, keeps only ATGC/atgc). The function updates the input sequence in place.
func Clean(s *fasta.Sequence) {
	d := s.Data()
	i := 0
	for _, c := range d {
		if c == 'A' || c == 'C' ||
			c == 'G' || c == 'T' ||
			c == 'a' || c == 'c' ||
			c == 'g' || c == 't' {
			d[i] = c
			i++
		}
	}
	d = d[:i]
	*s = *fasta.NewSequence(s.Header(), d)
}

// ReadAll reads all sequences from a file and returns a slice of Sequences.
func ReadAll(f *os.File) []*fasta.Sequence {
	sc := fasta.NewScanner(f)
	var s []*fasta.Sequence
	for sc.ScanSequence() {
		s = append(s, sc.Sequence())
	}
	f.Close()
	return s
}

// Concatenate accepts a slice of Sequences and a sentinel byte. It concatenates the slice into a single Sequence entry, where all headers and data are glued together. The concatenated headers and pieces of data are separated with the sentinel byte, if the latter is not zero.
func Concatenate(seqSlice []*fasta.Sequence,
	sentinel byte) (*fasta.Sequence, error) {
	var err error
	l := len(seqSlice)
	switch {
	case l > 1:
		h := []byte(seqSlice[0].Header())
		d := seqSlice[0].Data()
		for i := 1; i < l; i++ {
			if sentinel != 0 {
				h = append(h, sentinel)
				d = append(d, sentinel)
			}
			h = append(h, []byte(seqSlice[i].Header())...)
			d = append(d, seqSlice[i].Data()...)
		}
		cSeq := fasta.NewSequence(string(h), d)
		return cSeq, err
	case l == 1:
		return seqSlice[0], err
	default:
		err = fmt.Errorf("fastautils.Concatenate: " +
			"the input slice is empty\n")
		return nil, err
	}
}
