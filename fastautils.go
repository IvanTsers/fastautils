package fastautils

import (
	"fmt"
	"github.com/evolbioinf/fasta"
	"os"
	"regexp"
)

// Function Clean removes non-canonical nucleotides from a Sequence (that is, keeps only ATGC). The function updates the input sequence in place.
func Clean(s *fasta.Sequence) {
	d := s.Data()
	i := 0
	for _, c := range d {
		u := c &^ 0x20
		if isACGT[u] {
			d[i] = u
			i++
		}
	}
	d = d[:i]
	s.SetData(d)
}

var isACGT [256]bool

func init() {
	for _, c := range []byte("ACGT") {
		isACGT[c] = true
	}
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

// AddReverseComplement appends the reverse complement sequence to a fasta entry under the same header. The strands are separated with a hash (\#).
func AddReverseComplement(s *fasta.Sequence) {
	d := s.Data()
	var newD []byte
	rev := fasta.NewSequence("reverse", d)
	rev.ReverseComplement()
	newD = append(d, '#')
	newD = append(newD, rev.Data()...)
	s.SetData(newD)
}

// The function FindByHeader accepts a slice of sequences and returns a sub-slice of sequences, headers of which match the specified pattern.
func FindByHeader(ss []*fasta.Sequence, p string) []*fasta.Sequence {
	var res []*fasta.Sequence
	r := regexp.MustCompile(p)
	for _, s := range ss {
		h := s.Header()
		if r.MatchString(h) {
			res = append(res, s)
		}
	}
	return res
}
