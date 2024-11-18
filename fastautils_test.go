package fastautils

import (
	"github.com/evolbioinf/fasta"
	"testing"
)

func TestClean(t *testing.T) {
	s := "XXATATNGTnCactAploenTTg"
	w := "ATATGTCactATTg"
	seq := fasta.NewSequence("", []byte(s))
	Clean(seq)
	g := string(seq.Data())
	if g != w {
		t.Errorf("Clean() want:\n%s\nget:\n%s\n", w, g)
	}
}
