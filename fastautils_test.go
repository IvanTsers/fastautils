package fastautils

import (
	"github.com/evolbioinf/fasta"
	"os"
	"strconv"
	"testing"
)

func TestClean(t *testing.T) {
	s := "XXAtaTNGTnCactAploenTTg"
	w := "ATATGTCACTATTG"
	seq := fasta.NewSequence("", []byte(s))
	Clean(seq)
	g := string(seq.Data())
	if g != w {
		t.Errorf("Clean() want:\n%s\nget:\n%s\n", w, g)
	}
}
func TestReadAll(t *testing.T) {
	expectedLen := [9]int{0, 0, 0, 5, 70, 140, 700, 1000, 1000}
	for i := 1; i < 9; i++ {
		name := "./data/seq" + strconv.Itoa(i) +
			".fasta"
		f, err := os.Open(name)
		if err != nil {
			t.Errorf("couldn't open %q\n", name)
		}
		seqSlice := ReadAll(f)
		w := expectedLen[i-1]
		for entry, seq := range seqSlice {
			g := len(seq.Data())
			if g != w {
				t.Errorf("seq%d, entry %d - want: %d\nget: %d",
					i, entry, w, g)
			}
		}
	}
}
func TestConcatenate(t *testing.T) {
	wantDataLen := [9]int{0, 0, 0, 5, 2*70 + 1,
		2*140 + 1, 5*700 + 4, 5*1000 + 4}
	wantHeaders := [9]string{"",
		"",
		"seq3",
		"Rand_1; G/C=0.20",
		"Rand_1; G/C=0.41|Rand_2; G/C=0.54",
		"Rand_1; G/C=0.50|Rand_2; G/C=0.50",
		"Rand_1; G/C=0.50|Rand_2; G/C=0.50|Rand_3; " +
			"G/C=0.50|Rand_4; G/C=0.50|Rand_5; G/C=0.50",
		"Rand_1; G/C=0.50|Rand_2; G/C=0.50|Rand_3; " +
			"G/C=0.50|Rand_4; G/C=0.50|Rand_5; G/C=0.50",
		"Rand_1; G/C=0.50|Rand_2; G/C=0.50|Rand_3; " +
			"G/C=0.50|Rand_4; G/C=0.50|Rand_5; G/C=0.50"}
	for i := 1; i < 9; i++ {
		name := "./data/seq" + strconv.Itoa(i) +
			".fasta"
		f, err := os.Open(name)
		if err != nil {
			t.Errorf("couldn't open %q\n", name)
		}
		seqSlice := ReadAll(f)
		seq, err := Concatenate(seqSlice, '|')
		if err != nil {
			if i == 1 {
				continue
			}
		}
		wl := wantDataLen[i-1]
		gl := len(seq.Data())
		if gl != wl {
			t.Errorf("%s data:\nget:\n%d\nwant:\n%d\n", name, gl, wl)
		}
		wh := wantHeaders[i-1]
		gh := seq.Header()
		if gh != wh {
			t.Errorf("%s headers:\nget:\n%s\nwant:\n%s\n", name, gh, wh)
		}
	}
}
func TestFindByName(t *testing.T) {
	s := "AACACACACAC"
	seq1 := fasta.NewSequence("seq1_neighbor", []byte(s))
	s = "TGTGTGTGTG"
	seq2 := fasta.NewSequence("seq2_target", []byte(s))
	sequences := []*fasta.Sequence{seq1, seq2}

	get := FindByHeader(sequences, "target")
	want := seq2
	if get[0] != want {
		t.Errorf("get:\n%v\nwant:\n%v\n", get, want)
	}
}
