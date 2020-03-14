package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strings"
)

var p = fmt.Println

type fasta struct {
	id   string
	desc string
	seq  string
}

func build_fasta(header string, seq bytes.Buffer) (record fasta) {
	fields := strings.SplitN(header, " ", 2)

	if len(fields) > 1 {
		record.id = fields[0]
		record.desc = fields[1]
	} else {
		record.id = fields[0]
		record.desc = ""
	}

	record.seq = seq.String()

	return record
}

func parse(fastaFh io.Reader) chan fasta {

	outputChannel := make(chan fasta)

	scanner := bufio.NewScanner(fastaFh)
	// scanner.Split(bufio.ScanLines)
	header := ""
	var seq bytes.Buffer

	go func() {
		// Loop over the letters in inputString
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if len(line) == 0 {
				continue
			}

			// line := scanner.Text()

			if line[0] == '>' {
				// If we stored a previous identifier, get the DNA string and map to the
				// identifier and clear the string
				if header != "" {
					// outputChannel <- build_fasta(header, seq.String())
					outputChannel <- build_fasta(header, seq)
					// fmt.Println(record.id, len(record.seq))
					header = ""
					seq.Reset()
				}

				// Standard FASTA identifiers look like: ">id desc"
				header = line[1:]
			} else {
				// Append here since multi-line DNA strings are possible
				seq.WriteString(line)
			}

		}

		outputChannel <- build_fasta(header, seq)

		// Close the output channel, so anything that loops over it
		// will know that it is finished.
		close(outputChannel)
	}()

	return outputChannel
}

func codoncount(s, substr string) map[string]map[int][]int {
	var cys_map = map[string]map[int][]int{
		"forward": map[int][]int{},
	}
	cys_map["forward"][0] = make([]int, 0)
	cys_map["forward"][1] = make([]int, 0)
	cys_map["forward"][2] = make([]int, 0)
	// fmt.Println(s)
	original := s
	// n := 0
	nuc := -3
	for {
		i := strings.Index(s, substr)
		if i == -1 {
			return cys_map
			// return n
		}
		// n++
		nuc = nuc + (i + len(substr))
		if nuc+3 != len(original) {
			if original[nuc+3:nuc+3+1] == original[nuc+1:nuc+1+1] {
				// p("Whoopsy", nuc, original[nuc+3:nuc+3+1], original[nuc:nuc+1], original[nuc+1:nuc+1+1], original[nuc:nuc+3], original[nuc+1:nuc+4])
				// p("Whoopsy", nuc, original[nuc:nuc+3], original[nuc+2:nuc+3+2])
				cys_map["forward"][(nuc+2)%3] = append(cys_map["forward"][(nuc+2)%3], nuc+2)
			}
		}
		s = s[i+len(substr):]

		// p(nuc, nuc%3, original[nuc:nuc+3], original)
		cys_map["forward"][nuc%3] = append(cys_map["forward"][nuc%3], nuc)

		// if reading frame, add nuc (original position) to corresponding data structure
		// data structure should likely be a dict of the following structure:
		// cys_positions = {"forward": {0:[n1, n+1, ... nx]},1:[...],2:[...]}
		// reverse will have to be dealt with slightly differently
	}
}

func main() {

	fastaFh, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFh.Close()

	for record := range parse(fastaFh) {
		record_TGT := codoncount(record.seq, "TGT")
		record_TGC := codoncount(record.seq, "TGC")

		record_cys := record_TGT
		record_cys["forward"][0] = append(record_cys["forward"][0], record_TGC["forward"][0]...)
		sort.Ints(record_cys["forward"][0])
		record_cys["forward"][1] = append(record_cys["forward"][1], record_TGC["forward"][1]...)
		sort.Ints(record_cys["forward"][1])
		record_cys["forward"][2] = append(record_cys["forward"][2], record_TGC["forward"][2]...)
		sort.Ints(record_cys["forward"][2])

		p(record.id, len(record.seq),
			len(record_TGT["forward"][0]),
			len(record_TGT["forward"][1]),
			len(record_TGT["forward"][2]),
			len(record_TGC["forward"][0]),
			len(record_TGC["forward"][1]),
			len(record_TGC["forward"][2]))
		// p(record_cys)
		// record_forward_0 := append(record_TGT["forward"][0], record_TGC["forward"][0]...)
		// record_forward_1 := append(record_TGT["forward"][1], record_TGC["forward"][1]...)
		// record_forward_2 := append(record_TGT["forward"][2], record_TGC["forward"][2]...)
		// p(record_forward_0, record_forward_1, record_forward_2)
		// p(record_TGT["forward"][2])
	}
}
