package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"os/exec"
	"sort"
	"strconv"
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

func difList(d []int) []int {
	dif := []int{}
	for i := range d {
		if i+1 != len(d) {
			dif = append(dif, ((d[i+1]-d[i])/3)-1)
		} else {
			dif = append(dif, -1)
		}
	}
	return dif
}

func codoncount(s, substr string) map[int][]int {
	var cys_map = map[int][]int{}

	cys_map[0] = make([]int, 0)
	cys_map[1] = make([]int, 0)
	cys_map[2] = make([]int, 0)

	original := s
	// n := 0
	original_pos := -3

	for {
		i := strings.Index(s, substr)
		if i == -1 {
			return cys_map
			// return n
		}
		// n++
		original_pos = original_pos + (i + len(substr))
		// p("Now at:", original_pos, "Remaining:", len(original)-(original_pos+1))

		// if original_pos+3 != len(original) {
		if len(original)-(original_pos+1) > 3 {
			// p("Whoopsy", original_pos+2, original[original_pos+3:original_pos+3+1], original[original_pos+1:original_pos+1+1], original[original_pos:original_pos+3], original[original_pos:original_pos+3+2], original[original_pos+2:original_pos+3+2])
			if original[original_pos+3:original_pos+3+1] == original[original_pos+1:original_pos+1+1] {
				// p("Whoopsy", original_pos, original[original_pos+3:original_pos+3+1], original[original_pos:original_pos+1], original[original_pos+1:original_pos+1+1], original[original_pos:original_pos+3], original[original_pos+1:original_pos+4])
				// p("Whoopsy", original_pos, original[original_pos:original_pos+3], original[original_pos+2:original_pos+3+2])
				if original[original_pos:original_pos+3] == original[original_pos+2:original_pos+3+2] {
					// p("Whoopsy", original_pos+2, original[original_pos+3:original_pos+3+1], original[original_pos+1:original_pos+1+1], original[original_pos:original_pos+3], original[original_pos:original_pos+3+2], original[original_pos+2:original_pos+3+2])
					cys_map[(original_pos+2)%3] = append(cys_map[(original_pos+2)%3], original_pos+2)

				}
				// p("Whoopsy", original_pos+2, original[original_pos+3:original_pos+3+1], original[original_pos+1:original_pos+1+1], original[original_pos:original_pos+3], original[original_pos:original_pos+3+2], original[original_pos+2:original_pos+3+2])
				// cys_map[(original_pos+2)%3] = append(cys_map[(original_pos+2)%3], original_pos+2)
			}
		}

		// p(original_pos, original_pos%3, original[original_pos:original_pos+3], original)
		cys_map[original_pos%3] = append(cys_map[original_pos%3], original_pos)
		s = s[i+len(substr):]

		// if reading frame, add original_pos (original position) to corresponding data structure
		// data structure should likely be a dict of the following structure:
		// cys_positions = {"forward": {0:[n1, n+1, ... nx]},1:[...],2:[...]}
		// reverse will have to be dealt with slightly differently
	} // end for
}

func closeEnough(a []int, b [][2]int) bool {
	if len(a) != len(b) {
		return false
	}
	for pos, query := range a {
		if query < b[pos][0] || query > b[pos][1] {
			return false
		}
	}
	return true
}

func Equal(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i, v := range a {
		if v != b[i] { // instead of not equal, it can be <= b[i] + given_number, if b[i] > 0
			return false
		}
	}
	return true
}

func nearMatches(input [][2]int, query []int) [][2]int {
	// query := []int{6, 6, 1, 0, 5, 1, 6, 1, 11, 5, 10, 3, 3, -1}
	out_pairs := [][2]int{}
	// p("New out:", out_pairs, len(out_pairs))

	// p(Equal(input, query))
	// p(query)
	for i, value := range query {
		// p(len(query), query[i:], input[0], value, input[0] == value, len(input), len(query[i:]))
		// if input[0] == value {
		if value >= input[0][0] && value <= input[0][1] {
			// p(i, len(input), len(query[i:]))
			if len(query[i:]) >= len(input) {
				// p(i, len(input))
				// p("query[", i, ":", i, "+", len(input), "]")
				// p(query[i : i+len(input)])
				// p(Equal(input, query[i:i+len(input)]))
				if closeEnough(query[i:i+len(input)], input) {
					// p("Match at:", i, i+len(input))
					out_pairs = append(out_pairs, [2]int{i, i + len(input)})
				}
			}
		}
	} // end for
	return out_pairs
}

func getMatches(input, query []int) [][2]int {
	// query := []int{6, 6, 1, 0, 5, 1, 6, 1, 11, 5, 10, 3, 3, -1}
	out_pairs := [][2]int{}
	// p("New out:", out_pairs, len(out_pairs))

	// p(Equal(input, query))
	// p(query)
	for i, value := range query {
		// p(len(query), query[i:], input[0], value, input[0] == value, len(input), len(query[i:]))
		if input[0] == value {
			// p(i, len(input), len(query[i:]))
			if len(query[i:]) >= len(input) {
				// p(i, len(input))
				// p("query[", i, ":", i, "+", len(input), "]")
				// p(query[i : i+len(input)])
				// p(Equal(input, query[i:i+len(input)]))
				if Equal(input, query[i:i+len(input)]) {
					// p("Match at:", i, i+len(input))
					out_pairs = append(out_pairs, [2]int{i, i + len(input)})
				}
			}
		}
	} // end for
	return out_pairs
}

func reverseCompliment(codon string) string {
	// codon := "TGC"
	// p(codon, "A"[0])

	n := 0
	runes := make([]rune, len(codon))
	for _, r := range codon {

		// nuc := "G"
		switch string(r) {
		case "A":
			// fmt.Println("Switch from A to T")
			// r := []rune("ABC€")
			// fmt.Println(r) // [65 66 67 8364]
			runes[n] = []rune("T")[0]
		case "T":
			// fmt.Println("Switch from T to A")
			runes[n] = []rune("A")[0]
		case "C":
			// fmt.Println("Switch from C to G")
			runes[n] = []rune("G")[0]
		case "G":
			// fmt.Println("Switch from G to C")
			runes[n] = []rune("C")[0]
		default:
			// fmt.Println("Invalid nucleotide.")
			runes[n] = r

		}
		// p(i, r, (string(r)))

		n++
	}
	runes = runes[0:n]
	// p(runes)
	// Reverse
	for i := 0; i < n/2; i++ {
		runes[i], runes[n-1-i] = runes[n-1-i], runes[i]
		// p(string(runes[i]), i, "and", runes[n-1-i], "will now be", runes[n-1-i], "and", runes[i])
	}
	// Convert back to UTF-8.
	output := string(runes)
	// fmt.Println(output)
	return output
}

func reverse_pairs(numbers [][2]int) [][2]int {
	newNumbers := make([][2]int, len(numbers))
	for i, j := 0, len(numbers)-1; i <= j; i, j = i+1, j-1 {
		newNumbers[i], newNumbers[j] = numbers[j], numbers[i]
	}
	return newNumbers
}

func reverse(numbers []int) []int {
	newNumbers := make([]int, len(numbers))
	for i, j := 0, len(numbers)-1; i <= j; i, j = i+1, j-1 {
		newNumbers[i], newNumbers[j] = numbers[j], numbers[i]
	}
	return newNumbers
}

func translate(s string) string {
	// n := map[string]string{"foo": 1, "bar": 2}
	// if len(s%3) != 0 {
	// 	return "", error
	// }
	var aaSeq bytes.Buffer
	codon_map := map[string]string{
		"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
		"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
		"AAT": "N", "AAC": "N",
		"GAT": "D", "GAC": "D",
		"TGT": "C", "TGC": "C",
		"CAA": "Q", "CAG": "Q",
		"GAA": "E", "GAG": "E",
		"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
		"CAT": "H", "CAC": "H",
		"ATT": "I", "ATC": "I", "ATA": "I",
		"TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
		"AAA": "K", "AAG": "K",
		"ATG": "M",
		"TTT": "F", "TTC": "F",
		"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
		"TCA": "S", "TCT": "S", "TCC": "S", "TCG": "S", "AGT": "S", "AGC": "S",
		"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
		"TAA": "*", "TGA": "*", "TAG": "*",
		"TGG": "W",
		"TAT": "Y", "TAC": "Y",
		"GTA": "V", "GTT": "V", "GTC": "V", "GTG": "V"}

	// p(codon_map)
	strLength := len(s)
	splitedLength := int(math.Ceil(float64(strLength) / float64(3)))
	var start, stop int
	for i := 0; i < splitedLength; i += 1 {
		start = i * 3
		stop = start + 3
		if stop > strLength {
			stop = strLength
		}
		current_codon := s[start:stop]
		// p(current_codon)
		if aa, found := codon_map[current_codon]; found {
			aaSeq.WriteString(aa)
		} else {
			aaSeq.WriteString("X")
		}
		// aaSeq.WriteString(codon_map[current_codon])
	}

	return aaSeq.String()

}

func general2specific(general []int) string {
	var specific bytes.Buffer
	for _, num := range general {
		specific.WriteString("-C-" + strconv.Itoa(num))
	}
	return specific.String()[1:] + "-C"
}

func getIntrons(s, strand string) [][2]int {
	// var s_pos int
	var coordinates [][2]int
	// var donor, acceptor, intron string
	var donor, acceptor string
	if strand == "+" {
		donor = "GT"
		acceptor = "AG"
	} else if strand == "-" {
		donor = "CT"
		acceptor = "AC"
	}
	s_original := s
	len_s := len(s)
	// p(1, s)
	for {
		donor_pos := strings.Index(s, donor)
		orig_donor := (len_s - len(s)) + donor_pos
		if donor_pos != -1 {
			s = s[donor_pos:]
			// p(2, s)
			for {
				acceptor_pos := strings.Index(s, acceptor)
				// orig_acceptor := (len_s - len(s)) + acceptor_pos + 2
				orig_acceptor := (len_s - len(s)) + acceptor_pos + 1
				if acceptor_pos != -1 {
					coordinates = append(coordinates, [2]int{orig_donor, orig_acceptor})
					s = s[acceptor_pos+1:]
					// p(3, s)
				} else {
					s = s_original[orig_donor+2:]
					// p(4, s)
					break
				}
			}
		} else {
			return coordinates
		}
	}

} // end getIntrons

func get_start_exons(seq, strand string, introns [][2]int) map[[3]int]bool {
	var out_exons map[[3]int]bool
	out_exons = make(map[[3]int]bool)
	var possible_starts map[int][]int
	for intron := range introns {
		possible_exon := seq[:introns[intron][0]]
		if strand == "+" {
			possible_starts = codoncount(possible_exon, "ATG")
			// p("possible_starts", possible_starts)
		} else if strand == "-" {
			possible_starts = codoncount(possible_exon, "CAT")
		}
		if len(possible_starts) > 0 {
			for rf := range possible_starts {
				for start := range possible_starts[rf] {
					out_exons[[3]int{possible_starts[rf][start], introns[intron][0], introns[intron][1]}] = true
					// out_exons[[3]int{possible_starts[rf][start], introns[intron][0], introns[intron][1]}] = true
				}
			}
		}

	}
	return out_exons
}

func get_stop_exons(introns [][2]int, n int) map[[3]int]bool {
	var out_exons map[[3]int]bool
	out_exons = make(map[[3]int]bool)

	for intron := range introns {
		// out_exons = append(out_exons, [2]int{introns[intron][1] + 1, n - 1})
		out_exons[[3]int{introns[intron][0], introns[intron][1], n - 1}] = true
	}
	return out_exons
}

// func get_internal_exons(introns [][2]int) map[[4]int]bool {
// 	var out_exons map[[4]int]bool
// 	out_exons = make(map[[4]int]bool)
// 	// p(1, introns)
// 	for intron := range introns {
// 		// p(2, intron)
// 		exon_start := introns[intron][1]
// 		// p(3, exon_start)
// 		for next_intron := range introns[intron+1:] {
// 			// p(introns[intron], introns[intron+1:][next_intron])
// 			// p(4, next_intron)
// 			exon_stop := introns[intron+1:][next_intron][0]
// 			// p(5, exon_stop)
// 			// p(6, exon_stop > exon_start)
// 			if exon_stop > exon_start {
// 				intron_1 := introns[intron]
// 				intron_2 := introns[intron+1:][next_intron]
// 				out_exons[[4]int{intron_1[0], intron_1[1], intron_2[0], intron_2[1]}] = true
//
// 			}
// 		}
// 	}
// 	return out_exons
// }

func valid_orf(orf []int) bool {
	var is_valid bool
	var len_orf int
	if len(orf)%2 == 0 {
		if len(orf) == 2 {
			len_orf = (orf[1] + 1) - orf[0]
			// no introns
		} else {
			// first is start, final is stop, all internal pairs are introns
			// p("orf:", orf)
			// for i, value := range orf {
			for i := 0; i <= len(orf)-2; i += 2 {
				value := orf[i]
				// if value < len(orf)-1 {
				// if i < len(orf)-1 {
				if i == 0 || i == len(orf)-2 { // for start and stop
					len_orf += orf[i+1] - value
					// p(i, orf[i+1], value, len_orf)

				} else { // for all internal intron pairs
					len_orf += (orf[i+1] - 1) - value
					// p(i, (orf[i+1] - 1), value, len_orf)

				}

			}
			// add up length of all exons and set is_valid True if len%3 == 0
		}
		if len_orf%3 == 0 {
			is_valid = true
		}

		// no introns
	}
	return is_valid
}

func get_cds(starts, stops map[[3]int]bool, introns [][2]int, seq string) [][]int {
	var out_cds [][]int
	var start_cds, orf []int
	for start := range starts {
		for stop := range stops {
			if Equal(start[:][1:], stop[:][:len(stop)-1]) { // one intron
				// p(1, start, stop)
				orf = append(start[:], stop[:][len(stop)-1])
				// p(orf, orfSeq(seq, [][]int{orf}), valid_orf(orf))
				if valid_orf(orf) {
					out_cds = append(out_cds, orf)

				}
			} else if start[0] > stop[len(stop)-1] || start[1] > stop[0] || start[len(start)-1] > stop[0] {
				// p(2, start, stop)
				continue
			} else {
				// p(3, start, stop)
				if valid_orf(append(start[:], stop[:]...)) {
					out_cds = append(out_cds, append(start[:], stop[:]...))
				}
				for i, intron := range introns {

					if Equal(intron[:], start[:][1:]) || Equal(intron[:], stop[:][:len(stop)-1]) {
						// p(1, start, intron, stop)
						continue
					} else if intron[0] < start[len(start)-1] || intron[1] > stop[0] {
						// p(2, start, intron, stop)
						continue
					} else {
						// p(3, start, intron, stop)
						start_cds = append(start[:], intron[:]...)
						orf = append(start_cds, stop[:]...)
						// p(orf, orfSeq(seq, [][]int{orf}), valid_orf(orf))
						if valid_orf(orf) {
							out_cds = append(out_cds, orf)
						}
						for _, next_intron := range introns[i:] {
							// p(intron, next_intron)
							if Equal(next_intron[:], start_cds[:][len(start_cds)-2:]) || Equal(next_intron[:], stop[:][:len(stop)-1]) {
								continue
							} else if next_intron[0] < start_cds[len(start_cds)-1] || next_intron[1] > stop[0] {
								continue
							} else {
								start_cds = append(start_cds, next_intron[:]...)
								orf = append(start_cds, stop[:]...)
								// p(orf, orfSeq(seq, [][]int{orf}), valid_orf(orf))
								if valid_orf(orf) {
									out_cds = append(out_cds, orf)
								}
							}

						}
					}
				}
			}

		}
	}
	// 		if start[0] > stop[len(stop)-1] {
	// 			continue
	// 		} else if start[len(start)-1] > stop[0] {
	// 			continue
	// 		} else if Equal(start[1:][:], stop[:len(stop)-1][:]) { // one intron
	// 			start_cds = append(start[:], stop[len(stop)-1])
	// 			if valid_orf(start_cds) {
	// 				out_cds = append(out_cds, start_cds)
	// 				p(1, start, stop, start_cds, start_cds[:len(start_cds)-3])
	//
	// 			}
	// 			start_cds = start_cds[:len(start_cds)-3]
	// 			// if sum of length of the two exons are divisible by 3, return
	// 		} else {
	// 			if ((start[1]-start[0])+(stop[0]-start[len(start)-1]+1)+((stop[len(stop)-1]+1)-(stop[1]+1)))%3 == 0 { // two introns
	//
	// 			} else {
	// 				for intron := range introns {
	// 					if Equal(introns[intron][:], start[1:][:]) || Equal(introns[intron][:], stop[:len(stop)-1][:]) { // intron is already included in start or stop
	// 						if Equal(introns[intron][:], start[1:][:]) { // intron included in start
	// 							start_cds = append(start[:], stop[1:][:]...)
	// 						} else { // intron included in stop
	// 							start_cds = append(start[:len(start)-1][:], stop[:]...)
	// 						}
	// 					} else if introns[intron][0] > start[len(start)-1] {
	// 						start_cds = append(start[:], introns[intron][:]...)
	// 						start_cds = append(start_cds, stop[:]...)
	// 						if valid_orf(start_cds) {
	// 							out_cds = append(out_cds, start_cds)
	// 							p(2, start, stop, start_cds, start_cds[:len(start_cds)-3])
	//
	// 						}
	// 						start_cds = start_cds[:len(start_cds)-3]
	// 						for i := range introns[intron:] {
	// 							start_cds = append(start_cds, introns[intron:][i][:]...)
	// 							start_cds = append(start_cds, stop[:]...)
	// 							if valid_orf(start_cds) {
	// 								out_cds = append(out_cds, start_cds)
	// 								p(3, start, stop, start_cds, start_cds[:len(start_cds)-3])
	//
	// 							}
	// 							start_cds = start_cds[:len(start_cds)-3]
	//
	// 							// check if putative exon_length is divisible by 3
	// 							// afterwards, start_cds = append(start_cds,i[:]...)
	// 						}
	// 					}
	// 				}
	// 			}
	//
	// 		}
	//
	// 	}
	// }
	return out_cds
}

// func get_cds(starts, stops map[[2]int]bool, inners map[[4]int]bool) [][]int {
// 	var out_cds [][]int
// 	// p("starts", starts)
// 	for start := range starts {
// 		// p("start", start)
// 		start_cds := []int{start[0], start[1]}
// 		// p(0, "start_cds", start_cds)
// 		// start_cds = append(start_cds, start...)
// 		// p("stops", stops)
// 		for stop := range stops {
// 			// p("stop", stop)
// 			if start[1] == stop[0]-1 {
// 				if (stop[0]-start[0])%3 == 0 {
// 					out_cds = append(out_cds, []int{start[0], start[1], stop[0], stop[1]})
// 					// p(1, "out_cds", out_cds)
// 					// out_cds = append(out_cds, []int{stop[0],stop[1]})
// 					// return out_cds
// 					// break
//
// 					// return out_cds
// 				}
// 			} else if start[1] > stop[0] {
// 				// p(1, "evaluate", start[1], ">", stop[0])
// 				continue
// 			} else {
// 				// p("inners", inners)
// 				for inner := range inners {
// 					// p("inner", inner)
// 					// p(2, "evaluate", inner[0], start_cds[len(start_cds)-1]+1)
// 					if inner[0] == start_cds[len(start_cds)-1]+1 { // does inner plug into start
// 						// p(3, "evaluate", inner[len(start_cds)-1], stop[0]-1)
// 						if inner[len(inner)-1] == stop[0]-1 { // does inner plug into stop
//
// 							start_cds = append(start_cds, inner[:]...)
// 							len_exon := (start_cds[1] + 1) - start_cds[0]
// 							// p(1, "start_cds", start_cds, start_cds[2:], len_exon)
//
// 							for i := range start_cds[2:] {
// 								// p("i", i)
// 								// p(4, "evaluate", i+4, len(start_cds))
// 								if i+4 <= len(start_cds[2:]) {
// 									len_exon += (start_cds[2:][i+2]) - (start_cds[2:][i+1] + 1)
// 									// p("len_exon", start_cds[2:][i+2], start_cds[2:][i+1], (start_cds[2:][i+2]+1)-start_cds[2:][i+1], len_exon)
// 									i += 4
// 								}
// 							}
// 							// p(5, "evaluate", len_exon, stop[1], stop[0], (len_exon + (stop[1] - stop[0])))
// 							if (len_exon+((stop[1]+1)-stop[0]))%3 == 0 {
// 								start_cds = append(start_cds, stop[:]...)
// 								out_cds = append(out_cds, start_cds[:])
// 								// out_cds = append(out_cds, stop[:])
// 								// return out_cds
// 							}
// 						}
// 					} else { // might need elif here --> if start_cds[-2] == inner[0] or something like that
//
// 						start_cds = append(start_cds, inner[:]...)
// 						// p(2, "start_cds", start_cds)
// 					}
//
// 				} // end for inner
// 			}
//
// 		} // end for stop
// 	} // end for start
// 	return out_cds
// } // end get_cds

func get_orf(seq, strand string) [][]int {
	introns := getIntrons(seq, strand)
	p("introns", introns)
	start_exons := get_start_exons(seq, strand, introns)
	p("start_exons", start_exons)
	stop_exons := get_stop_exons(introns, len(seq))
	p("stop_exons", stop_exons)
	// internal_exons := get_internal_exons(introns)
	// p("internal_exons", internal_exons)

	cds := get_cds(start_exons, stop_exons, introns, seq)
	return cds
}
func orfSeq(seq string, orfs [][]int) []string {
	var seqList []string
	// p("seq", seq)
	for _, orf := range orfs {
		// p(coordinates)
		if len(orf)%2 == 0 {
			if len(orf) == 2 {
				// len_orf = (orf[1] + 1) - orf[0]
				seqList = append(seqList, seq[orf[0]:orf[1]+1])
				// no introns
			} else {
				// first is start, final is stop, all internal pairs are introns
				seq2write := ""

				// p("orf", orf)
				for i := 0; i <= len(orf)-2; i += 2 {
					value := orf[i]
					// for i, value := range orf {
					// if value < len(orf)-1 {
					if i == 0 || i == len(orf)-2 { // for start and stop
						// len_orf += orf[i+1] - value
						if i == 0 { // start
							// p(i, value, orf[i+1], seq[value:orf[i+1]])
							seq2write += seq[value:orf[i+1]]
						} else { // stop
							// p(i, value+1, orf[i+1], seq[value+1:orf[i+1]+1])
							seq2write += seq[value+1 : orf[i+1]+1]
						}

					} else { // for all internal intron pairs
						// len_orf += (orf[i+1] - 1) - value
						// p(i, value+1, orf[i+1], seq[value+1:orf[i+1]])
						seq2write += seq[value+1 : orf[i+1]]
					}

				}
				seqList = append(seqList, seq2write)
				// add up length of all exons and set is_valid True if len%3 == 0
			}
			// if len_orf%3 == 0 {
			// 	is_valid = true
			// }

			// no introns
		}
	}

	return seqList
}

// func possibleExons(introns [][2]int , s string) {
//
// }

// this goves the coordinates of all non-overlapping introns
//  I think what I need instead is also overlapping??
// meaning nXnGTnYnGTnZnAG could yield nXnGTnYnGTnZnAG or nYnGTnZnAG
// such that nYnGTnZn is ignored
// because it seems unlikely that zero GTs fall within an intronic sequence

func main() {
	var right_seq, left_seq string

	p(right_seq, left_seq)
	// toxin := []int{6, 6, 0, 4, 6}
	// toxin := []int{6, 6, 1, 0, 5, 1, 6, 1, 11, 5, 10, 3, 3}
	// p(toxin)
	// Open our jsonFile
	// jsonFile, err := os.Open("clucktig_C6_augment.json")
	jsonFile, err := os.Open("clucktig_general.json")
	// if we os.Open returns an error then handle it
	if err != nil {
		fmt.Println(err)
	}
	// fmt.Println("Successfully Opened clucktig.json")
	// defer the closing of our jsonFile so that we can parse it later on
	defer jsonFile.Close()

	byteValue, _ := ioutil.ReadAll(jsonFile)

	var patterns map[string][][2]int
	var list_of_cys []int
	json.Unmarshal([]byte(byteValue), &patterns)

	fastaFh, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFh.Close()

	for record := range parse(fastaFh) {
		record_TGT := codoncount(record.seq, "TGT")
		record_TGC := codoncount(record.seq, "TGC")
		record_ACA := codoncount(record.seq, "ACA")
		record_GCA := codoncount(record.seq, "GCA")

		var record_cys = map[string]map[int][]int{
			"forward": map[int][]int{},
			"reverse": map[int][]int{},
		}

		for reading_frame := 0; reading_frame <= 2; reading_frame++ {
			// record_cys["forward"][reading_frame] = make([]int, 0)
			record_cys["forward"][reading_frame] = append(record_TGT[reading_frame], record_TGC[reading_frame]...)
			// record_cys["reverse"][reading_frame] = make([]int, 0)
			record_cys["reverse"][reading_frame] = append(record_ACA[reading_frame], record_GCA[reading_frame]...)
			// p(record_cys)
			sort.Ints(record_cys["forward"][reading_frame])
			sort.Ints(record_cys["reverse"][reading_frame])

			for pattern, toxin := range patterns {
				var stop_found bool

				list_of_cys = difList(record_cys["forward"][reading_frame])
				toxin_matches := nearMatches(toxin, list_of_cys)
				if len(toxin_matches) > 0 {
					for _, match_pair := range toxin_matches {

						left := record_cys["forward"][reading_frame][match_pair[0]]
						right := record_cys["forward"][reading_frame][match_pair[1]]

						// p(record.id, "+", reading_frame, left, right, pattern, record.seq[left:right+3])

						// check for ORF below
						// if len(record.seq[:left]) >= 150 {
						// 	// p(record.seq[149:left])
						// 	// p("left:", len(record.seq[149:left]), record.seq[149:left])
						// 	// p(codoncount(record.seq[149:left], "ATG"))
						// 	// p(len(codoncount(record.seq[149:left], "ATG")[0]))
						// 	if len(codoncount(record.seq[149:left], "ATG")[0]) > 0 {
						// 		// p("Start found within 50 AA")
						// 		start_found = true
						// 	}
						// } else {
						// 	// p("frame", reading_frame, len(record.seq[reading_frame:left]))
						// 	if len(codoncount(record.seq[reading_frame:left], "ATG")[0]) > 0 {
						// 		// p("Start found between here and beginning")
						// 		start_found = true
						// 	}
						// } // end check start

						// p("Stop:", len(record.seq[right:]), record.seq[right:])
						if len(record.seq[right:]) >= 150 {
							if len(codoncount(record.seq[right+3:right+150], "TAA")[0]) > 0 || len(codoncount(record.seq[right+3:right+150], "TAG")[0]) > 0 || len(codoncount(record.seq[right+3:right+150], "TGA")[0]) > 0 {
								stop_found = true
								// p(record.seq[right+3 : right+150])
								// p("Stop found!")

							}
						} else {
							// p(len(record.seq[right:len(record.seq)]), record.seq[right:len(record.seq)])
							// p(len(record.seq[right+3:len(record.seq)]), record.seq[right+3:len(record.seq)])
							//
							// p(codoncount(record.seq[right+3:len(record.seq)], "TAA"))
							// p(codoncount(record.seq[right+3:len(record.seq)], "TAG"))
							// p(codoncount(record.seq[right+3:len(record.seq)], "TGA"))
							// p("Stop_found:", len(codoncount(record.seq[right+3:len(record.seq)], "TAA")) > 0 || len(codoncount(record.seq[right+3:len(record.seq)], "TAG")) > 0 || len(codoncount(record.seq[right+3:len(record.seq)], "TGA")) > 0)
							// right_seq := record.seq[right+3 : len(record.seq)]
							if len(codoncount(record.seq[right+3:len(record.seq)], "TAA")[0]) > 0 || len(codoncount(record.seq[right+3:len(record.seq)], "TAG")[0]) > 0 || len(codoncount(record.seq[right+3:len(record.seq)], "TGA")[0]) > 0 {
								// p("Stop found!!!")
								stop_found = true

							}

						} // end check stop
						internal_stops := strings.Count(translate(record.seq[left:right+3]), "*")
						p(record.id, "+", reading_frame, left, right, stop_found, internal_stops, pattern, general2specific(list_of_cys[match_pair[0]:match_pair[1]]), record.seq[left:right+3], translate(record.seq[left:right+3]))
						// p(match_pair, general2specific(list_of_cys[match_pair[0]:match_pair[1]]))
						// p("Start:", start_found, "Stop:", stop_found)
						// p("ORF:", start_found && stop_found)
					} // end for match_pair

				} // end forward

				toxin = reverse_pairs(toxin)
				list_of_cys = difList(record_cys["reverse"][reading_frame])
				toxin_matches = nearMatches(toxin, list_of_cys)
				// getMatches(toxin, difList(record_cys["reverse"][reading_frame]))
				if len(toxin_matches) > 0 {
					for _, match_pair := range toxin_matches {

						left := record_cys["reverse"][reading_frame][match_pair[0]]
						right := record_cys["reverse"][reading_frame][match_pair[1]]
						// p(record.id, "-", reading_frame, left, right, pattern, record.seq[left:right+3])

						// check for ORFs below
						// p("left:", left, "right:", right, len(record.seq))
						// p(record.seq[right+3:])
						// if len(record.seq[right+3:]) > 150 {
						// 	right_seq = record.seq[right+3 : right+150]
						// 	// p(len(right_seq), right_seq)
						// 	// p(codoncount(right_seq, "CAT"))
						// 	if len(codoncount(right_seq, "CAT")[0]) > 0 {
						// 		start_found = true
						// 	}
						//
						// } else {
						// 	right_seq = record.seq[right+3:]
						// 	if len(codoncount(right_seq, "CAT")[0]) > 0 {
						// 		start_found = true
						// 	}
						// }

						if len(record.seq[:left]) < 150 {
							left_seq = record.seq[reading_frame:left]

							if len(codoncount(left_seq, "TTA")[0]) > 0 || len(codoncount(left_seq, "CTA")[0]) > 0 || len(codoncount(left_seq, "TCA")[0]) > 0 {
								stop_found = true
							}
						} else {
							left_seq = record.seq[left-150 : left-3]

							if len(codoncount(left_seq, "TTA")[0]) > 0 || len(codoncount(left_seq, "CTA")[0]) > 0 || len(codoncount(left_seq, "TCA")[0]) > 0 {
								stop_found = true
							}

						}
						p(record.id, "-", reading_frame, left, right, stop_found, strings.Count(translate(reverseCompliment(record.seq[left:right+3])), "*"), pattern, general2specific(reverse(list_of_cys[match_pair[0]:match_pair[1]])), record.seq[left:right+3], translate(reverseCompliment(record.seq[left:right+3])))
						// p(record.id, "-", reading_frame, left, right, start_found && stop_found, strings.Count(translate(reverseCompliment(record.seq[left:right+3])), "*"), pattern, general2specific(reverse(list_of_cys[match_pair[0]:match_pair[1]])), record.seq[left:right+3], translate(reverseCompliment(record.seq[left:right+3])))

					}

				} // end reverse

			} // end toxin pattern for loop

		} // end reading_frame for loop

		// this is the end of the contig should now be able to evaluate data structure
		// if the data structure has enough ICKs past a certain threshold
		// look for gene models for each ICK

		//  once the gene models have been evaluated
		// I suppose they should just be written to a file
		// the next step would be to align them with mafft
		// I don't know if that can be done in GO...
		// IT CANN!!!!! https://golang.org/pkg/os/exec/

	} // end fasta iterater
	// p(getIntrons("aGTbccaaqww100000000000mnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmmnmnmnmnmnmnmnmmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmncdAGefgaGTbccaaqww100000000000000mmmmmmmmmmmmmmmmmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmmnmnmnmnmnmnmnmmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmnmncdAGefg", "+"))
	// test_seq0 := "AGnn1nGTnnn2nnGTnn3nnnAGnn4nnAG"
	// test_introns0 := getIntrons(test_seq0, "+")
	// p(test_introns0)
	// p(test_seq0[test_introns0[0][0]:test_introns0[0][1]])
	// p(1, get_stop_exons(test_introns0, len(test_seq0)))
	// p(2, get_stop_exons([][2]int{{2, 3}, {6, 7}}, 10))
	// test_seq := "NATG012GTNNNNNNNNNNAG...GTNNNNNNNNAG345TAA"
	// test_introns := getIntrons(test_seq, "+")
	// test_stops := get_stop_exons(test_introns, len(test_seq))
	//
	// test_internals := get_internal_exons(test_introns)
	// p(get_start_exons(test_seq, "+", getIntrons(test_seq, "+")))
	// p(3, test_introns)
	// p("get_internal_exons")
	// p(get_internal_exons(test_introns))

	// p(get_cds([][2]int{{0, 1}}, [][2]int{{8, 9}, {4, 9}}, [][4]int{{2, 3, 6, 7}}))
	// test_starts := get_start_exons(test_seq, "+", test_introns)
	// p("test_introns")
	// p(test_introns)
	// p("test_internals")
	// p(test_internals)
	// p("test_stops")
	// p(test_stops) // // NOTE: this reports duplicates some times
	// p("test_starts")
	// p(test_starts) // // NOTE: this reports duplicates some times

	// will need to make sure it doesn't exist in the list before it gets reported
	// perhaps worth considering if it should be put in a map??, the search  time shouldn't be rediculous though

	// p(get_orf(test_seq, "+"))
	// test_seq_1 := "ATGTAA"
	// p(get_orf(test_seq_1, "+"))
	// test_seqs := []string{
	// 	"ATGTAA",                                                                                   // no
	// 	"ATG...TAA",                                                                                // no
	// 	"ATGnnnNNNTAA",                                                                             // no
	// 	"ATGnnnGTxxxxxxxxxxxxxxxxxAGNNNTAA",                                                        // no, 1 intron
	// 	"ATGnnnGTxxxxxxxxxxxxxxxxxAGnnnGTxxxxxxxxxxxxxxxxxAGNNNTAA",                                // yes 2 introns
	// 	"NATG012GTNNNNNNNNNNAG...GTNNNNNNNNAG345TAA",                                               // yes
	// 	"NATGnnnGTxxxxxxxxxxxxxAGnnnGTxxxxxxxxxxxxAGNNNTAA",                                        // yes
	// 	"NATGnnnGTxxxxxxxxxGTxxxxxxxxxxxxxAGnnnGTxxxxxxxxxxxxAGxxxxAGnnnGTxxxxxxxxxxxxAGNNNTAA",    // no
	// 	"NATGnnnGTxxxxxxxxxxxxxAGnnnGTxxxxxxxxxxxxAGNNNGTxxxxxxxxxxxxxAGnnnGTxxxxxxxxxxxxAGTAA",    // no
	// 	"NATGnnnGTxxxxxxxxxxxxxAGnnnGTxxxxxxxxxxxxAGNNNGTxxxxxxxxxxxxxAGnnnGTxxxxxxxxxxxxAGNNNTAA", // no
	// 	"NATGnnnGTabcdefghijklmnopqrstuvwxyzAGnnnGTabcdefghijklmnopqrstuvwxyzAGNNNTAA"}             // yes
	// p(test_seqs)
	test_seq := "NATGaaaGT0000000000000AGbbbGT111111111111AGcccGT2222222222222AGdddGT333333333333AGeeeTAA"
	p("test_seq", test_seq)
	// test_introns := getIntrons(test_seq, "+")
	// p("test_introns", test_introns)
	// test_starts := get_start_exons(test_seq, "+", test_introns)
	// p("test_starts", test_starts)
	// test_stops := get_stop_exons(test_introns, len(test_seq))
	// p("test_stops", test_stops)
	// p("test_cds", get_cds(test_starts, test_stops, test_introns))
	test_orfs := get_orf(test_seq, "+")
	test_orf_seqs := orfSeq(test_seq, test_orfs)
	p("test_orfs", test_orfs)
	p("test_orf_seqs", test_orf_seqs)

	// p(test_seqs)
	// for _, test := range test_seqs {
	// 	p(test)
	// 	test_orfs := get_orf(test, "+")
	// 	// p(get_orf(test, "+"))
	// 	test_orf_seqs := orfSeq(test, test_orfs)
	// 	p(test_orfs)
	// 	p(test_orf_seqs)
	//
	// }
	// if true {
	// 	p("################ get_cds ################")
	// 	test_cds := get_cds(test_starts, test_stops, test_internals)
	//
	// 	p(test_cds)
	//
	// }

	// p(getIntrons("AGnn1nGTnnn2nnGTnn3nnnAGnn4nnAG", "+")) // <-- doesn't find GTnnn2nnGTnn3nnnAGnn4nnAG or GTnn3nnnAGnn4nnAG
	// [[4,20],[12,20],[4,27],[12,27]]

	// getIntrons("AGnn1nGTnnn2nnGTnn3nnnAGnn4nnAG", "+") <-- breaks
	//

	path, err := exec.LookPath("mafft")
	if err != nil {
		log.Fatal("installing mafft is in your future")
	}
	fmt.Printf("mafft is available at %s\n", path)
	cmd := exec.Command("sh", "-c", "mafft test.fa")
	stdoutStderr, err := cmd.CombinedOutput()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("%s\n", stdoutStderr)
	for record := range parse(bytes.NewReader(stdoutStderr)) {
		p(record.id)
	}

	// test_seq := "A"
	// p(tes)

}

// ATGaaaGT0000000000000AGbbbGT111111111111AGcccGT2222222222222AGdddGT333333333333AGeeeTAA
