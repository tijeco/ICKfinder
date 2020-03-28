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

// the below function can be modified to allow for partial matches
// just need to also figure out a way to modify the pattern string so that I can keep up with which matches are modified
// I also need to think about how parial matches of one pattern can cause it to look like another pattern
// I don't want a pattern to be searched or at least reported twice
// i.e [6,5,0,4,8] vs [6,5,0,3,8]
// perhaps it isn't possible to only search one pattern once
// if I store the results in a map, that will prevent them being reported Multiple times

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
	var aaSeq strings.Builder
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

	// var str strings.Builder

	// for i := 0; i < 1000; i++ {
	//     str.WriteString("a")
	// }

	// fmt.Println(str.String())

	return aaSeq.String()

}

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
				var start_found, stop_found bool

				toxin_matches := nearMatches(toxin, difList(record_cys["forward"][reading_frame]))
				if len(toxin_matches) > 0 {
					for _, match_pair := range toxin_matches {

						left := record_cys["forward"][reading_frame][match_pair[0]]
						right := record_cys["forward"][reading_frame][match_pair[1]]

						// p(record.id, "+", reading_frame, left, right, pattern, record.seq[left:right+3])
						// check for ORF below
						if len(record.seq[:left]) >= 150 {
							// p(record.seq[149:left])
							// p("left:", len(record.seq[149:left]), record.seq[149:left])
							// p(codoncount(record.seq[149:left], "ATG"))
							// p(len(codoncount(record.seq[149:left], "ATG")[0]))
							if len(codoncount(record.seq[149:left], "ATG")[0]) > 0 {
								// p("Start found within 50 AA")
								start_found = true
							}
						} else {
							// p("frame", reading_frame, len(record.seq[reading_frame:left]))
							if len(codoncount(record.seq[reading_frame:left], "ATG")[0]) > 0 {
								// p("Start found between here and beginning")
								start_found = true
							}
						} // end check start
						// p("Stop:", len(record.seq[right:]), record.seq[right:])
						if len(record.seq[right:]) >= 150 {
							// p(record.seq[right : right+150])
							// p(len(record.seq[right : right+150]))
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
						p(record.id, "+", reading_frame, left, right, start_found && stop_found, pattern, record.seq[left:right+3], translate(record.seq[left:right+3]))
						// p("Start:", start_found, "Stop:", stop_found)
						// p("ORF:", start_found && stop_found)
					} // end for match_pair

				} // end forward

				toxin = reverse_pairs(toxin)
				toxin_matches = nearMatches(toxin, difList(record_cys["reverse"][reading_frame]))
				// getMatches(toxin, difList(record_cys["reverse"][reading_frame]))
				if len(toxin_matches) > 0 {
					for _, match_pair := range toxin_matches {

						left := record_cys["reverse"][reading_frame][match_pair[0]]
						right := record_cys["reverse"][reading_frame][match_pair[1]]
						// p(record.id, "-", reading_frame, left, right, pattern, record.seq[left:right+3])
						// check for ORFs below
						// p("left:", left, "right:", right, len(record.seq))
						// p(record.seq[right+3:])
						if len(record.seq[right+3:]) > 150 {
							right_seq = record.seq[right+3 : right+150]
							// p(len(right_seq), right_seq)
							// p(codoncount(right_seq, "CAT"))
							if len(codoncount(right_seq, "CAT")[0]) > 0 {
								start_found = true
							}

						} else {
							right_seq = record.seq[right+3:]
							if len(codoncount(right_seq, "CAT")[0]) > 0 {
								start_found = true
							}
						}
						// p("Left:", left, record.seq[:left])
						if len(record.seq[:left]) < 150 {
							left_seq = record.seq[reading_frame:left]
							// p(len(left_seq), left_seq)
							// p(codoncount(left_seq, "TTA"))
							// p(codoncount(left_seq, "CTA"))
							// p(codoncount(left_seq, "TCA"))
							if len(codoncount(left_seq, "TTA")[0]) > 0 || len(codoncount(left_seq, "CTA")[0]) > 0 || len(codoncount(left_seq, "TCA")[0]) > 0 {
								stop_found = true
							}
						} else {
							left_seq = record.seq[left-150 : left-3]
							// p("Long:", len(left_seq), left_seq)
							// p(codoncount(left_seq, "TTA"))
							// p(codoncount(left_seq, "CTA"))
							// p(codoncount(left_seq, "TCA"))
							if len(codoncount(left_seq, "TTA")[0]) > 0 || len(codoncount(left_seq, "CTA")[0]) > 0 || len(codoncount(left_seq, "TCA")[0]) > 0 {
								stop_found = true
							}

						}
						p(record.id, "-", reading_frame, left, right, start_found && stop_found, pattern, record.seq[left:right+3], translate(reverseCompliment(record.seq[left:right+3])))
						// p("Start:", start_found, "Stop:", stop_found)
						// p("ORF:", start_found && stop_found)

					}

				} // end reverse

			} // end toxin pattern for loop

		} // end reading_frame for loop

	} // end fasta iterater

	// thingy := translate("ATNTAG")
	// p(thingy)

	// new_jsonFile, err := os.Open("clucktig_general.json")
	// // if we os.Open returns an error then handle it
	// if err != nil {
	// 	fmt.Println(err)
	// }
	// // fmt.Println("Successfully Opened clucktig.json")
	// // defer the closing of our new_jsonFile so that we can parse it later on
	// defer new_jsonFile.Close()
	//
	// new_byteValue, _ := ioutil.ReadAll(new_jsonFile)
	//
	// var new_patterns map[string][][2]int
	// json.Unmarshal([]byte(new_byteValue), &new_patterns)
	// p(new_patterns)
	// p(len(new_patterns["C-C-CC-C-C"]))
	// p(new_patterns["C-C-CC-C-C"][0])
	// p(new_patterns["C-C-CC-C-C"][0][0])
	// p(new_patterns["C-C-CC-C-C"][0][1])
	// test_num := 7
	// p(test_num >= new_patterns["C-C-CC-C-C"][0][0])
	// p(test_num <= new_patterns["C-C-CC-C-C"][0][1])
	// p(test_num >= new_patterns["C-C-CC-C-C"][0][0] && test_num <= new_patterns["C-C-CC-C-C"][0][1])
	// q_toxin := []int{6, 6, 1, 4, 6}
	// p(closeEnough(q_toxin, new_patterns["C-C-CC-C-C"]))
}
