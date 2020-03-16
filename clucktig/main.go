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

// May need to go back to the drawing board on this one, doesn't seem to work at all
func subset(first, second []int) [2]int {
	set := make(map[int]int)
	pos := make(map[int][]int)
	for i, value := range second {
		set[value] += 1
		pos[value] = append(pos[value], i)
	}
	// p(set)
	// p(pos)
	coordinates := [2]int{0, 0}
	first_found := true
	// p(coordinates, coordinates == [2]int{0, 0})
	// first_total := 0

	for _, value := range first {

		if count, found := set[value]; !found {
			// checks to see if value from first are in second
			// if first is a subset of second, then all values in first should be in second
			return [2]int{0, 0}
		} else if count < 1 {
			// sort of confused on this, but I believe it is if the count is less than 1
			// I believe this may be because in the next one it reduces cound
			return [2]int{0, 0}
		} else {
			if first_found {
				coordinates[0] = pos[value][len(pos[value])-count]
				first_found = false
			}

			// p(pos[value], value, count, len(pos[value])-count, pos[value][len(pos[value])-count], found)
			set[value] = count - 1
			coordinates[1] = len(first) + coordinates[0]
			// p("End:", len(first)-1+coordinates[0])
			// first_total += value * 3

		}

	}
	// p("End:", first_total)
	// p(coordinates)

	return coordinates
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
					cys_map["forward"][(original_pos+2)%3] = append(cys_map["forward"][(original_pos+2)%3], original_pos+2)

				}
				// p("Whoopsy", original_pos+2, original[original_pos+3:original_pos+3+1], original[original_pos+1:original_pos+1+1], original[original_pos:original_pos+3], original[original_pos:original_pos+3+2], original[original_pos+2:original_pos+3+2])
				// cys_map["forward"][(original_pos+2)%3] = append(cys_map["forward"][(original_pos+2)%3], original_pos+2)
			}
		}

		// p(original_pos, original_pos%3, original[original_pos:original_pos+3], original)
		cys_map["forward"][original_pos%3] = append(cys_map["forward"][original_pos%3], original_pos)
		s = s[i+len(substr):]

		// if reading frame, add original_pos (original position) to corresponding data structure
		// data structure should likely be a dict of the following structure:
		// cys_positions = {"forward": {0:[n1, n+1, ... nx]},1:[...],2:[...]}
		// reverse will have to be dealt with slightly differently
	}
}

func Equal(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i, v := range a {
		if v != b[i] {
			return false
		}
	}
	return true
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
func main() {
	toxin := []int{6, 6, 1, 0, 5, 1, 6, 1, 11, 5, 10, 3, 3}
	p(toxin)
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
		p(record.id, getMatches(toxin, difList(record_cys["forward"][2])))
		// p(record.id)
		// p(record.id, len(record.seq),
		// 	len(record_TGT["forward"][0]),
		// 	len(record_TGT["forward"][1]),
		// 	len(record_TGT["forward"][2]),
		// 	len(record_TGC["forward"][0]),
		// 	len(record_TGC["forward"][1]),
		// 	len(record_TGC["forward"][2]),
		// 	subset(toxin, difList(record_cys["forward"][0])),
		// 	subset(toxin, difList(record_cys["forward"][1])),
		// 	subset(toxin, difList(record_cys["forward"][2])))

		// p(record.id, difList(record_cys["forward"][0]),
		// 	difList(record_cys["forward"][1]),
		// 	difList(record_cys["forward"][2]))

		// p(record_cys["forward"][2])

		// if subset(toxin, difList(record_cys["forward"][2])) != [2]int{0, 0} {
		// 	p(toxin)
		// 	p("cys_pos:", record_cys["forward"][2])
		// 	p("cys_dist:", difList(record_cys["forward"][2]))
		// 	left := subset(toxin, difList(record_cys["forward"][2]))[0]
		// 	right := subset(toxin, difList(record_cys["forward"][2]))[1]
		// 	p("left_pos:", left, "right_pos:", right)
		// 	p("left_value:", record_cys["forward"][2][left], "right_value:", record_cys["forward"][2][right])
		// 	// p(record.id, record.seq[record_cys["forward"][2][left]:record_cys["forward"][2][right]+3])
		// }

		// p(record_cys["forward"][2][subset(toxin, difList(record_cys["forward"][2]))[0]])
		// p(record_cys)
		// p(difList(record_cys["forward"][0]))
		// p(difList(record_cys["forward"][1]))
		// p(difList(record_cys["forward"][2]))

		// p(subset(toxin, difList(record_cys["forward"][0])))
		// p(subset(toxin, difList(record_cys["forward"][1])))
		// p(toxin, "\n", difList(record_cys["forward"][2]))
		// p(subset(toxin, difList(record_cys["forward"][2])))
		// record_forward_0 := append(record_TGT["forward"][0], record_TGC["forward"][0]...)
		// record_forward_1 := append(record_TGT["forward"][1], record_TGC["forward"][1]...)
		// record_forward_2 := append(record_TGT["forward"][2], record_TGC["forward"][2]...)
		// p(record_forward_0, record_forward_1, record_forward_2)
		// p(record_TGT["forward"][2])

	}
	// true_thing := []int{6, 6, 1, 0, 5, 1, 6, 1, 11, 5, 10, 3, 3}
	// new_thing := []int{6, 6, 1, 0, 5, 1, 6, 1, 11, 5, 10, 3, 3, -1}
	// p(getMatches(toxin, new_thing))
}
