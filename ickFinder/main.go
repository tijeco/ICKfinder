package main

import (
	"fmt"
)

var p = fmt.Println

func runCommand(command string) (willRun bool) {
	return
}

func blastp(query, target string, evalue float64) (blastout string) {
	return
}

func hmmer(query, target string) (hmmerout string) {
	return
}

func combineBlastpHmmer(blastp, hmmer string, cys, maxAA int) (combinedOut string) {
	return
}

func header2seq(headers, sequences string) (seqOut string) {
	return
}

func signalP(pepSeq string) (signalpOut string) {
	return
}

func main() {
	var inPep, verifiedICK string
	blastResults := blastp(inPep, verifiedICK, 1e-3) // add err
	hmmerResults := hmmer(inPep, verifiedICK)        // add err
	bothResults := combineBlastpHmmer(blastResults, hmmerResults, 5, 200)
	blastHmmerSeqs := header2seq(bothResults, inPep)

	signalpResults := signalP(blastHmmerSeqs)

	p(signalpResults)

	// fmt.Println("Hello, playground")
}
