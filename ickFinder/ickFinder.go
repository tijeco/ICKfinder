package main

import (
	"bufio"
	"bytes"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"strings"
	"time"
)

var p = fmt.Println

var ickDB, queryPep, outDir, signalPath string
var help bool

func init() {

	flag.StringVar(&ickDB, "db", "", "fasta file with reference ICK")
	flag.StringVar(&queryPep, "pep", "", "fasta file with query peptides")
	flag.StringVar(&outDir, "out", "", "output directory")
	flag.StringVar(&outDir, "o", "", "output directory (shorthand)")
	flag.StringVar(&signalPath, "signalp", "", "full path to signalP binary")
	flag.BoolVar(&help, "help", false, "print usage")
	flag.BoolVar(&help, "h", false, "print usage (shorthand)")
	flag.Parse()

	someEmpty := ickDB == "" || queryPep == "" || outDir == ""
	// p(someEmpty)
	if someEmpty {
		flag.PrintDefaults()
		os.Exit(1)
	} else {
		outerr := os.MkdirAll(outDir, os.ModePerm)
		if outerr != nil {
			log.Fatal(outerr)
		}
		for _, file := range []string{ickDB, queryPep} {
			if fileExists(file) == false {
				log.Fatalln(file, "does not exist")
			}
		}

	}
	// p("checking if signalP has been requested")
	// log.Fatalln(signalPath == "")

}

type fasta struct {
	id   string
	desc string
	seq  string
}

func buildFasta(header string, seq bytes.Buffer) (record fasta) {
	fields := strings.SplitN(header, " ", 2)

	if len(fields) > 1 {
		record.id = fields[0]
		record.desc = fields[1]
	} else {
		record.id = fields[0]
		record.desc = ""
	}

	record.seq = strings.ToUpper(seq.String())

	return record
}

func parseFasta(fastaFh io.Reader) chan fasta {

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
					// outputChannel <- buildFasta(header, seq.String())
					outputChannel <- buildFasta(header, seq)
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

		outputChannel <- buildFasta(header, seq)

		// Close the output channel, so anything that loops over it
		// will know that it is finished.
		close(outputChannel)
	}()

	return outputChannel
}

// try using it to prevent further errors.
func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

func willRun(command string) (ran bool) {
	commandPath, err := exec.LookPath(command)
	p("Checking if " + command + " is installed")
	if err != nil {
		log.Fatal(command + " not found in $PATH, exiting now")
	} else {
		p(command+" is available at:\n", commandPath)
		ran = true
	}
	return
}

func blastp(query, target, outDir string, evalue float64) (blastTargets []string) {
	if willRun("makeblastdb") && willRun("blastp") {
		// makeblastdb -in spiderICK.fasta  -dbtype prot
		makeBlastDB := "makeblastdb -in " + target + " -dbtype prot"
		makeBlastDBcmd := exec.Command("sh", "-c", makeBlastDB)
		_, makeBlastDBerr := makeBlastDBcmd.Output()
		p("running:", makeBlastDB)
		if makeBlastDBerr != nil {
			p(":(")
			log.Fatal(makeBlastDBerr)
		} else {
			// p("WHOOP!", makeBlastDBstdout)
			// -query Stegotoxin.protein.fa -db spiderICK.fasta -outfmt 6
			blastP := "blastp -query " + query + " -db " + target + " -outfmt 6 -max_target_seqs 1 -evalue " + fmt.Sprintf("%f", evalue)
			blastPcmd := exec.Command("sh", "-c", blastP)
			blastPstdout, blastPerr := blastPcmd.Output()
			if blastPerr != nil {
				log.Fatal(blastPerr)
			} else {
				p("running", blastP)
				// p("WHOOP!", blastPstdout)
				blastTargets = targetFromBlast(string(blastPstdout))
				blastOutFile := outDir + "/blastpResults.txt"
				blastOutErr := ioutil.WriteFile(blastOutFile, blastPstdout, 0644)
				if blastOutErr != nil {
					log.Fatal(blastOutErr)
				}
			}

		}

	}
	return
}

func targetFromBlast(blastout string) (targets []string) {
	// p(strings.Split(blastout, "\n"))
	for _, row := range strings.Split(blastout, "\n") {
		columns := strings.Fields(row)
		// p(row)
		// p(len(columns))
		if len(columns) == 12 {
			targets = append(targets, columns[0])
		}
		// targets = append(targets, strings.Fields(row)[0])
	}
	return
}

func mafft(inFasta string) (alignment string) {
	if willRun("mafft") {
		mafftStr := "mafft --localpair --maxiterate 1000 " + inFasta
		mafftCmd := exec.Command("sh", "-c", mafftStr)
		p("running", mafftStr)
		mafftStdout, mafftErr := mafftCmd.Output()
		if mafftErr != nil {
			log.Fatal(mafftErr)
		} else {

			// p("WHOOP!", blastPstdout)
			alignment = string(mafftStdout)
		}
	}
	return
}

func hmmer(query, target, outDir string) (hmmerTargets []string) {
	if willRun("hmmbuild") && willRun("hmmsearch") {
		var msa []byte
		var msaOutErr error
		msaOutFile := outDir + "/ickDB.aln"
		if fileExists(msaOutFile) == false {
			msa = []byte(mafft(target))
			msaOutErr = ioutil.WriteFile(msaOutFile, msa, 0644)
		} else {
			p(msaOutFile, "exists")
		}

		if msaOutErr != nil {
			log.Fatal(msaOutErr)
		} else {
			hmmOutFile := msaOutFile + ".hmm"
			hmmbuild := "hmmbuild " + hmmOutFile + " " + msaOutFile
			hmmbuildCmd := exec.Command("sh", "-c", hmmbuild)
			p("running", hmmbuild)
			_, hmmbuildErr := hmmbuildCmd.Output()
			if hmmbuildErr != nil {
				log.Fatal(hmmbuildErr)
			} else {
				hmmsearch := "hmmsearch --notextw " + hmmOutFile + " " + query
				hmmsearchCmd := exec.Command("sh", "-c", hmmsearch)
				p("running", hmmsearch)
				hmmsearchOut, hmmsearchErr := hmmsearchCmd.Output()
				if hmmsearchErr != nil {
					log.Fatal(hmmsearchErr)
				} else {
					hmmsearchOutFile := outDir + "/hmmerResults.txt"
					hmmsearchOutErr := ioutil.WriteFile(hmmsearchOutFile, hmmsearchOut, 0644)
					if hmmsearchOutErr != nil {
						log.Fatal(hmmsearchOutErr)
					} else {
						hmmsearchOutStr := string(hmmsearchOut)
						hmmsearchRows := strings.Split(hmmsearchOutStr, "\n")[18:]
						// p(hmmsearchRows)
						for _, row := range hmmsearchRows {
							columns := strings.Fields(row)
							if len(columns) == 9 {
								hmmerTargets = append(hmmerTargets, columns[8])
							} else {
								break
							}
							// p(len(strings.Fields(row)), row)
						}

					}
				}

			}
		}

	}
	return
}

func combineBlastpHmmer(blastp, hmmer []string) map[string]bool {
	combinedOut := make(map[string]bool)
	bothLists := append(blastp, hmmer...)
	for _, h := range bothLists {
		// p("header", h)
		combinedOut[h] = true
	}
	// p("BLASTP:", len(blastp))
	// p("HMMER:", len(hmmer))
	return combinedOut
}

func header2seq(bothMap map[string]bool, seqFile string) map[string]string {
	seqOut := make(map[string]string)
	fastaFh, err := os.Open(seqFile)
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFh.Close()

	for record := range parseFasta(fastaFh) {
		if _, ok := bothMap[record.id]; ok {
			seqOut[record.id] = record.seq
			//do something here
		}

	}
	return seqOut
}

func writeSeqMap(seqIn map[string]string, outDir, outName string) string {
	outFile := outDir + "/" + outName + ".fa"
	out, outErr := os.Create(outFile)
	if outErr != nil {
		log.Fatal(outErr)
	}
	for header, sequence := range seqIn {
		_, headerErr := out.WriteString(">" + header + "\n")
		if headerErr != nil {
			out.Close()
			log.Fatal(headerErr)
		} else {
			_, seqErr := out.WriteString(sequence + "\n")
			if seqErr != nil {
				out.Close()
				log.Fatal(seqErr)
			}
		}
	}
	closeErr := out.Close()
	if closeErr != nil {
		log.Fatal(closeErr)
	}
	return outFile
}

func signalP(pepSeq map[string]string, signalPath, outDir string) (signalpOutStr string) {
	if willRun("signalp") {
		pepSeqFile := writeSeqMap(pepSeq, outDir, "preSignalP")
		signalPstr := signalPath + " -gff3 -prefix signalPout -fasta " + pepSeqFile
		signalPCmd := exec.Command("sh", "-c", signalPstr)
		p("running", signalPstr)
		signalPOut, signalPErr := signalPCmd.Output()
		if signalPErr != nil {
			log.Fatal(signalPErr)
		} else {
			signalpOutStr = string(signalPOut)
		}
	}

	return
}

func main() {
	// currentTime := time.Now()
	p(time.Now().Format("15:04:05 Mon Jan-02-2006"))
	// os.Exit(1)

	// var inPep, verifiedICK string

	blastResults := blastp(queryPep, ickDB, outDir, 1e-3) // add err
	// p("blastResults")
	// p(blastResults)
	// p(strings.Count(blastResults, "\n"))
	// blastTarget := targetFromBlast(blastResults)
	// p(blastTarget)

	hmmerResults := hmmer(queryPep, ickDB, outDir) // add err
	// p("hmmerResults")
	// p(hmmerResults)

	// p(len(intersect.Hash(blastResults, hmmerResults)), intersect.Hash(blastResults, hmmerResults))
	var bothResults map[string]bool
	bothResults = combineBlastpHmmer(blastResults, hmmerResults)
	// p("bothResults", bothResults)
	// p(len(bothResults))
	blastHmmerSeqs := header2seq(bothResults, queryPep) // maybe merge combineBlastpHmmer to go here and return map[string][string]

	if signalPath != "" {
		signalpResults := signalP(blastHmmerSeqs, signalPath, outDir)
		p(signalpResults)
		signalPfileName, signalPfileError := os.Open("signalPout.gff3")
		if signalPfileError != nil {
			log.Fatal(signalPfileError)
		}
		signalPfile := csv.NewReader(signalPfileName)
		signalPfile.Comma = '\t'
		signalPfile.Comment = '#'
		p(signalPfile)

	} else {
		p("Skipping signalP")
		writeSeqMap(blastHmmerSeqs, outDir, "noSignalP")
	}

	// p(willRun("mafft"))

	// pepSeqPtr := flag.String("pep", "", "query peptide sequence")
	// dbSeqPtr := flag.String("db", "", "reference peptide sequence")
	// outPtr := flag.String("out", "", "reference peptide sequence")

	// flag.Parse()
	// flag.PrintDefaults()

	// fmt.Println("Hello, playground")

	p(time.Now().Format("15:04:05 Mon Jan-02-2006"))
}
