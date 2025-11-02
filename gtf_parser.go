package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"log"
)

// parseGTFToLengths takes input a .gtf.gz file
// It calculate the valid length ( sum of exons) for every gene
// It returns map[gene_id_with_version] length (/Kilobases)
func parseGTFToLengths(gtfPath string) (map[string]float64, error) {
	file, err := os.Open(gtfPath)
	if err != nil {
		return nil, fmt.Errorf("cannot open GTF file %s: %w", gtfPath, err)
	}
	defer file.Close()

	gz, err := gzip.NewReader(file)
	if err != nil {
		return nil, fmt.Errorf("unable to create GTF reader: %w", err)
	}
	defer gz.Close()

	// GTF is **Tab-separated**
	bufferedGzipReader := bufio.NewReader(gz)
	reader := csv.NewReader(bufferedGzipReader)
	reader.Comma = '\t'
	reader.Comment = '#' // Lines starting with # are annotations
	reader.LazyQuotes = true

	// We will use a map to store the sum of exons for every gene
	// map[gene_id_version] -> total_base_pairs
	geneBasePairs := make(map[string]int)

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			// GTF may be not regular and formal, skip
			log.Printf("Warning: GTF parsing error,skip line: %v", err)
			continue
		}

		if len(record) < 9 {
			continue // incomplete line.
		}

		// We only focus **exon** length
		featureType := record[2]
		if featureType != "exon" {
			continue
		}

		// calculate
		start, err1 := strconv.Atoi(record[3])
		end, err2 := strconv.Atoi(record[4])
		if err1 != nil || err2 != nil {
			continue 
		}
		length := (end - start) + 1 
		// Bioinformatics coordinates are 1-based and closed.

		// Parsing column 9 (attributes)
		attributes, err := parseAttributes(record[8])
		if err != nil {
			continue
		}

		// Find "gene_id". 
		// The GTEx GCT file uses IDs with versions (such as ENSG... .15).
		geneID, ok := attributes["gene_id"]
		if !ok {
			continue // exon has no gene_id, discard.
		}

		geneBasePairs[geneID] += length
	}

	// transform "base pairs" (int) into "kilobases" (float64)
	geneKilobases := make(map[string]float64, len(geneBasePairs))
	for geneID, bp := range geneBasePairs {
		if bp > 0 {
			geneKilobases[geneID] = float64(bp) / 1000.0
		}
	}

	if len(geneKilobases) == 0 {
		return nil, errors.New("?")
	}

	return geneKilobases, nil
}

// parseAttributes
// gene_id "ENSG..."; transcript_id "ENST..." ; ... 
// Output: map["gene_id"] -> "ENSG..."
func parseAttributes(attrString string) (map[string]string, error) {
	attrs := make(map[string]string)
	
	// GTF arributes separate by "; " (its a ";" + " ")
	fields := strings.Split(strings.TrimSuffix(attrString, ";"), "; ")

	for _, field := range fields {
		parts := strings.SplitN(field, " ", 2)
		if len(parts) != 2 {
			continue
		}
		key := parts[0]
		// remove the ""
		value := strings.Trim(parts[1], "\"")
		attrs[key] = value
	}
	return attrs, nil
}