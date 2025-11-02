package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
)

// processGCTFile is to actually process the raw data.
// It filters and normalize the dataset with subroutines after.
func processGCTFile(
	gctPath string,
	geneLengthsKB map[string]float64,
	lowExprThreshold float64,
	lowVarPercentile float64,
) (
	finalMatrix [][]float64,
	finalGeneList []string,
	finalSampleList []string,
	err error,
) {

	// Pass 1: calculate "Per-Sample RPK Sum" (Used as the denominator of TPM)
	log.Println("  (GCT Pass 1/2) Calculating the TPM normalized factor...")
	perSampleRPKSum, sampleList, numSamples, err := gctPass1_CalculateRPKSums(gctPath, geneLengthsKB)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("GCT Pass 1 失败: %w", err)
	}
	log.Printf("  (GCT Pass 1/2) ...finished。 %d samples in the file。", numSamples)

	// Pass 2: Calculate TPM, perform Log2 conversion, and conduct two rounds of filtering

	log.Println("  (GCT Pass 2/2) Calculate TPM, perform Log2 conversion, and conduct two rounds of filtering...")
	
	finalMatrix, finalGeneList, err = gctPass2_FilterAndNormalize(
		gctPath,
		geneLengthsKB,
		perSampleRPKSum,
		numSamples,
		lowExprThreshold,
		lowVarPercentile,
	)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("GCT Pass 2 failed: %w", err)
	}
	
	log.Printf("  (GCT Pass 2/2) ...done")

	return finalMatrix, finalGeneList, sampleList, nil
}


// gctPass1_CalculateRPKSums realizes the first round of streaming read operation
func gctPass1_CalculateRPKSums(gctPath string, geneLengthsKB map[string]float64) (
	perSampleRPKSum []float64,
	sampleList []string,
	numSamples int,
	err error,
) {
	file, gz, reader, err := openGCTReader(gctPath)
	if err != nil {
		return nil, nil, 0, err
	}
	defer file.Close()
	defer gz.Close()

	// GCT File Format Processing

	// 1. Skip the first line (version "#1.2")
	if _, err := reader.Read(); err != nil {
		return nil, nil, 0, fmt.Errorf("1st line in GCT failed: %w", err)
	}
	// 2. Skip the second line (dimension "55584 578")
	if _, err := reader.Read(); err != nil {
		return nil, nil, 0, fmt.Errorf("2nd line in GCT failed: %w", err)
	}
	// 3. Read the 3rd line (headers)
	header, err := reader.Read()
	if err != nil {
		return nil, nil, 0, fmt.Errorf("3rd line in GCT failed: %w", err)
	}
	// GCT header format: [Name] [Description] [Sample1] [Sample2] ...
	if len(header) < 3 {
		return nil, nil, 0, errors.New("invalid GCT header format")
	}
	sampleList = header[2:]
	numSamples = len(sampleList)
	perSampleRPKSum = make([]float64, numSamples)

	// 4. Starting from the fourth line, process the data line by line.
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("Warning: Pass 1 skips a line of GCT: %v", err)
			continue
		}

		// GCT format: [gene_id_version] [gene_symbol] [count1] [count2] ...
		geneIDWithVersion := record[0] 
		
		lengthKB, ok := geneLengthsKB[geneIDWithVersion]
		if !ok || lengthKB == 0 {
			continue 
			// Gene length not found, it did not contribute to the total RPK value.
		}

		// Count the number of occurrences of this gene in all the samples
		for i := 0; i < numSamples; i++ {
			colIndex := i + 2 
			// +2 because the first two columns are "gene_id" and "description".
			count, err := strconv.ParseFloat(record[colIndex], 64)
			if err != nil {
				continue 
			}
			
			// RPK = Reads / Kilobase
			rpk := count / lengthKB 
			perSampleRPKSum[i] += rpk
		}
	}
	return perSampleRPKSum, sampleList, numSamples, nil
}


// gctPass2_FilterAndNormalize realizes second round of streaming read
// and two rounds of filtering
func gctPass2_FilterAndNormalize(
	gctPath string,
	geneLengthsKB map[string]float64,
	perSampleRPKSum []float64,
	numSamples int,
	lowExprThreshold float64,
	lowVarPercentile float64,
) ([][]float64, []string, error) {

	// These two slices are used to temporarily store the genes that have passed the "low expression" filter
	// We use the gene symbol (record[1]) as the human-readable ID
	var intermediateGenes []string 
	var intermediateData [][]float64

	file, gz, reader, err := openGCTReader(gctPath)
	if err != nil {
		return nil, nil, err
	}
	defer file.Close()
	defer gz.Close()

	// skip the first 3 lines of GCT
	_, _ = reader.Read()
	_, _ = reader.Read()
	_, _ = reader.Read()

	// Filter 1: filter the low expressions
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("Warning : Pass 2 skip a line of GCT: %v", err)
			continue
		}

		geneIDWithVersion := record[0]
		geneSymbol := record[1]
		
		lengthKB, ok := geneLengthsKB[geneIDWithVersion]
		if !ok || lengthKB == 0 {
			continue 
		}

		log2Values := make([]float64, numSamples)
		lowExprCount := 0

		for i := 0; i < numSamples; i++ {
			colIndex := i + 2
			count, _ := strconv.ParseFloat(record[colIndex], 64)
			
			// 1. RPK (Reads Per Kilobase)
			rpk := count / lengthKB

			// 2. TPM
			tpm := 0.0
			if perSampleRPKSum[i] > 0 {
				tpm = (rpk / perSampleRPKSum[i]) * 1_000_000
			}

			// 3. Log2 transformation
			log2Val := math.Log2(tpm + 1)
			log2Values[i] = log2Val

			// 4. Check low expression (log2(TPM+1) < 1)
			if log2Val < 1.0 { 
				lowExprCount++
			}
		}

		if float64(lowExprCount)/float64(numSamples) >= lowExprThreshold {
			continue // delete the low expression gene
		}

		// This gene has passed. Save it for variance filtering.
		// We save geneSymbol (e.g., "TP53") instead of "ENSG..."
		intermediateGenes = append(intermediateGenes, geneSymbol)
		intermediateData = append(intermediateData, log2Values)
	}

	if len(intermediateGenes) == 0 {
		return nil, nil, errors.New("no gene left after filtering low expression")
	}
	log.Printf("  (GCT Pass 2/2) ... %d genes passed the expression filtering。", len(intermediateGenes))



	// Filter 2: Low Variability Filtering

	log.Println("  (GCT Pass 2/2) Calculating variance and performing low variation filtering...")

	// 1. Calculate the variance of all genes
	type geneVar struct {
		geneSymbol string
		data   []float64
		v      float64
	}
	
	geneVariances := make([]geneVar, len(intermediateGenes))
	for i := 0; i < len(intermediateGenes); i++ {
		v := variance(intermediateData[i])
		geneVariances[i] = geneVar{
			geneSymbol: intermediateGenes[i],
			data:   intermediateData[i],
			v:      v,
		}
	}

	// 2. Sort by variance (in ascending order)
	sort.Slice(geneVariances, func(i, j int) bool {
		return geneVariances[i].v < geneVariances[j].v 
	})

	// 3. Find the 25th percentile index
	cutoffIndex := int(float64(len(geneVariances)) * lowVarPercentile)

	// 4. Construct the final matrix 
	// (retain only the genes with a frequency of > 25%)
	finalMatrix := make([][]float64, 0, len(geneVariances)-cutoffIndex)
	finalGeneList := make([]string, 0, len(geneVariances)-cutoffIndex)
	
	for i := cutoffIndex; i < len(geneVariances); i++ {
		finalGeneList = append(finalGeneList, geneVariances[i].geneSymbol)
		finalMatrix = append(finalMatrix, geneVariances[i].data)
	}

	return finalMatrix, finalGeneList, nil
}


// openGCTReader is used to open a gzip file.
func openGCTReader(gctPath string) (*os.File, *gzip.Reader, *csv.Reader, error) {
	file, err := os.Open(gctPath)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("failed to open the GCT file %s: %w", gctPath, err)
	}

	gz, err := gzip.NewReader(file)
	if err != nil {
		file.Close()
		return nil, nil, nil, fmt.Errorf("unable to create GCT gzip reader: %w", err)
	}


	bufferedGzipReader := bufio.NewReader(gz)
	// GCT is tab-separated

	reader := csv.NewReader(bufferedGzipReader) 
	reader.Comma = '\t'
	reader.LazyQuotes = true // GCT file may be not formal.
	reader.FieldsPerRecord = -1

	return file, gz, reader, nil
}