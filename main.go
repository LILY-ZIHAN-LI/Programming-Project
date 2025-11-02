package main

import (
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

const (
	// input1: GCT data path
	gctDataFile = "gene_reads_v10_thyroid.gct.gz"

	// input2: GTF annotation path (length of genes, for TPM)
	gtfAnnotationFile = "gencode.v44.annotation.gtf.gz"

	// output: matrix cleaned
	outputMatrixFile = "clean_thyroid_matrix.csv"

	// Filtering parameter
	// We delete a gene if it's expressed in 90% samples (based on log2(TPM+1) < 1)
	lowExpressionThreshold = 0.9

	// delete genes of low variance percentile
	lowVariancePercentile = 0.25
)

func main() {
	log.Println("Phase 1: Preprocessing the data (parsing & filtering)")
	// step 1: parsing GTF annotations (we need gene length for TPM)
	log.Println("Parsing GTF annotation...")
	// We need to perform **streaming parsing** of GTF because it becomes extremely large after decompression.
	// With parseGTFtoLengths, we get: 
	// map[gene_id_with_version] -> length_in_kilobases
	geneLengthsKB, err := parseGTFToLengths(gtfAnnotationFile)
	if err != nil {
		log.Fatalf("Failed: %v", err)
	}
	log.Printf("...succeed in parsing %d genes length\n", len(geneLengthsKB))

	
	// step 2: preprocessing GCT raw main counts
	
	log.Println("Preprocessing GCT raw counts")
	
	// With processGCTFile, we get a cleaned matrix.
	finalMatrix, finalGeneList, finalSampleList, err := processGCTFile(
		gctDataFile,
		geneLengthsKB,
		lowExpressionThreshold,
		lowVariancePercentile,
	)
	if err != nil {
		log.Fatalf("Failed: %v", err)
	}

	
	// Step 3: Write the matrix for following analyzing.
	
	log.Println("Generating the matrix:", outputMatrixFile)
	err = writeOutputCSV(outputMatrixFile, finalMatrix, finalGeneList, finalSampleList)
	if err != nil {
		log.Fatalf("Failed in writing the matrix: %v", err)
	}

}

// writeOutputCSV takes in filepath, processed matrix, genelist, samplelist input
// It saves the matrix in local environment.
func writeOutputCSV(filePath string, matrix [][]float64, geneList []string, sampleList []string) error {
	file, err := os.Create(filePath)
	if err != nil {
		return err
	}
	defer file.Close()

	// Write in sample id
	// add a gene_id column
	header := append([]string{"gene_id"}, sampleList...)
	_, err = fmt.Fprintln(file, strings.Join(header, ","))
	if err != nil {
		return err
	}

	// write in the gene data
	for i, geneID := range geneList {
		row := matrix[i]
		rowStr := make([]string, len(row)+1)
		rowStr[0] = geneID
		for j, val := range row {
			rowStr[j+1] = strconv.FormatFloat(val, 'f', 6, 64) 
		}
		_, err = fmt.Fprintln(file, strings.Join(rowStr, ","))
		if err != nil {
			
			log.Printf("Failed to write %s, %v", geneID, err)
		}
	}
	return nil
}