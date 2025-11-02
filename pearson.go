package main

import (
	"encoding/csv" 
	"fmt"
	"log"
	"math"
	"os" 
	"runtime"
	"strconv" 
	"sync"
)

// RunPhase2 starts the parallel calculation of the correlation matrix.
func RunPhase2(matrix [][]float64, geneList []string) ([][]float64, error) {
	numGenes := len(geneList)
	if numGenes == 0 {
		return nil, fmt.Errorf("matrix is empty")
	}
	numSamples := len(matrix[0])
	log.Printf("  (P2) Pre-calculating mean and stddev for %d genes...", numGenes)

	// 1. Pre-calculate Mean and StdDev for all genes .
	means := make([]float64, numGenes)
	stdDevs := make([]float64, numGenes)
	for i := 0; i < numGenes; i++ {
		means[i] = mean(matrix[i])                
		stdDevs[i] = math.Sqrt(variance(matrix[i])) 
	}
	log.Println("  (P2) ...Pre-calculation complete.")

	// 2. Setup worker pool.
	// We only calculate the upper triangle of the matrix (j > i).
	numJobs := numGenes * (numGenes - 1) / 2

	// jobs channel: stores pairs of indices to compare [i, j].
	jobs := make(chan [2]int, numJobs)
	// results channel: stores the result [i, j, correlation].
	results := make(chan [3]float64, numJobs)

	// Use all available CPU cores.
	numWorkers := runtime.NumCPU()
	log.Printf("  (P2) Starting %d workers for %d correlation jobs...", numWorkers, numJobs)

	var wg sync.WaitGroup

	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		// Start a worker goroutine.
		go correlationWorker(
			&wg,
			jobs,
			results,
			matrix,
			means,
			stdDevs,
			numSamples,
		)
	}

	// 3. Dispatch jobs.
	// Start a goroutine to feed all job pairs into the channel.
	go func() {
		for i := 0; i < numGenes; i++ {
			for j := i + 1; j < numGenes; j++ {
				jobs <- [2]int{i, j}
			}
		}
		// All jobs are sent, close the channel.
		close(jobs)
	}()

	// 4. Start result collector.
	// Initialize the final correlation matrix (using float64 for precision).
	corrMatrix := make([][]float64, numGenes)
	for i := range corrMatrix {
		corrMatrix[i] = make([]float64, numGenes)
	}

	// Start a goroutine to collect all results.
	var collectWg sync.WaitGroup
	collectWg.Add(1)
	go func() {
		defer collectWg.Done()
		
		// Wait for all results to come in
		for i := 0; i < numJobs; i++ {
			res := <-results
			idxI := int(res[0])
			idxJ := int(res[1])
			corr := res[2]

			// Fill the symmetric matrix
			corrMatrix[idxI][idxJ] = corr
			corrMatrix[idxJ][idxI] = corr
		}
	}()

	// 5. Wait for all goroutines to finish.
	wg.Wait()        // Wait for all workers to be done.
	close(results)   // Workers are done, close the results channel.
	collectWg.Wait() // Wait for the collector to finish.

	// Fill the diagonal (self-correlation is always 1).
	for i := 0; i < numGenes; i++ {
		corrMatrix[i][i] = 1.0
	}

	log.Println("  (P2) ...All correlation tasks complete!")
	return corrMatrix, nil
}

// correlationWorker calculates the Pearson correlation for jobs it receives.
func correlationWorker(
	wg *sync.WaitGroup,
	jobs <-chan [2]int,
	results chan<- [3]float64,
	matrix [][]float64,
	means, stdDevs []float64,
	numSamples int,
) {
	defer wg.Done()
	n := float64(numSamples)

	for job := range jobs {
		i := job[0]
		j := job[1]

		vA := matrix[i]
		vB := matrix[j]
		meanA := means[i]
		meanB := means[j]
		stdDevA := stdDevs[i]
		stdDevB := stdDevs[j]

		// Avoid division by zero if variance is zero (gene is constant).
		if stdDevA == 0 || stdDevB == 0 {
			results <- [3]float64{float64(i), float64(j), 0.0}
			continue
		}

		// Calculate covariance numerator.
		covariance := 0.0
		for k := 0; k < numSamples; k++ {
			covariance += (vA[k] - meanA) * (vB[k] - meanB)
		}

		// Pearson r = Cov(A, B) / (StdDev(A) * StdDev(B))
		// We use (covariance / n) for Cov(A, B)
		corr := (covariance / n) / (stdDevA * stdDevB)

		results <- [3]float64{float64(i), float64(j), corr}
	}
}

// writeCorrelationMatrix saves the final correlation matrix to a CSV file.
func writeCorrelationMatrix(
	filePath string,
	matrix [][]float64,
	geneList []string,
) error {
	// Create the output file
	file, err := os.Create(filePath)
	if err != nil {
		return fmt.Errorf("failed to create correlation file %s: %w", filePath, err)
	}
	defer file.Close()

	// Use a csv.Writer for efficient writing
	writer := csv.NewWriter(file)
	defer writer.Flush() // Ensure all buffered data is written

	numGenes := len(geneList)
	if numGenes == 0 {
		return fmt.Errorf("gene list is empty, nothing to write")
	}

	// 1. Write the header row
	// The header is: "gene_id", "GENE_1", "GENE_2", ...
	header := make([]string, numGenes+1)
	header[0] = "gene_id" // First column header
	copy(header[1:], geneList)

	if err := writer.Write(header); err != nil {
		return fmt.Errorf("failed to write correlation header: %w", err)
	}

	// 2. Write the matrix data row by row
	// Create a reusable slice to reduce memory allocations
	rowStr := make([]string, numGenes+1)

	for i := 0; i < numGenes; i++ {
		// The first column of each row is the gene name (row header)
		rowStr[0] = geneList[i]

		// Iterate through all correlation values in this row
		for j := 0; j < numGenes; j++ {
			// Convert the float64 correlation value to a string
			rowStr[j+1] = strconv.FormatFloat(matrix[i][j], 'f', 6, 64)
		}

		// Write the complete row to the file
		if err := writer.Write(rowStr); err != nil {
			// Log a warning but continue trying to write other rows
			log.Printf("warning: failed to write correlation row for gene %s: %v", geneList[i], err)
		}
	}

	log.Println("  (P2) successfully saved correlation matrix to:", filePath)
	return nil
}