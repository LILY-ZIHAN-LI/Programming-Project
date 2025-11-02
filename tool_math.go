package main

// variance calculate the variance of a slice
func variance(data []float64) float64 {
	if len(data) == 0 {
		return 0.0
	}
	m := mean(data)
	sumSq := 0.0
	for _, val := range data {
		sumSq += (val - m) * (val - m)
	}
	// WGCNA usually uses "population variance"
	return sumSq / float64(len(data)) 
}

// mean calculate the mean of a slice.
func mean(data []float64) float64 {
	if len(data) == 0 {
		return 0.0
	}
	sum := 0.0
	for _, val := range data {
		sum += val
	}
	return sum / float64(len(data))
}