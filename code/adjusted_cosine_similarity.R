#adjusted_cosine_similarity.R
adjusted_cosine_similarity <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  mean_A <- mean(A, na.rm = TRUE)
  mean_B <- rowMeans(B, na.rm = TRUE)  
  adjusted_A <- A - mean_A 
  adjusted_B <- sweep(B, 2, mean_B, FUN = "-")  
  similarity_scores <- numeric(nrow(B))  
  for (i in 1:nrow(B)) {
    numerator <- sum(adjusted_A * adjusted_B[i, ], na.rm = TRUE)
    denominator <- sqrt(sum(adjusted_A^2, na.rm = TRUE)) * sqrt(sum(adjusted_B[i, ]^2, na.rm = TRUE))
    similarity_scores[i] <- ifelse(denominator == 0, NA, numerator / denominator)
  }
  return(similarity_scores)
}