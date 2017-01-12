cv2 <- function(x) {
    # Computes cv2 for a vector
    out <- (sd(x)/mean(x))^2
}

