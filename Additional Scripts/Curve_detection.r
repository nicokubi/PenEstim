#Image_Detection

library(imager)

# 1. User Inputs: 
# This could be implemented in a UI or web interface
x_min <- as.numeric(readline(prompt="Enter x-axis minimum value: "))
x_max <- as.numeric(readline(prompt="Enter x-axis maximum value: "))
y_min <- as.numeric(readline(prompt="Enter y-axis minimum value: "))
y_max <- as.numeric(readline(prompt="Enter y-axis maximum value: "))

# 2. Image Pre-processing:
img <- load.image("/Users/nicolaskubista/Desktop/Screenshot 2023-11-02 at 16.52.21.png")
# Check the number of channels
num_channels <- dim(img)[3]

if (num_channels == 3) {
  # Remove the alpha channel
  img <- img[,,1:3]
} else if (num_channels == 1) {
  # Image is already grayscale
  img_gray <- img
} else if (num_channels == 3) {
  # Convert to grayscale
  img_gray <- grayscale(img)
}

# If it's not one of the above cases, you might want to inspect the image manually.

# Thresholding the image
img_bin <- img_gray > 0.5

# 3. Curve Detection:
edges <- cannyEdges(img_bin)

# Extract curve points
curve_points <- which(edges == 1, arr.ind = TRUE)

# 4. Data Extraction:
# Map pixel coordinates to x and y values
x_values <- with(curve_points, (ncol(img_bin) - col) / (ncol(img_bin) - 1) * (x_max - x_min) + x_min)
y_values <- with(curve_points, (nrow(img_bin) - row) / (nrow(img_bin) - 1) * (y_max - y_min) + y_min)

data_points <- data.frame(X=x_values, Y=y_values)
print(data_points)
