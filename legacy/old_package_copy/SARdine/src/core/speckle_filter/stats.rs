use ndarray::Array2;

/// Integral image structure for O(1) window statistics (sum, sumsq, count).
#[derive(Debug, Clone)]
pub struct IntegralImage {
    pub sum_table: Array2<f64>,
    pub sum_sq_table: Array2<f64>,
    pub count_table: Array2<u64>,
}

impl IntegralImage {
    pub fn new(image: &Array2<f32>) -> Self {
        let (h, w) = image.dim();
        let mut sum_table = Array2::zeros((h + 1, w + 1));
        let mut sum_sq_table = Array2::zeros((h + 1, w + 1));
        let mut count_table = Array2::zeros((h + 1, w + 1));

        for i in 1..=h {
            let mut row_sum = 0.0f64;
            let mut row_sum_sq = 0.0f64;
            let mut row_cnt = 0u64;
            for j in 1..=w {
                let pixel = image[[i - 1, j - 1]];
                if pixel.is_finite() && pixel >= 0.0 {
                    row_sum += pixel as f64;
                    row_sum_sq += (pixel as f64) * (pixel as f64);
                    row_cnt += 1;
                }
                sum_table[[i, j]] = sum_table[[i - 1, j]] + row_sum;
                sum_sq_table[[i, j]] = sum_sq_table[[i - 1, j]] + row_sum_sq;
                count_table[[i, j]] = count_table[[i - 1, j]] + row_cnt;
            }
        }

        Self {
            sum_table,
            sum_sq_table,
            count_table,
        }
    }

    /// Calculate window mean and variance in O(1).
    pub fn window_stats(&self, row1: usize, col1: usize, row2: usize, col2: usize) -> (f32, f32) {
        if row2 < row1 || col2 < col1 {
            return (0.0, 0.0);
        }
        let (h, w) = (self.sum_table.nrows() - 1, self.sum_table.ncols() - 1);
        if row2 >= h || col2 >= w {
            return (0.0, 0.0);
        }

        let sum = self.sum_table[[row2 + 1, col2 + 1]]
            - self.sum_table[[row1, col2 + 1]]
            - self.sum_table[[row2 + 1, col1]]
            + self.sum_table[[row1, col1]];

        let sum_sq = self.sum_sq_table[[row2 + 1, col2 + 1]]
            - self.sum_sq_table[[row1, col2 + 1]]
            - self.sum_sq_table[[row2 + 1, col1]]
            + self.sum_sq_table[[row1, col1]];

        let count = (self.count_table[[row2 + 1, col2 + 1]] as i128)
            - (self.count_table[[row1, col2 + 1]] as i128)
            - (self.count_table[[row2 + 1, col1]] as i128)
            + (self.count_table[[row1, col1]] as i128);
        let count = count.max(0) as u64;

        if count == 0 {
            return (0.0, 0.0);
        }

        let n = count as f64;
        let mean = sum / n;
        let variance = (sum_sq / n - mean * mean).max(0.0);
        (mean as f32, variance as f32)
    }
}
