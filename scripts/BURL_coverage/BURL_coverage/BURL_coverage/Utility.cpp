#include "Eigen/Dense"
#include "Utility.hpp"

#include <iomanip>
#include <iostream>
#include <utility>

std::pair<double, double> Utility::Bayesian::credibleInterval(Eigen::VectorXd v){
    std::vector<double> samples;
    for(int i = 0; i < v.size(); i++)
        samples.push_back(v(i));
    std::sort(samples.begin(), samples.end());
    
    size_t n = samples.size();
    
    // Calculate indices for 2.5th and 97.5th percentiles
    double lower_idx = (n - 1) * 0.025;
    double upper_idx = (n - 1) * 0.975;
    
    double lower_bound = samples[lower_idx];
    double upper_bound = samples[upper_idx];

    return {lower_bound, upper_bound};
}

std::pair<double, double> Utility::Bayesian::credibleIntervalBurnIn(Eigen::VectorXd v, double bi){
    std::vector<double> samples;
    for(int i = (int)(bi * v.size()); i < v.size(); i++)
        samples.push_back(v(i));

    std::sort(samples.begin(), samples.end());
    
    size_t n = samples.size();
    
    // Calculate indices for 2.5th and 97.5th percentiles
    double lower_idx = (n - 1) * 0.025;
    double upper_idx = (n - 1) * 0.975;
    
    double lower_bound = samples[lower_idx];
    double upper_bound = samples[upper_idx];

    return {lower_bound, upper_bound};
}

std::pair<double, double> Utility::Bayesian::hpdIntervalBurnIn(Eigen::VectorXd v, double bi){
// 1. Apply burn-in
    int burnin_idx = static_cast<int>(bi * v.size());
    std::vector<double> samples;
    samples.reserve(v.size() - burnin_idx);
    for(int i = burnin_idx; i < v.size(); i++){
        samples.push_back(v(i));
    }

    // 2. Sort samples
    std::sort(samples.begin(), samples.end());
    size_t n = samples.size();
    if(n == 0) throw std::runtime_error("No samples after burn-in!");

    // 3. Determine number of samples to include in the interval
    size_t interval_size = static_cast<size_t>(std::floor(0.95 * n));
    if(interval_size < 1) interval_size = 1;

    // 4. Slide a window of size interval_size to find the narrowest interval
    double min_width = samples[n-1] - samples[0];
    size_t min_idx = 0;
    for(size_t i = 0; i <= n - interval_size; i++){
        double width = samples[i + interval_size - 1] - samples[i];
        if(width < min_width){
            min_width = width;
            min_idx = i;
        }
    }

    // 5. Return the HPD interval
    double lower_bound = samples[min_idx];
    double upper_bound = samples[min_idx + interval_size - 1];

    return {lower_bound, upper_bound};
}

bool Utility::Bayesian::multimodalCredibleInterval(double x, Eigen::VectorXd v, double bi) {
    // Extract samples after burn-in
    std::vector<double> samples;
    for(int i = 0; i < v.size(); i++){
        if(i > (int)(bi * v.size()))
            samples.push_back(v(i));
    }
    std::sort(samples.begin(), samples.end());
    
    size_t n = samples.size();
    if(n == 0) return false;
    
    // Calculate total probability mass to exclude (5% = 0.05)
    double alpha = 0.05;
    int total_exclude = (int)(alpha * n);
    
    // Try all possible ways to exclude total_exclude samples
    // to find the shortest interval(s)
    double min_width = std::numeric_limits<double>::max();
    std::vector<std::pair<double, double>> best_intervals;
    
    // For each possible left tail exclusion from 0 to total_exclude
    for(int left_exclude = 0; left_exclude <= total_exclude; left_exclude++) {
        int right_exclude = total_exclude - left_exclude;
        
        // Check bounds
        if(left_exclude >= n || (n - 1 - right_exclude) < left_exclude) continue;
        
        double lower = samples[left_exclude];
        double upper = samples[n - 1 - right_exclude];
        double width = upper - lower;
        
        // If this is shorter than current best, replace
        if(width < min_width) {
            min_width = width;
            best_intervals.clear();
            best_intervals.push_back({lower, upper});
        }
        // If this is equal to current best, add to list (for multimodal case)
        else if(std::abs(width - min_width) < 1e-10) {
            best_intervals.push_back({lower, upper});
        }
    }
    
    // Check if x falls within any of the credible intervals
    for(const auto& interval : best_intervals) {
        if(x >= interval.first && x <= interval.second) {
            return true;
        }
    }
    
    // Additional check for true multimodal distributions
    // Look for gaps in the shortest interval that might indicate multiple modes
    if(best_intervals.size() == 1) {
        auto main_interval = best_intervals[0];
        
        // Find the density of samples within the main interval
        std::vector<double> interval_samples;
        for(double sample : samples) {
            if(sample >= main_interval.first && sample <= main_interval.second) {
                interval_samples.push_back(sample);
            }
        }
        
        // Look for significant gaps that might indicate multiple modes
        if(interval_samples.size() > 10) { // Need sufficient samples
            std::vector<double> gaps;
            for(size_t i = 1; i < interval_samples.size(); i++) {
                gaps.push_back(interval_samples[i] - interval_samples[i-1]);
            }
            std::sort(gaps.begin(), gaps.end());
            
            // If there are unusually large gaps, treat as multimodal
            double median_gap = gaps[gaps.size()/2];
            double large_gap_threshold = median_gap * 5.0; // Adjust as needed
            
            std::vector<std::pair<double, double>> subintervals;
            double current_start = interval_samples[0];
            
            for(size_t i = 1; i < interval_samples.size(); i++) {
                if(interval_samples[i] - interval_samples[i-1] > large_gap_threshold) {
                    // End current subinterval
                    subintervals.push_back({current_start, interval_samples[i-1]});
                    current_start = interval_samples[i];
                }
            }
            // Add final subinterval
            subintervals.push_back({current_start, interval_samples.back()});
            
            // If we found multiple subintervals, check each one
            if(subintervals.size() > 1) {
                for(const auto& subinterval : subintervals) {
                    if(x >= subinterval.first && x <= subinterval.second) {
                        return true;
                    }
                }
                return false; // x not in any subinterval
            }
        }
    }
    
    return false; // x not in credible interval
}

Eigen::MatrixXd Utility::Shapes::generateUnitCirclePoints(int n) {
    Eigen::MatrixXd points(n, 2);
    const double twoPi = 2.0 * M_PI;

    for (int i = 0; i < n; ++i) {
        double angle = twoPi * i / n;
        points(i, 0) = std::cos(angle);  // x coordinate
        points(i, 1) = std::sin(angle);  // y coordinate
    }

    return points;
}

Eigen::MatrixXd Utility::Shapes::generateUnitSpherePoints(int n) {
    Eigen::MatrixXd points(n, 3);
    const double goldenAngle = M_PI * (3.0 - std::sqrt(5.0));  // ~2.3999632

    for (int i = 0; i < n; ++i) {
        double y = 1.0 - (2.0 * i) / (n - 1);           // y goes from 1 to -1
        double radius = std::sqrt(1.0 - y * y);         // radius at y
        double theta = goldenAngle * i;                 // golden angle increment

        double x = radius * std::cos(theta);
        double z = radius * std::sin(theta);

        points(i, 0) = x;
        points(i, 1) = y;
        points(i, 2) = z;
    }

    return points;
}

Eigen::MatrixXd Utility::EigenUtils::vectorMatrix2Eigen(std::vector<std::vector<double>> v){
    Eigen::MatrixXd scratch = Eigen::MatrixXd::Zero(v.size(), v[0].size());
    for(int i = 0; i < v.size(); i++)
        for(int j = 0; j < v[0].size(); j++)
            scratch(i, j) = v[i][j];
    return scratch;
}

void Utility::EigenUtils::deleteRow(Eigen::MatrixXd& matrix, unsigned int rowToDelete)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToDelete < numRows )
        matrix.block(rowToDelete,0,numRows-rowToDelete,numCols) = matrix.block(rowToDelete+1,0,numRows-rowToDelete,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void Utility::EigenUtils::deleteColumn(Eigen::MatrixXd& matrix, unsigned int colToDelete)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToDelete < numCols )
        matrix.block(0,colToDelete,numRows,numCols-colToDelete) = matrix.block(0,colToDelete+1,numRows,numCols-colToDelete);

    matrix.conservativeResize(numRows,numCols);
}

Eigen::MatrixXd Utility::EigenUtils::kroneckerProduct(Eigen::MatrixXd& matrix0, Eigen::MatrixXd& matrix1){
    int rows0 = matrix0.rows();
    int cols0 = matrix0.cols();
    int rows1 = matrix1.rows();
    int cols1 = matrix1.cols();

    Eigen::MatrixXd result(rows0 * rows1, cols0 * cols1);

    for (int i = 0; i < rows0; ++i)
        for (int j = 0; j < cols0; ++j)
            result.block(i * rows1, j * cols1, rows1, cols1) = matrix0(i, j) * matrix1;

    return result;
}

void Utility::EigenUtils::printEigen(Eigen::MatrixXd e){
    std::cout << std::setprecision(3);
    for(int i = 0; i < e.rows(); i++){
//        std::cout <<"c(";
        for(int j = 0; j < e.cols(); j++){
//            if(j == e.cols()-1){
//                std::cout << e(i, j) << "),\t";
//            }else{
                std::cout << e(i, j) << ",\t";
//            }
        }
        std::cout << std::endl;
    }
}

void Utility::EigenUtils::printEigenR(Eigen::MatrixXd e){
    std::cout << std::setprecision(3);
    for(int i = 0; i < e.rows(); i++){
        std::cout <<"c(";
        for(int j = 0; j < e.cols(); j++){
            if(j == e.cols()-1){
                std::cout << e(i, j) << "),\t";
            }else{
                std::cout << e(i, j) << ",\t";
            }
        }
        std::cout << std::endl;
    }
}
