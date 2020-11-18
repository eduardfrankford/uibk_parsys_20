#include <ctype.h>
#include <string.h>
#include <chrono>
#include <iostream>
#include <vector>

typedef double value_t;

typedef value_t *Vector;

bool check_parameters(int argc, char *argv[]);

void compute_sequential(std::vector<value_t> &result, int height, int width, int depth, int timesteps, int heat_source_x, int heat_source_y, int heat_source_z, value_t min, value_t max);

void print_cube(const std::vector<value_t> &arr, int height, int width, int depth);

bool sequential_check = false;
bool print_result = false;
bool print_progression = false;

int main(int argc, char **argv) {
    if (!check_parameters(argc, argv)) {
        std::cout << "expected input: ./[file_name] timesteps height width depth [optional parameters]" << std::endl;
        return EXIT_FAILURE;
    }

    // ---------- setup ----------
    const value_t max = 273 + 60;
    const value_t min = 273;
    const int timesteps = atoi(argv[1]);
    const int height = atoi(argv[2]);
    const int width = atoi(argv[3]);
    const int depth = atoi(argv[4]);
    const int heat_source_x = width / 2;
    const int heat_source_y = height / 2;
    const int heat_source_z = depth / 2;
    bool verification_success = true;
    std::vector<value_t> stencil(height * width * depth);

    // initialize stencil
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            for (size_t k = 0; k < depth; k++) {
                if (i == heat_source_y && j == heat_source_x && k == heat_source_z)
                    stencil[i * width * depth + j * depth + k] = max;
                else
                    stencil[i * width * depth + j * depth + k] = min;
            }
        }
    }

    std::cout << "Computing heat-distribution" << std::endl;
    std::cout << "Height: " << height << " Width: " << width << " Depth: " << depth << " Timesteps: " << timesteps << std::endl;

    // ---------- compute ----------
    auto start = std::chrono::high_resolution_clock::now();
    compute_sequential(stencil, height, width, depth, timesteps, heat_source_x, heat_source_y, heat_source_z, min, max);
    auto end = std::chrono::high_resolution_clock::now();

    // ---------- check ----------
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            for (size_t k = 0; k < depth; k++) {
                value_t temp = stencil[i * width + j];
                if (!(273 <= temp && temp <= 273 + 60)) {
                    std::cout << "Error at position (" << i << ", " << j << "), value: " << stencil[i * width + j] << std::endl;
                    verification_success = false;
                    break;
                }
            }
            if (!verification_success)
                break;
        }
        if (!verification_success)
            break;
    }

    if (print_result && verification_success) {
        std::cout << "Result stencil:" << std::endl;
        print_cube(stencil, height, width, depth);
    }

    if (verification_success)
        std::cout << "Verification OK" << std::endl;
    else
        std::cout << "Verification FAILED" << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time: " << duration.count() / 1e6 << "sec" << std::endl;

    return (verification_success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

bool check_parameters(int argc, char *argv[]) {
    if (argc < 5) {
        return false;
    }

    bool is_numeric = true;
    for (size_t i = 1; i < 5; i++) {
        for (size_t j = 0; j < strlen(argv[i]); j++)
            if (!isdigit(argv[i][j]))
                is_numeric = false;
    }

    if (argc >= 6) {
        for (size_t i = 5; i < argc; i++) {
            if (!strcmp(argv[i], "-s"))
                sequential_check = true;
            else if (!strcmp(argv[i], "-p"))
                print_result = true;
            else if (!strcmp(argv[i], "-pr"))
                print_progression = true;
        }
    }

    return is_numeric;
}

void compute_sequential(std::vector<value_t> &result, int height, int width, int depth, int timesteps, int heat_source_x, int heat_source_y, int heat_source_z, value_t min, value_t max) {
    std::vector<value_t> temp_vec(height * width * depth);

    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            for (size_t k = 0; k < depth; k++) {
                if (i == heat_source_y && j == heat_source_x && k == heat_source_z)
                    result[i * width * depth + j * depth + k] = max;
                else
                    result[i * width * depth + j * depth + k] = min;
            }
        }
    }

    for (size_t t = 0; t < timesteps; t++) {
        for (size_t i = 0; i < height; i++) {
            for (size_t j = 0; j < width; j++) {
                for (size_t k = 0; k < depth; k++) {
                    // heat source stays at the same temperature
                    if (j == heat_source_x && i == heat_source_y && k == heat_source_z) {
                        temp_vec[i * width * depth + j * depth + k] = result[i * width * depth + j * depth + k];
                        continue;
                    }

                    // get current temperature at (i,j,k)
                    value_t tc = result[i * width * depth + j * depth + k];

                    // get temperatures left/right and up/down
                    value_t tl = (j != 0) ? result[i * width * depth + (j - 1) * depth + k] : tc;
                    value_t tr = (j != width - 1) ? result[i * width * depth + (j + 1) * depth + k] : tc;
                    value_t tu = (i != 0) ? result[(i - 1) * width * depth + j * depth + k] : tc;
                    value_t td = (i != height - 1) ? result[(i + 1) * width * depth + j * depth + k] : tc;
                    value_t tf = (k != 0) ? result[i * width * depth + j * depth + (k - 1)] : tc;
                    value_t tb = (k != depth - 1) ? result[i * width * depth + j * depth + (k + 1)] : tc;

                    // update temperature at current point
                    temp_vec[i * width * depth + j * depth + k] = tc + 0.1f * (tl + tr + tu + td + tf + tb + (-6.0f * tc));
                }
            }
        }

        // swap buffer (just pointers, not content)
        result.swap(temp_vec);
        if (print_progression && (t % 100 == 0 || t == timesteps - 1))
            std::cout << "\rProgression: " << (int)(((float)t) / (timesteps - 1) * 100) << "%" << std::flush;
    }
    if (print_progression)
        std::cout << std::endl;
}

void print_cube(const std::vector<value_t> &arr, int height, int width, int depth) {
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            for (size_t k = 0; k < depth; k++) {
                if (k == depth)
                    std::cout << arr[i * width * depth + j * depth + k];
                else
                    std::cout << arr[i * width * depth + j * depth + k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
