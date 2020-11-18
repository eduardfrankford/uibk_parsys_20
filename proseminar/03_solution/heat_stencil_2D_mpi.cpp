#include <ctype.h>
#include <mpi.h>
#include <string.h>
#include <chrono>
#include <iostream>
#include <vector>

typedef double value_t;

typedef value_t *Vector;

Vector createVector(int height, int width);

void releaseVector(Vector m);

bool check_parameters(int argc, char *argv[]);

void compute_sequential(std::vector<value_t> &result, int height, int width, int timesteps, int heat_source_x, int heat_source_y, value_t min, value_t max);

void printTemperature(const Vector r, int N, int M, value_t min, value_t max);

void print_matrix(const Vector arr, int height, int width);

bool sequential_check = false;
bool print_result = false;
bool print_progression = false;
bool draw_inter = false;
bool split_height = true;

int main(int argc, char **argv) {
    int rank, num_Procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_Procs);

    if (!check_parameters(argc, argv)) {
        if (rank == 0)
            std::cout << "expected input: ./[file_name] timesteps height width [optional parameters]" << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // ---------- setup ----------
    const value_t max = 273 + 60;
    const value_t min = 273;
    const int timesteps = atoi(argv[1]);
    const int height = atoi(argv[2]);
    const int width = atoi(argv[3]);
    const int heat_source_x = width / 2;
    const int heat_source_y = height / 2;
    bool verification_success = true;
    int subrange_height = height / num_Procs;
    int chunk = subrange_height * width;

    if (height % num_Procs != 0) {
        std::cout << "No possible tiling, the height be divisible by the number of total cores." << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int dim[1] = {num_Procs};
    int period[1] = {false};
    int coord[1];
    int global_pos, global_x, global_y;
    Vector stencil = createVector(height, width);
    Vector temp_stencil = createVector(height, width);
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dim, period, true, &comm);
    MPI_Cart_coords(comm, rank, 2, coord);

    // initialize stencil
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            if (i == heat_source_y && j == heat_source_x)
                stencil[i * width + j] = 273 + 60;
            else
                stencil[i * width + j] = 273;
        }
    }

    MPI_Bcast(stencil, height * width, MPI_DOUBLE, 0, comm);

    if (rank == 0) {
        std::cout << "Computing heat-distribution" << std::endl;
        std::cout << "Height: " << height << " Width: " << width << " Timesteps: " << timesteps << std::endl;
    }

    // ---------- compute ----------
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t t = 0; t < timesteps; t++) {
        // Send and recv data
        if ((coord[0] % 2) == 0) {
            // send upper side
            if (coord[0] != 0)
                MPI_Send(&stencil[coord[0] * chunk], width, MPI_DOUBLE, coord[0] - 1, 0, comm);
            // send lower side
            if (coord[0] != num_Procs - 1)
                MPI_Send(&stencil[(coord[0] + 1) * chunk - width], width, MPI_DOUBLE, coord[0] + 1, 0, comm);

            // receive lower side
            if (coord[0] != num_Procs - 1)
                MPI_Recv(&stencil[(coord[0] + 1) * chunk], width, MPI_DOUBLE, coord[0] + 1, 0, comm, MPI_STATUS_IGNORE);
            // receive upper side
            if (coord[0] != 0)
                MPI_Recv(&stencil[coord[0] * chunk - width], width, MPI_DOUBLE, coord[0] - 1, 0, comm, MPI_STATUS_IGNORE);

        } else {
            // receive lower side
            if (coord[0] != num_Procs - 1)
                MPI_Recv(&stencil[(coord[0] + 1) * chunk], width, MPI_DOUBLE, coord[0] + 1, 0, comm, MPI_STATUS_IGNORE);
            // receive upper side
            if (coord[0] != 0)
                MPI_Recv(&stencil[coord[0] * chunk - width], width, MPI_DOUBLE, coord[0] - 1, 0, comm, MPI_STATUS_IGNORE);

            // send upper side
            if (coord[0] != 0)
                MPI_Send(&stencil[coord[0] * chunk], width, MPI_DOUBLE, coord[0] - 1, 0, comm);
            // send lower side
            if (coord[0] != num_Procs - 1)
                MPI_Send(&stencil[(coord[0] + 1) * chunk - width], width, MPI_DOUBLE, coord[0] + 1, 0, comm);
        }

        // calculate
        for (size_t i = 0; i < subrange_height; i++) {
            for (size_t j = 0; j < width; j++) {
                global_x = split_height ? j : coord[0] * width + j;
                global_y = split_height ? coord[0] * subrange_height + i : i;
                global_pos = global_y * width + global_x;

                if (heat_source_x == global_x && heat_source_y == global_y) {
                    temp_stencil[global_pos] = stencil[global_pos];
                    continue;
                }

                // get current temperature at (i,j)
                value_t tc = stencil[global_pos];

                // get temperatures left/right and up/down
                value_t tl = (global_x != 0) ? stencil[global_pos - 1] : tc;
                value_t tr = (global_x != width - 1) ? stencil[global_pos + 1] : tc;
                value_t tu = (global_y != 0) ? stencil[global_pos - width] : tc;
                value_t td = (global_y != height - 1) ? stencil[global_pos + width] : tc;

                // update temperature at current point
                temp_stencil[global_pos] = tc + 0.2 * (tl + tr + tu + td + (-4.0f * tc));
            }
        }

        // draw
        if (draw_inter && (t % 1000 == 0 || t == timesteps - 1)) {
            MPI_Gather(&stencil[coord[0] * chunk], chunk, MPI_DOUBLE, stencil, chunk, MPI_DOUBLE, 0, comm);
            if (rank == 0) {
                std::cout << "Step t=" << t << ":\t" << std::endl;
                printTemperature(stencil, height, width, min, max);
                std::cout << std::endl;
            }
        }

        // swap matrices (just pointers, not content)
        Vector H = stencil;
        stencil = temp_stencil;
        temp_stencil = H;
    }
    MPI_Gather(&stencil[coord[0] * chunk], chunk, MPI_DOUBLE, stencil, chunk, MPI_DOUBLE, 0, comm);
    auto end = std::chrono::high_resolution_clock::now();

    // ---------- check ----------
    if (rank == 0) {
        if (sequential_check) {
            std::vector<value_t> sequential(height * width);
            std::cout << "Computing sequentially ..." << std::endl;
            for (size_t i = 0; i < height; i++) {
                for (size_t j = 0; j < width; j++) {
                    if (i == heat_source_y && j == heat_source_x)
                        sequential[i * width + j] = max;
                    else
                        sequential[i * width + j] = min;
                }
            }

            auto start_seq = std::chrono::high_resolution_clock::now();
            compute_sequential(sequential, height, width, timesteps, heat_source_x, heat_source_y, min, max);
            auto end_seq = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < height; i++) {
                for (size_t j = 0; j < width; j++) {
                    value_t diff = stencil[i * width + j] - sequential[i * width + j];
                    if (diff > 10E-3 || diff < -10E-3) {
                        std::cout << "Error at position (" << i << ", " << j << "), expected value " << sequential[i * width + j] << ", actual value " << stencil[i * width + j] << std::endl;
                        verification_success = false;
                        break;
                    }
                }
                if (!verification_success)
                    break;
            }
            auto duration_seq = std::chrono::duration_cast<std::chrono::microseconds>(end_seq - start_seq);
            std::cout << "Time sequential: " << duration_seq.count() / 1e6 << "sec" << std::endl;

        } else {
            for (size_t i = 0; i < height; i++) {
                for (size_t j = 0; j < width; j++) {
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
        }

        if (print_result /*&& verification_success*/) {
            std::cout << "Result stencil:" << std::endl;
            print_matrix(stencil, height, width);
        }
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time parallel: " << duration.count() / 1e6 << "sec" << std::endl;

        if (verification_success)
            std::cout << "Verification OK" << std::endl;
        else
            std::cout << "Verification FAILED" << std::endl;
    }

    // ---------- cleanup ----------
    releaseVector(stencil);
    releaseVector(temp_stencil);
    MPI_Finalize();
    return (verification_success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

bool check_parameters(int argc, char *argv[]) {
    if (argc < 4) {
        return false;
    }

    bool is_numeric = true;
    for (size_t i = 1; i < 4; i++) {
        for (size_t j = 0; j < strlen(argv[i]); j++)
            if (!isdigit(argv[i][j]))
                is_numeric = false;
    }

    if (argc >= 5) {
        for (size_t i = 4; i < argc; i++) {
            if (!strcmp(argv[i], "-s"))
                sequential_check = true;
            else if (!strcmp(argv[i], "-p"))
                print_result = true;
            else if (!strcmp(argv[i], "-d"))
                draw_inter = true;
            else if (!strcmp(argv[i], "-pr"))
                print_progression = true;
        }
    }

    return is_numeric;
}

void compute_sequential(std::vector<value_t> &result, int height, int width, int timesteps, int heat_source_x, int heat_source_y, value_t min, value_t max) {
    std::vector<value_t> temp_vec(height * width);

    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            if (i == heat_source_y && j == heat_source_x)
                result[i * width + j] = max;
            else
                result[i * width + j] = min;
        }
    }

    for (size_t t = 0; t < timesteps; t++) {
        for (size_t i = 0; i < height; i++) {
            for (size_t j = 0; j < width; j++) {
                // heat source stays at the same temperature
                if (j == heat_source_x && i == heat_source_y) {
                    temp_vec[i * width + j] = result[i * width + j];
                    continue;
                }

                // get current temperature at (i,j)
                value_t tc = result[i * width + j];

                // get temperatures left/right and up/down
                value_t tl = (j != 0) ? result[i * width + (j - 1)] : tc;
                value_t tr = (j != width - 1) ? result[i * width + (j + 1)] : tc;
                value_t tu = (i != 0) ? result[(i - 1) * width + j] : tc;
                value_t td = (i != height - 1) ? result[(i + 1) * width + j] : tc;

                // update temperature at current point
                temp_vec[i * width + j] = tc + 0.2 * (tl + tr + tu + td + (-4.0f * tc));
            }
        }

        // swap buffers (just pointers, not content)
        result.swap(temp_vec);
        if (print_progression && (t % 100 == 0 || t == timesteps - 1))
            std::cout << "\rProgression: " << (int)(((float)t) / (timesteps - 1) * 100) << "%" << std::flush;
    }
    if (print_progression)
        std::cout << std::endl;
}

void printTemperature(const Vector r, int N, int M, value_t min, value_t max) {
    const char *colors = " .-:=+*#%@";
    const int numColors = 10;

    // set the 'render' resolution
    int H = 30;
    int W = 50;

    // step size in each dimension
    int sH = N / H;
    int sW = M / W;

    // upper wall
    for (int i = 0; i < W + 2; i++) {
        std::cout << "X";
    }
    std::cout << std::endl;

    // room
    for (int i = 0; i < H; i++) {
        // left wall
        std::cout << "X";

        // actual room
        for (int j = 0; j < W; j++) {
            // get max temperature in this tile
            value_t max_t = 0;
            for (int x = sH * i; x < sH * i + sH; x++) {
                for (int y = sW * j; y < sW * j + sW; y++) {
                    max_t = (max_t < r[x * M + y]) ? r[x * M + y] : max_t;
                }
            }
            value_t temp = max_t;

            // pick the 'color'
            int c = ((temp - min) / (max - min)) * numColors;
            c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

            // print the average temperature
            std::cout << colors[c];
        }

        // right wall
        std::cout << "X" << std::endl;
    }

    // lower wall
    for (int i = 0; i < W + 2; i++) {
        std::cout << "X";
    }
    std::cout << std::endl;
}

void print_matrix(const Vector arr, int height, int width) {
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            if (j == width)
                std::cout << arr[i * width + j];
            else
                std::cout << arr[i * width + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

Vector createVector(int height, int width) {
    // create data and index vector
    return (value_t *)malloc(sizeof(value_t) * height * width);
}

void releaseVector(Vector m) { free(m); }