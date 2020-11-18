#include <ctype.h>
#include <mpi.h>
#include <string.h>
#include <chrono>
#include <iostream>
#include <vector>

typedef double value_t;

typedef value_t *Vector;

Vector createVector(int height, int width, int depth);

void releaseVector(Vector m);

bool check_parameters(int argc, char *argv[]);

void compute_sequential(std::vector<value_t> &result, int height, int width, int depth, int timesteps, int heat_source_x, int heat_source_y, int heat_source_z, value_t min, value_t max);

void print_cube(const Vector arr, int height, int width, int depth);

void send_recv_height_split(Vector stecil, int *coord, int chunk, int subrange_width, int num_Procs, MPI_Comm comm);

bool sequential_check = false;
bool print_result = false;
bool print_progression = false;

int main(int argc, char **argv) {
    int rank, num_Procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_Procs);

    if (!check_parameters(argc, argv)) {
        if (rank == 0)
            std::cout << "expected input: ./[file_name] timesteps height width depth [optional parameters]" << std::endl;
        MPI_Finalize();
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
    int subrange_height = height / num_Procs;
    const int slice = depth * width;
    const int chunk = subrange_height * width * depth;

    if (height % num_Procs != 0) {
        std::cout << "No possible tiling, the height be divisible by the number of total cores." << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int dim[1] = {num_Procs};
    int period[1] = {false};
    int coord[1];
    int global_pos, global_x, global_y, global_z;
    Vector stencil = createVector(height, width, depth);
    Vector temp_stencil = createVector(height, width, depth);
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dim, period, true, &comm);
    MPI_Cart_coords(comm, rank, 2, coord);

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

    MPI_Bcast(stencil, height * width * depth, MPI_DOUBLE, 0, comm);

    if (rank == 0) {
        std::cout << "Computing heat-distribution" << std::endl;
        std::cout << "Height: " << height << " Width: " << width << " Depth: " << depth << " Timesteps: " << timesteps << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t t = 0; t < timesteps; t++) {
        // Send and recv data
        if ((coord[0] % 2) == 0) {
            // send upper side
            if (coord[0] != 0)
                MPI_Send(&stencil[coord[0] * chunk], slice, MPI_DOUBLE, coord[0] - 1, 0, comm);

            // send lower side
            if (coord[0] != num_Procs - 1)
                MPI_Send(&stencil[(coord[0] + 1) * chunk - slice], slice, MPI_DOUBLE, coord[0] + 1, 0, comm);

            // receive lower side
            if (coord[0] != num_Procs - 1)
                MPI_Recv(&stencil[(coord[0] + 1) * chunk], slice, MPI_DOUBLE, coord[0] + 1, 0, comm, MPI_STATUS_IGNORE);

            // receive upper side
            if (coord[0] != 0)
                MPI_Recv(&stencil[coord[0] * chunk - slice], slice, MPI_DOUBLE, coord[0] - 1, 0, comm, MPI_STATUS_IGNORE);

        } else {
            // receive lower side
            if (coord[0] != num_Procs - 1)
                MPI_Recv(&stencil[(coord[0] + 1) * chunk], slice, MPI_DOUBLE, coord[0] + 1, 0, comm, MPI_STATUS_IGNORE);

            // receive upper side
            if (coord[0] != 0)
                MPI_Recv(&stencil[coord[0] * chunk - slice], slice, MPI_DOUBLE, coord[0] - 1, 0, comm, MPI_STATUS_IGNORE);

            // send upper side
            if (coord[0] != 0)
                MPI_Send(&stencil[coord[0] * chunk], slice, MPI_DOUBLE, coord[0] - 1, 0, comm);

            // send lower side
            if (coord[0] != num_Procs - 1)
                MPI_Send(&stencil[(coord[0] + 1) * chunk - slice], slice, MPI_DOUBLE, coord[0] + 1, 0, comm);
        }

        // calculate
        for (size_t i = 0; i < subrange_height; i++) {
            for (size_t j = 0; j < width; j++) {
                for (size_t k = 0; k < depth; k++) {
                    global_x = j;
                    global_y = coord[0] * subrange_height + i;
                    global_z = k;
                    global_pos = global_y * width * depth + global_x * depth + global_z;

                    if (heat_source_x == global_x && heat_source_y == global_y && heat_source_z == global_z) {
                        temp_stencil[global_pos] = stencil[global_pos];
                        continue;
                    }

                    // get current temperature at (i,j)
                    value_t tc = stencil[global_pos];

                    // get temperatures left/right and up/down
                    value_t tl = (global_x != 0) ? stencil[global_pos - depth] : tc;
                    value_t tr = (global_x != width - 1) ? stencil[global_pos + depth] : tc;
                    value_t tu = (global_y != 0) ? stencil[global_pos - slice] : tc;
                    value_t td = (global_y != height - 1) ? stencil[global_pos + slice] : tc;
                    value_t tf = (global_z != 0) ? stencil[global_pos - 1] : tc;
                    value_t tb = (global_z != depth - 1) ? stencil[global_pos + 1] : tc;

                    // update temperature at current point
                    temp_stencil[global_pos] = tc + 0.1f * (tl + tr + tu + td + tf + tb + (-6.0f * tc));
                }
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
            std::vector<value_t> sequential(height * width * depth);
            std::cout << "Computing sequentially ..." << std::endl;
            for (size_t i = 0; i < height; i++) {
                for (size_t j = 0; j < width; j++) {
                    for (size_t k = 0; k < depth; k++) {
                        if (i == heat_source_y && j == heat_source_x && k == heat_source_z)
                            sequential[i * width * depth + j * depth + k] = max;
                        else
                            sequential[i * width * depth + j * depth + k] = min;
                    }
                }
            }

            auto start_seq = std::chrono::high_resolution_clock::now();
            compute_sequential(sequential, height, width, depth, timesteps, heat_source_x, heat_source_y, heat_source_z, min, max);
            auto end_seq = std::chrono::high_resolution_clock::now();

            for (size_t i = 0; i < height; i++) {
                for (size_t j = 0; j < width; j++) {
                    for (size_t k = 0; k < depth; k++) {
                        value_t diff = stencil[i * width * depth + j * depth + k] - sequential[i * width * depth + j * depth + k];
                        if (diff > 10E-3 || diff < -10E-3) {
                            std::cout << "Error at position (" << i << ", " << j << "), expected value " << sequential[i * width * depth + j * depth + k] << ", actual value " << stencil[i * width * depth + j * depth + k] << std::endl;
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

            auto duration_seq = std::chrono::duration_cast<std::chrono::microseconds>(end_seq - start_seq);
            std::cout << "Time sequential: " << duration_seq.count() / 1e6 << "sec" << std::endl;

        } else {
            for (size_t i = 0; i < height; i++) {
                for (size_t j = 0; j < width; j++) {
                    for (size_t k = 0; k < depth; k++) {
                        value_t temp = stencil[i * width * depth + j * depth + k];
                        if (!(temp >= min && temp <= max)) {
                            std::cout << "Error at position (" << i << ", " << j << "),"
                                      << "value: " << stencil[i * width * depth + j * depth + k] << std::endl;
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
        }

        if (print_result /*&& verification_success*/) {
            std::cout << "Result stencil:" << std::endl;
            print_cube(stencil, height, width, depth);
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

void send_recv_height_split(Vector stencil, int *coord, int chunk, int subrange_width, int num_Procs, MPI_Comm comm) {
    if ((coord[0] % 2) == 0) {
        // send upper side
        if (coord[0] != 0)
            MPI_Send(&stencil[coord[0] * chunk], subrange_width, MPI_DOUBLE, coord[0] - 1, 0, comm);
        // send lower side
        if (coord[0] != num_Procs - 1)
            MPI_Send(&stencil[(coord[0] + 1) * chunk - subrange_width], subrange_width, MPI_DOUBLE, coord[0] + 1, 0, comm);

        // receive lower side
        if (coord[0] != num_Procs - 1)
            MPI_Recv(&stencil[(coord[0] + 1) * chunk], subrange_width, MPI_DOUBLE, coord[0] + 1, 0, comm, MPI_STATUS_IGNORE);
        // receive upper side
        if (coord[0] != 0)
            MPI_Recv(&stencil[coord[0] * chunk - subrange_width], subrange_width, MPI_DOUBLE, coord[0] - 1, 0, comm, MPI_STATUS_IGNORE);

    } else {
        // receive lower side
        if (coord[0] != num_Procs - 1)
            MPI_Recv(&stencil[(coord[0] + 1) * chunk], subrange_width, MPI_DOUBLE, coord[0] + 1, 0, comm, MPI_STATUS_IGNORE);
        // receive upper side
        if (coord[0] != 0)
            MPI_Recv(&stencil[coord[0] * chunk - subrange_width], subrange_width, MPI_DOUBLE, coord[0] - 1, 0, comm, MPI_STATUS_IGNORE);

        // send upper side
        if (coord[0] != 0)
            MPI_Send(&stencil[coord[0] * chunk], subrange_width, MPI_DOUBLE, coord[0] - 1, 0, comm);
        // send lower side
        if (coord[0] != num_Procs - 1)
            MPI_Send(&stencil[(coord[0] + 1) * chunk - subrange_width], subrange_width, MPI_DOUBLE, coord[0] + 1, 0, comm);
    }
}

void print_cube(const Vector arr, int height, int width, int depth) {
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

Vector createVector(int height, int width, int depth) {
    // create data and index vector
    return (value_t *)malloc(sizeof(value_t) * height * width * depth);
}

void releaseVector(Vector m) { free(m); }