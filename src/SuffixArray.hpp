// #include <cstdint>

#include <atomic>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>

#include "parlay/parallel.h"
// #include <divsufsort64.h>


typedef __uint128_t hash_t;
typedef uint64_t index_t;

class SuffixArray {
    private:
		//enum to specify merge type
		enum MergeType {max_parallel, binary_parallel, serial};
		// struct to store Suffix information
		struct Suffix{
			hash_t hash;
		};

		// genome string as char array
        char* genome;
		// length of genome string
        index_t genome_length;
		// array to store suffix objects
        // Suffix* suffixes;
        Suffix** suffixes;
		// defines block size - i.e. how many characters to use as a hash for the suffix
        index_t block_size;
		// stores actual suffix array with indices to the corresponding suffix in the genome
		index_t* suffix_array;

		MergeType merge_type;

		// to keep track of the number of comparisons we do in sorting
		std::atomic<uint64_t> num_extra_compares;
		std::atomic<uint64_t> num_extra_compares_over_2;
		std::atomic<uint64_t> num_extra_compares_over_4;
		std::atomic<uint64_t> num_extra_compares_over_8;
		std::atomic<uint64_t> num_extra_compares_over_16;
		std::atomic<uint64_t> num_extra_compares_over_32;
		std::atomic<uint64_t> num_extra_compares_over_64;

		bool keep_track_of_extra_comps;
		


		hash_t calc_suffix_hash(index_t genome_index) {
			hash_t hash = 0;
			char* suffix = &this->genome[genome_index] ;
			for (hash_t i = 0; i < this->block_size; i++) {
				hash *= 5;
				if (genome_index + i < genome_length){
					switch (suffix[i]){
					case 'A':
						hash += 1;
						break;
					case 'C':
						hash += 2;
						break;

					case 'G':
						hash += 3;
						break;

					case 'T':
						hash += 4;
						break; 

					case '$':
						break;
					}

				}
				// else {
				// 	hash += 1;
				// }
			}
			return hash;
		}

        // reads fasta file, and stores it as the genome
		void read_fasta(char* file_location) {
			// genome = new std::string();
			std::ifstream fasta_file;
			fasta_file.open(file_location, std::ios::in);
			if (!fasta_file.is_open()) {
				printf("Couldn't open fasta file\n");
				throw std::invalid_argument("couldn't open FASTA file");
			} else {

				// parsing fasta, reading into a temp c++ string, then copying into a char*
				std::string line, description;
				this->genome_length = 0; // even for an empty genome, we need space to fit the $ and \0 chars
				this->genome = nullptr;
				while(getline(fasta_file, line)) {
					if(!line.empty()) {
						if(line[0] == '>') { // description line
							// only reading in one genome
							if(description.empty()) {
								description += line;
							} else {
								break;
							}
						} else {
							// reallocating for the line
							this->genome = (char*) realloc(this->genome, this->genome_length + 2 + line.length()); // since we need to leave space for the last two chars
							//copying read into the genome string
							line.copy(this->genome + this->genome_length, line.length());
							// need to update the genome length
							this->genome_length += line.length();
						}
					}
					
				}
				this->genome[this->genome_length++] = '$';
				this->genome[this->genome_length] = '\0';
				
			}

		}

		// makes every letter uppercase in genome
		void genome_to_upper() {
			// for(index_t i = 0; i < this->genome_length - 1; i++) {
			// 	this->genome[i] = toupper(this->genome[i]);
			// }
			parlay::parallel_for(0, this->genome_length - 1, [&] (index_t i) {
                this->genome[i] = toupper(this->genome[i]);
            });
		}

		// Initializes array of Suffix structs and suffix array
        void set_up_suffixes() {
           	// assumes done after reading genome
            this->suffixes = (Suffix**) calloc(block_size, sizeof(Suffix*));

			
			for (index_t i = 0; i < block_size; i++) {
				
				hash_t block_length = genome_length / block_size;
				
				if (i < genome_length % block_size) {
					block_length++;
				}

				suffixes[i] = (Suffix*) calloc(block_length, sizeof(Suffix));
				
				// for(index_t j = 0; j < block_length; j++) {
				// 	index_t genome_index = j * block_size + i;
				// 	suffixes[i][j].hash = calc_suffix_hash(genome_index);
				// }
				parlay::parallel_for(0, block_length, [&] (index_t j) {
                    index_t genome_index = j * block_size + i;
                    suffixes[i][j].hash = calc_suffix_hash(genome_index);
                });
			}
			// this->suffixes = (Suffix*) calloc(genome_length, sizeof(Suffix));

			// set up suffix_array
			this->suffix_array = (index_t*) calloc(genome_length, sizeof(index_t));
			// #pragma omp parallel for
			// for(index_t genome_i = 0; genome_i < genome_length; genome_i++) {
			// 	suffix_array[genome_i] = genome_i;
			// }
			parlay::parallel_for(0, genome_length, [&] (index_t genome_i) {
				suffix_array[genome_i] = genome_i;
			});

        }

        void sort_suffixes() {

			index_t* tmp_arr = (index_t*) calloc(genome_length, sizeof(index_t));
			// call sort
			parallel_merge_sort(suffix_array, genome_length, tmp_arr, 0);
			// very_slow_sort();
        }

		void parallel_merge_sort(index_t* array, index_t length, index_t* tmp_arr, uint depth) {
			// base case
			if(length < 2) {
				for(index_t i = 0; i < length; i++) {
                    tmp_arr[i] = array[i];
                }
				return;
			}

			// if(true) {
			// if(depth > 0) {
			// 	// serial_merge_sort(array, length, tmp_arr, depth);
			// 	use_external_sort(offset, length, tmp_arr_start, depth);
			// 	return;
			// }

			// if(depth < 5) {
			// 	for(uint i = 0; i < depth; i++) {
   //                  std::cout << "\t";
   //              }
			// 	std::cout << "starting sorting " << length << " indices." << std::endl;
			// }

			index_t left_length, right_length;
			index_t* left_arr, *right_arr;
			index_t* left_tmp_arr, *right_tmp_arr;

			left_length = length / 2;
			right_length = length - left_length;

			left_tmp_arr = tmp_arr;
			right_tmp_arr = tmp_arr + left_length;

			left_arr = array;
			right_arr = array + left_length;
			parlay::par_do(
				[&] {
					// call recursive sort
					parallel_merge_sort(left_arr, left_length, left_tmp_arr, depth + 1);
				}, [&] {
					// call recursive sort
					parallel_merge_sort(right_arr, right_length, right_tmp_arr, depth + 1);
				}
			);

			// merge left and right halves
			merge(left_tmp_arr, left_length, right_tmp_arr, right_length, array);

			// if(depth < 5) {
			// 	std::string prefix = "";
			// 	for(uint i = 0; i < depth; i++) {
			// 		prefix.append("\t");
			// 	}
			// 	std::cout << prefix << "finished sorting " << length << " indices." << std::endl;
			// }
		}

		// void use_external_sort(index_t offset, index_t length, index_t* tmp_arr, uint depth) {
		// 	if(divsufsort64( (unsigned char *) genome + offset, suffix_array + offset, length) != 0) {
		// 		throw std::invalid_argument("Couldn't build suffix array");
		// 	}
		// 	parlay::parallel_for(0, length, [&] (index_t i) {
		//      tmp_arr[i + offset] = suffix_array[i + offset];
		//  });
		// }

		void serial_merge_sort(index_t* array, index_t length, index_t* tmp_arr, uint depth) {
			// base case
			if(length < 3) {
				for(index_t i = 0; i < length; i++) {
                    tmp_arr[i] = array[i];
                }
				return;
			}

			// if(depth > 3) {
			// 	parallel_merge_sort(array, length, tmp_arr, depth);
			// 	return;
			// }

			if(length > 10000000) {
				for(uint i = 0; i < depth; i++) {
					std::cout << "\t";
				}
				std::cout << "started sorting " << length << " indices." << std::endl;
			}

			index_t left_length, right_length;
			index_t* left_arr, *right_arr;
			index_t* left_tmp_arr, *right_tmp_arr;

			left_length = length / 2;
			right_length = length - left_length;

			left_tmp_arr = tmp_arr;
			right_tmp_arr = tmp_arr + left_length;

			left_arr = array;
			right_arr = array + left_length;


			// call recursive sort
			serial_merge_sort(left_arr, left_length, left_tmp_arr, depth + 1);

			// call recursive sort
			serial_merge_sort(right_arr, right_length, right_tmp_arr, depth + 1);

			// merge left and right halves
			serial_merge(left_tmp_arr, left_length, right_tmp_arr, right_length, array);

			if(length > 10000000) {
				for(uint i = 0; i < depth; i++) {
					std::cout << "\t";
				}
				std::cout << "finished sorting " << length << " indices." << std::endl;
			}

		}

		index_t binary_search(index_t* arr, index_t length, index_t key) {
			index_t low = 0;
			if(length < 1) {
				return low;
			}
			index_t high = length;
			index_t mid;
			while(low < high) {
				mid = (low + high) / 2;
				// if(arr[mid] == key) {
				// 	return mid;
				// }
				if(suffix_less_than(arr[mid], key)) {
					low = mid + 1;
				} else {
					high = mid;
				}
			}
			return low;
		}

		void max_parallel_merge(index_t* left_arr, index_t left_length, index_t* right_arr, index_t right_length, index_t* result_arr) {
			// we'll start with the left
			parlay::par_do([&] {
				// find global index for suffix in left
				parlay::parallel_for(0, left_length, [&] (index_t i) {
					index_t key = left_arr[i];
					index_t index = binary_search(right_arr, right_length, key);
					// insert key into result array
					result_arr[i + index] = key;
                });
			}, [&] {
				// find global index for suffic in right
				parlay::parallel_for(0, right_length, [&] (index_t i) {
					index_t key = right_arr[i];
					index_t index = binary_search(left_arr, left_length, key);
					result_arr[i + index] = key;
                });
			});
		}

		void binary_parallel_merge(index_t* left_arr, index_t left_length, index_t* right_arr, index_t right_length, index_t* result_arr) {

			parlay::par_do([&] {
				index_t offset = 0;
				for(index_t i = 0; i < left_length; i++) {
					index_t key = left_arr[i];
					index_t index = binary_search(right_arr + offset, right_length - offset, key);
					// insert key into result array
					result_arr[i + index + offset] = key;
					offset += index;
				}
			}, [&] {
				index_t offset = 0;
				for(index_t i = 0; i < right_length; i++) {
					index_t key = right_arr[i];
					index_t index = binary_search(left_arr + offset, left_length - offset, key);
					result_arr[i + index + offset] = key;
					offset += index;
				}
			});
		}

		void serial_merge(index_t* left_arr, index_t left_length, index_t* right_arr, index_t right_length, index_t* result_arr) {
			
			index_t left_i = 0, right_i = 0, result_i = 0;

			while(left_i < left_length && right_i < right_length) {
				if(suffix_less_than(right_arr[right_i], left_arr[left_i])) {
					// right is smaller merge that
					result_arr[result_i++] = right_arr[right_i++];
				} else {
					result_arr[result_i++] = left_arr[left_i++];
				}
			}
			// when the above loop finishes either the left or right array may be 
			// not traversed but not both
			while(left_i < left_length) {
				result_arr[result_i++] = left_arr[left_i++];
			}
			while(right_i < right_length) {
				result_arr[result_i++] = right_arr[right_i++];
			}

		}

		void merge(index_t* left_arr, index_t left_length, index_t* right_arr, index_t right_length, index_t* result_arr) {
			if (merge_type == MergeType::max_parallel) {
                max_parallel_merge(left_arr, left_length, right_arr, right_length, result_arr);
            } else if (merge_type == MergeType::binary_parallel) {
                binary_parallel_merge(left_arr, left_length, right_arr, right_length, result_arr);
            } else {
                serial_merge(left_arr, left_length, right_arr, right_length, result_arr);
            }


			// copy result array back into left array
			index_t total_length = left_length + right_length;

			parlay::parallel_for(0, total_length, [&] (index_t i) {
				left_arr[i] = result_arr[i];
			});
        }
		// N^2 sequential sort
		void very_slow_sort(index_t* array, index_t length) {

			// sort suffixes
			for (index_t i = 0; i < length; i++) {
				
				// find index of min suffix
				index_t min_index = i;
				for (index_t j = i + 1; j < length; j++) {

					index_t suffix_index, curr_min_index; 
					suffix_index = array[j];
					curr_min_index = array[min_index];
					if (suffix_less_than(suffix_index, curr_min_index)) {
						min_index = j;
					}
				}

				// swap least index into current
				index_t temp = array[i];
				array[i] = array[min_index];
				array[min_index] = temp;
			}
		}

        int64_t compare_suffixes(uint64_t genome_index_0, uint64_t genome_index_1) {
			uint64_t s_i_0, s_i_1;
			if (genome_index_0 == genome_index_1) {
				return 0;
			}
			Suffix* s0 = get_suffix(genome_index_0);
			Suffix* s1 = get_suffix(genome_index_1);
			int64_t difference = 0;
			uint64_t num_misses = 0;
			while(difference == 0) {
				difference = s0->hash;
				difference -= s1->hash;
				num_misses++;

				s0 = s0 + 1;
				s1 = s1 + 1;
			}
			num_extra_compares += num_misses - 1;

			return difference;
        }

        bool suffix_less_than(index_t genome_index_0, index_t genome_index_1) {
			if (genome_index_0 == genome_index_1) {
				return 0;
			}
			Suffix* s0;
			Suffix* s1;
			// int64_t difference = 0;
			// bool equal;
			uint64_t num_misses = 0;
			do {
				if(genome_index_0 > genome_length) {
					return  genome_index_0 > genome_index_1;
				}
				s0 = get_suffix(genome_index_0); 
				s1 = get_suffix(genome_index_1);
				// difference = s0->hash;
				// difference -= s1->hash;
				// equal = (s0->hash == s1->hash);
				num_misses++;

				// s0 = s0 + 1;
				// s1 = s1 + 1;
				genome_index_0 += block_size;
				genome_index_1 += block_size;
			} while (s0->hash == s1->hash);
			
			
			// if(num_misses > 1) {
			// 	num_extra_compares += (num_misses - 1);
			// }
			// if(num_misses > 2) {
			// 	num_extra_compares_over_2 += 1;
			// }
			// if(num_misses > 4) {
			// 	num_extra_compares_over_4 += 1;
			// }
			// if(num_misses > 8) {
			// 	num_extra_compares_over_8 += 1;
			// }
			// if(num_misses > 16) {
			// 	num_extra_compares_over_16 += 1;
			// }
			// if(num_misses > 32) {
			// 	num_extra_compares_over_32 += 1;
			// }
			// if(num_misses > 64) {
			// 	num_extra_compares_over_64 += 1;
			// }

			return s0->hash < s1->hash;
        }

        Suffix* get_suffix(index_t genome_index) {

			index_t block_length = genome_length / block_size;
				
			index_t num_over =  genome_length % block_size;
			
			index_t i = genome_index % block_size;

			if(i < num_over) {
				num_over = i;
			}
			
			index_t j = genome_index / block_size;

            Suffix* suffix = suffixes[i] + j;
            // Suffix* suffix = suffixes + (i * block_length + num_over + j);

			// num_over * (block_length + 1) + (i - num_over) * (block_length)
			// = num_over * block_length + num_over + i * block_length - num_over * block_length
			// = num_over + i * block_length

            return suffix;
        }

        SuffixArray(char* fasta_file, uint32_t block_size){
			merge_type = MergeType::max_parallel;
			keep_track_of_extra_comps = false;
            num_extra_compares = 0;
            num_extra_compares_over_2 = 0;
            num_extra_compares_over_4 = 0;
            num_extra_compares_over_8 = 0;
			num_extra_compares_over_16 = 0;
			num_extra_compares_over_32 = 0;
			num_extra_compares_over_64 = 0;
			this->block_size = block_size;
			printf("Starting to read genome\n");
            this->read_fasta(fasta_file);
			printf("Genome lengthh is: %d\n", genome_length);
			if(parlay::num_workers() < std::log2(genome_length)) {
			    merge_type = MergeType::serial;
				std::cout << "Using Serial Merge" << std::endl;
			} else {
				std::cout << "Using Max Parallel Merge" << std::endl;
			}

			uint64_t suffix_sizes = sizeof(Suffix) * genome_length;
			uint64_t sa_size = sizeof(index_t) * genome_length;
			uint64_t total_size = genome_length + suffix_sizes + (sa_size * 2);
			std::cout << "Total estimated memory consuption: " << total_size << " Bytes" << std::endl;
			
			printf("Starting genome to upper\n");
			this->genome_to_upper();

			printf("Starting setup suffixes\n");
			auto setup_start_time = std::chrono::high_resolution_clock::now();
            this->set_up_suffixes();
			auto setup_end_time = std::chrono::high_resolution_clock::now();
			auto set_up_time = std::chrono::duration_cast<std::chrono::milliseconds>(setup_end_time - setup_start_time);
			std::cout << "Set up took " << set_up_time.count() << " (ms)" << std::endl;

			printf("Starting to sort suffixes\n");
			auto sort_start_time = std::chrono::high_resolution_clock::now();
			this->sort_suffixes();
			auto sort_end_time = std::chrono::high_resolution_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(sort_end_time - sort_start_time);
			printf("finished sorting\n");
			printf("That took %llu comparisons\n", num_extra_compares.load());
			printf("That took %llu comparisons\n", num_extra_compares_over_2.load());
			printf("That took %llu comparisons\n", num_extra_compares_over_4.load());
			printf("That took %llu comparisons\n", num_extra_compares_over_8.load());
			printf("That took %llu comparisons\n", num_extra_compares_over_16.load());
			printf("That took %llu comparisons\n", num_extra_compares_over_32.load());
			printf("That took %llu comparisons\n", num_extra_compares_over_64.load());
			// printf("\t%d (ns)\n", time);
			std::cout << time.count() << " (ms)" << std::endl;


        }
	
	public:
		static SuffixArray* from_fasta(std::string& fasta_file, u_int32_t block_size) {
			char* fasta_loc = (char*) fasta_file.data();
			return new SuffixArray(fasta_loc, block_size);
		}

		void print_sa() {
			char* file_location = "influenza-A-out.txt";
			FILE *fptr;
			// Open a file in writing mode
			fptr = fopen(file_location, "w");

			// Write some text to the file
			// fprintf(fptr, "Some text");
			for(uint64_t i = 0; i < genome_length; i++){
				fprintf(fptr, "%llu\n", suffix_array[i]);
				// printf("%llu\n", suffix_array[i]);

			}

			// Close the file
			fclose(fptr);
		}
};