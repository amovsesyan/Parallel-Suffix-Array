
#include "SuffixArray.hpp"


int main() {
    std::string file_location = "../input/fasta_files/influenza_A.fna";
    // std::string file_location = "../input/fasta_files/salmonella.fa";
    uint32_t block_size = 27;
    SuffixArray* sa = SuffixArray::from_fasta(file_location, block_size);
    sa->print_sa();
}
