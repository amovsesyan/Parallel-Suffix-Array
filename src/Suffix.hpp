class Suffix {
	private:
		char* suffix;
		uint32_t block_size;

		void calculate_hash() {
			hash = 0;
            # pragma omp parallel for
			for (int i = 0; i < this->block_size; i++) {
				hash *= 4;    
				switch (suffix[i]){
				case 'A':
					// value stays the same (+= 0)
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

				default:
					break;
				}
			}
		}
   
	public:
		uint32_t hash;
		Suffix(char* suffix, uint32_t block_size) {
			this->suffix = suffix;
			this->block_size = block_size;
            this->hash = 0;
            this->calculate_hash();
		}


		
};