#include <cstdint>
#include <array>
//#include "export/use_all.inc"
#define USE_JK_BD_02
#define USE_JK_BD_06

using namespace std;

struct DrumEnv
{
	int8_t exponent, magnitude0, magnitude1;
	uint8_t falloff0, falloff1;
};

struct DrumFilter
{
	int8_t exponent, b1, b2, a1, a2, a3;
};

struct Drum
{
	DrumEnv treble_env;
	DrumFilter treble_filter;
	DrumEnv bass_env;
	std::array<int8_t, 4> bass_palette;
	uint16_t bass_len;
	const uint8_t* bass_data;
};

namespace drums
{
	enum Drums
	{
		#include "export/enums.inc"
	};

	#include "export/datas.inc"
	Drum drums[] = {
		#include "export/structs.inc"
		/*
		{
			{ 6, 69, 1, 255, 255 }, // treble envelope
			{ 6, 10, -19, 64, -114, 57 }, // treble filter
			{ 7, 106, 1, 249, 255 }, // bass envelope
			{ -28, 8, -8, -73 }, // bass palette
			clap_808_data // bass data
		}
		*/
	};
}


class AdpcmDecoder
{
	const uint8_t* ptr;
	const uint8_t* end;
	array<int8_t, 4> palette;

	uint8_t word;
	uint8_t subidx;
	int8_t getIdx()
	{
		if (--subidx) {
			word >>= 2;
		} else if (ptr!=end) {
			subidx = 4;
			ptr++;
			word = *ptr;
		}
		return word & 3;
	}

public:
	void trigger(const uint8_t* data, uint16_t len, const array<int8_t, 4>& palette)
	{
		this->ptr = data;
		this->end = data+len;
		this->palette = palette;
	}

	void next()
	{
		idx = getIdx();
		pal = palette[idx];
        prediction = recon_1 + recon_1 - recon_2;
        if (recon_1 < 0)
            recon = wrap8(prediction - palette_val);
        else
            recon = wrap8(prediction + palette_val);
        end
        data_out(i) = recon;
        recon_2 = recon_1;
        recon_1 = recon;
    end
		
	}

	int8_t get()
	{
		return val;
	}

};

class DrumDecoder
{
	uint8_t bass_data_idx;
	uint8_t bass_data_byte;
	const uint8_t* bass_data_ptr;
	const uint8_t* bass_data_end;
	array<int8_t, 4> bass_palette;
	int8_t bass_val;
public:

	int8_t bass_data_iter()
	{
	}

	void trigger(const Drum& drum)
	{
		this->bass_palette = drum.bass_palette;

		this->bass_data_ptr = drum.bass_data;
		this->bass_data_end = this->bass_data_ptr + drums.bass_len;
		this->bass_data_idx = this->bass_data;

		// state
		this->bass_val = 0;
	}

	void decodeBlock(int16_t* dest)
	{
		if (bass_data_pos != bass_data_end) {
			if (--bass_data_subidx == 0) {

			}
		}
			

		
	}
};

void write_wav(const char* filename, int16_t* data, int num_samples, int sample_rate)
{
}


int main(int argc, char* argv[])
{
	static const int sample_rate = 22050;
	static const int block_size = 8;
	static const int num_blocks = sample_rate * 2 / block_size;
	std::array<int16_t, num_blocks*block_size> buffer;

	DrumDecoder decoder;
	for (const auto iter : drums::drums) 
	{
		decoder.trigger(iter);
		for (int i=0; i<num_blocks; i++) {
			decoder.decodeBlock(&buffer[i*block_size]);
		}
	}
}