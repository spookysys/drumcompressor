#include <cstdint>
#include <array>
//#include "export/use_all.inc"
#define USE_JK_BD_02
#define USE_JK_BD_06

using namespace std;

struct DrumEnv
{
	int8_t exponent;
	int8_t magnitude0, magnitude1;
	uint8_t falloff0, falloff1;
};

struct DrumFilter
{
	int8_t a2, a3, b1, b2, b3;
};

struct AdpcmSample
{
	std::array<int8_t, 4> palette;
	uint16_t len;
	const uint8_t* data;
};

struct Drum
{
	DrumEnv treble_env;
	DrumFilter treble_filter;
	DrumEnv bass_env;
	AdpcmSample bass_sample;
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
	};
}


class AdpcmDecoder
{
	array<int8_t, 4> palette;
	const uint8_t* ptr;
	const uint8_t* end;
	uint8_t subidx;
	uint8_t word;

	uint8_t getIdx()
	{
		uint8_t ret = word & 3;
		word >>= 2;
		subidx--;
		if (subidx==0) {
			subidx = 4;
			word = ((ptr+1)!=end) ? *ptr : 0;
			ptr++;
		}
		return ret;
	}

	int8_t recon1;
	int8_t recon2;
public:
	void reset()
	{
		ptr = end;
	}

	void trigger(const AdpcmSample& sample)
	{
		this->palette = sample.palette;
		this->ptr = sample.data + 1;
		this->end = sample.data + sample.len + 1;
		if (isActive()) {
			this->word = *ptr;
			this->subidx = 4;
		}
		this->recon1 = 0;
		this->recon2 = 0;
	}

	bool isActive()
	{
		return ptr!=end;
	}	

	int8_t get()
	{
		uint8_t idx = getIdx();
		int8_t palval = palette[idx];
		int8_t prediction = recon_1 + recon_1 - recon_2;
		int8_t recon = prediction; 
		if (recon_1 < 0) {
			recon -= palette_val;
		} else {
			recon += palette_val;
		}
		this->recon_2 = this->recon_1;
		this->recon_1 = this->recon;
		return recon;
	}

};


class DrumDecoder
{
	AdpcmDecoder adpcmDecoder;

public:

	void trigger(const Drum& drum)
	{
		adpcmDecoder.trigger(drum.bass_sample);
		adpcmDecoder.next();
	}

	void decodeBlock(int16_t* dest)
	{
		if (adpcmDecoder.isActive()) {
			int8_t bass_last = this->bass_val;
			this->bass_val = adpcmDecoder.get();

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