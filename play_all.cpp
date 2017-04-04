#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <math.h>
#include <array>
#include <algorithm>
#include "export/use_all.inc"
//#define USE_JK_BD_02
//#define USE_JK_BD_06

static const int block_size = 8;
static const int sample_rate = 22050;

using namespace std;

struct DrumExpEnv
{
	uint16_t amplitude;
	uint16_t falloff;
};

struct DrumExpFilter
{
	int8_t a2, a3, b1, b2, b3;
};

struct AdpcmSample
{
	std::array<int8_t, 4> palette;
	uint16_t len;
	const uint8_t* data;
};

struct DrumExp
{
	DrumExpEnv treble_env;
	DrumExpFilter treble_filter;
	AdpcmSample bass_sample;
};


namespace drums
{
	enum Drums
	{
		#include "export/enums.inc"
	};

	#include "export/datas.inc"
	DrumExp drums[] = {
		#include "export/structs.inc"
	};
	const char* names[] = {
		#include "export/names.inc"
	};
	const int num = sizeof(drums)/sizeof(*drums);
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

	int8_t recon_1;
	int8_t recon_2;
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
			this->word = *sample.data;
			this->subidx = 4;
		}
		this->recon_1 = 0;
		this->recon_2 = 0;
	}

	bool isActive()
	{
		return ptr!=end;
	}	

	int8_t get()
	{
		uint8_t idx = getIdx();
		int8_t palette_val = palette[idx];
		int8_t prediction = recon_1 + recon_1 - recon_2;
		int8_t recon = prediction; 
		if (recon_1 < 0) {
			recon -= palette_val;
		} else {
			recon += palette_val;
		}
		this->recon_2 = this->recon_1;
		this->recon_1 = recon;
		return recon;
	}

};


class DrumDecoder
{
	AdpcmDecoder adpcmDecoder;

	int8_t bass_val;
	int16_t bass_filter;

	uint16_t treble_amplitude;
	uint16_t treble_falloff;
public:

	void trigger(const DrumExp& drum)
	{
		adpcmDecoder.trigger(drum.bass_sample);
		bass_val = adpcmDecoder.get();
		bass_filter = 0;
		treble_amplitude = drum.treble_env.amplitude;
		treble_falloff = drum.treble_env.falloff;
	}

	void decodeBlock(int16_t* dest)
	{
		for (int i=0; i<block_size; i++) {
			dest[i] = 0;
		}

		// bass
		if (adpcmDecoder.isActive()) {
			int8_t bass_first = this->bass_val;
			this->bass_val = adpcmDecoder.get();
			int8_t bass_last = this->bass_val;

			int16_t bass_inter = int16_t(bass_first) << 8;
			int16_t bass_delta = int16_t(bass_last - bass_first) << (8-3);
			for (int i=0; i<block_size; i++) {
				bass_filter = (bass_filter>>1) + (bass_inter>>1);
				bass_inter += bass_delta;
				dest[i] += bass_filter;
			}
		}

		// treble
		{
			uint16_t amplitude = this->treble_amplitude >> 6;
			this->treble_amplitude = (uint32_t(this->treble_amplitude) * uint32_t(uint16_t(-1 - this->treble_falloff))) >> 16;

			for (int i=0; i<block_size; i++) {
				int32_t noiz = int16_t(rand() | (rand() << 8));
				noiz = int32_t(noiz * amplitude) >> 8;
				dest[i] += noiz;
			}			
		}
	}
};

struct WavHeader
{
	char chunk_id_0 = 'R';
	char chunk_id_1 = 'I';
	char chunk_id_2 = 'F';
	char chunk_id_3 = 'F';
	uint32_t chunk_size = 0;
	char format_0 = 'W';
	char format_1 = 'A';
	char format_2 = 'V';
	char format_3 = 'E';
	char subchunk1_id_0 = 'f';
	char subchunk1_id_1 = 'm';
	char subchunk1_id_2 = 't';
	char subchunk1_id_3 = 0x20;
	uint32_t subchunk1_size = 16;
	uint16_t audio_format = 1; // 1 = PCM
	uint16_t num_channels = 1;
	uint32_t sample_rate = 22050; 
	uint32_t bytes_per_second = 2*22050;
	uint16_t bytes_per_frame = 2;
	uint16_t bits_per_sample = 16;
	char subchunk2_id_0 = 'd';
	char subchunk2_id_1 = 'a';
	char subchunk2_id_2 = 't';
	char subchunk2_id_3 = 'a';
	uint32_t subchunk2_size = 0;
};


void write_wav(const char* filename, int16_t* data, int num_samples)
{
	const int wavheader_size = sizeof(WavHeader);
	static WavHeader wav_header;
	wav_header.chunk_size = 2*num_samples + 36;//sizeof(WavHeader)-8;
	wav_header.subchunk2_size = 2*num_samples - 8;
	FILE* fid = fopen(filename, "wb");
	fwrite(&wav_header, 1, sizeof(WavHeader), fid);
	fwrite(data, 2, num_samples, fid);
	fclose(fid);
}


int main(int argc, char* argv[])
{
	static const int num_blocks = (sample_rate * 3) / block_size;
	std::array<int16_t, num_blocks*block_size> buffer;

	DrumDecoder decoder;
	for (int drum_i=0; drum_i<drums::num; drum_i++) {
		const DrumExp& drum = drums::drums[drum_i];
		const char* name = drums::names[drum_i];
		decoder.trigger(drum);
		for (int block_i=0; block_i<num_blocks; block_i++) {
			std::array<int16_t, block_size> block_buffer;
			decoder.decodeBlock(block_buffer.data());
			for (int i=0; i<block_size; i++) {
				int16_t val = block_buffer[i];
				val = val * (32760. / 32768.);
				buffer[block_i*block_size + i] = val;
			}
		}
		{
			static char filename[256];
			char* d = filename;
			const char* s = name;
			do {
				*d++ = *s++;
			} while (*s);
			d[0] = '.';
			d[1] = 'w';
			d[2] = 'a';
			d[3] = 'v';
			d[4] = '\0';
			write_wav(filename, buffer.data(), buffer.size());
		}
	}
}