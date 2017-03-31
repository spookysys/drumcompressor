#include <cstdint>
#include <array>
#include <fcntl.h>
#include "export/use_all.inc"
//#define USE_JK_BD_02
//#define USE_JK_BD_06

static const int block_size = 8;
static const int sample_rate = 22050;
static const bool adpcm_dithering = false;

using namespace std;

struct DrumExpEnv
{
	int8_t exponent;
	int8_t magnitude0, magnitude1;
	uint8_t falloff0, falloff1;
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
	DrumExpEnv bass_env;
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

struct ActiveEnv
{
	int16_t magnitude0, magnitude1;
	uint8_t falloff0, falloff1;

	static const uint8_t precision = 4;
	void trigger(const DrumExpEnv& op)
	{
		uint8_t mag_shl = precision + 8 - op.exponent;
		if (mag_shl < 0) mag_shl = 0;
		this->magnitude0 = int16_t(op.magnitude0) << mag_shl;
		this->magnitude1 = int16_t(op.magnitude1) << mag_shl;
		this->falloff0 = op.falloff0;
		this->falloff1 = op.falloff1;
	}

	int16_t get()
	{
		int16_t ret = (this->magnitude0 + this->magnitude1) >> precision;
		const uint8_t decay_booster = 2;
		const uint16_t bias = 1 << (7+decay_booster);
		this->magnitude0 -= (int32_t(this->magnitude0) * falloff0) >> (8+decay_booster);
		this->magnitude1 -= (int32_t(this->magnitude1) * falloff1) >> (8+decay_booster);
		return ret;
	}
};

class DrumDecoder
{
	AdpcmDecoder adpcmDecoder;
	int16_t bass_val;
	
	ActiveEnv bass_env;
	
	int16_t get_bass_val()
	{
		int16_t adpcm = adpcmDecoder.get();
		int16_t volume = bass_env.get();
		int32_t ret = adpcm * volume;
		return ret;
	}
public:

	void trigger(const DrumExp& drum)
	{
		adpcmDecoder.trigger(drum.bass_sample);
		bass_env.trigger(drum.bass_env);
		bass_val = get_bass_val();
	}

	void decodeBlock(int16_t* dest)
	{
		if (adpcmDecoder.isActive()) {
			int16_t prev_bass_val = this->bass_val;
			this->bass_val = get_bass_val();

			for (int i=0; i<block_size; i++) {
				int16_t v0 = prev_bass_val >> 4;
				int16_t v1 = this->bass_val >> 4;
				int8_t dither_val = 0;
				if (adpcm_dithering) {
					dither_val = ((rand() % 3) - 1) << 4;
				}
				dest[i] = ((v0 * (block_size-i) + v1 * i + dither_val) / block_size) >> 4;
			}
		} else {
			for (int i=0; i<block_size; i++) {
				dest[i] = 0;
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
	uint32_t bytes_per_second = 22050;
	uint16_t bytes_per_frame = 1;
	uint16_t bits_per_sample = 8;
	char subchunk2_id_0 = 'd';
	char subchunk2_id_1 = 'a';
	char subchunk2_id_2 = 't';
	char subchunk2_id_3 = 'a';
	uint32_t subchunk2_size = 0;
};


void write_wav(const char* filename, int8_t* data, int num_samples)
{
	const int wavheader_size = sizeof(WavHeader);
	static WavHeader wav_header;
	wav_header.chunk_size = num_samples + 36;//sizeof(WavHeader)-8;
	wav_header.subchunk2_size = num_samples;
	int fid = open(filename, O_WRONLY|O_CREAT|O_TRUNC);
	write(fid, &wav_header, sizeof(WavHeader));
	write(fid, data, num_samples);
	close(fid);
}


int main(int argc, char* argv[])
{
	static const int num_blocks = (sample_rate * 2) / block_size;
	std::array<int8_t, num_blocks*block_size> buffer;

	DrumDecoder decoder;
	for (int drum_i=0; drum_i<drums::num; drum_i++) {
		const DrumExp& drum = drums::drums[drum_i];
		const char* name = drums::names[drum_i];
		decoder.trigger(drum);
		for (int block_i=0; block_i<num_blocks; block_i++) {
			std::array<int16_t, block_size> block_buffer;
			decoder.decodeBlock(block_buffer.data());
			std::array<int8_t, block_size> block_buffer_clipped;
			for (int i=0; i<block_size; i++) {
				int16_t val = block_buffer[i];
				if (val > 127) val = 127;
				if (val < -128) val = -128;
				val = val + 128;
				block_buffer_clipped[i] = val;
			}
			for (int i=0; i<block_size; i++) {
				buffer[block_i*block_size + i] = block_buffer_clipped[i];
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