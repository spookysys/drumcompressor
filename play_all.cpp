#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <math.h>
#include <array>
#include "export/use_all.inc"
//#define USE_JK_BD_02
//#define USE_JK_BD_06

static const int block_size = 8;
static const int sample_rate = 22050;
static const int vol_shr = 2;


using namespace std;

struct AdpcmSample
{
	uint16_t len;
	const uint8_t* data;
};



namespace drums
{
	struct DrumEnv
	{
		uint16_t increment;
		uint16_t initial;
	};

	struct DrumFilter
	{
		int8_t a2, a3, b1, b2, b3;
	};

	struct Drum
	{
		DrumEnv treble_env;
		DrumFilter treble_filter;
		AdpcmSample bass_sample;
	};

	enum Drums
	{
		#include "export/enums.inc"
	};

	#include "export/datas.inc"
	const Drum drums[] = {
		#include "export/structs.inc"
	};
	const char* names[] = {
		#include "export/names.inc"
	};
	const int num = sizeof(drums)/sizeof(*drums);
}


class AdpcmDecoder
{
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
	int8_t modifier;
public:
	void reset()
	{
		ptr = end;
	}

	void trigger(const AdpcmSample& sample)
	{
		if (sample.data) {
			this->ptr = sample.data + 1;
			this->end = sample.data + sample.len + 1;
		} else {
			this->ptr = nullptr;
			this->end = nullptr;
		}
		if (isActive()) {
			this->word = *sample.data;
			this->subidx = 4;
		}
		this->recon_1 = 0;
		this->recon_2 = 0;
		this->modifier = 0;
	}

	bool isActive()
	{
		return ptr!=end;
	}	

	int8_t get()
	{
		uint8_t idx = getIdx();
		switch (idx) {
			case 0: break;
			case 1: 
				modifier = -modifier;
				break;
			case 2: 
				modifier = modifier ? (modifier<<1) : 1;
				break;
			case 3:
				modifier = modifier ? (modifier>>1) : -1;
				break;
		}
		int8_t prediction = recon_1 + recon_1 - recon_2;
		int8_t recon = prediction + modifier;
		this->recon_2 = this->recon_1;
		this->recon_1 = recon;
		return recon;
	}

};

class Filter
{
	float a2, a3, b1, b2, b3;
	float xn1, xn2;
	float yn1, yn2;
public:
	void init(const drums::DrumFilter& op) 
	{
		a2 = op.a2 / 64.f;
		a3 = op.a3 / 128.f;
		b1 = op.b1 / 256.f;
		b2 = op.b2 / 256.f;
		b3 = op.b3 / 128.f;
		xn1 = 0;
		yn1 = 0;
		xn2 = 0;
		yn2 = 0;
	}
	int16_t get(int16_t op)
	{
		float xn = op;

		float yn = b1*xn + b2*xn1 + b3*xn2 - a2*yn1 - a3*yn2;

		this->xn2 = this->xn1;
		this->yn2 = this->yn1;
		this->xn1 = xn;
		this->yn1 = yn;

		return yn;
	}
};

class DrumDecoder
{
	AdpcmDecoder adpcmDecoder;

	int8_t bass_val;
	int16_t bass_filtered;

	int32_t treble_amplitude;
	uint16_t treble_increment;

	Filter treble_filter;
public:

	void trigger(const drums::Drum& drum)
	{
		adpcmDecoder.trigger(drum.bass_sample);
		if (adpcmDecoder.isActive()) bass_val = adpcmDecoder.get();
		bass_filtered = 0;
		treble_amplitude = uint32_t(drum.treble_env.initial) << 8;
		treble_increment = drum.treble_env.increment;
		treble_filter.init(drum.treble_filter);
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
				bass_filtered = ((bass_filtered+(rand()&1))>>1) + ((bass_inter+(rand()&1))>>1);
				bass_inter += bass_delta;
				dest[i] += bass_filtered >> vol_shr;
			}
		}

		// treble
		if (this->treble_amplitude > 0) {
			uint16_t amplitude = this->treble_amplitude >> 8;

			for (int i=0; i<block_size; i++) {
				int32_t noiz = int16_t(rand() | (rand() << 8));
				noiz = int32_t(noiz * amplitude) >> 8;
				noiz >>= vol_shr;
				int16_t val = treble_filter.get(noiz);
				dest[i] += val;
			}			

			this->treble_amplitude -= this->treble_increment;
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
	uint32_t bytes_per_second = 0;
	uint16_t bytes_per_frame = 0;
	uint16_t bits_per_sample = 0;
	char subchunk2_id_0 = 'd';
	char subchunk2_id_1 = 'a';
	char subchunk2_id_2 = 't';
	char subchunk2_id_3 = 'a';
	uint32_t subchunk2_size = 0;
};


static WavHeader wav_header;
void write_wav(const char* filename, int16_t* data, int num_samples, bool only_8=false)
{
	const int bytes_per_frame = only_8 ? 1 : 2;
	wav_header.chunk_size = bytes_per_frame * num_samples + 36;
	wav_header.subchunk2_size = bytes_per_frame * num_samples - 8;
	wav_header.bytes_per_second = bytes_per_frame * wav_header.sample_rate;
	wav_header.bytes_per_frame = bytes_per_frame;
	wav_header.bits_per_sample = bytes_per_frame * 8;
	FILE* fid = fopen(filename, "wb");
	fwrite(&wav_header, 1, sizeof(WavHeader), fid);
	if (only_8) {
		int8_t* wdata = (int8_t*)malloc(num_samples);
		for (int i=0; i<num_samples; i++) {
			int16_t val = data[i];
			int limit = 1<<(8-vol_shr);
			int dither = rand() % limit;
			val += dither;
			val >>= (8-vol_shr);
			val += 128;
			if (val<0) val = 0;
			if (val>255) val = 255;
			wdata[i] = val;
		}
		fwrite(wdata, 1, num_samples, fid);
		free(wdata);
	} else {
		fwrite(data, 2, num_samples, fid);
	}
	fclose(fid);
}

void write_dat(const char* filename, const drums::Drum& op)
{
	uint16_t bass_data_ptr_placeholder = 0;
	FILE* fid = fopen(filename, "wb");
	fwrite(&op.treble_env, sizeof(op.treble_env), 1, fid);
	fwrite(&op.treble_filter, sizeof(op.treble_filter), 1, fid);
	fwrite(&op.bass_sample.len, sizeof(op.bass_sample.len), 1, fid);
	fwrite(&bass_data_ptr_placeholder, sizeof(bass_data_ptr_placeholder), 1, fid);
	fwrite(op.bass_sample.data, 1, op.bass_sample.len, fid);
	fclose(fid);
}


int main(int argc, char* argv[])
{
	static const int num_blocks = (sample_rate * 3) / block_size;
	std::array<int16_t, num_blocks*block_size> buffer;

	DrumDecoder decoder;
	for (int drum_i=0; drum_i<drums::num; drum_i++) {
		const drums::Drum& drum = drums::drums[drum_i];
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
			write_wav(filename, buffer.data(), num_blocks*block_size, true);
			d[0] = '.';
			d[1] = 'd';
			d[2] = 'a';
			d[3] = 't';
			d[4] = '\0';
			write_dat(filename, drum);
		}
	}
}

#ifdef CRINKLE
#define WINDOWS_LEAN_AND_MEAN
#include <windows.h>
extern "C" void __stdcall mainCRTStartup()
{
	main(0, nullptr);
	ExitProcess(0);
}
#endif
