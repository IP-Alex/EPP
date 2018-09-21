#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <string>
#include <tiff.h>
#include <tiffio.h>
#include "custom_types.h"
#include "xml_aux.h"
#include "dct8x8_block.h"
#include <vector>
#include "gettimeofday.h"
#include "config.h"
#include <immintrin.h>
#include "opencl_helper.h"
  
using namespace std;

#define print false

void loadImage(int number, string path, Image** photo) {
    string filename;
    TIFFRGBAImage img;
    char emsg[1024];
    
	filename = path + to_string(number) + ".tiff";
    TIFF* tif = TIFFOpen(filename.c_str(), "r");
    if(tif==NULL) fprintf(stderr,"Failed opening image: %s\n", filename);
    if (!(TIFFRGBAImageBegin(&img, tif, 0, emsg))) TIFFError(filename.c_str(), emsg);
     
    uint32 w, h;
    size_t npixels;
    uint32* raster;
     
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    npixels = w * h;
    raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
     
    TIFFReadRGBAImage(tif, w, h, raster, 0);
     
    if(*photo==NULL) 
	*photo = new Image((int)w, (int)h, FULLSIZE);
     
    //Matlab and LibTIFF store the image diferently.
    //Necessary to mirror the image horizonatly to be consistent
	/*
		IMPROVEMENT: Read in one CL at a time for each array. 
	*/
    for (int j=0; j<(int)w; j++) {
        for (int i=0; i<(int)h; i++) {  
            // The inversion is ON PURPOSE
            (*photo)->rc->data[(h-1-i)*w+j] = (float)TIFFGetR(raster[i*w+j]);
            (*photo)->gc->data[(h-1-i)*w+j] = (float)TIFFGetG(raster[i*w+j]);
            (*photo)->bc->data[(h-1-i)*w+j] = (float)TIFFGetB(raster[i*w+j]);
        }
    }
     
    _TIFFfree(raster);
    TIFFRGBAImageEnd(&img);
    TIFFClose(tif);
     
}
#define gpu 1
#define normal 2
#define mode normal
#define _CORES 1

// in_r, in_g, in_b, out_r, out_g, out_b
void convertRGBtoYCbCr(Image* in, Image* out) {

    int width = in->width;
    int height = in->height;

		/*
			IMPROVEMENT: Flip x and y loop. 
			Look for conflicts between R,G and B CLs.
			SIMD: RAKT på!
		*/
#ifdef SIMD
		const __m256 Y299 = _mm256_set1_ps(0.299);
		const __m256 Y587 = _mm256_set1_ps(0.587);
		const __m256 Y113 = _mm256_set1_ps(0.113);

		const __m256 Cb168736 = _mm256_set1_ps(0.168736);
		const __m256 Cb331264 = _mm256_set1_ps(0.331264);
		const __m256 Cb5 = _mm256_set1_ps(0.5);

		const __m256 Cr5 = _mm256_set1_ps(0.5);
		const __m256 Cr418688 = _mm256_set1_ps(0.418688);
		const __m256 Cr081312 = _mm256_set1_ps(0.081312);

		const __m256 _128 = _mm256_set1_ps(128);

		for(int y=0; y<width; y+=8) {





			//float R = in->rc->data[x*width+y];
			//float G = in->gc->data[x*width+y];
			//float B = in->bc->data[x*width+y];
			//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
			//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
			//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
			__m256 red, green, blue, resultY, resultY2, resultY3, resultCb, resultCb2, resultCb3, resultCr, resultCr2, resultCr3;

			for (int x = 0; x<height; x++) {
				
				red = _mm256_load_ps(&in->rc->data[x*width + y]);
				green = _mm256_load_ps(&in->gc->data[x*width + y]);
				blue = _mm256_load_ps(&in->bc->data[x*width + y]);

				//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
				resultY = _mm256_mul_ps(red, Y299);
				resultY2 = _mm256_mul_ps(green, Y587);
				resultY3 = _mm256_mul_ps(blue, Y113);
				resultY = _mm256_add_ps(resultY, resultY2);
				resultY = _mm256_add_ps(resultY, resultY3);


				//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
				resultCb = _mm256_mul_ps(red, Cb168736);
				resultCb2 = _mm256_mul_ps(green, Cb331264);
				resultCb3 = _mm256_mul_ps(blue, Cb5);
				resultCb = _mm256_sub_ps(_128, resultCb);
				resultCb = _mm256_sub_ps(resultCb, resultCb2);
				resultCb = _mm256_add_ps(resultCb, resultCb3);
			




				//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
				resultCr = _mm256_mul_ps(red, Cr5);
				resultCr2 = _mm256_mul_ps(green, Cr418688);
				resultCr3 = _mm256_mul_ps(blue, Cr081312);
				resultCr = _mm256_sub_ps(resultCr, resultCr2);
				resultCr = _mm256_sub_ps(resultCr, resultCr3);
				resultCr = _mm256_add_ps(resultCr, _128);


				/*
				r[1,2,3,4,5,6,7,8]
				*
				const[1..8]
				=
				r'[1..8]
				+
				g'[1..8]
				mm_mul()
				mm_add(R1*const, )
				*/

				_mm256_store_ps(&out->rc->data[x*width + y], resultY);
				_mm256_store_ps(&out->gc->data[x*width + y], resultCb);
				_mm256_store_ps(&out->bc->data[x*width + y], resultCr);
#endif

#ifdef PSIMD
				const __m256 Y299 = _mm256_set1_ps(0.299);
				const __m256 Y587 = _mm256_set1_ps(0.587);
				const __m256 Y113 = _mm256_set1_ps(0.113);

				const __m256 Cb168736 = _mm256_set1_ps(0.168736);
				const __m256 Cb331264 = _mm256_set1_ps(0.331264);
				const __m256 Cb5 = _mm256_set1_ps(0.5);

				const __m256 Cr5 = _mm256_set1_ps(0.5);
				const __m256 Cr418688 = _mm256_set1_ps(0.418688);
				const __m256 Cr081312 = _mm256_set1_ps(0.081312);

				const __m256 _128 = _mm256_set1_ps(128);

				#pragma omp parallel num_threads(_CORES) for
				for (int y = 0; y<width; y += 8) {





					//float R = in->rc->data[x*width+y];
					//float G = in->gc->data[x*width+y];
					//float B = in->bc->data[x*width+y];
					//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
					//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
					//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
					__m256 red, green, blue, resultY, resultY2, resultY3, resultCb, resultCb2, resultCb3, resultCr, resultCr2, resultCr3;

					for (int x = 0; x<height; x++) {

						red = _mm256_load_ps(&in->rc->data[x*width + y]);
						green = _mm256_load_ps(&in->gc->data[x*width + y]);
						blue = _mm256_load_ps(&in->bc->data[x*width + y]);

						//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
						resultY = _mm256_mul_ps(red, Y299);
						resultY2 = _mm256_mul_ps(green, Y587);
						resultY3 = _mm256_mul_ps(blue, Y113);
						resultY = _mm256_add_ps(resultY, resultY2);
						resultY = _mm256_add_ps(resultY, resultY3);


						//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
						resultCb = _mm256_mul_ps(red, Cb168736);
						resultCb2 = _mm256_mul_ps(green, Cb331264);
						resultCb3 = _mm256_mul_ps(blue, Cb5);
						resultCb = _mm256_sub_ps(_128, resultCb);
						resultCb = _mm256_sub_ps(resultCb, resultCb2);
						resultCb = _mm256_add_ps(resultCb, resultCb3);





						//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
						resultCr = _mm256_mul_ps(red, Cr5);
						resultCr2 = _mm256_mul_ps(green, Cr418688);
						resultCr3 = _mm256_mul_ps(blue, Cr081312);
						resultCr = _mm256_sub_ps(resultCr, resultCr2);
						resultCr = _mm256_sub_ps(resultCr, resultCr3);
						resultCr = _mm256_add_ps(resultCr, _128);


						/*
						r[1,2,3,4,5,6,7,8]
						*
						const[1..8]
						=
						r'[1..8]
						+
						g'[1..8]
						mm_mul()
						mm_add(R1*const, )
						*/

						_mm256_store_ps(&out->rc->data[x*width + y], resultY);
						_mm256_store_ps(&out->gc->data[x*width + y], resultCb);
						_mm256_store_ps(&out->bc->data[x*width + y], resultCr);
#endif

#ifdef PCS
						const __m256 Y299 = _mm256_set1_ps(0.299);
						const __m256 Y587 = _mm256_set1_ps(0.587);
						const __m256 Y113 = _mm256_set1_ps(0.113);

						const __m256 Cb168736 = _mm256_set1_ps(0.168736);
						const __m256 Cb331264 = _mm256_set1_ps(0.331264);
						const __m256 Cb5 = _mm256_set1_ps(0.5);

						const __m256 Cr5 = _mm256_set1_ps(0.5);
						const __m256 Cr418688 = _mm256_set1_ps(0.418688);
						const __m256 Cr081312 = _mm256_set1_ps(0.081312);

						const __m256 _128 = _mm256_set1_ps(128);

						#pragma omp parallel num_threads(_CORES) for
						for (int x = 0; x<height; x++) {
							//float R = in->rc->data[x*width+y];
							//float G = in->gc->data[x*width+y];
							//float B = in->bc->data[x*width+y];
							//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
							//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
							//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
							__m256 red, green, blue, resultY, resultY2, resultY3, resultCb, resultCb2, resultCb3, resultCr, resultCr2, resultCr3;

							for (int y = 0; y<width; y+=8) {

								red = _mm256_load_ps(&in->rc->data[x*width + y]);
								green = _mm256_load_ps(&in->gc->data[x*width + y]);
								blue = _mm256_load_ps(&in->bc->data[x*width + y]);

								//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
								resultY = _mm256_mul_ps(red, Y299);
								resultY2 = _mm256_mul_ps(green, Y587);
								resultY3 = _mm256_mul_ps(blue, Y113);
								resultY = _mm256_add_ps(resultY, resultY2);
								resultY = _mm256_add_ps(resultY, resultY3);


								//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
								resultCb = _mm256_mul_ps(red, Cb168736);
								resultCb2 = _mm256_mul_ps(green, Cb331264);
								resultCb3 = _mm256_mul_ps(blue, Cb5);
								resultCb = _mm256_sub_ps(_128, resultCb);
								resultCb = _mm256_sub_ps(resultCb, resultCb2);
								resultCb = _mm256_add_ps(resultCb, resultCb3);





								//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
								resultCr = _mm256_mul_ps(red, Cr5);
								resultCr2 = _mm256_mul_ps(green, Cr418688);
								resultCr3 = _mm256_mul_ps(blue, Cr081312);
								resultCr = _mm256_sub_ps(resultCr, resultCr2);
								resultCr = _mm256_sub_ps(resultCr, resultCr3);
								resultCr = _mm256_add_ps(resultCr, _128);


								/*
								r[1,2,3,4,5,6,7,8]
								*
								const[1..8]
								=
								r'[1..8]
								+
								g'[1..8]
								mm_mul()
								mm_add(R1*const, )
								*/

								_mm256_store_ps(&out->rc->data[x*width + y], resultY);
								_mm256_store_ps(&out->gc->data[x*width + y], resultCb);
								_mm256_store_ps(&out->bc->data[x*width + y], resultCr);
#endif

#ifdef CS
								const __m256 Y299 = _mm256_set1_ps(0.299);
								const __m256 Y587 = _mm256_set1_ps(0.587);
								const __m256 Y113 = _mm256_set1_ps(0.113);

								const __m256 Cb168736 = _mm256_set1_ps(0.168736);
								const __m256 Cb331264 = _mm256_set1_ps(0.331264);
								const __m256 Cb5 = _mm256_set1_ps(0.5);

								const __m256 Cr5 = _mm256_set1_ps(0.5);
								const __m256 Cr418688 = _mm256_set1_ps(0.418688);
								const __m256 Cr081312 = _mm256_set1_ps(0.081312);

								const __m256 _128 = _mm256_set1_ps(128);

								for (int x = 0; x<height; x++) {
									//float R = in->rc->data[x*width+y];
									//float G = in->gc->data[x*width+y];
									//float B = in->bc->data[x*width+y];
									//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
									//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
									//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
									__m256 red, green, blue, resultY, resultY2, resultY3, resultCb, resultCb2, resultCb3, resultCr, resultCr2, resultCr3;

									for (int y = 0; y<width; y += 8) {

										red = _mm256_load_ps(&in->rc->data[x*width + y]);
										green = _mm256_load_ps(&in->gc->data[x*width + y]);
										blue = _mm256_load_ps(&in->bc->data[x*width + y]);

										//float Y = 0+((float)0.299*R)+((float)0.587*G)+((float)0.113*B);
										resultY = _mm256_mul_ps(red, Y299);
										resultY2 = _mm256_mul_ps(green, Y587);
										resultY3 = _mm256_mul_ps(blue, Y113);
										resultY = _mm256_add_ps(resultY, resultY2);
										resultY = _mm256_add_ps(resultY, resultY3);


										//float Cb = 128-((float)0.168736*R)-((float)0.331264*G)+((float)0.5*B);
										resultCb = _mm256_mul_ps(red, Cb168736);
										resultCb2 = _mm256_mul_ps(green, Cb331264);
										resultCb3 = _mm256_mul_ps(blue, Cb5);
										resultCb = _mm256_sub_ps(_128, resultCb);
										resultCb = _mm256_sub_ps(resultCb, resultCb2);
										resultCb = _mm256_add_ps(resultCb, resultCb3);





										//float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
										resultCr = _mm256_mul_ps(red, Cr5);
										resultCr2 = _mm256_mul_ps(green, Cr418688);
										resultCr3 = _mm256_mul_ps(blue, Cr081312);
										resultCr = _mm256_sub_ps(resultCr, resultCr2);
										resultCr = _mm256_sub_ps(resultCr, resultCr3);
										resultCr = _mm256_add_ps(resultCr, _128);


										/*
										r[1,2,3,4,5,6,7,8]
										*
										const[1..8]
										=
										r'[1..8]
										+
										g'[1..8]
										mm_mul()
										mm_add(R1*const, )
										*/

										_mm256_store_ps(&out->rc->data[x*width + y], resultY);
										_mm256_store_ps(&out->gc->data[x*width + y], resultCb);
										_mm256_store_ps(&out->bc->data[x*width + y], resultCr);
#endif

#ifdef openMP
#pragma omp parallel num_threads(1) for
	for (int y = 0; y<width; y++) {
		for (int x = 0; x<height; x++) {
			float R = in->rc->data[x*width + y];
			float G = in->gc->data[x*width + y];
			float B = in->bc->data[x*width + y];
			float Y = 0 + ((float)0.299*R) + ((float)0.587*G) + ((float)0.113*B);
			float Cb = 128 - ((float)0.168736*R) - ((float)0.331264*G) + ((float)0.5*B);
			float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
			out->rc->data[x*width + y] = Y;
			out->gc->data[x*width + y] = Cb;
			out->bc->data[x*width + y] = Cr;
#endif

#ifdef PC
#pragma omp parallel num_threads(_CORES) for
			
			for (int x = 0; x<height; x++) {
				for (int y = 0; y<width; y++) {
					float R = in->rc->data[x*width + y];
					float G = in->gc->data[x*width + y];
					float B = in->bc->data[x*width + y];
					float Y = 0 + ((float)0.299*R) + ((float)0.587*G) + ((float)0.113*B);
					float Cb = 128 - ((float)0.168736*R) - ((float)0.331264*G) + ((float)0.5*B);
					float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
					out->rc->data[x*width + y] = Y;
					out->gc->data[x*width + y] = Cb;
					out->bc->data[x*width + y] = Cr;
#endif

#if mode == normal

		for (int y = 0; y<width; y++) {
			for (int x = 0; x<height; x++) {
				float R = in->rc->data[x*width + y];
				float G = in->gc->data[x*width + y];
				float B = in->bc->data[x*width + y];
				float Y = 0 + ((float)0.299*R) + ((float)0.587*G) + ((float)0.113*B);
				float Cb = 128 - ((float)0.168736*R) - ((float)0.331264*G) + ((float)0.5*B);
				float Cr = 128 + ((float)0.5*R) - ((float)0.418688*G) - ((float)0.081312*B);
				out->rc->data[x*width + y] = Y;
				out->gc->data[x*width + y] = Cb;
				out->bc->data[x*width + y] = Cr;
			}
		}
#endif
}

Channel* lowPass(Channel* in, Channel* out){
    // Applies a simple 3-tap low-pass filter in the X- and Y- dimensions.
    // E.g., blur
    // weights for neighboring pixels
    float a=0.25;
    float b=0.5;
    float c=0.25;

	int width = in->width; 
	int height = in->height;
     
    //out = in; TODO Is this necessary?
	for(int i=0; i<width*height; i++) out->data[i] =in->data[i];
     
    // In X
    for (int y=1; y<(width-1); y++) {
        for (int x=1; x<(height-1); x++) {
            out->data[x*width+y] = a*in->data[(x-1)*width+y]+b*in->data[x*width+y]+c*in->data[(x+1)*width+y];
        }
    }
    // In Y
    for (int y=1; y<(width-1); y++) {
        for (int x=1; x<(height-1); x++) {
            out->data[x*width+y] = a*out->data[x*width+(y-1)]+b*out->data[x*width+y]+c*out->data[x*width+(y+1)];
        }
    }
     
    return out;
}


std::vector<mVector>* motionVectorSearch(Frame* source, Frame* match, int width, int height) {
    std::vector<mVector> *motion_vectors = new std::vector<mVector>(); // empty list of ints

    float Y_weight = 0.5;
    float Cr_weight = 0.25;
    float Cb_weight = 0.25;
     
    //Window size is how much on each side of the block we search
    int window_size = 16;
    int block_size = 16;
     
    //How far from the edge we can go since we don't special case the edges
    int inset = (int) max((float)window_size, (float)block_size);
    int iter=0;
     
    for (int my=inset; my<height-(inset+window_size)+1; my+=block_size) {
      for (int mx=inset; mx<width-(inset+window_size)+1; mx+=block_size) {
         
            float best_match_sad = 1e10;
            int best_match_location[2] = {0, 0};
             
            for(int sy=my-window_size; sy<my+window_size; sy++) {
                for(int sx=mx-window_size; sx<mx+window_size; sx++) {        
                    float current_match_sad = 0;
                    // Do the SAD
                    for (int y=0; y<block_size; y++) {
                        for (int x=0; x<block_size; x++) {                
                            int match_x = mx+x;
                            int match_y = my+y;
                            int search_x = sx+x;
                            int search_y = sy+y;
                            float diff_Y = abs(match->Y->data[match_x*width+match_y] - source->Y->data[search_x*width+search_y]);
                            float diff_Cb = abs(match->Cb->data[match_x*width+match_y] - source->Cb->data[search_x*width+search_y]);
                            float diff_Cr = abs(match->Cr->data[match_x*width+match_y] - source->Cr->data[search_x*width+search_y]);
                             
                            float diff_total = Y_weight*diff_Y + Cb_weight*diff_Cb + Cr_weight*diff_Cr;
                            current_match_sad = current_match_sad + diff_total;

                        }
                    } //end SAD
                     
                    if (current_match_sad <= best_match_sad){
                        best_match_sad = current_match_sad;
                        best_match_location[0] = sx-mx;
                        best_match_location[1] = sy-my;
                    }        
                }
            }
             
            mVector v;
            v.a=best_match_location[0];
            v.b=best_match_location[1];
            motion_vectors->push_back(v);
 
        }
    }
     
    return motion_vectors;
}


Frame* computeDelta(Frame* i_frame_ycbcr, Frame* p_frame_ycbcr, std::vector<mVector>* motion_vectors){
    Frame *delta = new Frame(p_frame_ycbcr);
 
    int width = i_frame_ycbcr->width;
    int height = i_frame_ycbcr->height;
    int window_size = 16;
    int block_size = 16;
    // How far from the edge we can go since we don't special case the edges
    int inset = (int) max((float) window_size, (float)block_size);
     
    int current_block = 0;
    for(int my=inset; my<width-(inset+window_size)+1; my+=block_size) {
        for(int mx=inset; mx<height-(inset+window_size)+1; mx+=block_size) {
            int vector[2];
            vector[0]=(int)motion_vectors->at(current_block).a;
            vector[1]=(int)motion_vectors->at(current_block).b;
             
            // copy the block
                for(int y=0; y<block_size; y++) {
                    for(int x=0; x<block_size; x++) {
 
                    int src_x = mx+vector[0]+x;
                    int src_y = my+vector[1]+y;
                    int dst_x = mx+x;
                    int dst_y = my+y;
                    delta->Y->data[dst_x*width+dst_y] = delta->Y->data[dst_x*width+dst_y] - i_frame_ycbcr->Y->data[src_x*width+src_y];
                    delta->Cb->data[dst_x*width+dst_y] = delta->Cb->data[dst_x*width+dst_y] - i_frame_ycbcr->Cb->data[src_x*width+src_y];
                    delta->Cr->data[dst_x*width+dst_y] = delta->Cr->data[dst_x*width+dst_y] - i_frame_ycbcr->Cr->data[src_x*width+src_y];
                }
            }
 
            current_block = current_block + 1;
        }
    }
    return delta;
}
 

Channel* downSample(Channel* in){
	int width = in->width;
	int height = in->height;
	int w2=width/2;
	int h2=height/2;

	Channel* out = new Channel((width/2),(height/2));

	/*
		IMPROVENT: Flip x and y loop.
		SIMD: Ta fyra, skicka över.
	*/
#if mode == normal || mode == gpu
	for(int y2=0,y=0; y2<w2; y2++) {
		for (int x2 = 0, x = 0; x2<h2; x2++) {
			out->data[x2*w2+y2]= in->data[x*width+y];
            x+=2;
        }
        y+=2;
    }
#endif
#ifdef PSIMD

#pragma omp parallel num_threads(_CORES) for
	for (int y2 = 0, y = 0; y2<w2; y2++) {
		for (int x2 = 0, x = 0; x2<h2; x2++) {
			out->data[x2*w2 + y2] = in->data[x*width + y];
			x += 2;
		}
		y += 2;
	}
#endif
#ifdef PC
#pragma omp parallel num_threads(_CORES) for
	for (int x2 = 0, x = 0; x2<h2; x2++) {
		for (int y2 = 0, y = 0; y2<w2; y2++) {
			out->data[x2*w2 + y2] = in->data[x*width + y];
			y += 2;
		}
		x += 2;
	}
#endif
#ifdef CS
	for (int x2 = 0, x = 0; x2<h2; x2++) {
		for (int y2 = 0, y = 0; y2<w2; y2++) {
			out->data[x2*w2 + y2] = in->data[x*width + y];
			y += 2;
		}
		x += 2;
	}
#endif
#ifdef PCS
#pragma omp parallel num_threads(_CORES) for
	for (int x2 = 0, x = 0; x2<h2; x2++) {
		for (int y2 = 0, y = 0; y2<w2; y2++) {
			out->data[x2*w2 + y2] = in->data[x*width + y];
			y += 2;
		}
		x += 2;
	}
#endif

    return out;
}


void dct8x8(Channel* in, Channel* out){
	int width = in->width; 
	int height = in -> height;

    // 8x8 block dct on each block
    for(int i=0; i<width*height; i++) {
        in->data[i] -= 128;
        out->data[i] = 0; //zeros
    }
 
    for(int y=0; y<width; y+=8) {
        for(int x=0; x<height; x+=8) {
            dct8x8_block(&(in->data[x*width+y]),&(out->data[x*width+y]), width);
        }
    }
}


void round_block(float* in, float* out, int stride){
    float quantMatrix[8][8] ={
        {16, 11, 10, 16,  24,  40,  51,  61},
        {12, 12, 14, 19,  26,  58,  60,  55},
        {14, 13, 16, 24,  40,  57,  69,  56},
        {14, 17, 22, 29,  51,  87,  80,  62},
        {18, 22, 37, 56,  68, 109, 103,  77},
        {24, 35, 55, 64,  81, 104, 113,  92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99},
    };
 
    for(int y=0; y<8; y++) {
        for(int x=0; x<8; x++) { 
            quantMatrix[x][y] = ceil(quantMatrix[x][y]/QUALITY);
            out[x*stride+y] = (float)round(in[x*stride+y]/quantMatrix[x][y]);
        }
    }
}


void quant8x8(Channel* in, Channel* out) {
	int width = in->width;
	int height = in->height;

    for(int i=0; i<width*height; i++) {
        out->data[i]=0; //zeros
    }
    
    for (int y=0; y<width; y+=8) {
        for (int x=0; x<height; x+=8) {       
            round_block(&(in->data[x*width+y]), &(out->data[x*width+y]), width);
        }
    }
}


void dcDiff(Channel* in, Channel* out) {
	int width = in->width;
	int height = in->height;

    int number_of_dc = width*height/64;
	double* dc_values_transposed = new double[number_of_dc];
    double* dc_values = new double[number_of_dc];
 
    int iter = 0;
    for(int j=0; j<width; j+=8){
        for(int i=0; i<height; i+=8) {
            dc_values_transposed[iter] = in->data[i*width+j];
            dc_values[iter] = in->data[i*width+j];
            iter++;
        }
    }
 
    int new_w = (int) max((float)(width/8), 1);
    int new_h = (int) max((float)(height/8), 1);
     
    out->data[0] = (float)dc_values[0];
  
    double prev = 0.;
    iter = 0;
    for (int j=0; j<new_w; j++) {
        for (int i=0; i<new_h; i++) {
            out->data[iter]= (float)(dc_values[i*new_w+j] - prev);
            prev = dc_values[i*new_w+j];
            iter++;
        }
    }
	delete dc_values_transposed;
	delete dc_values;

}

void cpyBlock(float* in, float* out, int blocksize, int stride) {
    for (int j=0; j<blocksize; j++) {
        for (int i=0; i<blocksize; i++) {
            out[i*blocksize+j] = in[i*stride+j];
        }
    }
}


void zigZagOrder(Channel* in, Channel* ordered) {
	int width = in->width;
	int height = in->height;
    int zigZagIndex[64]={0,1,8,16,9,2,3,10,17,24,32,25,18,11,4,5,12,19,26,33,40,
        48,41,34,27,20,13,6,7,14,21,28,35,42,49,56,57,50,43,36,29,22,15,23,30,37,
        44,51,58,59,52,45,38,31,39,46,53,60,61,54,47,55,62,63};
     
    int blockNumber=0;
    float _block[MPEG_CONSTANT];
 
    for(int x=0; x<height; x+=8) {
        for(int y=0; y<width; y+=8) {
             cpyBlock(&(in->data[x*width+y]), _block, 8, width); //block = in(x:x+7,y:y+7);
            //Put the coefficients in zig-zag order
            float zigZagOrdered[MPEG_CONSTANT] = { 0 };
            for (int index=0; index < MPEG_CONSTANT; index++){
                zigZagOrdered[index] = _block[zigZagIndex[index]];
            }
            for (int i=0; i<MPEG_CONSTANT; i++) 
                ordered->data[blockNumber*MPEG_CONSTANT+i] = zigZagOrdered[i];
            blockNumber++;
        }
    }
}


void encode8x8(Channel* ordered, SMatrix* encoded){
	int width = encoded->height;
	int height = encoded->width;
    int num_blocks = height;
    
	for(int i=0; i<num_blocks; i++) {
		std::string block_encode[MPEG_CONSTANT];
		for (int j=0; j<MPEG_CONSTANT; j++) {
            block_encode[j]="\0"; //necessary to initialize every string position to empty string
        }

		double* block = new double[width];
        for(int y=0; y<width; y++) block[y] = ordered->data[i*width+y];
        int num_coeff = MPEG_CONSTANT; //width
        int encoded_index = 0;
        int in_zero_run = 0;
        int zero_count = 0;
 
        // Skip DC coefficient
        for(int c=1; c<num_coeff; c++){
            double coeff = block[c];
            if (coeff == 0){
                if (in_zero_run == 0){
                    zero_count = 0;
                    in_zero_run = 1;
                }
                zero_count = zero_count + 1;
            }
            else {
                if (in_zero_run == 1){
                    in_zero_run = 0;
					block_encode[encoded_index] = "Z" + std::to_string(zero_count);
                    encoded_index = encoded_index+1;
                }
				block_encode[encoded_index] = std::to_string((int)coeff);
                encoded_index = encoded_index+1;
            }
        }
 
        // If we were in a zero run at the end attach it as well.    
        if (in_zero_run == 1) {
            if (zero_count > 1) {
				block_encode[encoded_index] = "Z" + std::to_string(zero_count);
            } else {
				block_encode[encoded_index] = "0";
            }
        }
 
 
        for(int it=0; it < MPEG_CONSTANT; it++) {
			if (block_encode[it].length() > 0) 
				encoded->data[i*width+it] = new std::string(block_encode[it]);
            else
                it = MPEG_CONSTANT;
        }
		delete block;
    }
}

void gpu_convertRGBtoYCbCr(
	Image* in, Image* out,
	cl_mem* in_r_buffer, cl_mem* in_g_buffer, cl_mem* in_b_buffer, cl_mem* out_r_buffer, cl_mem* out_g_buffer, cl_mem* out_b_buffer,
	cl_context* context, cl_command_queue* queue, cl_kernel* kernel, size_t* global_dimension, cl_mem* thread_workload, cl_mem* image_size)
{

//Start data movement tiemr
#if RUN_MODE == TIME_GPU
	struct timeval starttime, endtime;
	double runtime[3] = { 0 };
	gettimeofday(&starttime, NULL);
#endif

	// write to in buffers
	clEnqueueWriteBuffer(*queue, *in_r_buffer, CL_FALSE, 0, SIZE_BYTES, in->rc->data, 0, NULL, NULL);
	clEnqueueWriteBuffer(*queue, *in_g_buffer, CL_FALSE, 0, SIZE_BYTES, in->gc->data, 0, NULL, NULL);
	clEnqueueWriteBuffer(*queue, *in_b_buffer, CL_FALSE, 0, SIZE_BYTES, in->bc->data, 0, NULL, NULL);
	
//Wait for the write to the kernel and time the duration
#if RUN_MODE == TIME_GPU
	clFinish(*queue);//Wait for kernel activity to finish

	gettimeofday(&endtime, NULL);
	runtime[0] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec) / 1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec) / 1000.0f; //in ms
#endif

	// make in args
	clSetKernelArg(*kernel, 0, sizeof(*in_r_buffer), in_r_buffer);
	clSetKernelArg(*kernel, 1, sizeof(*in_g_buffer), in_g_buffer);
	clSetKernelArg(*kernel, 2, sizeof(*in_b_buffer), in_b_buffer);

	// make out args
	clSetKernelArg(*kernel, 3, sizeof(*out_r_buffer), out_r_buffer);
	clSetKernelArg(*kernel, 4, sizeof(*out_g_buffer), out_g_buffer);
	clSetKernelArg(*kernel, 5, sizeof(*out_b_buffer), out_b_buffer);

	// Pixels per thread arg
	clSetKernelArg(*kernel, 6, sizeof(*thread_workload), thread_workload);

	//Image size arg
	clSetKernelArg(*kernel, 7, sizeof(*image_size), image_size);

#if RUN_MODE == TIME_GPU
	gettimeofday(&starttime, NULL);
#endif
	// enqueue kernel
	clEnqueueNDRangeKernel(*queue, *kernel, 1, NULL, global_dimension, NULL, 0, NULL, NULL);


#if RUN_MODE == TIME_GPU
	clFinish(*queue);
	gettimeofday(&endtime, NULL);
	runtime[1] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec) / 1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec) / 1000.0f; //in m

	//Start timer for write back
	gettimeofday(&starttime, NULL);
#endif

	// read output buffers
	clEnqueueReadBuffer(*queue, *out_r_buffer, CL_FALSE, 0, SIZE_BYTES, out->rc->data, 0, NULL, NULL);
	clEnqueueReadBuffer(*queue, *out_g_buffer, CL_FALSE, 0, SIZE_BYTES, out->gc->data, 0, NULL, NULL);
	clEnqueueReadBuffer(*queue, *out_b_buffer, CL_FALSE, 0, SIZE_BYTES, out->bc->data, 0, NULL, NULL);

	clFinish(*queue);


#if RUN_MODE == TIME_GPU
	gettimeofday(&endtime, NULL);
	runtime[2] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec) / 1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec) / 1000.0f; //in ms
	printf("GPU Data Movement %lf: \n", runtime[0]+runtime[2]);
	printf("GPU Kernel execution %lf: \n", runtime[1]);
#endif
}

int encode() {
    int end_frame = int(N_FRAMES);
    int i_frame_frequency = int(I_FRAME_FREQ);
	struct timeval starttime, endtime;
	double runtime[10] = {0};
    
    // Hardcoded paths
    string image_path =  "..\\..\\inputs\\" + string(image_name) + "\\" + image_name + ".";
	string stream_path = "..\\..\\outputs\\stream_c_" + string(image_name) + ".xml";

    xmlDocPtr stream = NULL;
 
    Image* frame_rgb = NULL;
    Image* previous_frame_rgb = NULL;
    Frame* previous_frame_lowpassed = NULL;
 
    loadImage(0, image_path, &frame_rgb);
 
    int width = frame_rgb->width;
    int height = frame_rgb->height;
    int npixels = width*height;
 
	delete frame_rgb;

	createStatsFile();
    stream = create_xml_stream(width, height, QUALITY, WINDOW_SIZE, BLOCK_SIZE);
    vector<mVector>* motion_vectors = NULL;

#if RUN_MODE == TIME_GPU || RUN_MODE == TIME_RUNTIME
	gettimeofday(&starttime, NULL);
#endif


	//opencl
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	cl_kernel kernel;
	cl_program program_cl;

	setup_cl(&device, &context, &queue);
	compile_kernel(&device, &context, &kernel, &program_cl);

	//Create in and out array
	float* in_r;
	float* in_g;
	float* in_b;

	float* out_r;
	float* out_g;
	float* out_b;

	in_r = (float*)malloc(SIZE_BYTES);
	in_g = (float*)malloc(SIZE_BYTES);
	in_b = (float*)malloc(SIZE_BYTES);
	
	out_r = (float*)malloc(SIZE_BYTES);
	out_g = (float*)malloc(SIZE_BYTES);
	out_b = (float*)malloc(SIZE_BYTES);

	
	//Set dimension size after requested threads
	int dim_size;// = SIZE * SIZE;
	int requested_threads = 32 * GPU_CORES;
	SIZE*SIZE > requested_threads ? dim_size = requested_threads : dim_size = SIZE * SIZE;
	size_t global_dimension[] = { dim_size, 0, 0 };

	//Assign workload per thread 
	int t_workload = max(ceil(SIZE * SIZE / requested_threads),1);

	// kernel init
	cl_mem in_r_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, SIZE_BYTES, NULL, NULL);
	cl_mem in_g_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, SIZE_BYTES, NULL, NULL);
	cl_mem in_b_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, SIZE_BYTES, NULL, NULL);

	cl_mem out_r_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, SIZE_BYTES, NULL, NULL);
	cl_mem out_g_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, SIZE_BYTES, NULL, NULL);
	cl_mem out_b_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, SIZE_BYTES, NULL, NULL);

	cl_mem thread_workload = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);
	clEnqueueWriteBuffer(queue, thread_workload, CL_FALSE, 0, sizeof(int), &t_workload, 0, NULL, NULL);

	int isize = SIZE * SIZE;
	cl_mem image_size = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, NULL);
	clEnqueueWriteBuffer(queue, image_size, CL_FALSE, 0, sizeof(int), &isize, 0, NULL, NULL);
	


#if RUN_MODE == TIME_GPU || RUN_MODE == TIME_RUNTIME
	gettimeofday(&endtime, NULL);
	runtime[0] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec) / 1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec) / 1000.0f; //in ms
	printf("GPU Initialization %lf: \n", runtime[0]);
#endif

    for (int frame_number = 0 ; frame_number < end_frame ; frame_number++) {
		frame_rgb = NULL;

        loadImage(frame_number, image_path, &frame_rgb);
		
        //  Convert to YCbCr
		print("Convert to YCbCr...");
        
		Image* frame_ycbcr = new Image(width, height, FULLSIZE);

		// start timer
		gettimeofday(&starttime, NULL);


#if RUN_MODE == TIME_RUNTIME
		gettimeofday(&starttime, NULL);
#endif

#define gpu 1
#define cpu 2
#define mode gpu

#if mode == gpu
		gpu_convertRGBtoYCbCr(
			frame_rgb, frame_ycbcr,
			&in_r_buffer, &in_g_buffer, &in_b_buffer, &out_r_buffer, &out_g_buffer, &out_b_buffer,
			&context, &queue, &kernel, global_dimension, &thread_workload, &image_size);
#elif mode == cpu
		convertRGBtoYCbCr(frame_rgb, frame_ycbcr);
#endif


#if RUN_MODE == TIME_RUNTIME
		gettimeofday(&endtime, NULL);
		runtime[0] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec) / 1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec) / 1000.0f; //in ms
		printf("GPU TOTAL: %lf\n", runtime[0]);
#endif
		
		// stop timer
		gettimeofday(&endtime, NULL);
		runtime[0] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms
		
        
		dump_image(frame_ycbcr, "frame_ycbcr", frame_number);
		delete frame_rgb;
 
        // We low pass filter Cb and Cr channesl
        print("Low pass filter..."); 

		gettimeofday(&starttime, NULL);
		Channel* frame_blur_cb = new Channel(width, height);
        Channel* frame_blur_cr = new Channel(width, height);
		Frame *frame_lowpassed = new Frame(width, height, FULLSIZE);
		
		lowPass(frame_ycbcr->gc, frame_blur_cb);
		lowPass(frame_ycbcr->bc, frame_blur_cr);

		//If we are just copying data around before deleting the object then maybe the object is unessecary. 
		frame_lowpassed->Y->copy(frame_ycbcr->rc);
		frame_lowpassed->Cb->copy(frame_blur_cb);
        frame_lowpassed->Cr->copy(frame_blur_cr);
		gettimeofday(&endtime, NULL);
		runtime[1] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms   
        
		dump_frame(frame_lowpassed, "frame_ycbcr_lowpass", frame_number);
		delete frame_ycbcr; 
		delete frame_blur_cb; 
		delete frame_blur_cr;

        Frame *frame_lowpassed_final = NULL;
 
		
        if (frame_number % i_frame_frequency != 0) { 
            // We have a P frame 
            // Note that in the first iteration we don't enter this branch!
            
			//Compute the motion vectors
			print("Motion Vector Search...");

			gettimeofday(&starttime, NULL);
            motion_vectors = motionVectorSearch(previous_frame_lowpassed, frame_lowpassed, frame_lowpassed->width, frame_lowpassed->height);
            gettimeofday(&endtime, NULL);
			runtime[2] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms 

			print("Compute Delta...");
			gettimeofday(&starttime, NULL);
            frame_lowpassed_final = computeDelta(previous_frame_lowpassed, frame_lowpassed, motion_vectors);
			gettimeofday(&endtime, NULL);
			runtime[3] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms 

        } else {
            // We have a I frame 
            motion_vectors = NULL;
            frame_lowpassed_final = new Frame(frame_lowpassed);
        }
		delete frame_lowpassed; frame_lowpassed=NULL;

		if (frame_number > 0) delete previous_frame_lowpassed;
        previous_frame_lowpassed = new Frame(frame_lowpassed_final); 
		
 
        // Downsample the difference
		print("Downsample...");
		
		gettimeofday(&starttime, NULL);
        Frame* frame_downsampled = new Frame(width, height, DOWNSAMPLE);
 		
        // We don't touch the Y frame
		frame_downsampled->Y->copy(frame_lowpassed_final->Y);
        Channel* frame_downsampled_cb = downSample(frame_lowpassed_final->Cb);
		frame_downsampled->Cb->copy(frame_downsampled_cb);       
        Channel* frame_downsampled_cr = downSample(frame_lowpassed_final->Cr);
		frame_downsampled->Cr->copy(frame_downsampled_cr);
        gettimeofday(&endtime, NULL);
		runtime[4] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms 

        dump_frame(frame_downsampled, "frame_downsampled", frame_number);
		delete frame_lowpassed_final; 
		delete frame_downsampled_cb; 
		delete frame_downsampled_cr; 

        // Convert to frequency domain
		print("Convert to frequency domain...");
        
		gettimeofday(&starttime, NULL);
		Frame* frame_dct = new Frame(width, height, DOWNSAMPLE);
		
        dct8x8(frame_downsampled->Y, frame_dct->Y);
        dct8x8(frame_downsampled->Cb, frame_dct->Cb);
        dct8x8(frame_downsampled->Cr, frame_dct->Cr);
        gettimeofday(&endtime, NULL);
		runtime[5] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms 

        dump_frame(frame_dct, "frame_dct", frame_number);
		delete frame_downsampled;
		
        //Quantize the data
		print("Quantize...");
		
		gettimeofday(&starttime, NULL);
        Frame* frame_quant = new Frame(width, height, DOWNSAMPLE);

        quant8x8(frame_dct->Y, frame_quant->Y);
		quant8x8(frame_dct->Cb, frame_quant->Cb);
		quant8x8(frame_dct->Cr, frame_quant->Cr);
        gettimeofday(&endtime, NULL);
		runtime[6] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms      
        
		dump_frame(frame_quant, "frame_quant", frame_number);
		delete frame_dct;

        //Extract the DC components and compute the differences
		print("Compute DC differences...");
        
		gettimeofday(&starttime, NULL);  
		Frame* frame_dc_diff = new Frame(1, (width/8)*(height/8), DCDIFF); //dealocate later
		   
        dcDiff(frame_quant->Y, frame_dc_diff->Y);
        dcDiff(frame_quant->Cb, frame_dc_diff->Cb);
        dcDiff(frame_quant->Cr, frame_dc_diff->Cr);
		gettimeofday(&endtime, NULL);
		runtime[7] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms      
         
        dump_dc_diff(frame_dc_diff, "frame_dc_diff", frame_number);

		// Zig-zag order for zero-counting
		print("Zig-zag order...");
        gettimeofday(&starttime, NULL);
		
		Frame* frame_zigzag = new Frame(MPEG_CONSTANT, width*height/MPEG_CONSTANT, ZIGZAG);
		
        zigZagOrder(frame_quant->Y, frame_zigzag->Y);
        zigZagOrder(frame_quant->Cb, frame_zigzag->Cb);
        zigZagOrder(frame_quant->Cr, frame_zigzag->Cr);
		gettimeofday(&endtime, NULL);
		runtime[8] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms  
        
		dump_zigzag(frame_zigzag, "frame_zigzag", frame_number);
		delete frame_quant;
 
        // Encode coefficients
		print("Encode coefficients...");
        
		gettimeofday(&starttime, NULL);
		FrameEncode* frame_encode = new FrameEncode(width, height, MPEG_CONSTANT);
 
        encode8x8(frame_zigzag->Y, frame_encode->Y);
        encode8x8(frame_zigzag->Cb, frame_encode->Cb);
        encode8x8(frame_zigzag->Cr, frame_encode->Cr);
		gettimeofday(&endtime, NULL);
		runtime[9] = double(endtime.tv_sec)*1000.0f + double(endtime.tv_usec)/1000.0f - double(starttime.tv_sec)*1000.0f - double(starttime.tv_usec)/1000.0f; //in ms  
       
		delete frame_zigzag;
	       
        stream_frame(stream, frame_number, motion_vectors, frame_number-1, frame_dc_diff, frame_encode);
        write_stream(stream_path, stream);

		delete frame_dc_diff;
		delete frame_encode;
 
        if (motion_vectors != NULL) {
            free(motion_vectors);
			motion_vectors = NULL;
        }

		writestats(frame_number, frame_number % i_frame_frequency, runtime);

    }
	
	free(in_r);
	free(in_g);
	free(in_b);

	free(out_r);
	free(out_g);
	free(out_b);

	cleanup(&kernel, &program_cl);
	closeStats();
	/* Uncoment to prevent visual studio output window from closing */
	system("pause");
    
	return 0;
}
 
int main(int args, char** argv) {
    encode();
    return 0;
}