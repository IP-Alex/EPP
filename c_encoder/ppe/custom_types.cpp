#include "custom_types.h"

#include <stdlib.h>
#include <vector>
#include <math.h>
#include "config.h"


float max(float a, float b) { 
    if(a > b) return a;
    return b;
}

/*
float round(float number){
    return number < 0.0 ? (float)ceil(number - 0.5) : (float)floor(number + 0.5);
}
*/

Channel::Channel(int _width, int _height){
	width = _width;
	height = _height;
	data = new float[_width*_height];
}

Channel::Channel(Channel* in){
	width = in->width;
	height = in->height;
	
	int npixels = in->width*in->height;
	data = new float[npixels];
	
	for(int i=0; i<npixels; i++)
		data[i] = in->data[i];
}

Channel::Channel(Channel& in) {
	width = in.width;
	height = in.height;
	int npixels = in.width*in.height;

	data = new float[npixels];

#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i<npixels; i++)
		this->data[i] = in.data[i];
}

Channel::Channel(Channel&& rhs) {
	width = rhs.width;
	height = rhs.height;
	data = rhs.data;
	rhs.data = nullptr;
	//printf("Move Channel\n");
}


Channel& Channel::operator=(Channel& rhs) {
	if (this != &rhs) {
		Channel temp(rhs);
		std::swap(*this, temp);
		//*this = std::move(temp);
	}
	return *this;
}

Channel& Channel::operator=(Channel&& rhs) {
	if (this != &rhs) {
		width = rhs.width;
		height = rhs.height;
		data = rhs.data;
		rhs.data = nullptr;
	}

	//printf("Move Assign Channel\n");
	return *this;
}


Channel::~Channel(){
	if(data != nullptr)
		delete[] data;
}


void Channel::copy(Channel* ch){
	int npixels = ch->width * ch->height;
	
	this->width = ch->width;
	this->height = ch->height;
	
#pragma omp parallel for num_threads(NUM_THREADS)
	for(int i=0; i<npixels; i++)
		this->data[i] = ch->data[i];
}

Image::Image(int _w, int _h, int _type){
	width = _w;
	height = _h;
	type = _type;

	rc = new Channel(_w, _h);
	if (_type == DOWNSAMPLE) {
		_w = _w/2;
		_h = _h/2;
	}
	gc = new Channel(_w, _h);
	bc = new Channel(_w, _h);
}

Image::~Image(){
	delete this->rc;
	delete this->gc;
	delete this->bc;
	this->rc=NULL;
	this->gc=NULL;
	this->bc=NULL;
}

Frame::Frame(int _w, int _h, int _type){
	width = _w;
	height = _h;
	type = _type;

	Y = new Channel(_w, _h);
	if (_type == DOWNSAMPLE) {
		_w = _w/2;
		_h = _h/2;
	}
	if (_type == DCDIFF){
		_w = (int)max(float(_w/4), 1.);
		_h = (int)max(float(_h/4), 1.);
	}
	if (_type == ZIGZAG){
		_h = _h/4;
	}

	Cb = new Channel(_w, _h);
	Cr = new Channel(_w, _h);
}

Frame::Frame(Frame* in){
	width = in->width;
	height = in->height;
	type = in->type;

	Y = new Channel(in->Y);
	Cb = new Channel(in->Cb);
	Cr = new Channel(in->Cr);

}


Frame::Frame(Frame& in) {
	width = in.width;
	height = in.height;
	type = in.type;

	Y = new Channel(in.Y);
	Cb = new Channel(in.Cb);
	Cr = new Channel(in.Cr);
}


Frame::Frame(Frame&& in) {
	width = in.width;
	height = in.height;
	type = in.type;

	Y = in.Y;
	Cb = in.Cb;
	Cr = in.Cr;

	in.Y = nullptr;
	in.Cb = nullptr;
	in.Cr = nullptr;

	//printf("Move Frame\n");
}

Frame& Frame::operator=(Frame& rhs) {
	if (this != &rhs) {
		Frame temp(rhs);
		//*this = std::move(temp);
		std::swap(*this, rhs);
	}
	return *this;
}

Frame& Frame::operator=(Frame&& rhs) {
	if (this != &rhs) {
		width = rhs.width;
		height = rhs.height;
		type = rhs.type;

		Y = rhs.Y;
		Cb = rhs.Cb;
		Cr = rhs.Cr;

		rhs.Y = nullptr;
		rhs.Cb = nullptr;
		rhs.Cr = nullptr;
	}

	//printf("Move Assign Frame\n");
	return *this;
}

Frame::~Frame(){
	if (Y != nullptr) {
		delete Y;
		Y = nullptr;
	}
	if (Cb != nullptr)
	{
		delete Cb;
		Cb = nullptr;
	}
	if (Cr != nullptr) {
		delete Cr;
		Cr = nullptr;
	}
}


SMatrix::SMatrix(int _width, int _height){
	width = _width;
	height = _height;
	//data = (std::string**)malloc(_width*_height*sizeof(std::string*));
	data = new std::string*[_width*_height];
	for(int i=0; i<_width*_height; i++) data[i] = NULL;
}

SMatrix::~SMatrix(){
	int n = width * height;
	for(int i=0; i<n; i++) delete data[i];
	delete data;
}


FrameEncode::FrameEncode(int _w, int _h, int _mpg){
	width = _w*_h/MPEG_CONSTANT;
	height = MPEG_CONSTANT;

	Y = new SMatrix(_w*_h/MPEG_CONSTANT, MPEG_CONSTANT);
	Cb = new SMatrix((_w/2)*(_h/2)/MPEG_CONSTANT, MPEG_CONSTANT);
	Cr = new SMatrix((_w/2)*(_h/2)/MPEG_CONSTANT, MPEG_CONSTANT);
}


FrameEncode::~FrameEncode(){
	delete Y;
	delete Cb;
	delete Cr;
	Y=NULL;
	Cb=NULL;
	Cr=NULL;
}