#ifndef config_h
#define config_h

#define FULLSIZE 0
#define DOWNSAMPLE 1
#define ZIGZAG 2
#define DCDIFF 3
 
#define QUALITY 1
#define WINDOW_SIZE 16
#define BLOCK_SIZE 16

#define MPEG_CONSTANT 64

#define NUM_THREADS 4


/*
--Mem leak stuff
--For more details, see https://docs.microsoft.com/en-us/visualstudio/debugger/finding-memory-leaks-using-the-crt-library?view=vs-2017
*/
#define _CRTDBG_MAP_ALLOC  
#include <crtdbg.h>  

#ifdef _DEBUG
	#define DBG_NEW new// ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#else
	#define DBG_NEW new
#endif
/*-------------------------------------------------------------------------*/

#include "test_setup.h"

#endif